#include "uns/io/MeshIO.h"

#include <sstream>
#include <stdexcept>

#include "uns/io/BinContainer.h"

namespace uns::io
{
namespace
{
void put_inc(BinWriter &w, const std::string &n, const Incidence &t)
{
    w.put(n + ".ptr", t.ptr);
    w.put(n + ".idx", t.idx);
    w.put(n + ".sgn", t.sgn);
}
void get_inc(const BinReader &r, const std::string &n, Incidence &t)
{
    t.ptr = r.get<std::int64_t>(n + ".ptr");
    t.idx = r.get<Id>(n + ".idx");
    t.sgn = r.get<std::int8_t>(n + ".sgn");
}
void put_lists(BinWriter &w, const std::string &n, const CommLists &l)
{
    w.put(n + ".peers", l.peers);
    w.put(n + ".ptr", l.ptr);
    w.put(n + ".idx", l.idx);
}
void get_lists(const BinReader &r, const std::string &n, CommLists &l)
{
    l.peers = r.get<int>(n + ".peers");
    l.ptr = r.get<std::int64_t>(n + ".ptr");
    l.idx = r.get<Id>(n + ".idx");
}
} // namespace

void write_local_mesh(const std::string &path, const LocalMesh &lm)
{
    const CellComplex &cx = lm.mesh.cx;
    BinWriter w;
    std::vector<std::int64_t> meta = {1, lm.rank, lm.nranks};
    for (int d = 0; d < kNumDims; ++d) meta.push_back(cx.count(static_cast<Dim>(d)));
    for (int d = 0; d < kNumDims; ++d) meta.push_back(cx.ent[d].n_owned);
    meta.push_back(static_cast<std::int64_t>(topology_hash(cx)));
    w.put("meta", meta);
    put_inc(w, "e2n", cx.e2n);
    put_inc(w, "f2e", cx.f2e);
    put_inc(w, "c2f", cx.c2f);
    for (int d = 0; d < kNumDims; ++d)
    {
        const std::string p = "ent" + std::to_string(d);
        w.put(p + ".gid", cx.ent[d].gid);
        w.put(p + ".owner", cx.ent[d].owner);
        w.put(p + ".alias", cx.ent[d].alias_sign);
    }
    w.put("face_tag", cx.face_tag);
    w.put("cell_material", cx.cell_material);
    std::vector<double> xyz;
    xyz.reserve(lm.mesh.x.size() * 3);
    for (const Vec3 &p : lm.mesh.x) { xyz.push_back(p.x); xyz.push_back(p.y); xyz.push_back(p.z); }
    w.put("x", xyz, 3);
    std::string names;
    for (const auto &s : lm.mesh.tag_names) names += s + "\n";
    w.put("tag_names", std::vector<std::int8_t>(names.begin(), names.end()));
    for (int d = 0; d < kNumDims; ++d)
    {
        put_lists(w, "plan.send" + std::to_string(d), lm.plan.send[d]);
        put_lists(w, "plan.recv" + std::to_string(d), lm.plan.recv[d]);
    }
    put_lists(w, "plan.reduce_send", lm.plan.reduce_send);
    put_lists(w, "plan.reduce_recv", lm.plan.reduce_recv);
    w.write(path);
}

LocalMesh read_local_mesh(const std::string &path)
{
    BinReader r(path);
    LocalMesh lm;
    CellComplex &cx = lm.mesh.cx;
    const auto meta = r.get<std::int64_t>("meta");
    if (meta.size() != 12 || meta[0] != 1) throw std::runtime_error("unsupported complex file version: " + path);
    lm.rank = static_cast<int>(meta[1]);
    lm.nranks = static_cast<int>(meta[2]);
    get_inc(r, "e2n", cx.e2n);
    get_inc(r, "f2e", cx.f2e);
    get_inc(r, "c2f", cx.c2f);
    for (int d = 0; d < kNumDims; ++d)
    {
        const std::string p = "ent" + std::to_string(d);
        cx.ent[d].gid = r.get<GId>(p + ".gid");
        cx.ent[d].owner = r.get<int>(p + ".owner");
        cx.ent[d].alias_sign = r.get<std::int8_t>(p + ".alias");
        cx.ent[d].n_owned = static_cast<Id>(meta[7 + d]);
        if (static_cast<std::int64_t>(cx.ent[d].gid.size()) != meta[3 + d]) throw std::runtime_error("corrupt complex file: " + path);
    }
    cx.face_tag = r.get<int>("face_tag");
    cx.cell_material = r.get<int>("cell_material");
    cx.build_up();
    if (static_cast<std::int64_t>(topology_hash(cx)) != meta[11]) throw std::runtime_error("topology hash mismatch: " + path);
    const auto xyz = r.get<double>("x");
    lm.mesh.x.resize(xyz.size() / 3);
    for (std::size_t i = 0; i < lm.mesh.x.size(); ++i) lm.mesh.x[i] = {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
    const auto nm = r.get<std::int8_t>("tag_names");
    std::istringstream is(std::string(nm.begin(), nm.end()));
    for (std::string s; std::getline(is, s);) lm.mesh.tag_names.push_back(s);
    for (int d = 0; d < kNumDims; ++d)
    {
        get_lists(r, "plan.send" + std::to_string(d), lm.plan.send[d]);
        get_lists(r, "plan.recv" + std::to_string(d), lm.plan.recv[d]);
    }
    get_lists(r, "plan.reduce_send", lm.plan.reduce_send);
    get_lists(r, "plan.reduce_recv", lm.plan.reduce_recv);
    return lm;
}

void write_fields(const std::string &path, std::uint64_t topology_hash, const std::vector<FieldRecord> &fields)
{
    BinWriter w;
    w.put("hash", std::vector<std::int64_t>{static_cast<std::int64_t>(topology_hash)});
    std::string index;
    for (const auto &f : fields)
    {
        index += f.name + "\n";
        w.put("f." + f.name + ".info", std::vector<std::int32_t>{static_cast<std::int32_t>(f.dim), f.ncomp});
        w.put("f." + f.name, f.values, static_cast<std::uint32_t>(f.ncomp));
    }
    w.put("index", std::vector<std::int8_t>(index.begin(), index.end()));
    w.write(path);
}

std::vector<FieldRecord> read_fields(const std::string &path, std::uint64_t *topology_hash)
{
    BinReader r(path);
    if (topology_hash) *topology_hash = static_cast<std::uint64_t>(r.get<std::int64_t>("hash").at(0));
    const auto idx = r.get<std::int8_t>("index");
    std::istringstream is(std::string(idx.begin(), idx.end()));
    std::vector<FieldRecord> out;
    for (std::string s; std::getline(is, s);)
    {
        FieldRecord f;
        f.name = s;
        const auto info = r.get<std::int32_t>("f." + s + ".info");
        f.dim = static_cast<Dim>(info.at(0));
        f.ncomp = info.at(1);
        f.values = r.get<double>("f." + s);
        out.push_back(std::move(f));
    }
    return out;
}
} // namespace uns::io
