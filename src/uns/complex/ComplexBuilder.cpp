#include "uns/complex/ComplexBuilder.h"

#include <algorithm>
#include <stdexcept>

namespace uns
{
std::size_t ComplexBuilder::VecHash::operator()(const std::vector<Id> &v) const
{
    std::size_t h = 1469598103934665603ull;
    for (Id x : v)
    {
        h ^= static_cast<std::size_t>(static_cast<std::uint32_t>(x)) + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
    }
    return h;
}

ComplexBuilder::ComplexBuilder(Id n_nodes) : n_nodes_(n_nodes) {}

bool ComplexBuilder::canonicalize(std::vector<Id> &loop)
{
    const std::size_t m = loop.size();
    const std::size_t s = std::min_element(loop.begin(), loop.end()) - loop.begin();
    std::rotate(loop.begin(), loop.begin() + s, loop.end());
    if (loop[1] > loop[m - 1])
    {
        std::reverse(loop.begin() + 1, loop.end());
        return true;
    }
    return false;
}

Id ComplexBuilder::edge_id_(Id a, Id b)
{
    const Id lo = std::min(a, b), hi = std::max(a, b);
    const std::uint64_t key = (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint32_t>(hi);
    auto it = edge_map_.find(key);
    if (it != edge_map_.end()) return it->second;
    const Id id = static_cast<Id>(edges_.size());
    edges_.emplace_back(lo, hi);
    edge_map_.emplace(key, id);
    return id;
}

Id ComplexBuilder::add_cell(const std::vector<std::vector<Id>> &outward_loops, int material)
{
    if (outward_loops.size() < 4) throw std::runtime_error("ComplexBuilder: cell with fewer than 4 faces");
    for (const auto &raw : outward_loops)
    {
        if (raw.size() < 3) throw std::runtime_error("ComplexBuilder: face with fewer than 3 nodes");
        std::vector<Id> loop = raw;
        std::vector<Id> key = raw;
        std::sort(key.begin(), key.end());
        const bool flipped = canonicalize(loop);
        auto it = face_map_.find(key);
        Id fid;
        if (it == face_map_.end())
        {
            fid = static_cast<Id>(face_map_.size());
            face_map_.emplace(std::move(key), fid);
            const std::size_t m = loop.size();
            for (std::size_t i = 0; i < m; ++i)
            {
                const Id a = loop[i], b = loop[(i + 1) % m];
                f2e_.push(edge_id_(a, b), a < b ? +1 : -1);
            }
            f2e_.close_row();
        }
        else
            fid = it->second;
        c2f_.push(fid, flipped ? -1 : +1);
    }
    c2f_.close_row();
    cell_material_.push_back(material);
    return static_cast<Id>(cell_material_.size()) - 1;
}

Id ComplexBuilder::add_hex(const Id v[8], int material)
{
    static const int L[6][4] = {{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}};
    std::vector<std::vector<Id>> loops(6);
    for (int f = 0; f < 6; ++f)
        for (int k = 0; k < 4; ++k) loops[f].push_back(v[L[f][k]]);
    return add_cell(loops, material);
}

void ComplexBuilder::tag_face(const std::vector<Id> &nodes, int tag)
{
    std::vector<Id> key = nodes;
    std::sort(key.begin(), key.end());
    pending_tags_[key] = tag;
}

Id ComplexBuilder::face_id(const std::vector<Id> &nodes) const
{
    std::vector<Id> key = nodes;
    std::sort(key.begin(), key.end());
    auto it = face_map_.find(key);
    return it == face_map_.end() ? -1 : it->second;
}

CellComplex ComplexBuilder::finish()
{
    CellComplex cx;
    const Id ne = static_cast<Id>(edges_.size());
    for (Id e = 0; e < ne; ++e)
    {
        cx.e2n.push(edges_[e].first, -1);
        cx.e2n.push(edges_[e].second, +1);
        cx.e2n.close_row();
    }
    cx.f2e = std::move(f2e_);
    cx.c2f = std::move(c2f_);
    cx.cell_material = std::move(cell_material_);
    const Id nf = cx.f2e.rows(), nc = cx.c2f.rows();
    const Id counts[4] = {n_nodes_, ne, nf, nc};
    for (int d = 0; d < kNumDims; ++d)
    {
        EntityInfo &ei = cx.ent[d];
        ei.gid.resize(counts[d]);
        for (Id i = 0; i < counts[d]; ++i) ei.gid[i] = i;
        ei.owner.assign(counts[d], 0);
        ei.alias_sign.assign(counts[d], 1);
        ei.n_owned = counts[d];
    }
    cx.face_tag.assign(nf, 0);
    for (const auto &kv : face_map_)
    {
        auto it = pending_tags_.find(kv.first);
        if (it != pending_tags_.end()) cx.face_tag[kv.second] = it->second;
    }
    cx.build_up();
    return cx;
}
} // namespace uns
