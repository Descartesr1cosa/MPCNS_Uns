#include "uns/io/BinContainer.h"

#include <bit>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace uns::io
{
static_assert(std::endian::native == std::endian::little, "UNSBIN is little endian");
namespace
{
const char kMagic[8] = {'U', 'N', 'S', 'B', 'I', 'N', 0, 1};
std::size_t elem_size(DType t)
{
    switch (t)
    {
    case DType::I8: return 1;
    case DType::I32: case DType::F32: return 4;
    default: return 8;
    }
}
struct File
{
    std::FILE *f;
    File(const std::string &p, const char *mode) : f(std::fopen(p.c_str(), mode))
    {
        if (!f) throw std::runtime_error("cannot open '" + p + "' (" + mode + ")");
    }
    ~File() { if (f) std::fclose(f); }
};
void wr(std::FILE *f, const void *p, std::size_t n)
{
    if (n && std::fwrite(p, 1, n, f) != n) throw std::runtime_error("write failed");
}
void rd(std::FILE *f, void *p, std::size_t n)
{
    if (n && std::fread(p, 1, n, f) != n) throw std::runtime_error("read failed / truncated file");
}
} // namespace

void BinWriter::put_raw(const std::string &name, DType t, std::uint32_t cols, std::uint64_t rows, const void *p, std::size_t nbytes)
{
    Record r;
    r.type = t; r.cols = cols; r.rows = rows;
    r.bytes.assign(static_cast<const char *>(p), static_cast<const char *>(p) + nbytes);
    recs_.emplace_back(name, std::move(r));
}

void BinWriter::write(const std::string &path) const
{
    File fl(path, "wb");
    wr(fl.f, kMagic, 8);
    const std::uint64_t n = recs_.size();
    wr(fl.f, &n, 8);
    static const char zeros[8] = {0};
    for (const auto &kv : recs_)
    {
        const std::uint32_t nl = static_cast<std::uint32_t>(kv.first.size());
        wr(fl.f, &nl, 4);
        wr(fl.f, kv.first.data(), nl);
        const std::uint8_t t[4] = {static_cast<std::uint8_t>(kv.second.type), 0, 0, 0};
        wr(fl.f, t, 4);
        wr(fl.f, &kv.second.cols, 4);
        wr(fl.f, &kv.second.rows, 8);
        wr(fl.f, kv.second.bytes.data(), kv.second.bytes.size());
        wr(fl.f, zeros, (8 - kv.second.bytes.size() % 8) % 8);
    }
}

BinReader::BinReader(const std::string &path)
{
    File fl(path, "rb");
    char magic[8];
    rd(fl.f, magic, 8);
    if (std::memcmp(magic, kMagic, 8) != 0) throw std::runtime_error("'" + path + "' is not an UNSBIN file");
    std::uint64_t n;
    rd(fl.f, &n, 8);
    for (std::uint64_t i = 0; i < n; ++i)
    {
        std::uint32_t nl;
        rd(fl.f, &nl, 4);
        std::string name(nl, '\0');
        rd(fl.f, name.data(), nl);
        std::uint8_t t[4];
        rd(fl.f, t, 4);
        Record r;
        r.type = static_cast<DType>(t[0]);
        rd(fl.f, &r.cols, 4);
        rd(fl.f, &r.rows, 8);
        r.bytes.resize(r.rows * r.cols * elem_size(r.type));
        rd(fl.f, r.bytes.data(), r.bytes.size());
        char pad[8];
        rd(fl.f, pad, (8 - r.bytes.size() % 8) % 8);
        recs_.emplace(std::move(name), std::move(r));
    }
}

const Record &BinReader::find_(const std::string &name, DType t) const
{
    auto it = recs_.find(name);
    if (it == recs_.end()) throw std::runtime_error("UNSBIN: missing record '" + name + "'");
    if (it->second.type != t) throw std::runtime_error("UNSBIN: record '" + name + "' has a different data type");
    return it->second;
}

std::vector<std::string> BinReader::names() const
{
    std::vector<std::string> v;
    for (const auto &kv : recs_) v.push_back(kv.first);
    return v;
}
} // namespace uns::io
