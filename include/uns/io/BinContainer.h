#pragma once
// Minimal self-describing binary container ("UNSBIN").  Used for meshes, per-rank
// sub-complexes, plans and entity-numbered field data (restart / coarse solution).
//
// Layout (little endian):  magic[8] | u64 n_records | records...
// record:  u32 name_len | name | u8 dtype | u8[3] pad | u32 cols | u64 rows | payload (padded to 8 bytes)
#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace uns::io
{
enum class DType : std::uint8_t { I8 = 0, I32 = 1, I64 = 2, F32 = 3, F64 = 4 };

template <class T> struct DTypeOf;
template <> struct DTypeOf<std::int8_t>  { static constexpr DType v = DType::I8; };
template <> struct DTypeOf<std::int32_t> { static constexpr DType v = DType::I32; };
template <> struct DTypeOf<std::int64_t> { static constexpr DType v = DType::I64; };
template <> struct DTypeOf<float>        { static constexpr DType v = DType::F32; };
template <> struct DTypeOf<double>       { static constexpr DType v = DType::F64; };

struct Record
{
    DType type = DType::F64;
    std::uint32_t cols = 1;
    std::uint64_t rows = 0;
    std::vector<char> bytes;
};

class BinWriter
{
public:
    template <class T> void put(const std::string &name, const std::vector<T> &v, std::uint32_t cols = 1)
    {
        put_raw(name, DTypeOf<T>::v, cols, v.size() / cols, v.data(), v.size() * sizeof(T));
    }
    void write(const std::string &path) const; // throws on I/O error
private:
    void put_raw(const std::string &name, DType t, std::uint32_t cols, std::uint64_t rows, const void *p, std::size_t nbytes);
    std::vector<std::pair<std::string, Record>> recs_;
};

class BinReader
{
public:
    explicit BinReader(const std::string &path); // throws on I/O/format error
    bool has(const std::string &name) const { return recs_.count(name) != 0; }
    template <class T> std::vector<T> get(const std::string &name, std::uint32_t *cols = nullptr) const
    {
        const Record &r = find_(name, DTypeOf<T>::v);
        if (cols) *cols = r.cols;
        std::vector<T> v(r.bytes.size() / sizeof(T));
        std::copy(r.bytes.begin(), r.bytes.end(), reinterpret_cast<char *>(v.data()));
        return v;
    }
    std::vector<std::string> names() const;
private:
    const Record &find_(const std::string &name, DType t) const;
    std::map<std::string, Record> recs_;
};
} // namespace uns::io
