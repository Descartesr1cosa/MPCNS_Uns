#pragma once
// Internal: lattice frames of the coarse hexahedra and canonical node keys.
//
// A coarse hex gets a right-handed local frame (node q at bits (q&1, q>>1&1, q>>2&1)).
// A lattice point (a,b,c) in [0,N]^3 of that frame is identified by the coarse entity it lies on
// (node / edge / face / cell interior) and its coordinates in that entity's *canonical* frame, so
// points shared by several coarse cells get identical keys without comparing coordinates.
#include <array>
#include <cstdint>
#include <map>
#include <unordered_map>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns::refine_detail
{
struct NodeKey
{
    std::uint64_t a = 0, b = 0;
    bool operator==(const NodeKey &o) const { return a == o.a && b == o.b; }
};
struct NodeKeyHash
{
    std::size_t operator()(const NodeKey &k) const
    {
        std::uint64_t h = k.a * 0x9E3779B97F4A7C15ull;
        h ^= (k.b + 0x7F4A7C15ull + (h << 6) + (h >> 2)) * 0xBF58476D1CE4E5B9ull;
        return static_cast<std::size_t>(h ^ (h >> 29));
    }
};

struct SideFrame
{
    Id face = -1;
    int sa = 0, ta = 1; // free-axis index (0=u,1=v) along which canonical s / t grow
    int ou = 0, ov = 0; // lattice corner (0/1) of the canonical origin
};

struct CellFrame
{
    Id node[8];
    SideFrame side[6]; // side = 2*axis + (0: low, 1: high)
};

class Lattice
{
public:
    Lattice(const Mesh &coarse, int max_level, bool spherical_snap);
    int N() const { return N_; }
    int max_level() const { return Lmax_; }
    const CellFrame &frame(Id c) const { return fr_[c]; }
    NodeKey key(Id cell, int a, int b, int c) const;
    Vec3 position(Id cell, int a, int b, int c) const;
    // canonical in-plane coordinates of lattice point (cu,cv) on a side of the cell
    void side_st(Id cell, int side, int cu, int cv, int &s, int &t) const;
private:
    const Mesh &m_;
    int Lmax_, N_;
    bool snap_;
    std::vector<CellFrame> fr_;
    std::unordered_map<std::uint64_t, Id> edge_of_;
};
} // namespace uns::refine_detail
