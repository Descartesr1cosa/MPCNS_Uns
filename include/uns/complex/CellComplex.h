#pragma once
// Signed cell complex  cell -> face -> edge -> node.
//
//   c2f : cell  -> faces  (sign = +1 if the stored face orientation is outward)
//   f2e : face  -> edges  (ordered closed loop; sign = traversal vs. edge orientation)
//   e2n : edge  -> 2 nodes(-1 start, +1 end)
//
// d0, d1, d2 of the discrete exterior calculus are these tables applied as
//   (d0 phi)_e = sum_n e2n[e,n] phi_n,  (d1 a)_f = sum_e f2e[f,e] a_e,
//   (d2 b)_c   = sum_f c2f[c,f] b_f.
// d1*d0 = 0 and d2*d1 = 0 hold as integer identities of these tables only.
#include <array>
#include <string>
#include <vector>

#include "uns/complex/Incidence.h"

namespace uns
{
// Per-dimension bookkeeping.  Local numbering is [owned | ghost].
struct EntityInfo
{
    std::vector<GId> gid;
    std::vector<int> owner;               // owning rank
    std::vector<std::int8_t> alias_sign;  // orientation of this copy relative to the owner's
    Id n_owned = 0;
};

class CellComplex
{
public:
    Incidence e2n, f2e, c2f; // downward (primary)
    Incidence n2e, e2f, f2c; // upward (derived by build_up)
    std::array<EntityInfo, kNumDims> ent;

    std::vector<int> face_tag;      // boundary condition id, 0 = none/interior
    std::vector<int> cell_material; // uns::Material

    Id count(Dim d) const { return static_cast<Id>(ent[di(d)].gid.size()); }
    Id n_nodes() const { return count(Dim::Node); }
    Id n_edges() const { return count(Dim::Edge); }
    Id n_faces() const { return count(Dim::Face); }
    Id n_cells() const { return count(Dim::Cell); }

    void build_up();

    // Node ids along the stored loop of a face, and edge end points.
    void edge_nodes(Id e, Id &start, Id &end) const;
    void face_loop(Id f, std::vector<Id> &nodes) const;
    // Outward loop of face f as seen from cell c (loop reversed if sign = -1).
    void face_loop_outward(Id c, Id f, int sign, std::vector<Id> &nodes) const;
    void cell_nodes(Id c, std::vector<Id> &nodes) const; // unique, sorted
};

struct Mesh
{
    CellComplex cx;
    std::vector<Vec3> x; // node coordinates
    std::vector<std::string> tag_names; // index = face tag id (0 reserved)
};
} // namespace uns

namespace uns
{
// FNV-1a hash of the incidence tables (ties solution files to a complex).
std::uint64_t topology_hash(const CellComplex &cx);
} // namespace uns
