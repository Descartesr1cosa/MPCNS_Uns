#pragma once
// Builds a CellComplex from cells given as lists of outward-oriented face loops.
// Cell type is only relevant here; everything downstream is type-agnostic.
#include <unordered_map>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns
{
class ComplexBuilder
{
public:
    explicit ComplexBuilder(Id n_nodes);

    // Generic polyhedron: faces are node loops, counter-clockwise seen from outside.
    Id add_cell(const std::vector<std::vector<Id>> &outward_loops, int material = kFluid);
    // Hexahedron, node i at (i&1, (i>>1)&1, (i>>2)&1) of a right-handed cell.
    Id add_hex(const Id v[8], int material = kFluid);
    // Tag a (boundary) face given any ordering of its nodes; applied by finish().
    void tag_face(const std::vector<Id> &nodes, int tag);

    // Face id (creation order == final numbering) of the face with this node set, -1 if absent.
    Id face_id(const std::vector<Id> &nodes) const;
    Id n_faces() const { return static_cast<Id>(face_map_.size()); }

    CellComplex finish();

    // Rotate/reverse a loop to the canonical form (start at the smallest node,
    // second node smaller than the last).  Returns true if reversed.
    static bool canonicalize(std::vector<Id> &loop);

private:
    struct VecHash { std::size_t operator()(const std::vector<Id> &v) const; };
    Id n_nodes_;
    std::unordered_map<std::uint64_t, Id> edge_map_;
    std::unordered_map<std::vector<Id>, Id, VecHash> face_map_;
    std::vector<std::pair<Id, Id>> edges_; // start < end
    Incidence f2e_, c2f_;
    std::vector<int> cell_material_;
    std::unordered_map<std::vector<Id>, int, VecHash> pending_tags_;
    Id edge_id_(Id a, Id b);
};
} // namespace uns
