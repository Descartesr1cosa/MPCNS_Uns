# Conventions (orientation, numbering, geometry)

## Complex
- Entities of every dimension are numbered contiguously; local numbering is `[owned | ghost]`.
- `e2n`: edge → nodes, sign −1 = start, +1 = end.  `f2e`: face → edges in loop order; sign +1 if the
  loop traverses the edge from start to end.  `c2f`: cell → faces; sign +1 if the stored face loop is
  counter-clockwise seen from outside (outward normal by the right-hand rule).
- Canonical orientation (derived from global node ids, identical on every rank without communication):
  edge start = smaller node id; face loop starts at the smallest node and its second node is smaller than
  its last node.  Hence alias signs are +1 for meshes produced by this code; the tables still carry them so
  imported/periodic meshes work (tests flip random orientations).
- `d0 d1 d2` of the DEC are the signed tables applied as sums; `d1·d0 = 0` and `d2·d1 = 0` are integer
  identities checked by `count_nonzero_products`.

## Hex cell (builder input)
Node `i` sits at `(i&1, (i>>1)&1, (i>>2)&1)` of a right-handed cell (positive volume).  Face loops (outward):
`x-:[0,4,6,2] x+:[1,3,7,5] y-:[0,1,5,4] y+:[2,6,7,3] z-:[0,2,3,1] z+:[4,5,7,6]`.

## Ownership
`owner(entity) = min rank over the ranks of all cells incident to it`.  Sub-complex = owned cells + one
vertex-neighbour layer.  Hence every edge/face of an owned cell sees all of its incident cells locally.

## Geometry (kept separate from topology)
- Face centre = area-weighted centroid of the triangle fan; face area vector from triangles (c_f, v_i, v_i+1).
  Warped faces are therefore well defined, the cell surface is watertight and the volume
  `V = Σ_tri a_tri·(c_f − x_c)/3` does not depend on the apex `x_c`.
- Cell centre = tet-volume-weighted centroid.  `w_{e,c}` = share of cell c in the dual face of edge e.
- Hodge stars are diagonal and strictly positive for a valid mesh: `⋆1 = |A_dual|/|e|`,
  `⋆2 = (Σ_c |c_f − x_c|)/|A_f|`, `⋆0` = nodal dual volume, `⋆3 = 1/V`.
