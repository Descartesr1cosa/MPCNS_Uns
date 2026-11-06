# File formats

## UNSBIN (`*.unsb`) — mesh, sub-complex, plans, entity fields
Self-describing little-endian container: `magic | n_records | {name, dtype, cols, rows, payload}`.
Complex files hold the incidence CSRs (`e2n/f2e/c2f` `.ptr/.idx/.sgn`), per-dimension `gid/owner/alias`,
`x`, `face_tag`, `cell_material`, `tag_names` and the communication plan (`plan.send<d>`, `plan.recv<d>`,
`plan.reduce_send/recv`), plus a topology hash that is verified on read.  Field files (`write_fields`)
store entity-numbered arrays and the hash of the complex they belong to (restart and coarse solution).
Any language can read it (see `tools/check_vtu.py` for a numpy-style parser of the VTU twin).

## VTU / PVTU / PVD — visualisation (ParaView, VisIt)
VTK XML unstructured grid, binary appended raw data, `header_type=UInt64`.  Every cell is a
`VTK_POLYHEDRON` (type 42) with outward face loops, so hexahedra, hanging-node polyhedra and later cell types
need no special handling.  One `.vtu` per rank (owned cells only), `.pvtu` master, `.pvd` for time series.
`tools/check_vtu.py` recomputes polyhedron volumes from the file and compares with the written `volume` array.

### Tecplot (binary .plt)
Tecplot cannot import VTK polyhedral cells (VTK_POLYHEDRON, type 42; face-based connectivity), which is what the solvers write ("missing
connectivity").  Use the converter `tools/vtu2tec/vtu2tec.py` (Tecplot binary only, format `#!TDV112`), controlled by a yml file instead of
command-line options: edit `tools/vtu2tec/config.yml` (`input_dir`, `output_dir`, `name`, `step`, `vars`, `double`), then
`python3 tools/vtu2tec/vtu2tec.py [other.yml]`.
It merges the rank pieces of a step (duplicate interface nodes are merged by exact coordinates) and writes a zone `hex` (`FEBRICK`) with all
hexahedral cells (every unrefined Mars / Mercury mesh -> a single zone) and, for refined meshes, a zone `poly` (`FETETRAHEDRON`) in which each
hanging-node polyhedron is cut into tetrahedra about its centroid and face centres (cell data copied to the sub-tetrahedra).  Cell data stay
cell-centred, vectors become `name_x/_y/_z`, values are 32-bit unless `double: true`, `SOLUTIONTIME` comes from `solution.pvd`; `step: all` writes
one file per step (`NAME_<step>.plt`, `STRANDID` set so Tecplot animates them).
Checked without Tecplot: the files parse with the independent pure-python package `tecio` (zones, variable locations, connectivity); all bricks
have positive Jacobian, all tetrahedra positive volume, total volume of a refined test mesh within 0.04 % of the solver's volumes.  Loading in a
real Tecplot installation has not been tried here.

## Coarse / restart solution (`*_fields.unsb`)
UNSBIN written by `io::write_fields`, entity-numbered in the *global* numbering of the complex named by the
stored topology hash (the refinement checks it).  Record names used by `mpcns_pre task=refine`:
`B` (face, 1: flux in the canonical face orientation), `rho` and `p` (cell, 1), `B_cell` and `J_cell` (cell, 3),
and any other cell field (`ncomp` values per cell, prolonged conservatively).  The fine file has the same layout
in the global numbering of the refined complex.
