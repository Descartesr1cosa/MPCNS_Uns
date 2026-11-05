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

### Tecplot
Not verified here.  Tecplot's native `.szplt` needs the TecIO library (not installed); loading polyhedral
VTU in Tecplot 360 depends on its version.  Plan (M8): add a Tecplot-friendly writer after testing against
a real Tecplot installation — candidates: `.szplt` via TecIO (supports polyhedra), or EnSight Gold binary
(nfaced), plus an option that decomposes polyhedra into tetrahedra for tools without polyhedron support.
