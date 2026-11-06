# vtu2tec — VTU (ParaView) -> Tecplot binary

Tecplot cannot read the solvers' VTK polyhedral cells.  `vtu2tec.py` merges the rank pieces of an output step and writes a Tecplot binary
`.plt` (TDV112): hexahedral meshes -> one `FEBRICK` zone; refined meshes with hanging-node polyhedra -> an extra `FETETRAHEDRON` zone
(polyhedra cut into tetrahedra, cell data copied).  Cell data stay cell-centred, vectors become `name_x/_y/_z`.

1. edit `config.yml` (input/output directory, file name, step, variables)
2. `python3 tools/vtu2tec/vtu2tec.py`            (or `python3 tools/vtu2tec/vtu2tec.py my_other.yml`)

Needs numpy (PyYAML optional: a minimal built-in parser is used without it).  See docs/FileFormats.md for the checks that were done.
