# boundary/  (planned, M5)
Boundary conditions are driven by face tags (`CellComplex::face_tag`) and cell materials
(fluid / solid).  Semantics follow `legacy/.../Z4_Mercury/src/boundary`; the cubed-sphere
pole treatment disappears (variable-valence edges are ordinary edges here).
