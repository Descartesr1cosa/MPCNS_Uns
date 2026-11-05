#include "uns/generate/HexMesh.h"

#include "uns/complex/ComplexBuilder.h"

namespace uns
{
Mesh build_mesh(const HexMesh &h)
{
    ComplexBuilder b(static_cast<Id>(h.x.size()));
    for (std::size_t c = 0; c < h.hex.size(); ++c)
        b.add_hex(h.hex[c].data(), h.material.empty() ? kFluid : h.material[c]);
    for (const auto &t : h.tagged) b.tag_face(std::vector<Id>(t.nodes.begin(), t.nodes.end()), t.tag);
    Mesh m;
    m.cx = b.finish();
    m.x = h.x;
    m.tag_names = h.tag_names;
    return m;
}
} // namespace uns
