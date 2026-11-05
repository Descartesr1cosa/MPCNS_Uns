#pragma once
// Source-neutral hexahedral mesh description produced by the generators and the
// readers (box, cubed-sphere shell, legacy plot3d, gmsh, ...).  It is turned into
// a Mesh (cell complex + coordinates) by build_mesh().
#include <array>
#include <string>
#include <vector>

#include "uns/complex/CellComplex.h"

namespace uns
{
struct HexMesh
{
    std::vector<Vec3> x;
    std::vector<std::array<Id, 8>> hex; // node i at (i&1,(i>>1)&1,(i>>2)&1), right-handed
    std::vector<int> material;          // per hex, uns::Material
    struct TaggedFace { std::array<Id, 4> nodes; int tag; };
    std::vector<TaggedFace> tagged;
    std::vector<std::string> tag_names; // index = tag id, [0] unused
};

Mesh build_mesh(const HexMesh &h);

// Cartesian box [lo,hi] with nx*ny*nz cells.  Boundary tags 1..6 = xmin,xmax,ymin,ymax,zmin,zmax.
struct BoxSpec
{
    int nx = 8, ny = 8, nz = 8;
    Vec3 lo{0, 0, 0}, hi{1, 1, 1};
    bool periodic_x = false; // reserved
};
HexMesh make_box(const BoxSpec &s);

// Cubed-sphere spherical shell r in [r_in, r_out]; layers 0..n_solid-1 are Solid.
// Tags: 1 = inner surface, 2 = outer surface.  ncube cells per panel edge.
struct CubedSphereSpec
{
    int ncube = 8;
    int nr = 8;               // total radial layers (solid + fluid)
    int n_solid = 0;
    double r_in = 1.0, r_out = 10.0;
    double r_stretch = 1.0;   // geometric ratio of successive radial layer thicknesses

    // Planet + ellipsoidal domain (Mercury setup).  If r_surf > 0 the first n_solid layers are uniform between r_in and
    // r_surf, the fluid layers grow geometrically (ratio r_stretch) from r_surf to the outer surface.  If `ellipsoid` is set the
    // outer surface is the ellipsoid with bounding box [ell_lo, ell_hi] (hit by the ray from the origin) instead of the sphere r_out.
    double r_surf = 0.0;
    bool ellipsoid = false;
    Vec3 ell_lo{-12, -10, -10}, ell_hi{4, 10, 10};
};
HexMesh make_cubed_sphere_shell(const CubedSphereSpec &s);
} // namespace uns
