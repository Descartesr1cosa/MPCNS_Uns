#pragma once
// Basic scalar/vector types shared by every uns module.
#include <array>
#include <cmath>
#include <cstdint>

namespace uns
{
using Id = std::int32_t;  // local entity index
using GId = std::int64_t; // global entity index

enum class Dim : int
{
    Node = 0,
    Edge = 1,
    Face = 2,
    Cell = 3
};
constexpr int kNumDims = 4;
constexpr int di(Dim d) { return static_cast<int>(d); }

struct Vec3
{
    double x = 0, y = 0, z = 0;
    Vec3() = default;
    Vec3(double a, double b, double c) : x(a), y(b), z(c) {}
    Vec3 &operator+=(const Vec3 &o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3 &operator-=(const Vec3 &o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vec3 &operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};
inline Vec3 operator+(Vec3 a, const Vec3 &b) { return a += b; }
inline Vec3 operator-(Vec3 a, const Vec3 &b) { return a -= b; }
inline Vec3 operator*(Vec3 a, double s) { return a *= s; }
inline Vec3 operator*(double s, Vec3 a) { return a *= s; }
inline double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3 &a, const Vec3 &b)
{
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double norm(const Vec3 &a) { return std::sqrt(dot(a, a)); }

// Material of a cell: decides the cost model and the physics that is solved.
enum Material : int
{
    kFluid = 0, // fluid + electromagnetics
    kSolid = 1  // electromagnetics only (conducting shell / planet)
};
} // namespace uns
