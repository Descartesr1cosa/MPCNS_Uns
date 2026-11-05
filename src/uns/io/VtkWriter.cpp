#include "uns/io/VtkWriter.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace uns::io
{
namespace
{
// Collects typed arrays and writes them as <DataArray> elements whose payloads
// live in one <AppendedData encoding="raw"> block (UInt64 length headers).
class Appended
{
public:
    template <class T> std::string add(const char *vtk_type, const std::string &name, int ncomp, const std::vector<T> &v)
    {
        std::ostringstream s;
        s << "<DataArray type=\"" << vtk_type << "\" Name=\"" << name << "\" NumberOfComponents=\"" << ncomp
          << "\" format=\"appended\" offset=\"" << bytes_ << "\"/>";
        const std::uint64_t n = v.size() * sizeof(T);
        blocks_.emplace_back(reinterpret_cast<const char *>(v.data()), n);
        bytes_ += 8 + n;
        return s.str();
    }
    void write(std::ostream &os) const
    {
        os << "<AppendedData encoding=\"raw\">\n_";
        for (const auto &b : blocks_)
        {
            os.write(reinterpret_cast<const char *>(&b.second), 8);
            os.write(b.first, static_cast<std::streamsize>(b.second));
        }
        os << "\n</AppendedData>\n";
    }
private:
    std::uint64_t bytes_ = 0;
    std::vector<std::pair<const char *, std::uint64_t>> blocks_;
};
} // namespace

void write_vtu(const std::string &path, const Mesh &m, Id n_out, const std::vector<VtkField> &cell_data,
               const std::vector<VtkField> &point_data)
{
    const CellComplex &cx = m.cx;
    std::vector<double> pts;
    pts.reserve(m.x.size() * 3);
    for (const Vec3 &p : m.x) { pts.push_back(p.x); pts.push_back(p.y); pts.push_back(p.z); }
    std::vector<std::int64_t> conn, offs, faces, foffs;
    std::vector<std::uint8_t> types(n_out, 42);
    std::vector<Id> loop, cn;
    for (Id c = 0; c < n_out; ++c)
    {
        cx.cell_nodes(c, cn);
        for (Id v : cn) conn.push_back(v);
        offs.push_back(static_cast<std::int64_t>(conn.size()));
        faces.push_back(cx.c2f.count(c));
        for (auto kf = cx.c2f.begin(c); kf < cx.c2f.end(c); ++kf)
        {
            cx.face_loop_outward(c, cx.c2f.idx[kf], cx.c2f.sgn[kf], loop);
            faces.push_back(static_cast<std::int64_t>(loop.size()));
            for (Id v : loop) faces.push_back(v);
        }
        foffs.push_back(static_cast<std::int64_t>(faces.size()));
    }
    Appended ap;
    std::ostringstream hdr;
    hdr << "<?xml version=\"1.0\"?>\n<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
        << "<UnstructuredGrid>\n<Piece NumberOfPoints=\"" << m.x.size() << "\" NumberOfCells=\"" << n_out << "\">\n";
    hdr << "<Points>" << ap.add("Float64", "Points", 3, pts) << "</Points>\n<Cells>\n";
    hdr << ap.add("Int64", "connectivity", 1, conn) << "\n" << ap.add("Int64", "offsets", 1, offs) << "\n"
        << ap.add("UInt8", "types", 1, types) << "\n" << ap.add("Int64", "faces", 1, faces) << "\n"
        << ap.add("Int64", "faceoffsets", 1, foffs) << "\n</Cells>\n";
    auto emit = [&](const char *tag, const std::vector<VtkField> &fs) {
        hdr << "<" << tag << ">\n";
        for (const auto &f : fs) hdr << ap.add("Float64", f.name, f.ncomp, f.values) << "\n";
        hdr << "</" << tag << ">\n";
    };
    emit("CellData", cell_data);
    emit("PointData", point_data);
    hdr << "</Piece>\n</UnstructuredGrid>\n";
    std::ofstream os(path, std::ios::binary);
    if (!os) throw std::runtime_error("cannot write " + path);
    os << hdr.str();
    ap.write(os);
    os << "</VTKFile>\n";
}

void write_pvtu(const std::string &path, const std::vector<std::string> &pieces, const std::vector<VtkField> &cell_data,
                const std::vector<VtkField> &point_data)
{
    std::ofstream os(path);
    if (!os) throw std::runtime_error("cannot write " + path);
    os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"PUnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
       << "<PUnstructuredGrid GhostLevel=\"0\">\n<PPoints><PDataArray type=\"Float64\" NumberOfComponents=\"3\"/></PPoints>\n";
    auto emit = [&](const char *tag, const std::vector<VtkField> &fs) {
        os << "<" << tag << ">\n";
        for (const auto &f : fs) os << "<PDataArray type=\"Float64\" Name=\"" << f.name << "\" NumberOfComponents=\"" << f.ncomp << "\"/>\n";
        os << "</" << tag << ">\n";
    };
    emit("PCellData", cell_data);
    emit("PPointData", point_data);
    for (const auto &p : pieces) os << "<Piece Source=\"" << p << "\"/>\n";
    os << "</PUnstructuredGrid>\n</VTKFile>\n";
}

void write_pvd(const std::string &path, const std::vector<std::pair<double, std::string>> &entries)
{
    std::ofstream os(path);
    if (!os) throw std::runtime_error("cannot write " + path);
    os << "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n";
    for (const auto &e : entries) os << "<DataSet timestep=\"" << e.first << "\" part=\"0\" file=\"" << e.second << "\"/>\n";
    os << "</Collection>\n</VTKFile>\n";
}
} // namespace uns::io
