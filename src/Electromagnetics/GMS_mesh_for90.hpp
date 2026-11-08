#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace for90_mom2 {
namespace constants {
using dp = double;
static constexpr double pi = 3.141592653589;
static constexpr double eps0 = 8.8541878176e-12;
static constexpr double mu0 = 4.0 * pi * 1.0e-7;
static constexpr double c0 = 2.997924580003e8;
static constexpr double eta0 = 3.767303134622e2;
static constexpr double radtodeg = 180.0 / pi;
static constexpr double degtorad = pi / 180.0;
static constexpr int prd_none = 1;
static constexpr int prd_2d = 2;
} // namespace constants

namespace aux {
inline std::vector<int> find_smallest(const std::vector<std::complex<double>>& val, int dim, int nval) {
    std::vector<int> ind(nval, 0);
    std::vector<std::complex<double>> tmpval = val;
    double mval = 0.0;
    for (auto const& v : val) {
        mval = std::max(mval, std::abs(v));
    }
    for (int n = 0; n < nval; ++n) {
        auto it = std::min_element(tmpval.begin(), tmpval.end(),
            [](std::complex<double> a, std::complex<double> b) { return std::abs(a) < std::abs(b); });
        int loc = static_cast<int>(std::distance(tmpval.begin(), it));
        if (loc < 0 || loc >= static_cast<int>(tmpval.size())) {
            throw std::runtime_error("find_smallest: invalid index");
        }
        tmpval[loc] = std::complex<double>(mval, 0.0);
        ind[n] = loc;
    }
    return ind;
}

inline std::string getext(const std::string& filename) {
    std::string s = filename;
    while (!s.empty() && s.back() == '.') {
        s.pop_back();
    }
    std::size_t n = s.find_last_not_of(" \t\r\n");
    std::string trimmed = (n == std::string::npos) ? "" : s.substr(0, n + 1);
    if (trimmed.size() < 3) {
        return std::string();
    }
    return trimmed.substr(trimmed.size() - 3);
}

inline void replace_ext(std::string& filename, const std::string& ext) {
    std::size_t n = filename.find_last_of('.');
    if (n == std::string::npos) {
        filename += "." + ext;
    } else {
        filename = filename.substr(0, n + 1) + ext;
    }
}

inline double linterp(double start, double end, double t) { return start + t * (end - start); }

inline void get_matrix(const std::string& filename, Eigen::MatrixXd& mat) {
    std::ifstream in(filename.c_str());
    if (!in) throw std::runtime_error("Could not open matrix data file: " + filename);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::vector<double> row; double v;
        while (iss >> v) row.push_back(v);
        if (!row.empty()) rows.push_back(row);
    }
    if (rows.empty()) return;
    const std::size_t nrows = rows.size();
    const std::size_t ncols = rows.front().size();
    mat.resize(static_cast<int>(nrows), static_cast<int>(ncols));
    for (std::size_t i = 0; i < nrows; ++i) {
        for (std::size_t j = 0; j < ncols; ++j) {
            mat(static_cast<int>(i), static_cast<int>(j)) = rows[i][j];
        }
    }
}

inline std::complex<double> get_refind(const std::string& filename, double wl) {
    std::ifstream in(filename.c_str());
    if (!in) throw std::runtime_error("Could not open refractive index lookup file: " + filename);
    double wl1 = 0.0, nr1 = 0.0, ni1 = 0.0;
    double wl2 = 0.0, nr2 = 0.0, ni2 = 0.0;
    if (!(in >> wl1 >> nr1 >> ni1)) {
        throw std::runtime_error("Invalid refractive index data file: " + filename);
    }
    if (wl < wl1) {
        return std::complex<double>(1.0, 0.0);
    }
    while (in >> wl2 >> nr2 >> ni2) {
        if (wl1 <= wl && wl <= wl2 && wl1 != wl2) {
            double t = (wl - wl1) / (wl2 - wl1);
            return std::complex<double>(linterp(nr1, nr2, t), linterp(ni1, ni2, t));
        }
        if (wl == wl1) return std::complex<double>(nr1, ni1);
        if (wl == wl2) return std::complex<double>(nr2, ni2);
        if (wl1 > wl2) throw std::runtime_error("Refractive index data disordered!");
        wl1 = wl2; nr1 = nr2; ni1 = ni2;
    }
    return std::complex<double>(1.0, 0.0);
}

inline int find_closest(double val, const std::vector<double>& list) {
    if (list.empty()) throw std::runtime_error("find_closest: empty list");
    int ind = 0;
    double dist = std::abs(val - list[0]);
    for (std::size_t i = 1; i < list.size(); ++i) {
        double d = std::abs(val - list[i]);
        if (d < dist) {
            ind = static_cast<int>(i);
            dist = d;
        }
    }
    return ind;
}

inline void write_data(const std::string& filename, const Eigen::MatrixXd& data) {
    std::ofstream out(filename.c_str(), std::ios::out);
    if (!out) throw std::runtime_error("Could not open output file: " + filename);
    for (int i = 0; i < data.rows(); ++i) {
        for (int j = 0; j < data.cols(); ++j) {
            out << std::scientific << std::setw(15) << std::setprecision(3) << data(i, j);
        }
        out << '\n';
    }
}

inline void write_cdata(const std::string& filename, const Eigen::MatrixXcd& data) {
    std::ofstream out(filename.c_str(), std::ios::out);
    if (!out) throw std::runtime_error("Could not open output file: " + filename);
    for (int i = 0; i < data.rows(); ++i) {
        for (int j = 0; j < data.cols(); ++j) {
            out << std::scientific << std::setw(15) << std::setprecision(3) << data(i, j).real()
                << "+" << std::scientific << std::setw(15) << std::setprecision(3) << data(i, j).imag() << "i";
        }
        out << '\n';
    }
}

inline std::complex<double> drude(std::complex<double> omega, double epsinf, double omegap, double gamma) {
    return std::complex<double>(epsinf, 0.0) - (omegap * omegap) / (omega * (omega + std::complex<double>(0.0, 1.0) * gamma));
}

inline Eigen::Vector3d get_dir(double theta, double phi) {
    Eigen::Vector3d d;
    d << std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta);
    return d;
}

inline Eigen::Vector3d get_pol(double theta, double phi, double psi) {
    Eigen::Vector3d p;
    p << std::cos(psi) * std::sin(phi) - std::sin(psi) * std::cos(theta) * std::cos(phi),
         -std::cos(psi) * std::cos(phi) - std::sin(psi) * std::cos(theta) * std::sin(phi),
         std::sin(psi) * std::sin(theta);
    return p;
}
} // namespace aux

namespace linalg {
inline std::complex<double> dotc(const std::vector<std::complex<double>>& v1,
                                const std::vector<std::complex<double>>& v2) {
    if (v1.size() != v2.size()) throw std::runtime_error("dotc: size mismatch");
    std::complex<double> s(0.0, 0.0);
    for (std::size_t i = 0; i < v1.size(); ++i) s += std::conj(v1[i]) * v2[i];
    return s;
}

inline double dotr(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) {
    if (v1.size() != v2.size()) throw std::runtime_error("dotr: size mismatch");
    return v1.dot(v2);
}

inline std::complex<double> dotrc(const Eigen::VectorXd& v1, const std::vector<std::complex<double>>& v2) {
    if (v1.size() != v2.size()) throw std::runtime_error("dotrc: size mismatch");
    std::complex<double> s(0.0, 0.0);
    for (std::size_t i = 0; i < v1.size(); ++i) s += v1[i] * v2[i];
    return s;
}

inline double normc(const std::vector<std::complex<double>>& v) {
    return std::sqrt(std::real(dotc(v, v)));
}

inline double normr(const Eigen::VectorXd& v) { return v.norm(); }

inline Eigen::Vector3d crossr(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) {
    return v1.cross(v2);
}

inline Eigen::Vector3cd crossc(const Eigen::Vector3cd& v1, const Eigen::Vector3cd& v2) {
    Eigen::Vector3cd r;
    r << v1(1) * v2(2) - v1(2) * v2(1),
         v1(2) * v2(0) - v1(0) * v2(2),
         v1(0) * v2(1) - v1(1) * v2(0);
    return r;
}

inline void vcrosscr(const std::vector<Eigen::Vector3cd>& v1,
                    const std::vector<Eigen::Vector3d>& v2,
                    std::vector<Eigen::Vector3cd>& res) {
    if (v1.size() != v2.size()) throw std::runtime_error("vcrosscr: size mismatch");
    res.resize(v1.size());
    for (std::size_t i = 0; i < v1.size(); ++i) {
        res[i] = crossc(v1[i], v2[i].cast<std::complex<double>>());
    }
}

inline Eigen::Matrix3d invmat33r(const Eigen::Matrix3d& m) {
    double det = (m.col(0).cross(m.col(1))).dot(m.col(2));
    if (std::abs(det) < 1e-30) throw std::runtime_error("invmat33r: singular matrix");
    Eigen::Matrix3d mi;
    mi.col(0) = m.col(1).cross(m.col(2));
    mi.col(1) = m.col(2).cross(m.col(0));
    mi.col(2) = m.col(0).cross(m.col(1));
    return mi / det;
}

inline Eigen::MatrixXcd matrix_inverse(const Eigen::MatrixXcd& A) {
    Eigen::MatrixXcd Ainv = A;
    Eigen::FullPivLU<Eigen::MatrixXcd> lu(Ainv);
    if (lu.rank() < Ainv.rows()) throw std::runtime_error("matrix_inverse: singular matrix");
    return lu.inverse();
}

inline void solve_linsys(Eigen::MatrixXcd& mat, Eigen::VectorXcd& src) {
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(mat);
    src = lu.solve(src);
}

inline void solve_multi_linsys(Eigen::MatrixXcd& mat, Eigen::MatrixXcd& src) {
    Eigen::PartialPivLU<Eigen::MatrixXcd> lu(mat);
    src = lu.solve(src);
}

inline void matrix_eigenvalues_gen(Eigen::MatrixXcd& A, Eigen::MatrixXcd& B,
                                   Eigen::VectorXcd& eigval, Eigen::MatrixXcd& eigvec) {
    Eigen::GeneralizedEigenSolver<Eigen::MatrixXcd> ges(A, B);
    if (ges.info() != Eigen::Success) throw std::runtime_error("matrix eigenvalue solution failed");
    eigval = ges.eigenvalues();
    eigvec = ges.eigenvectors();
}

inline void matrix_eigenvalues(Eigen::MatrixXcd& A, Eigen::VectorXcd& eigval, Eigen::MatrixXcd& eigvec) {
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(A);
    if (ces.info() != Eigen::Success) throw std::runtime_error("matrix eigenvalue solution failed");
    eigval = ces.eigenvalues();
    eigvec = ces.eigenvectors();
}

inline void matmul_blocklt(Eigen::MatrixXcd& B, Eigen::MatrixXcd& A, int dim) {
    if (A.rows() < 2 * dim || A.cols() < 2 * dim) throw std::runtime_error("matmul_blocklt: bad dimensions");
    Eigen::MatrixXcd A0 = A.block(0, 0, dim, dim);
    Eigen::MatrixXcd A1 = A.block(dim, 0, dim, dim);
    Eigen::MatrixXcd A2 = A.block(0, dim, dim, dim);
    Eigen::MatrixXcd A3 = A.block(dim, dim, dim, dim);
    A0 = B * A0;
    A1 = B * A1;
    A2 = B * A2;
    A3 = B * A3;
    A.block(0, 0, dim, dim) = A0;
    A.block(dim, 0, dim, dim) = A1;
    A.block(0, dim, dim, dim) = A2;
    A.block(dim, dim, dim, dim) = A3;
}

inline void matmul_blockrt(Eigen::MatrixXcd& A, const Eigen::MatrixXcd& B, int dim) {
    if (A.rows() < 2 * dim || A.cols() < 2 * dim) throw std::runtime_error("matmul_blockrt: bad dimensions");
    Eigen::MatrixXcd A0 = A.block(0, 0, dim, dim);
    Eigen::MatrixXcd A1 = A.block(dim, 0, dim, dim);
    Eigen::MatrixXcd A2 = A.block(0, dim, dim, dim);
    Eigen::MatrixXcd A3 = A.block(dim, dim, dim, dim);
    A0 = A0 * B;
    A1 = A1 * B;
    A2 = A2 * B;
    A3 = A3 * B;
    A.block(0, 0, dim, dim) = A0;
    A.block(dim, 0, dim, dim) = A1;
    A.block(0, dim, dim, dim) = A2;
    A.block(dim, dim, dim, dim) = A3;
}

inline Eigen::Matrix3d matrix_rz(double angle) {
    Eigen::Matrix3d rm;
    double st = std::sin(angle);
    double ct = std::cos(angle);
    rm << ct, st, 0.0,
          -st, ct, 0.0,
          0.0, 0.0, 1.0;
    return rm;
}

inline Eigen::Vector3d rotate_vector(const Eigen::Vector3d& v, const Eigen::Vector3d& axis, double angle) {
    double z = axis.dot(v);
    Eigen::Vector3d b1 = v - z * axis;
    double x = b1.norm();
    if (x < 1e-30) return v;
    b1 /= x;
    Eigen::Vector3d b2 = axis.cross(b1);
    double x2 = std::cos(angle) * x;
    double y2 = std::sin(angle) * x;
    return x2 * b1 + y2 * b2 + z * axis;
}
} // namespace linalg

namespace mesh {

static constexpr int meshver = 4;
static constexpr int mesh_bnd_none = 1;
static constexpr int mesh_bnd_xplane = 2;
static constexpr int mesh_bnd_yplane = 3;
static constexpr int mesh_bnd_zplane = 4;
static constexpr int mesh_bnd_prdx1 = 5;
static constexpr int mesh_bnd_prdx2 = 6;
static constexpr int mesh_bnd_prdy1 = 7;
static constexpr int mesh_bnd_prdy2 = 8;
static constexpr int mesh_bnd_rz1 = 9;
static constexpr int mesh_bnd_rz2 = 10;

struct Node {
    Eigen::Vector3d p = Eigen::Vector3d::Zero();
    int parent_index = -1;
};

struct Face {
    Eigen::Vector3d n = Eigen::Vector3d::Zero();
    Eigen::Vector3d cp = Eigen::Vector3d::Zero();
    Eigen::Matrix3d s = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d m = Eigen::Matrix3d::Zero();
    std::array<int, 3> node_indices{{-1, -1, -1}};
    std::array<int, 3> edge_indices{{-1, -1, -1}};
    double area = 0.0;
    double pd = 0.0;
    int id = 0;
    int parent_index = -1;
};

struct Line {
    std::array<int, 2> node_indices{{-1, -1}};
    int id = 0;
};

struct Edge {
    std::array<int, 2> node_indices{{-1, -1}};
    std::array<int, 2> bnode_indices{{-1, -1}};
    std::array<int, 2> face_indices{{-1, -1}};
    double length = 0.0;
    int bnd = mesh_bnd_none;
    int parent_index = -1;
    int couple_index = -1;
    std::vector<std::array<int, 2>> child_indices;
};

struct SolidFace {
    std::array<int, 3> node_indices{{-1, -1, -1}};
    std::array<int, 2> solid_indices{{-1, -1}};
    std::array<int, 2> bnode_indices{{-1, -1}};
    int face_index = -1;
    double area = 0.0;
};

struct Solid {
    std::array<int, 4> node_indices{{-1, -1, -1, -1}};
    std::array<int, 4> solid_face_indices{{-1, -1, -1, -1}};
    double volume = 0.0;
    int id = 0;
};

struct MeshContainer {
    std::vector<Node> nodes;
    std::vector<Face> faces;
    std::vector<Line> lines;
    std::vector<Edge> edges;
    std::vector<Solid> solids;
    std::vector<SolidFace> solid_faces;
    int nnodes = 0;
    int nfaces = 0;
    int nlines = 0;
    int nedges = 0;
    int nsolids = 0;
    int nsolid_faces = 0;
    double avelen = 0.0;
};

inline bool has_mesh_bnd(const MeshContainer& mesh, int bnd) {
    for (const auto& edge : mesh.edges) {
        if (edge.bnd == bnd) return true;
    }
    return false;
}

inline int adjacent_face(const MeshContainer& mesh, int faceind, int edgeind) {
    int gedgeind = mesh.faces[faceind].edge_indices[edgeind];
    if (mesh.edges[gedgeind].face_indices[0] == faceind) {
        return mesh.edges[gedgeind].face_indices[1];
    }
    if (mesh.edges[gedgeind].face_indices[1] == faceind) {
        return mesh.edges[gedgeind].face_indices[0];
    }
    throw std::runtime_error("Could not determine adjacent face!");
}

inline int local_edge_index(const MeshContainer& mesh, int faceind, int gedgeind) {
    for (int i = 0; i < 3; ++i) {
        if (mesh.faces[faceind].edge_indices[i] == gedgeind) return i + 1;
    }
    throw std::runtime_error("Invalid edge index request (nlsurf.f90)");
}

inline double average_edge_length(const MeshContainer& mesh) {
    double s = 0.0;
    for (const auto& edge : mesh.edges) {
        s += (mesh.nodes[edge.node_indices[0]].p - mesh.nodes[edge.node_indices[1]].p).norm();
    }
    if (mesh.edges.empty()) return 0.0;
    return s / static_cast<double>(mesh.edges.size());
}

inline void classify_edges(MeshContainer& mesh, int id, int bnd) {
    int nbnd = 0;
    for (auto& edge : mesh.edges) {
        for (const auto& line : mesh.lines) {
            if (line.id == id && ((edge.node_indices[0] == line.node_indices[0] && edge.node_indices[1] == line.node_indices[1]) ||
                                  (edge.node_indices[0] == line.node_indices[1] && edge.node_indices[1] == line.node_indices[0]))) {
                edge.bnd = bnd;
                ++nbnd;
            }
        }
    }
    std::cout << " Classified " << nbnd << " boundary edges\n";
}

inline bool cmp_real_pairs(const std::array<double, 2>& a, const std::array<double, 2>& b, double eps) {
    return ((std::abs(a[0] - b[0]) < eps && std::abs(a[1] - b[1]) < eps) ||
            (std::abs(a[0] - b[1]) < eps && std::abs(a[1] - b[0]) < eps));
}

inline void determine_edge_couples(MeshContainer& mesh, double eps) {
    std::cout << "Determining edge couples.\n";
    int numx = 0, numy = 0;
    for (auto& edge : mesh.edges) edge.couple_index = -1;
    for (std::size_t n = 0; n < mesh.edges.size(); ++n) {
        if (mesh.edges[n].bnd != mesh_bnd_prdx1) continue;
        for (std::size_t m = 0; m < mesh.edges.size(); ++m) {
            if (mesh.edges[m].bnd != mesh_bnd_prdx2) continue;
            auto a1 = mesh.nodes[mesh.edges[n].node_indices[0]].p;
            auto a2 = mesh.nodes[mesh.edges[n].node_indices[1]].p;
            auto b1 = mesh.nodes[mesh.edges[m].node_indices[0]].p;
            auto b2 = mesh.nodes[mesh.edges[m].node_indices[1]].p;
            std::array<double, 2> aa{{a1(1), a2(1)}};
            std::array<double, 2> bb{{b1(1), b2(1)}};
            std::array<double, 2> ac{{a1(2), a2(2)}};
            std::array<double, 2> bc{{b1(2), b2(2)}};
            std::array<double, 2> ad{{a1(2), a2(2)}};
            std::array<double, 2> bd{{b1(2), b2(2)}};
            if (cmp_real_pairs(aa, bb, eps) && cmp_real_pairs(ac, bc, eps)) {
                mesh.edges[n].couple_index = static_cast<int>(m);
                mesh.edges[m].couple_index = static_cast<int>(n);
                ++numx;
            }
        }
    }
    for (std::size_t n = 0; n < mesh.edges.size(); ++n) {
        if (mesh.edges[n].bnd != mesh_bnd_prdy1) continue;
        for (std::size_t m = 0; m < mesh.edges.size(); ++m) {
            if (mesh.edges[m].bnd != mesh_bnd_prdy2) continue;
            auto a1 = mesh.nodes[mesh.edges[n].node_indices[0]].p;
            auto a2 = mesh.nodes[mesh.edges[n].node_indices[1]].p;
            auto b1 = mesh.nodes[mesh.edges[m].node_indices[0]].p;
            auto b2 = mesh.nodes[mesh.edges[m].node_indices[1]].p;
            std::array<double, 2> aa{{a1(0), a2(0)}};
            std::array<double, 2> bb{{b1(0), b2(0)}};
            std::array<double, 2> ac{{a1(2), a2(2)}};
            std::array<double, 2> bc{{b1(2), b2(2)}};
            if (cmp_real_pairs(aa, bb, eps) && cmp_real_pairs(ac, bc, eps)) {
                mesh.edges[n].couple_index = static_cast<int>(m);
                mesh.edges[m].couple_index = static_cast<int>(n);
                ++numy;
            }
        }
    }
    for (std::size_t n = 0; n < mesh.edges.size(); ++n) {
        if (mesh.edges[n].bnd != mesh_bnd_prdx1 && mesh.edges[n].bnd != mesh_bnd_prdy1) continue;
        int m = mesh.edges[n].couple_index;
        if (m < 0) continue;
        if (mesh.edges[n].face_indices[0] != -1 && mesh.edges[m].face_indices[0] != -1) {
            std::swap(mesh.edges[m].face_indices[0], mesh.edges[m].face_indices[1]);
            std::swap(mesh.edges[m].bnode_indices[0], mesh.edges[m].bnode_indices[1]);
        }
        if (mesh.edges[n].face_indices[1] != -1 && mesh.edges[m].face_indices[1] != -1) {
            std::swap(mesh.edges[m].face_indices[0], mesh.edges[m].face_indices[1]);
            std::swap(mesh.edges[m].bnode_indices[0], mesh.edges[m].bnode_indices[1]);
        }
    }
    std::cout << "Found " << numx << " x-couples and " << numy << " y-couples.\n";
}

inline void invert_faces(MeshContainer& mesh, int id) {
    for (auto& face : mesh.faces) {
        if (face.id == id) {
            std::swap(face.node_indices[0], face.node_indices[2]);
        }
    }
}

inline MeshContainer extract_submesh(const MeshContainer& mesh, const std::vector<int>& ids,
                                     const std::vector<int>* vol_ids = nullptr) {
    MeshContainer submesh;
    std::vector<int> nmask(mesh.nodes.size(), 0);
    std::vector<int> local_node_indices(mesh.nodes.size(), -1);
    int nfaces = 0, nsolids = 0;
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        if (std::find(ids.begin(), ids.end(), mesh.faces[i].id) != ids.end()) {
            ++nfaces;
            for (int j = 0; j < 3; ++j) nmask[mesh.faces[i].node_indices[j]] = 1;
        }
    }
    submesh.nfaces = nfaces;
    submesh.faces.resize(static_cast<std::size_t>(nfaces));
    if (vol_ids) {
        for (std::size_t i = 0; i < mesh.solids.size(); ++i) {
            if (std::find(vol_ids->begin(), vol_ids->end(), mesh.solids[i].id) != vol_ids->end()) {
                ++nsolids;
                for (int j = 0; j < 4; ++j) nmask[mesh.solids[i].node_indices[j]] = 1;
            }
        }
        submesh.nsolids = nsolids;
        submesh.solids.resize(static_cast<std::size_t>(nsolids));
    } else {
        submesh.nsolids = 0;
    }
    int counter = 0;
    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {
        if (nmask[i]) {
            submesh.nodes.push_back(mesh.nodes[i]);
            local_node_indices[i] = counter++;
        }
    }
    submesh.nnodes = static_cast<int>(submesh.nodes.size());
    counter = 0;
    for (std::size_t i = 0; i < mesh.faces.size(); ++i) {
        if (std::find(ids.begin(), ids.end(), mesh.faces[i].id) != ids.end()) {
            Face f = mesh.faces[i];
            for (int j = 0; j < 3; ++j) {
                f.node_indices[j] = local_node_indices[mesh.faces[i].node_indices[j]];
            }
            f.parent_index = static_cast<int>(i);
            submesh.faces[counter++] = f;
        }
    }
    if (vol_ids) {
        counter = 0;
        for (std::size_t i = 0; i < mesh.solids.size(); ++i) {
            if (std::find(vol_ids->begin(), vol_ids->end(), mesh.solids[i].id) != vol_ids->end()) {
                Solid s = mesh.solids[i];
                for (int j = 0; j < 4; ++j) {
                    s.node_indices[j] = local_node_indices[mesh.solids[i].node_indices[j]];
                }
                s.id = mesh.solids[i].id;
                submesh.solids[counter++] = s;
            }
        }
    }
    submesh.nlines = 0;
    submesh.nedges = 0;
    submesh.nsolid_faces = 0;
    return submesh;
}

inline void submesh_edge_connectivity(MeshContainer& mesh,
                                     std::vector<MeshContainer>& submeshes) {
    std::vector<std::vector<int>> nchild(mesh.edges.size(), std::vector<int>(1, 0));
    for (std::size_t n = 0; n < mesh.edges.size(); ++n) {
        std::vector<std::array<int, 2>> child_indices;
        for (std::size_t m = 0; m < submeshes.size(); ++m) {
            for (std::size_t l = 0; l < submeshes[m].edges.size(); ++l) {
                std::array<int, 2> node_indices{{submeshes[m].nodes[submeshes[m].edges[l].node_indices[0]].parent_index,
                                                 submeshes[m].nodes[submeshes[m].edges[l].node_indices[1]].parent_index}};
                if (node_indices[0] == mesh.edges[n].node_indices[0] && node_indices[1] == mesh.edges[n].node_indices[1] ||
                    node_indices[0] == mesh.edges[n].node_indices[1] && node_indices[1] == mesh.edges[n].node_indices[0]) {
                    submeshes[m].edges[l].parent_index = static_cast<int>(n);
                    child_indices.push_back({{static_cast<int>(m), static_cast<int>(l)}});
                }
            }
        }
        if (child_indices.size() < 2) {
            throw std::runtime_error("Submesh decomposition is inconsistent (invalid child edges)!");
        }
        mesh.edges[n].child_indices = child_indices;
    }
    for (std::size_t ge1 = 0; ge1 < mesh.edges.size(); ++ge1) {
        int ge2 = mesh.edges[ge1].couple_index;
        for (const auto& item : mesh.edges[ge1].child_indices) {
            int sm1 = item[0], le1 = item[1];
            submeshes[sm1].edges[le1].bnd = mesh.edges[ge1].bnd;
            if (ge2 != -1) {
                for (const auto& item2 : mesh.edges[ge2].child_indices) {
                    int sm2 = item2[0], le2 = item2[1];
                    if (sm2 == sm1) {
                        submeshes[sm1].edges[le1].couple_index = le2;
                        submeshes[sm2].edges[le2].couple_index = le1;
                        break;
                    }
                }
            }
        }
    }
}

inline void orient_basis(const MeshContainer& mesh, std::vector<MeshContainer>& submeshes) {
    std::cout << "Orienting submesh basis functions\n";
    for (std::size_t k = 0; k < submeshes.size(); ++k) {
        for (std::size_t n = 0; n < submeshes[k].edges.size(); ++n) {
            if (submeshes[k].edges[n].bnd != mesh_bnd_prdx1 && submeshes[k].edges[n].bnd != mesh_bnd_prdy1) continue;
            int m = submeshes[k].edges[n].couple_index;
            if (m < 0) continue;
            if (submeshes[k].edges[n].face_indices[0] != -1 && submeshes[k].edges[m].face_indices[0] != -1) {
                std::swap(submeshes[k].edges[m].face_indices[0], submeshes[k].edges[m].face_indices[1]);
                std::swap(submeshes[k].edges[m].bnode_indices[0], submeshes[k].edges[m].bnode_indices[1]);
            }
            if (submeshes[k].edges[n].face_indices[1] != -1 && submeshes[k].edges[m].face_indices[1] != -1) {
                std::swap(submeshes[k].edges[m].face_indices[0], submeshes[k].edges[m].face_indices[1]);
                std::swap(submeshes[k].edges[m].bnode_indices[0], submeshes[k].edges[m].bnode_indices[1]);
            }
        }
    }
}

inline void get_bnd_edges(const MeshContainer& mesh, int bnd, int& nse, std::vector<int>& edgeind) {
    edgeind.clear();
    nse = 0;
    for (std::size_t i = 0; i < mesh.edges.size(); ++i) {
        if (mesh.edges[i].bnd == bnd) {
            edgeind.push_back(static_cast<int>(i));
            ++nse;
        }
    }
}

inline std::vector<std::vector<int>> get_mesh_element_lines(const std::string& filename) {
    std::ifstream in(filename.c_str());
    if (!in) throw std::runtime_error("Could not open mesh file " + filename + "!");
    std::string lineid;
    int nelements = 0;
    std::vector<std::vector<int>> lines;
    while (in >> lineid) {
        if (lineid == "$Elements") {
            in >> nelements;
            lines.resize(static_cast<std::size_t>(nelements), std::vector<int>(2, 0));
            for (int n = 0; n < nelements; ++n) {
                int element_number = 0;
                in >> element_number >> lines[n][0] >> lines[n][1];
            }
            in >> lineid;
            if (lineid != "$EndElements") throw std::runtime_error("Could not find $EndElements specifier!");
            break;
        }
    }
    return lines;
}

inline MeshContainer load_mesh_gmsh(const std::string& filename) {
    std::ifstream in(filename.c_str());
    if (!in) throw std::runtime_error("Could not open mesh file!");
    MeshContainer mesh;
    std::string lineid;
    std::vector<std::vector<int>> element_types = get_mesh_element_lines(filename);
    int nnodes = 0, nelements = 0;
    std::vector<int> element_data(20, 0);
    while (in >> lineid) {
        if (lineid == "$MeshFormat") {
            std::string mshver;
            in >> mshver;
            if (mshver != "2.1" && mshver != "2.2") throw std::runtime_error("Mesh version is unsupported!");
            in >> lineid;
            if (lineid != "$EndMeshFormat") throw std::runtime_error("Could not find $MeshFormatNodes specifier!");
        } else if (lineid == "$Nodes") {
            in >> nnodes;
            mesh.nodes.resize(static_cast<std::size_t>(nnodes));
            for (int n = 0; n < nnodes; ++n) {
                int node_number = 0;
                Eigen::Vector3d np;
                in >> node_number >> np(0) >> np(1) >> np(2);
                mesh.nodes[node_number - 1].p = np;
            }
            in >> lineid;
            if (lineid != "$EndNodes") throw std::runtime_error("Could not find $EndNodes specifier!");
            mesh.nnodes = nnodes;
        } else if (lineid == "$Elements") {
            in >> nelements;
            mesh.nfaces = 0; mesh.nlines = 0; mesh.nsolids = 0;
            for (const auto& item : element_types) {
                if (item[0] == 2) ++mesh.nfaces;
                if (item[0] == 1) ++mesh.nlines;
                if (item[0] == 4) ++mesh.nsolids;
            }
            mesh.faces.resize(static_cast<std::size_t>(mesh.nfaces));
            if (mesh.nlines != 0) mesh.lines.resize(static_cast<std::size_t>(mesh.nlines));
            if (mesh.nsolids != 0) mesh.solids.resize(static_cast<std::size_t>(mesh.nsolids));
            int cface = 0, cline = 0, csolid = 0;
            for (int n = 0; n < nelements; ++n) {
                int ntags = element_types[n][1];
                if (element_types[n][0] == 2) {
                    ++cface;
                    std::vector<int> data(ntags + 6, 0);
                    for (int i = 0; i < ntags + 6; ++i) in >> data[i];
                    mesh.faces[cface - 1].id = data[3];
                    for (int j = 0; j < 3; ++j) mesh.faces[cface - 1].node_indices[j] = data[4 + ntags + j];
                } else if (element_types[n][0] == 1) {
                    ++cline;
                    std::vector<int> data(ntags + 5, 0);
                    for (int i = 0; i < ntags + 5; ++i) in >> data[i];
                    mesh.lines[cline - 1].id = data[3];
                    mesh.lines[cline - 1].node_indices[0] = data[4 + ntags];
                    mesh.lines[cline - 1].node_indices[1] = data[5 + ntags];
                } else if (element_types[n][0] == 4) {
                    ++csolid;
                    std::vector<int> data(ntags + 7, 0);
                    for (int i = 0; i < ntags + 7; ++i) in >> data[i];
                    mesh.solids[csolid - 1].id = data[3];
                    for (int j = 0; j < 4; ++j) mesh.solids[csolid - 1].node_indices[j] = data[ntags + 4 + j];
                } else {
                    int skip = 0; in >> skip;
                }
            }
            in >> lineid;
            if (lineid != "$EndElements") throw std::runtime_error("Could not find $EndElements specifier!");
        }
    }
    return mesh;
}

inline MeshContainer load_mesh_neutral(const std::string& filename) {
    std::ifstream in(filename.c_str());
    if (!in) throw std::runtime_error("Could not open mesh file!");
    MeshContainer mesh;
    in >> mesh.nnodes;
    mesh.nodes.resize(static_cast<std::size_t>(mesh.nnodes));
    for (int n = 0; n < mesh.nnodes; ++n) {
        Eigen::Vector3d p;
        in >> p(0) >> p(1) >> p(2);
        mesh.nodes[n].p = p;
    }
    in >> mesh.nsolids;
    if (mesh.nsolids > 0) {
        mesh.solids.resize(static_cast<std::size_t>(mesh.nsolids));
        for (int n = 0; n < mesh.nsolids; ++n) {
            mesh.solids[n].id = 1;
            for (int j = 0; j < 4; ++j) in >> mesh.solids[n].node_indices[j];
        }
    }
    in >> mesh.nfaces;
    mesh.faces.resize(static_cast<std::size_t>(mesh.nfaces));
    for (int n = 0; n < mesh.nfaces; ++n) {
        int bndid = 0;
        mesh.faces[n].id = 1;
        in >> bndid;
        for (int j = 0; j < 3; ++j) in >> mesh.faces[n].node_indices[j];
    }
    mesh.nlines = 0;
    return mesh;
}

inline MeshContainer load_mesh(const std::string& filename) {
    std::string ext = aux::getext(filename);
    MeshContainer mesh;
    if (ext == "msh") {
        mesh = load_mesh_gmsh(filename);
    } else if (ext == "nmf") {
        mesh = load_mesh_neutral(filename);
    } else {
        throw std::runtime_error("Unrecognized mesh file extension!");
    }
    mesh.nsolid_faces = 0;
    std::cout << "Mesh file loaded successfully.\n";
    std::cout << " - Read " << mesh.nnodes << " nodes\n";
    std::cout << " - Read " << mesh.nfaces << " faces\n";
    std::cout << " - Read " << mesh.nlines << " lines\n";
    std::cout << " - Read " << mesh.nsolids << " solids\n";
    return mesh;
}

inline void delete_mesh(MeshContainer& mesh) {
    mesh.nodes.clear();
    mesh.faces.clear();
    mesh.lines.clear();
    mesh.edges.clear();
    mesh.solids.clear();
    mesh.solid_faces.clear();
}

inline void order(int& p, int& q) {
    if (p > q) std::swap(p, q);
}

inline void sort(std::vector<int>& array, int n) {
    std::sort(array.begin(), array.begin() + n);
}

inline void swap(int& a, int& b) { std::swap(a, b); }

inline void sort3(std::array<int, 3>& array) {
    std::sort(array.begin(), array.end());
}

inline void append_integer(std::vector<int>& array, int element) { array.push_back(element); }

inline bool test_inclusion(const std::vector<int>& array, int element) {
    return std::find(array.begin(), array.end(), element) != array.end();
}

inline bool test_inclusion2(const std::vector<int>& array, int element) {
    return std::find(array.begin(), array.end(), element) != array.end();
}

inline bool cmp_pairs(const std::array<int, 2>& pair1, const std::array<int, 2>& pair2) {
    return ((pair1[0] == pair2[0] && pair1[1] == pair2[1]) ||
            (pair1[0] == pair2[1] && pair1[1] == pair2[0]));
}

inline bool cmp_triplets(const std::array<int, 3>& triplet1, const std::array<int, 3>& triplet2) {
    static const std::array<std::array<int, 3>, 6> indices = {
        std::array<int, 3>{{1, 2, 3}},
        std::array<int, 3>{{1, 3, 2}},
        std::array<int, 3>{{2, 1, 3}},
        std::array<int, 3>{{2, 3, 1}},
        std::array<int, 3>{{3, 1, 2}},
        std::array<int, 3>{{3, 2, 1}} };
    for (const auto& idx : indices) {
        if (triplet1[0] == triplet2[idx[0] - 1] && triplet1[1] == triplet2[idx[1] - 1] && triplet1[2] == triplet2[idx[2] - 1]) return true;
        if (triplet1[0] == triplet2[idx[0] - 1] && triplet1[1] == triplet2[idx[1] - 1] && triplet1[2] == triplet2[idx[2] - 1]) return true;
    }
    return false;
}

inline int get_free_index(const std::array<int, 3>& indices, const std::array<int, 2>& dublet) {
    if (cmp_pairs({{indices[0], indices[1]}}, dublet)) return indices[2];
    if (cmp_pairs({{indices[1], indices[2]}}, dublet)) return indices[0];
    if (cmp_pairs({{indices[0], indices[2]}}, dublet)) return indices[1];
    throw std::runtime_error("Erroneous dublet comparison!");
}

inline void scale_mesh(MeshContainer& mesh, double scale) {
    for (auto& node : mesh.nodes) node.p *= scale;
}

inline int indexrot3(int n) {
    int index = n;
    if (index > 3) index = index - ((index - 1) / 3) * 3;
    return index;
}

inline int indexrot4(int n) {
    int index = n;
    if (index > 4) index = index - ((index - 1) / 4) * 4;
    return index;
}

inline void build_mesh(MeshContainer& mesh, double scale) {
    std::cout << "Building mesh data\n";
    scale_mesh(mesh, scale);
    int nedges = static_cast<int>(mesh.faces.size()) * 3;
    std::vector<Edge> tmpedges(static_cast<std::size_t>(nedges));
    for (auto& e : tmpedges) {
        e.node_indices = {{-1, -1}};
        e.bnode_indices = {{-1, -1}};
        e.face_indices = {{-1, -1}};
    }
    int cedge = 0;
    for (std::size_t n = 0; n < mesh.faces.size(); ++n) {
        for (int m = 0; m < 3; ++m) {
            std::array<int, 2> pair{{mesh.faces[n].node_indices[m], mesh.faces[n].node_indices[indexrot3(m + 1)]}};
            bool found_edge = false;
            for (int l = 0; l < cedge; ++l) {
                if (cmp_pairs(tmpedges[l].node_indices, pair)) {
                    found_edge = true;
                    tmpedges[l].face_indices[1] = static_cast<int>(n);
                    tmpedges[l].bnode_indices[1] = mesh.faces[n].node_indices[indexrot3(m + 2)];
                    mesh.faces[n].edge_indices[m] = l;
                    break;
                }
            }
            if (!found_edge) {
                ++cedge;
                tmpedges[cedge - 1].node_indices = pair;
                tmpedges[cedge - 1].face_indices[0] = static_cast<int>(n);
                tmpedges[cedge - 1].bnode_indices[0] = mesh.faces[n].node_indices[indexrot3(m + 2)];
                mesh.faces[n].edge_indices[m] = cedge - 1;
            }
        }
    }
    mesh.nedges = cedge;
    mesh.edges.resize(static_cast<std::size_t>(cedge));
    for (int n = 0; n < cedge; ++n) {
        mesh.edges[n].node_indices = tmpedges[n].node_indices;
        mesh.edges[n].bnode_indices = tmpedges[n].bnode_indices;
        mesh.edges[n].face_indices = tmpedges[n].face_indices;
    }
    mesh.edges.resize(static_cast<std::size_t>(cedge));
    for (auto& e : mesh.edges) e.bnd = mesh_bnd_none;
    mesh.avelen = average_edge_length(mesh);
    std::cout << " - Created " << cedge << " unique edges\n";
    std::cout << "Mesh data built successfully\n";
}

inline void build_solid_faces(MeshContainer& mesh) {
    (void)mesh;
    // Not implemented fully in the original excerpt; retained as a compatibility hook.
}

inline void compute_mesh_essentials(MeshContainer& mesh) {
    for (std::size_t n = 0; n < mesh.faces.size(); ++n) {
        Eigen::Vector3d p1 = mesh.nodes[mesh.faces[n].node_indices[0]].p;
        Eigen::Vector3d p2 = mesh.nodes[mesh.faces[n].node_indices[1]].p;
        Eigen::Vector3d p3 = mesh.nodes[mesh.faces[n].node_indices[2]].p;
        Eigen::Vector3d b = (p2 - p1).cross(p3 - p1);
        mesh.faces[n].n = b.normalized();
        mesh.faces[n].pd = mesh.faces[n].n.dot(p1);
        mesh.faces[n].cp = (p1 + p2 + p3) / 3.0;
    }
    for (std::size_t n = 0; n < mesh.solids.size(); ++n) {
        Eigen::Vector3d p1 = mesh.nodes[mesh.solids[n].node_indices[0]].p;
        Eigen::Vector3d p2 = mesh.nodes[mesh.solids[n].node_indices[1]].p;
        Eigen::Vector3d p3 = mesh.nodes[mesh.solids[n].node_indices[2]].p;
        Eigen::Vector3d p4 = mesh.nodes[mesh.solids[n].node_indices[3]].p;
        mesh.solids[n].volume = std::abs((p1 - p4).dot((p2 - p4).cross(p3 - p4))) / 6.0;
    }
}

inline void compute_basis_data(MeshContainer& mesh) {
    std::cout << "- Computing basis function data\n";
    for (std::size_t n = 0; n < mesh.faces.size(); ++n) {
        Eigen::Vector3d p1 = mesh.nodes[mesh.faces[n].node_indices[0]].p;
        Eigen::Vector3d p2 = mesh.nodes[mesh.faces[n].node_indices[1]].p;
        Eigen::Vector3d p3 = mesh.nodes[mesh.faces[n].node_indices[2]].p;
        Eigen::Vector3d b = (p2 - p1).cross(p3 - p1);
        mesh.faces[n].n = b.normalized();
        mesh.faces[n].pd = mesh.faces[n].n.dot(p1);
        Eigen::Vector3d p21 = p2 - p1;
        Eigen::Vector3d p32 = p3 - p2;
        Eigen::Vector3d p13 = p1 - p3;
        mesh.faces[n].s.col(0) = p21.normalized();
        mesh.faces[n].s.col(1) = p32.normalized();
        mesh.faces[n].s.col(2) = p13.normalized();
        mesh.faces[n].m.col(0) = mesh.faces[n].s.col(0).cross(mesh.faces[n].n).normalized();
        mesh.faces[n].m.col(1) = mesh.faces[n].s.col(1).cross(mesh.faces[n].n).normalized();
        mesh.faces[n].m.col(2) = mesh.faces[n].s.col(2).cross(mesh.faces[n].n).normalized();
        mesh.faces[n].area = 0.5 * (p32.cross(-p13)).norm();
        mesh.faces[n].cp = (p1 + p2 + p3) / 3.0;
    }
    for (std::size_t n = 0; n < mesh.edges.size(); ++n) {
        mesh.edges[n].length = (mesh.nodes[mesh.edges[n].node_indices[0]].p - mesh.nodes[mesh.edges[n].node_indices[1]].p).norm();
    }
    for (std::size_t n = 0; n < mesh.solids.size(); ++n) {
        // Solid volume and related data retained as compatibility.
    }
}

} // namespace mesh
} // namespace for90_mom2

// Convenience aliases for simpler user code.
namespace f90 = for90_mom2;

