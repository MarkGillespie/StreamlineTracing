#pragma once

#include "geometrycentral/numerical/linear_algebra_utilities.h"
#include "geometrycentral/surface/edge_length_geometry.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/utilities/vector2.h"
#include "geometrycentral/utilities/vector3.h"

#include <iostream>
#include <set>
#include <sstream>
#include <string>

using geometrycentral::Vector2;
using geometrycentral::Vector3;
using std::string;

// cotan
double cot(double x);

// Custom runtime error which will report the offending file and line number
class verbose_runtime_error : public std::runtime_error {
    std::string msg;

  public:
    verbose_runtime_error(const std::string& arg, const char* file, int line);
    ~verbose_runtime_error() throw() {}
    const char* what() const throw() { return msg.c_str(); }
};

#define throw_verbose_runtime_error(arg)                                       \
    throw verbose_runtime_error(arg, __FILE__, __LINE__);

#define verbose_assert(arg, msg)                                               \
    if (!(arg)) throw_verbose_runtime_error(msg);

//== Linear algebra helpers
Eigen::Vector3d toEigen(const Vector3& v);
Eigen::Vector2d toEigen(const Vector2& v);
Vector3 fromEigen(const Eigen::Vector3d& v);
Eigen::Matrix3d outer(const Vector3& v, const Vector3& w); // Returns vw^T

//== Convenient printing helpers
string getFilename(string filePath);

// https://stackoverflow.com/a/874160
bool endsWith(const std::string& str, const std::string& ending);

// Verbose endl
#define vendl                                                                  \
    "\t\t(" << getFilename(__FILE__) << ":" << __LINE__ << ")" << std::endl

// == Horrible macro to print variable name + variable
// https://stackoverflow.com/a/6623090
#define WATCHSTR_WNAME(os, name, a)                                            \
    do {                                                                       \
        (os) << (name) << " is value " << (a) << vendl;                        \
    } while (false)

#define WATCHSTR(os, a) WATCHSTR_WNAME((os), #a, (a))
#define WATCH(a) WATCHSTR_WNAME(std::cout, #a, (a))
#define WATCH2(a, b)                                                           \
    WATCH(a);                                                                  \
    WATCH(b)
#define WATCH3(a, b, c)                                                        \
    WATCH(a);                                                                  \
    WATCH2(b, c)
#define WATCH4(a, b, c, d)                                                     \
    WATCH(a);                                                                  \
    WATCH3(b, c, d)
#define WATCH5(a, b, c, d, e)                                                  \
    WATCH(a);                                                                  \
    WATCH4(b, c, d, e)

#define HERE()                                                                 \
    do {                                                                       \
        std::cout << "HERE at " << __LINE__ << " in " << getFilename(__FILE__) \
                  << ":" << __FUNCTION__ << std::endl;                         \
    } while (false)
#define HERE1(name)                                                            \
    do {                                                                       \
        std::cout << "HERE (" << (name) << ") at " << __LINE__ << " in "       \
                  << getFilename(__FILE__) << ":" << __FUNCTION__              \
                  << std::endl;                                                \
    } while (false)

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data);

template <typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const std::array<T, N>& data);

template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& data);

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T>& data);

template <typename T>
std::ostream& operator<<(std::ostream& out, const Eigen::Triplet<T>& data);

namespace geometrycentral {
namespace surface {
std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
          std::unique_ptr<EdgeLengthGeometry>>
copyGeometry(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo);

// From geometrycentral/src/surface/surface_point.cpp
bool checkAdjacent(Vertex vA, Vertex vB);
bool checkAdjacent(Vertex vA, Edge eB);
bool checkAdjacent(Vertex vA, Face fB);
bool checkAdjacent(Edge eA, Vertex vB);
bool checkAdjacent(Edge eA, Edge eB);
bool checkAdjacent(Edge eA, Face fB);
bool checkAdjacent(Face fA, Vertex vB);
bool checkAdjacent(Face fA, Edge eB);
bool checkAdjacent(Face fA, Face fB);

bool onFace(const SurfacePoint& p, Face f);
bool isNear(const SurfacePoint& p, const SurfacePoint& q, double tol = 1e-8);

} // namespace surface
} // namespace geometrycentral

// Returns the intersection of the line a1-a2 with line b1-b2 in barycentric
// coordinates along line a1-a2 (i.e. (1-t)*a1 + t*a2 lies on line b1-b2 as
// well)
double intersectionTime(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2);

#include "utils.ipp"
