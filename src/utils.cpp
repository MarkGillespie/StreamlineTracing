#include "utils.h"

double cot(double x) { return cos(x) / sin(x); }

verbose_runtime_error::verbose_runtime_error(const std::string& arg,
                                             const char* file, int line)
    : std::runtime_error(arg) {

    std::string niceName = getFilename(file);

    std::ostringstream o;
    o << arg << " at " << niceName << ":" << line;
    msg = o.str();
}

Eigen::Vector3d toEigen(const Vector3& v) {
    Eigen::Vector3d ret;
    ret << v.x, v.y, v.z;
    return ret;
}
Eigen::Vector2d toEigen(const Vector2& v) {
    Eigen::Vector2d ret;
    ret << v.x, v.y;
    return ret;
}

Vector3 fromEigen(const Eigen::Vector3d& v) {
    return Vector3{v[0], v[1], v[2]};
}

// Returns vw^T
Eigen::Matrix3d outer(const Vector3& v, const Vector3& w) {
    return toEigen(v) * toEigen(w).transpose();
}


string getFilename(string filePath) {

    // stolen from polyscope/utilities.cpp
    size_t startInd = 0;
    for (std::string sep : {"/", "\\"}) {
        size_t pos = filePath.rfind(sep);
        if (pos != std::string::npos) {
            startInd = std::max(startInd, pos + 1);
        }
    }

    return filePath.substr(startInd, filePath.size());
}

// https://stackoverflow.com/a/874160
bool endsWith(const std::string& str, const std::string& ending) {
    return (str.length() >= ending.length()) &&
           (0 == str.compare(str.length() - ending.length(), ending.length(),
                             ending));
}

namespace geometrycentral {
namespace surface {
std::pair<std::unique_ptr<ManifoldSurfaceMesh>,
          std::unique_ptr<EdgeLengthGeometry>>
copyGeometry(ManifoldSurfaceMesh& mesh, IntrinsicGeometryInterface& geo) {
    std::unique_ptr<ManifoldSurfaceMesh> meshCopy = mesh.copy();
    geo.requireEdgeLengths();
    EdgeData<double> lengthsCopy = geo.edgeLengths.reinterpretTo(*meshCopy);
    std::unique_ptr<EdgeLengthGeometry> geoCopy;
    geoCopy.reset(new EdgeLengthGeometry(*meshCopy, lengthsCopy));
    return std::make_pair(std::move(meshCopy), std::move(geoCopy));
}

bool checkAdjacent(Vertex vA, Vertex vB) {
    bool adjacent = false;
    for (Vertex vN : vA.adjacentVertices()) {
        if (vN == vB) adjacent = true;
    }
    return adjacent;
}

bool checkAdjacent(Vertex vA, Edge eB) {
    bool adjacent = false;
    for (Halfedge he : vA.outgoingHalfedges()) {
        if (eB == he.next().edge()) adjacent = true;
        if (eB == he.edge()) adjacent = true;
    }
    return adjacent;
}

bool checkAdjacent(Vertex vA, Face fB) {
    bool adjacent = false;
    for (Face fN : vA.adjacentFaces()) {
        if (fN == fB) adjacent = true;
    }
    return adjacent;
}

bool checkAdjacent(Edge eA, Vertex vB) { return checkAdjacent(vB, eA); }

bool checkAdjacent(Edge eA, Edge eB) {
    bool adjacent = false;

    // Must have a shared face
    adjacent |= (eA.halfedge().face() == eB.halfedge().face());
    adjacent |= (eA.halfedge().twin().face() == eB.halfedge().face());
    adjacent |= (eA.halfedge().face() == eB.halfedge().twin().face());
    adjacent |= (eA.halfedge().twin().face() == eB.halfedge().twin().face());

    return adjacent;
}

bool checkAdjacent(Edge eA, Face fB) {
    bool adjacent = false;
    for (Edge eN : fB.adjacentEdges()) {
        if (eN == eA) adjacent = true;
    }
    return adjacent;
}


bool checkAdjacent(Face fA, Vertex vB) { return checkAdjacent(vB, fA); }
bool checkAdjacent(Face fA, Edge eB) { return checkAdjacent(eB, fA); }
bool checkAdjacent(Face fA, Face fB) { return fA == fB; }

bool onFace(const SurfacePoint& p, Face f) {
    if (p == SurfacePoint() || f == Face()) return false;
    switch (p.type) {
    case SurfacePointType::Vertex:
        return checkAdjacent(p.vertex, f);
    case SurfacePointType::Edge:
        return checkAdjacent(p.edge, f);
    case SurfacePointType::Face:
        return checkAdjacent(p.face, f);
    }
}

bool isNear(const SurfacePoint& p, const SurfacePoint& q, double tol) {
    if (p.type != q.type) return false;
    if ((p.type == SurfacePointType::Vertex && p.vertex == q.vertex) ||
        (p.type == SurfacePointType::Edge && p.edge == q.edge &&
         abs(p.tEdge - q.tEdge) < tol) ||
        (p.type == SurfacePointType::Face && p.face == q.face &&
         (p.faceCoords - q.faceCoords).norm() < tol)) {
        return true;
    } else {
        return false;
    }
}

} // namespace surface
} // namespace geometrycentral

// Returns the intersection of the line a1-a2 with line b1-b2 in barycentric
// coordinates along line a1-a2 (i.e. (1-t)*a1 + t*a2 lies on line b1-b2 as
// well)
double intersectionTime(Vector2 a1, Vector2 a2, Vector2 b1, Vector2 b2) {
    // TODO: are these checks time-consuming?
    if ((a1 - b1).norm() < 1e-12 || (a1 - b2).norm() < 1e-12) {
        return 1;
    } else if ((a2 - b1).norm() < 1e-12 || (a2 - b2).norm() < 1e-12) {
        return 0;
    } else if ((a1 - a2).norm() < 1e-12 && (b1 - b2).norm() < 1e-12) {
        // If a and b coincided, we would have caught that in a previous case.
        // So the "lines" don't intersect, and we just return a big number
        return 50;
    } else if ((b1 - b2).norm() < 1e-12) {
        Vector2 dA = (a2 - a1).normalize();
        return 1 - dot(b1 - a1, dA) / dot(a2 - a1, dA);
    } else if ((a1 - a2).norm() < 1e-12) {
        Vector2 dB = (b2 - b1).normalize();
        double t   = dot(a1 - b1, dB) / dot(b2 - b1, dB);

        if (0 <= t && t <= 1) {
            Vector2 projAB = b1 + t * dB;
            if ((projAB - a1).norm() < 1e-12) {
                return 0.5;
            } else {
                return 50;
            }
        } else {
            return 50;
        }

    } else {
        double m1 = a2.x - a1.x;
        double m2 = b2.x - b1.x;
        double m3 = a2.y - a1.y;
        double m4 = b2.y - b1.y;
        double vx = b1.x - a1.x;
        double vy = b1.y - a1.y;
        double t  = (m2 * vy - m4 * vx) / (m2 * m3 - m1 * m4);

        if (t != t) {
            std::cerr << "intersectionTime returned NaN" << vendl;
            std::cerr << "a1: " << a1 << vendl;
            std::cerr << "a2: " << a2 << vendl;
            std::cerr << "b1: " << b1 << vendl;
            std::cerr << "b2: " << b2 << vendl;
            std::cerr << "det: " << m2 * m3 - m1 * m4 << vendl;
        }

        return t;
    }
}
