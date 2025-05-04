#include "streamlines.h"

#include "geometrycentral/surface/barycentric_coordinate_helpers.h" // permuteBarycentric...

namespace geometrycentral {
namespace surface {

const SvgCurveOptions defaultCurveOptions;

void draw_mesh_curves_to_svg(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const CornerData<Vector2>& uv,
    const std::vector<std::vector<SurfacePoint>>& curves, std::string name,
    SvgCurveOptions opt) {

    if (!opt.curveColorFunction) {
        opt.curveColorFunction = [&](size_t iCurve,
                                     size_t nCurves) -> std::string {
            return opt.curveColor;
        };
    }

    geom.requireEdgeLengths();
    auto uvLen = [&](Halfedge ij) -> double {
        return (uv[ij.next().corner()] - uv[ij.corner()]).norm();
    };

    auto lengthScaling = [&](Face ijk, Vector3 b) -> double {
        Halfedge ij = ijk.halfedge(), jk = ijk.halfedge().next();
        Halfedge ki = ijk.halfedge().next().next();
        if (b[0] < 1e-8) { // point lies on jk
            return uvLen(jk) / geom.edgeLengths[jk.edge()];
        } else if (b[1] < 1e-8) { // point lies on ki
            return uvLen(ki) / geom.edgeLengths[ki.edge()];
        } else if (b[2] < 1e-8) { // point lies on ij
            return uvLen(ij) / geom.edgeLengths[ij.edge()];
        } else { // otherwise just return average
            return (uvLen(ij) / geom.edgeLengths[ij.edge()] +
                    uvLen(jk) / geom.edgeLengths[jk.edge()] +
                    uvLen(ki) / geom.edgeLengths[ki.edge()]) /
                   3;
        }
    };

    auto meanScaling = [&](const std::vector<SurfacePoint>& curve) {
        double totalScaling = 0, nP = 0;
        for (size_t iP = 0; iP + 1 < curve.size(); iP++) {
            SurfacePoint pi = curve[iP], pj = curve[iP + 1];
            Face f = sharedFace(pi, pj);
            totalScaling += lengthScaling(f, pi.inFace(f).faceCoords);
            totalScaling += lengthScaling(f, pj.inFace(f).faceCoords);
            nP += 2.;
        }
        return totalScaling / nP;
    };

    auto flip         = [](Vector2 v) -> Vector2 { return {v.x, 1 - v.y}; };
    const double& dim = opt.imageSize;
    std::fstream out;
    out.open(name, std::ios::out | std::ios::trunc);
    if (out.is_open()) {
        out << "<svg viewBox=\"0 0 " << dim << " " << dim << "\">" << std::endl;
        if (!opt.backgroundColor.empty() && opt.backgroundColor != "None") {
            out << "<polygon fill=\"" << opt.backgroundColor
                << "\" stroke=\"none\" points=\"0,0 " << dim << ",0 " << dim
                << "," << dim << " 0," << dim << "\"/>" << std::endl;
        }

        out << "<g>" << std::endl;

        for (size_t iCurve = 0; iCurve < curves.size(); iCurve++) {
            double scale =
                dim * opt.worldspaceStrokeWidth * meanScaling(curves[iCurve]);
            Vector2 qPrev = dim * Vector2{100, 100};
            for (size_t iP = 0; iP + 1 < curves[iCurve].size(); iP++) {
                SurfacePoint pi = curves[iCurve][iP],
                             pj = curves[iCurve][iP + 1];
                Face f          = sharedFace(pi, pj);
                Vector2 qi =
                    dim * flip(interpUV(f, pi.inFace(f).faceCoords, uv));
                Vector2 qj =
                    dim * flip(interpUV(f, pj.inFace(f).faceCoords, uv));
                double d = (qj - qi).norm();

                if ((qi - qPrev).norm() < dim * 1e-5) { // continuing old curve
                    out << qPrev.x << "," << qPrev.y << " ";
                } else { // start a new curve
                    if (iP > 0) out << "\"/>" << std::endl;
                    out << "<polyline fill=\"none\" stroke=\""
                        << opt.curveColorFunction(iCurve, curves.size())
                        << "\" stroke-width=\"" << scale << "\" points=\"";
                    out << qi.x << "," << qi.y << " ";
                    out << qj.x << "," << qj.y << " ";
                }
                qPrev = qj;
            }
            out << "\"/>" << std::endl;
        }
        out << "</g>" << std::endl;
        out << "</svg>" << std::endl;
        out.close();
        std::cout << "File " << name << " written successfully." << std::endl;
    } else {
        std::cout << "Could not save svg '" << name << "'!" << std::endl;
    }
    geom.unrequireEdgeLengths();
}

std::array<Vector2, 3>
vertexCoordinatesInTriangle(IntrinsicGeometryInterface& geom, Face face) {
    return {Vector2{0., 0.}, geom.halfedgeVectorsInFace[face.halfedge()],
            -geom.halfedgeVectorsInFace[face.halfedge().next().next()]};
}

// Return type from tracing subroutines, slightly modified from TraceSubResult
// in trace_geodesic.cpp
struct StreamlineSubResult {
    bool terminated;    // Did the trace end?
    double traceLength; // Length of segment within current face

    // One of the two sets of values will be defined:

    // If the trace continues (terminated == false)
    Halfedge crossHe; // halfedge we hit (crossHe.twin().face() is new face)
    double tCross;    // t along crossHe
    Vector2 traceDir; // direction followed in current face

    // If the trace ends (terminated == true)
    SurfacePoint endPoint; // ending location
};

// General form for tracing barycentrically within a face
// Assumes that approriate projects have already been performed such that
// startPoint and vectors are valid (inside triangle and pointing in the right
// direction) Returns halfedge hit, time of intersection along halfedge,
// direction followed (which may not be vecCartesian for n-vector fields), and
// length of segment in face
//
// This function is tightly coupled with the routines which call it. They
// prepare the values startPoint, vecBary, and vecCartesian, ensuring that those
// values satisify basic properties (essentially that the trace points in a
// valid direction).

// if nDirections>1, treats vecCartesian as one of the directions in an n-vector
// field, and searches for the direction closest to oldDirCartesian
// std::tuple<Halfedge, double, Vector2, double>
inline StreamlineSubResult traceStreamlineInFace(
    IntrinsicGeometryInterface& geom, Face face, Vector3 startPoint,
    Vector2 vecCartesian, double remainingLength, size_t nDirections = 1,
    Vector2 oldDirCartesian = Vector2{1, 0}, bool TRACE_PRINT = false) {

    if (nDirections > 1) { // find closest direction
        Vector2 vClosest = vecCartesian;
        double dClosest  = (oldDirCartesian - vecCartesian).norm2();
        for (size_t iRot = 1; iRot < nDirections; iRot++) {
            Vector2 rot = Vector2::fromAngle(2. * M_PI / (double)nDirections *
                                             (double)iRot);
            double d    = (oldDirCartesian - rot * vecCartesian).norm2();
            if (d < dClosest) {
                dClosest = d;
                vClosest = rot * vecCartesian;
            }
        }
        vecCartesian = vClosest;
    }

    // TODO: fix GC headers so we can use the GC version
    auto cartesianVectorToBarycentric =
        [](const std::array<Vector2, 3>& vertCoords,
           Vector2 faceVec) -> Vector3 {
        // Build matrix for linear transform problem
        // (last constraint comes from chosing the displacement vector with sum
        // = 0)
        Eigen::Matrix3d A;
        Eigen::Vector3d rhs;
        const std::array<Vector2, 3>& c = vertCoords; // short name
        A << c[0].x, c[1].x, c[2].x, c[0].y, c[1].y, c[2].y, 1., 1., 1.;
        rhs << faceVec.x, faceVec.y, 0.;

        // Solve
        Eigen::Vector3d result = A.colPivHouseholderQr().solve(rhs);
        Vector3 resultBary{result(0), result(1), result(2)};

        resultBary = normalizeBarycentricDisplacement(resultBary);

        return resultBary;
    };

    // Gather values
    double speed = vecCartesian.norm();
    double tMax  = remainingLength / speed;
    std::array<Vector2, 3> vertexCoords =
        vertexCoordinatesInTriangle(geom, face);
    Vector3 triangleLengths{
        geom.edgeLengths[face.halfedge().edge()],
        geom.edgeLengths[face.halfedge().next().edge()],
        geom.edgeLengths[face.halfedge().next().next().edge()]};
    Vector3 vecBary = cartesianVectorToBarycentric(vertexCoords, vecCartesian);

    if (sum(startPoint) < 0.5) {
        if (TRACE_PRINT) {
            std::cout << "  bad bary point: " << startPoint << std::endl;
        }
        throw std::runtime_error("bad bary point");
    }

    if (TRACE_PRINT) {
        std::cout << "  general trace in face: " << std::endl;
        std::cout << "  face: " << face << " startPoint " << startPoint
                  << " vecBary = " << vecBary << " vecCartesian "
                  << vecCartesian << std::endl;
    }

    // Find first hit along opposite edges
    double tRay      = std::numeric_limits<double>::infinity();
    Halfedge crossHe = Halfedge();
    int iOppVertEnd  = -777;
    Halfedge currHe  = face.halfedge();
    for (int i = 0; i < 3; i++) {
        currHe = currHe.next(); // always opposite the i'th vertex

        if (vecBary[i] >= 0) continue;

        // Check the crossing
        double tRayThis = -startPoint[i] / vecBary[i];

        if (TRACE_PRINT) {
            std::cout << "    considering intersection:" << std::endl;
            std::cout << "      vecBary[i]: " << vecBary[i] << std::endl;
            std::cout << "      startPoint[i]: " << startPoint[i] << std::endl;
            std::cout << "      tRayThis: " << tRayThis << std::endl;
        }

        if (tRayThis < tRay) {
            // This is the new closest intersection
            tRay        = tRayThis;
            crossHe     = currHe;
            iOppVertEnd = i;
        }
    }

    if (TRACE_PRINT) {
        std::cout << "    selected intersection:" << std::endl;
        std::cout << "      crossHe: " << crossHe << std::endl;
        std::cout << "      tRay: " << tRay << std::endl;
        std::cout << "      iOppVertEnd: " << iOppVertEnd << std::endl;
    }

    if (crossHe == Halfedge()) {
        // throw std::logic_error("no halfedge intersection was selected,
        // precondition problem?");
        if (TRACE_PRINT) {
            std::cout << "    PROBLEM PROBLEM NO INTERSECTION:" << std::endl;
        }

        // End immediately
        StreamlineSubResult result;
        result.terminated = true;
        result.endPoint   = SurfacePoint(face, startPoint);
        return result;
    }

    // If the hit point happens after tMax, we terminate inside the triangle
    if (tRay >= tMax) {
        StreamlineSubResult result;
        result.terminated  = true;
        result.endPoint    = SurfacePoint(face, startPoint + tMax * vecBary);
        result.traceLength = tMax * speed;
        return result;
    }

    // Otherwise, compute some useful info about the endpoint
    Vector3 endPointOnEdge = startPoint + tRay * vecBary;
    double tCross          = endPointOnEdge[(iOppVertEnd + 2) % 3] /
                    (endPointOnEdge[(iOppVertEnd + 1) % 3] +
                     endPointOnEdge[(iOppVertEnd + 2) % 3]);
    if (TRACE_PRINT) {
        std::cout << "    end point on edge: " << endPointOnEdge << std::endl;
        std::cout << "    tCross raw: " << tCross << std::endl;
    }
    tCross = clamp(tCross, 0., 1.);

    StreamlineSubResult result;
    result.terminated  = false;
    result.crossHe     = crossHe;
    result.tCross      = tCross;
    result.traceDir    = vecCartesian;
    result.traceLength = tRay * speed;
    return result;
}

const TraceStreamlineOptions defaultTraceStreamlineOptions;
std::vector<SurfacePoint>
traceStreamline(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                SurfacePoint pStart, const FaceData<Vector2>& vector_field,
                size_t nSym, TraceStreamlineOptions opt) {
    std::vector<SurfacePoint> forwardStreamline{pStart}, reverseStreamline;

    // always start inside a face
    pStart = pStart.inSomeFace();

    if (opt.nVisits) (*opt.nVisits)[pStart.face]++;

    if (nSym > 1) geom.requireTransportVectorsAcrossHalfedge();

    Vector2 vStart = vector_field[pStart.face];
    if (nSym > 1) {
        size_t iRot = randomIndex(nSym);
        Vector2 rot =
            Vector2::fromAngle(2. * M_PI / (double)nSym * (double)iRot);
        vStart = rot * vStart;
    }

    if (opt.verbosity > 0) {
        std::cout << "=== tracing streamline from " << pStart
                  << " in direction " << vStart << std::endl;
    }

    double forwardLength = 0;
    double forwardLimit  = opt.maxLen / 2.;
    Vector2 vCurr        = vStart;

    bool printSteps            = opt.verbosity > 1;
    StreamlineSubResult result = traceStreamlineInFace(
        geom, pStart.face, pStart.faceCoords, vCurr,
        forwardLimit - forwardLength, nSym, vCurr, printSteps);
    forwardLength += result.traceLength;

    // if trace already terminated, grab endpoint. Then heCurr will equal
    // Halfedge(), so the remaining forward code will be skipped
    if (result.terminated) forwardStreamline.push_back(result.endPoint);

    Halfedge heCurr = result.crossHe;
    double tCurr    = result.tCross;
    if (nSym > 1 && heCurr != Halfedge())
        vCurr = geom.transportVectorsAcrossHalfedge[heCurr] * vCurr;

    if (opt.verbosity > 0)
        std::cout << "forwardLength: " << forwardLength << std::endl;

    auto should_visit_face = [&](Halfedge he, Vector2 v) {
        if (heCurr == Halfedge() || heCurr.edge().isBoundary()) return false;

        return !opt.nVisits || (*opt.nVisits)[heCurr.twin().face()] <
                                   (*opt.maxVisits)[heCurr.twin().face()];
    };

    while (should_visit_face(heCurr, vCurr) &&
           2 * forwardStreamline.size() < opt.maxSegments &&
           forwardLength < forwardLimit) {
        //=== step across heCurr to enter heCurr.twin().face()
        if (opt.nVisits) (*opt.nVisits)[heCurr.twin().face()]++;
        forwardStreamline.push_back(SurfacePoint(heCurr, tCurr));
        Vector3 bary =
            forwardStreamline.back().inFace(heCurr.twin().face()).faceCoords;
        result = traceStreamlineInFace(geom, heCurr.twin().face(), bary,
                                       vector_field[heCurr.twin().face()],
                                       forwardLimit - forwardLength, nSym,
                                       vCurr, printSteps);
        heCurr = result.crossHe, tCurr = result.tCross;
        forwardLength += result.traceLength;
        if (result.terminated) forwardStreamline.push_back(result.endPoint);

        if (nSym > 1 && heCurr != Halfedge())
            vCurr = geom.transportVectorsAcrossHalfedge[heCurr] * vCurr;
        if (opt.verbosity > 0)
            std::cout << "forwardLength: " << forwardLength << std::endl;
    }

    vCurr                 = -vStart;
    double backwardLength = 0;
    double backwardLimit  = opt.maxLen - forwardLength;
    result = traceStreamlineInFace(geom, pStart.face, pStart.faceCoords, vCurr,
                                   backwardLimit - backwardLength, nSym, vCurr,
                                   printSteps);
    heCurr = result.crossHe, tCurr = result.tCross;
    backwardLength += result.traceLength;

    // if trace already terminated, grab endpoint. Then heCurr will equal
    // Halfedge(), so the remaining backward code will be skipped
    if (result.terminated) reverseStreamline.push_back(result.endPoint);


    if (nSym > 1 && heCurr != Halfedge())
        vCurr = geom.transportVectorsAcrossHalfedge[heCurr] * vCurr;
    while (should_visit_face(heCurr, vCurr) &&
           forwardStreamline.size() + reverseStreamline.size() <
               opt.maxSegments &&
           backwardLength < backwardLimit) {
        //=== step across heCurr to enter heCurr.twin().face()
        if (opt.nVisits) (*opt.nVisits)[heCurr.twin().face()]++;
        reverseStreamline.push_back(SurfacePoint(heCurr, tCurr));
        Vector3 bary =
            reverseStreamline.back().inFace(heCurr.twin().face()).faceCoords;
        result = traceStreamlineInFace(geom, heCurr.twin().face(), bary,
                                       -vector_field[heCurr.twin().face()],
                                       backwardLimit - backwardLength, nSym,
                                       vCurr, printSteps);
        heCurr = result.crossHe, tCurr = result.tCross;
        backwardLength += result.traceLength;
        if (result.terminated) reverseStreamline.push_back(result.endPoint);

        if (nSym > 1 && heCurr != Halfedge())
            vCurr = geom.transportVectorsAcrossHalfedge[heCurr] * vCurr;
    }

    std::vector<SurfacePoint> streamline;
    for (int iP = reverseStreamline.size() - 1; iP >= 0; iP--) {
        streamline.push_back(reverseStreamline[iP]);
    }
    for (int iP = 0; iP < int(forwardStreamline.size()); iP++) {
        streamline.push_back(forwardStreamline[iP]);
    }

    if (nSym > 1) geom.unrequireTransportVectorsAcrossHalfedge();

    return streamline;
}

Vector2 interpUV(Face f, Vector3 b, const CornerData<Vector2>& uv) {
    return b[0] * uv[f.halfedge().corner()] +
           b[1] * uv[f.halfedge().next().corner()] +
           b[2] * uv[f.halfedge().next().next().corner()];
}

} // namespace surface
} // namespace geometrycentral
