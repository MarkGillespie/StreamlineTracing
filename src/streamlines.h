#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/surface_point.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"

#include "utils.h"

namespace geometrycentral {
namespace surface {

struct SvgCurveOptions {
    double worldspaceStrokeWidth = .05;
    std::string backgroundColor  = "#fff";
    std::string curveColor       = "#000000";
    std::function<std::string(size_t, size_t)> curveColorFunction;
    double imageSize = 500;
};
extern const SvgCurveOptions defaultSvgCurveOptions;

void draw_mesh_curves_to_svg(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const CornerData<Vector2>& uv,
    const std::vector<std::vector<SurfacePoint>>& curves, std::string name,
    SvgCurveOptions opt = defaultSvgCurveOptions);

struct TraceStreamlineOptions {
    size_t maxSegments             = 10;
    double maxLen                  = std::numeric_limits<double>::infinity();
    FaceData<int>* nVisits         = nullptr;
    const FaceData<int>* maxVisits = nullptr;
    int verbosity                  = 0;
};
extern const TraceStreamlineOptions defaultTraceStreamlineOptions;

std::vector<SurfacePoint>
traceStreamline(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                SurfacePoint pStart, const FaceData<Vector2>& field,
                size_t nSym                = 1,
                TraceStreamlineOptions opt = defaultTraceStreamlineOptions);

// helpers
Vector2 interpUV(Face f, Vector3 b, const CornerData<Vector2>& uv);
} // namespace surface
} // namespace geometrycentral
