#include "geometrycentral/surface/direction_fields.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/poisson_disk_sampler.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include "utils.h"

#include "streamlines.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;

bool loadedParam = false;
CornerData<Vector2> uv;

// field data
FaceData<Vector2> faceField;
int nSym              = 1;
bool saveField        = false;
std::string fieldName = "field.faceField";

// tracing options
int maxSegments = 150;
float maxLen    = 1.;

// svg options
bool saveSVG = false, saveOBJ = false;
float worldspaceStrokeWidth = .05;
float svgForeground[4]      = {0.f, 0.f, 0.f, 1.f};
float svgBackground[4]      = {1.f, 1.f, 1.f, 1.f};
float imageSize             = 500;
std::string svgName         = "streamlines.svg";
std::string objName         = "streamlines.obj";

// seed sampling options
float minDist = .75; // minimum distance r between samples

// visualization data
bool viz = false;
polyscope::SurfaceMesh* psMesh;

bool hasTextureCoordinates(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string firstToken;
        if (iss >> firstToken) { // Check if each line starts with "vt"
            if (firstToken == "vt") return true;
        }
    }

    return false;
}

polyscope::SurfaceFaceTangentVectorQuantity*
drawFaceField(polyscope::SurfaceMesh& psMesh, std::string name,
              ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
              const FaceData<Vector2>& faceField) {
    geom.requireFaceTangentBasis();
    FaceData<Vector3> fX(mesh), fY(mesh); // record face tangent bases
    for (Face f : mesh.faces()) {
        fX[f] = geom.faceTangentBasis[f][0];
        fY[f] = geom.faceTangentBasis[f][1];
    }
    geom.unrequireFaceTangentBasis();

    return psMesh.addFaceTangentVectorQuantity(name, faceField, fX, fY);
}

polyscope::CurveNetwork*
drawCurves(std::string name, ManifoldSurfaceMesh& mesh,
           VertexPositionGeometry& geom,
           const std::vector<std::vector<SurfacePoint>>& curves) {
    std::vector<Vector3> pts;
    std::vector<std::array<size_t, 2>> segs;
    std::vector<double> randomColors;
    for (const std::vector<SurfacePoint>& curve : curves) {
        double color = unitRand();
        for (size_t iP = 0; iP < curve.size(); iP++) {
            if (iP + 1 < curve.size())
                segs.push_back({pts.size(), pts.size() + 1});
            pts.push_back(curve[iP].interpolate(geom.vertexPositions));
            randomColors.push_back(color);
        }
    }
    auto psCurves = polyscope::registerCurveNetwork(name, pts, segs);
    psCurves->addNodeScalarQuantity("color", randomColors)->setEnabled(true);
    return psCurves;
}

void writeFaceField(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                    const FaceData<Vector2>& field, std::string name) {
    std::fstream out;
    out.open(name, std::ios::out | std::ios::trunc);
    if (out.is_open()) {
        geom.requireFaceTangentBasis();
        for (Face f : mesh.faces()) {
            Vector3 fX = geom.faceTangentBasis[f][0];
            Vector3 fY = geom.faceTangentBasis[f][1];
            Vector3 v  = field[f].x * fX + field[f].y * fY;
            out << v.x << " " << v.y << " " << v.z << std::endl;
        }
        geom.unrequireFaceTangentBasis();
        out.close();
        std::cout << "File " << name << " written successfully." << std::endl;
    } else {
        std::cout << "Could not save svg '" << name << "'!" << std::endl;
    }
}

void writeCurvesToOBJ(ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
                      const std::vector<std::vector<SurfacePoint>>& curves,
                      std::string name) {
    std::fstream out;
    out.open(name, std::ios::out | std::ios::trunc);
    if (out.is_open()) {
        // write vertices
        for (const std::vector<SurfacePoint>& curve : curves) {
            for (const SurfacePoint& p : curve) {
                Vector3 pos = p.interpolate(geom.vertexPositions);
                out << "v " << pos.x << " " << pos.y << " " << pos.z
                    << std::endl;
            }
        }
        // write lines
        size_t iStart = 1;
        for (const std::vector<SurfacePoint>& curve : curves) {
            out << "l";
            for (size_t iP = 0; iP < curve.size(); iP++)
                out << " " << (iP + iStart);
            out << std::endl;
            iStart += curve.size();
        }
        out.close();
        std::cout << "File " << name << " written successfully." << std::endl;
    } else {
        std::cout << "Could not save svg '" << name << "'!" << std::endl;
    }
}

std::vector<std::vector<SurfacePoint>> traceManyStreamlines(
    ManifoldSurfaceMesh& mesh, VertexPositionGeometry& geom,
    const FaceData<Vector2>& field, size_t nSym = 1,
    TraceStreamlineOptions opt = defaultTraceStreamlineOptions) {
    std::vector<std::vector<SurfacePoint>> streamlines;

    // generate Poisson disk samples to trace from
    PoissonDiskSampler sampler(mesh, geom);
    PoissonDiskOptions pdOpt;
    pdOpt.minDist = minDist;
    for (const SurfacePoint& p : sampler.sample(pdOpt)) {
        streamlines.push_back(traceStreamline(mesh, geom, p, field, nSym, opt));
    }
    return streamlines;
}


polyscope::SurfaceFaceTangentVectorQuantity* generateFaceField() {
    faceField = computeSmoothestFaceDirectionField(*geom, nSym);
    for (Face f : mesh->faces())
        faceField[f] = faceField[f].pow(1. / (double)nSym);

    if (saveField) writeFaceField(*mesh, *geom, faceField, fieldName);

    if (viz) {
        return drawFaceField(*psMesh, "smooth direction field", *mesh, *geom,
                             faceField);
    } else {
        return nullptr;
    }
}

std::string rgbaToHex(const float color[4]) {
    // Convert float values (0.0 - 1.0) to integer values (0 - 255)
    int r = static_cast<int>(std::round(color[0] * 255.0f));
    int g = static_cast<int>(std::round(color[1] * 255.0f));
    int b = static_cast<int>(std::round(color[2] * 255.0f));
    int a = static_cast<int>(std::round(color[3] * 255.0f));

    // Clamp values to ensure they're within valid range
    r = std::max(0, std::min(255, r));
    g = std::max(0, std::min(255, g));
    b = std::max(0, std::min(255, b));
    a = std::max(0, std::min(255, a));

    // Create hex string
    std::stringstream ss;
    ss << "#" << std::hex << std::uppercase << std::setfill('0');
    ss << std::setw(2) << r;
    ss << std::setw(2) << g;
    ss << std::setw(2) << b;

    // Include alpha channel if it's not fully opaque
    if (a < 255) {
        ss << std::setw(2) << a;
    }

    return ss.str();
}

bool hexToRGBA(const std::string& hex, float rgba[4]) {
    std::string hexColor = hex;

    // Remove '#' if present
    if (hexColor[0] == '#') {
        hexColor = hexColor.substr(1);
    }

    // Check if the hex string is valid
    if (hexColor.length() != 6 && hexColor.length() != 8) {
        return false;
    }

    // Parse RGB values
    unsigned int r, g, b, a = 255;

    try {
        std::stringstream ss;
        ss << std::hex << hexColor.substr(0, 2);
        ss >> r;
        ss.clear();
        ss << std::hex << hexColor.substr(2, 2);
        ss >> g;
        ss.clear();
        ss << std::hex << hexColor.substr(4, 2);
        ss >> b;

        // Parse alpha if present
        if (hexColor.length() == 8) {
            ss.clear();
            ss << std::hex << hexColor.substr(6, 2);
            ss >> a;
        }

        // Convert to float values between 0 and 1
        rgba[0] = r / 255.0f;
        rgba[1] = g / 255.0f;
        rgba[2] = b / 255.0f;
        rgba[3] = a / 255.0f;

        return true;
    } catch (...) {
        return false;
    }
}

FaceData<Vector2> readFaceField(const std::string& filename,
                                ManifoldSurfaceMesh& mesh,
                                VertexPositionGeometry& geom) {
    std::ifstream in(filename);
    if (in.is_open()) {
        FaceData<Vector2> faceField(mesh);
        geom.requireFaceTangentBasis();
        std::string line;
        size_t iF = 0;
        while (std::getline(in, line)) {
            // Skip empty lines and comment lines starting with #
            if (line.empty() || line[0] == '#') continue;

            std::istringstream iss(line);
            double x, y, z;
            iss >> x >> y >> z;
            Vector3 v    = {x, y, z}; // extrinsic vector
            Face f       = mesh.face(iF);
            Vector3 fX   = geom.faceTangentBasis[f][0];
            Vector3 fY   = geom.faceTangentBasis[f][1];
            faceField[f] = {dot(fX, v), dot(fY, v)};
            iF++;
        }

        in.close();
        geom.requireFaceTangentBasis();
        return faceField;
    } else {
        throw std::runtime_error("Could not open file: " + filename);
    }
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {
    ImGui::InputInt("nSym", &nSym);
    if (ImGui::Button("Generate n-Vector Field")) {
        generateFaceField();
    }
    if (ImGui::Button("Trace many streamlines")) {
        if (faceField.size() == 0) generateFaceField();
        TraceStreamlineOptions opt = defaultTraceStreamlineOptions;
        opt.maxSegments            = maxSegments;
        opt.maxLen                 = maxLen;
        std::vector<std::vector<SurfacePoint>> streamlines =
            traceManyStreamlines(*mesh, *geom, faceField, size_t(nSym), opt);
        drawCurves("streamlines", *mesh, *geom, streamlines);
        if (saveOBJ) {
            writeCurvesToOBJ(*mesh, *geom, streamlines, objName);
        }
        if (saveSVG) {
            if (!loadedParam) {
                std::cout << "ERROR: cannot write streamlines to SVG because "
                             "mesh does not have UV coordinates"
                          << std::endl;
            } else {
                SvgCurveOptions svgOpt;
                svgOpt.worldspaceStrokeWidth = worldspaceStrokeWidth;
                svgOpt.backgroundColor       = rgbaToHex(svgBackground);
                svgOpt.curveColor            = rgbaToHex(svgForeground);
                svgOpt.imageSize             = imageSize;
                draw_mesh_curves_to_svg(*mesh, *geom, uv, streamlines, svgName,
                                        svgOpt);
            }
        }
    }
    ImGui::Checkbox("Save to SVG", &saveSVG);
    ImGui::SameLine();
    ImGui::Checkbox("Save to OBJ", &saveOBJ);
    ImGui::SameLine();
    ImGui::Checkbox("Save field", &saveField);

    if (ImGui::TreeNode("Streamline options")) {
        ImGui::InputFloat("Max streamline length", &maxLen);
        ImGui::InputInt("Max streamline segments", &maxSegments);
        ImGui::InputFloat("Seed dist", &minDist);
        ImGui::TreePop();
    }

    if (saveSVG) {
        if (ImGui::TreeNode("SVG options")) {
            ImGui::InputFloat("Worldspace stroke width",
                              &worldspaceStrokeWidth);
            ImGui::ColorEdit4("Curve color", svgForeground);
            ImGui::ColorEdit4("Background color", svgBackground);
            ImGui::InputFloat("Image size", &imageSize);
            ImGui::TreePop();
        }
    }
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Streamline tracer");
    args::Positional<std::string> inputFilenameArg(
        parser, "mesh", "Mesh to be processed.", args::Options::Required);

    args::Flag vizArg(parser, "viz", "Run program with GUI", {'v', "viz"});
    args::ValueFlag<std::string> svgNameArg(
        parser, "svg-name", "Save streamline SVG to this file", {"svg-name"});
    args::ValueFlag<std::string> objNameArg(
        parser, "obj-name", "Save streamline obj to this file", {"obj-name"});

    // Field options
    args::Group fieldGroup(parser, "Field options:");
    args::ValueFlag<string> inputFieldArg(
        fieldGroup, "input-field",
        "Read vector field from this file. (If none is provided, a smooth "
        "field will be generated)",
        {"input-field"});
    args::ValueFlag<string> outputFieldArg(
        fieldGroup, "output-field", "Write generated vector field to this file",
        {"output-field"});
    args::ValueFlag<int> nSymArg(fieldGroup, "nSym", "symmetry of vector field",
                                 {"nSym"});

    // Tracing options
    args::Group tracingGroup(parser, "Tracing options:");
    args::ValueFlag<int> maxSegmentsArg(tracingGroup, "segments",
                                        "Maximum number of segments",
                                        {"max-segments"});
    args::ValueFlag<float> maxLenArg(tracingGroup, "length",
                                     "Maximum segment length", {"max-len"});
    args::ValueFlag<float> minDistArg(tracingGroup, "distance",
                                      "Minimum distance between samples",
                                      {"min-dist"});

    // SVG options
    args::Group svgGroup(parser, "SVG options:");
    args::ValueFlag<float> strokeWidthArg(
        svgGroup, "width", "Worldspace stroke width", {"stroke-width"});
    args::ValueFlag<std::string> foregroundArg(
        svgGroup, "color",
        "SVG foreground color (hex format: #RRGGBB or #RRGGBBAA)",
        {"foreground"});
    args::ValueFlag<std::string> backgroundArg(
        svgGroup, "color",
        "SVG background color (hex format: #RRGGBB or #RRGGBBAA)",
        {"background"});
    args::ValueFlag<float> imageSizeArg(svgGroup, "size", "Image size",
                                        {"image-size"});
    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = args::get(inputFilenameArg); // argument was required

    // Set values from arguments
    if (nSymArg) nSym = args::get(nSymArg);
    if (maxSegmentsArg) maxSegments = args::get(maxSegmentsArg);
    if (maxLenArg) maxLen = args::get(maxLenArg);
    if (strokeWidthArg) worldspaceStrokeWidth = args::get(strokeWidthArg);
    if (imageSizeArg) imageSize = args::get(imageSizeArg);
    if (minDistArg) minDist = args::get(minDistArg);
    if (vizArg) viz = true;
    if (svgNameArg) {
        svgName = args::get(svgNameArg);
        saveSVG = true;
    }
    if (objNameArg) {
        objName = args::get(objNameArg);
        saveOBJ = true;
    }
    if (outputFieldArg) {
        fieldName = args::get(outputFieldArg);
        saveField = true;
    }

    if (foregroundArg) {
        std::string hexColor = args::get(foregroundArg);
        if (!hexToRGBA(hexColor, svgForeground)) {
            std::cout << "Error: Invalid foreground color format. Use #RRGGBB "
                         "or #RRGGBBAA"
                      << std::endl;
            return 1;
        }
    }

    if (backgroundArg) {
        std::string hexColor = args::get(backgroundArg);
        if (!hexToRGBA(hexColor, svgBackground)) {
            std::cout << "Error: Invalid background color format. Use #RRGGBB "
                         "or #RRGGBBAA"
                      << std::endl;
            return 1;
        }
    }

    if (viz) {
        // Initialize polyscope
        polyscope::init();

        // Set the callback function
        polyscope::state::userCallback = myCallback;
    }

    // Load mesh
    if (endsWith(filename, ".obj")) {
        // if input is obj, try to load texture coords
        if (hasTextureCoordinates(filename)) {
            std::unique_ptr<CornerData<Vector2>> uvPtr;
            std::tie(mesh, geom, uvPtr) =
                readParameterizedManifoldSurfaceMesh(filename);
            uv          = *uvPtr; // copy param to plain CornerData
            loadedParam = true;
        } else {
            std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);
            loadedParam          = false;
        }
    } else {
        std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);
        loadedParam          = false;
    }

    if (loadedParam) {
        std::cout << "  loaded mesh with uv coordinates" << std::endl;
    } else {
        std::cout << "  loaded mesh, but did not find any uv coordinates"
                  << std::endl;
    }

    // if present, load input field
    if (inputFieldArg) {
        faceField = readFaceField(args::get(inputFieldArg), *mesh, *geom);
    }

    if (viz) {
        // Register the mesh with polyscope
        psMesh = polyscope::registerSurfaceMesh(
            polyscope::guessNiceNameFromPath(filename), geom->vertexPositions,
            mesh->getFaceVertexList(), polyscopePermutations(*mesh));

        if (inputFieldArg) {
            drawFaceField(*psMesh, "smooth direction field", *mesh, *geom,
                          faceField);

            TraceStreamlineOptions opt = defaultTraceStreamlineOptions;
            opt.maxSegments            = maxSegments;
            opt.maxLen                 = maxLen;
            std::vector<std::vector<SurfacePoint>> streamlines =
                traceManyStreamlines(*mesh, *geom, faceField, size_t(nSym),
                                     opt);
            drawCurves("streamlines", *mesh, *geom, streamlines);
            if (saveOBJ) {
                writeCurvesToOBJ(*mesh, *geom, streamlines, objName);
            }
            if (saveSVG) {
                if (!loadedParam) {
                    std::cout
                        << "ERROR: cannot write streamlines to SVG because "
                           "mesh does not have UV coordinates"
                        << std::endl;
                } else {
                    SvgCurveOptions svgOpt;
                    svgOpt.worldspaceStrokeWidth = worldspaceStrokeWidth;
                    svgOpt.backgroundColor       = rgbaToHex(svgBackground);
                    svgOpt.curveColor            = rgbaToHex(svgForeground);
                    svgOpt.imageSize             = imageSize;
                    draw_mesh_curves_to_svg(*mesh, *geom, uv, streamlines,
                                            svgName, svgOpt);
                }
            }
        }

        // Give control to the polyscope gui
        polyscope::show();
    } else {
        if (faceField.size() == 0) generateFaceField();
        TraceStreamlineOptions opt = defaultTraceStreamlineOptions;
        opt.maxSegments            = maxSegments;
        opt.maxLen                 = maxLen;
        std::vector<std::vector<SurfacePoint>> streamlines =
            traceManyStreamlines(*mesh, *geom, faceField, size_t(nSym), opt);
        if (saveOBJ) {
            writeCurvesToOBJ(*mesh, *geom, streamlines, objName);
        }
        if (saveSVG) {
            if (!loadedParam) {
                std::cout << "ERROR: cannot write streamlines to SVG because "
                             "mesh does not have UV coordinates"
                          << std::endl;
            } else {
                SvgCurveOptions svgOpt;
                svgOpt.worldspaceStrokeWidth = worldspaceStrokeWidth;
                svgOpt.backgroundColor       = rgbaToHex(svgBackground);
                svgOpt.curveColor            = rgbaToHex(svgForeground);
                svgOpt.imageSize             = imageSize;
                draw_mesh_curves_to_svg(*mesh, *geom, uv, streamlines, svgName,
                                        svgOpt);
            }
        }
    }

    return EXIT_SUCCESS;
}
