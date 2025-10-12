#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>

// MapLibre coordinate system constants (from your integration)
constexpr int16_t EXTENT = 8192;  // From util::EXTENT in MapLibre

// Simulate MapLibre's GeometryCoordinate type
struct GeometryCoordinate {
    int16_t x, y;
    GeometryCoordinate(int16_t x_, int16_t y_) : x(x_), y(y_) {}
};

// MapLibre-style triangulation function that replicates your integration exactly
struct TriangulationResult {
    std::vector<uint32_t> indices;
    std::vector<std::pair<double, double>> vertices;
    bool usedMeshcut;
    
    TriangulationResult(std::vector<uint32_t> idx, bool meshcut = false) 
        : indices(std::move(idx)), usedMeshcut(meshcut) {}
    TriangulationResult(std::vector<uint32_t> idx, std::vector<std::pair<double, double>> verts, bool meshcut = true) 
        : indices(std::move(idx)), vertices(std::move(verts)), usedMeshcut(meshcut) {}
};

TriangulationResult maplibreTriangulate(const std::vector<double>& polygon, bool useMeshcut = true, uint32_t gridSize = 32) {
    if (useMeshcut) {
        std::cout << "=== MeshCut Triangulation (Replicating MapLibre Integration) ===\n";
        
        // Replicate your exact MapLibre meshcut integration
        meshcut::MeshCutOptions options;
        
        // Configure grid EXACTLY like in your integration
        const double extentDouble = static_cast<double>(EXTENT);  // 8192.0
        options.gridWidth = gridSize;                             // 32
        options.gridHeight = gridSize;                            // 32 
        options.cellSize = (2.0 * extentDouble) / gridSize;      // (2 * 8192)/32 = 512.0
        options.gridOriginX = -extentDouble;                      // -8192
        options.gridOriginY = -extentDouble;                      // -8192
        options.diagonalNE = true;
        
        std::cout << "Grid config: " << options.gridWidth << "x" << options.gridHeight << " cells\n";
        std::cout << "Cell size: " << options.cellSize << " units\n";
        std::cout << "Grid origin: (" << options.gridOriginX << "," << options.gridOriginY << ")\n";
        std::cout << "Grid covers: " << options.gridOriginX << " to " << (options.gridOriginX + options.gridWidth * options.cellSize) << "\n";
        
        // Log input polygon bounds
        double inputMinX = 999999, inputMinY = 999999, inputMaxX = -999999, inputMaxY = -999999;
        for (size_t i = 0; i < polygon.size(); i += 2) {
            double x = polygon[i];
            double y = polygon[i + 1];
            inputMinX = std::min(inputMinX, x);
            inputMaxX = std::max(inputMaxX, x);
            inputMinY = std::min(inputMinY, y);
            inputMaxY = std::max(inputMaxY, y);
        }
        
        std::cout << "Input polygon bounds: (" << inputMinX << "," << inputMinY << ") to (" << inputMaxX << "," << inputMaxY << ")\n";
        
        try {
            auto meshResult = meshcut::meshcut_full(polygon, {}, options);
            
            if (!meshResult.indices.empty() && !meshResult.vertices.empty()) {
                size_t originalVertexCount = polygon.size() / 2;
                size_t meshVertexCount = meshResult.vertices.size() / 2;
                
                std::cout << "MeshCut result: " << originalVertexCount << " -> " << meshVertexCount << " vertices\n";
                std::cout << "Triangles: " << meshResult.indices.size() / 3 << "\n";
                
                // Check coordinate ranges in meshcut output
                double minX = 999999, minY = 999999, maxX = -999999, maxY = -999999;
                size_t overflowCount = 0;
                for (size_t i = 0; i < meshResult.vertices.size(); i += 2) {
                    double x = meshResult.vertices[i];
                    double y = meshResult.vertices[i + 1];
                    minX = std::min(minX, x);
                    maxX = std::max(maxX, x);
                    minY = std::min(minY, y);
                    maxY = std::max(maxY, y);
                    
                    // Check for int16_t overflow (CRITICAL for MapLibre)
                    if (x < std::numeric_limits<int16_t>::min() || x > std::numeric_limits<int16_t>::max() ||
                        y < std::numeric_limits<int16_t>::min() || y > std::numeric_limits<int16_t>::max()) {
                        overflowCount++;
                    }
                }
                
                std::cout << "Output vertex range: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\n";
                std::cout << "⚠️  Coordinate overflow: " << overflowCount << " vertices exceed int16_t range\\n";
                
                // Convert to coordinate pairs
                std::vector<std::pair<double, double>> allVertices;
                allVertices.reserve(meshVertexCount);
                
                for (size_t i = 0; i < meshResult.vertices.size(); i += 2) {
                    allVertices.emplace_back(meshResult.vertices[i], meshResult.vertices[i + 1]);
                }
                
                return TriangulationResult(std::move(meshResult.indices), std::move(allVertices), true);
            }
        } catch (const std::exception& e) {
            std::cout << "MeshCut failed: " << e.what() << ", falling back to earcut\n";
        }
    }
    
    // Earcut fallback
    std::cout << "=== Earcut Triangulation ===\n";
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({});
    for (size_t i = 0; i < polygon.size(); i += 2) {
        earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
    }
    auto earcutIndices = mapbox::earcut<uint32_t>(earcutPolygon);
    
    std::cout << "Earcut result: " << polygon.size() / 2 << " vertices, " << earcutIndices.size() / 3 << " triangles\n";
    
    return TriangulationResult(std::move(earcutIndices), false);
}

// SVG visualization for edge artifacts
class EdgeArtifactVisualizer {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double width, height, scale;
    
public:
    EdgeArtifactVisualizer(const std::string& filename, const std::vector<double>& polygon) {
        // Calculate bounds
        minX = maxX = polygon[0];
        minY = maxY = polygon[1];
        for (size_t i = 0; i < polygon.size(); i += 2) {
            minX = std::min(minX, polygon[i]);
            maxX = std::max(maxX, polygon[i]);
            minY = std::min(minY, polygon[i+1]);
            maxY = std::max(maxY, polygon[i+1]);
        }
        
        // Expand to show tile boundaries
        minX = std::min(minX, 0.0);
        maxX = std::max(maxX, 4096.0);
        minY = std::min(minY, 0.0);
        maxY = std::max(maxY, 4096.0);
        
        // Add some padding
        double padX = (maxX - minX) * 0.05;
        double padY = (maxY - minY) * 0.05;
        minX -= padX; maxX += padX;
        minY -= padY; maxY += padY;
        
        width = maxX - minX;
        height = maxY - minY;
        scale = 600.0 / std::max(width, height);
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"" << (width * scale * 2 + 150) << "\" ";  // Side by side
        file << "height=\"" << (height * scale + 200) << "\" ";
        file << "viewBox=\"0 0 " << (width * scale * 2 + 150) << " " << (height * scale + 200) << "\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    }
    
    ~EdgeArtifactVisualizer() {
        file << "</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale;
        double ty = (maxY - y) * scale;
        if (rightPanel) tx += width * scale + 75;  // Move to right panel
        return {tx + 50, ty + 100};  // Add margins
    }
    
    void drawTitle(const std::string& title, bool rightPanel = false) {
        double titleX = rightPanel ? width * scale + 75 + (width * scale / 2) : (width * scale / 2);
        file << "<text x=\"" << (titleX + 50) << "\" y=\"50\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"18\" font-weight=\"bold\">" << title << "</text>\n";
    }
    
    void drawTileBoundaries(bool rightPanel = false) {
        // Draw tile boundaries at 0, 4096
        file << "<g stroke=\"#ff0000\" stroke-width=\"2\" stroke-dasharray=\"5,5\" fill=\"none\" opacity=\"0.8\">\n";
        
        // Vertical boundaries at x=0 and x=4096
        auto [x_left_top, y_left_top] = transform(0, maxY, rightPanel);
        auto [x_left_bottom, y_left_bottom] = transform(0, minY, rightPanel);
        auto [x_right_top, y_right_top] = transform(4096, maxY, rightPanel);
        auto [x_right_bottom, y_right_bottom] = transform(4096, minY, rightPanel);
        
        file << "<line x1=\"" << x_left_top << "\" y1=\"" << y_left_top << "\" x2=\"" << x_left_bottom << "\" y2=\"" << y_left_bottom << "\"/>\n";
        file << "<line x1=\"" << x_right_top << "\" y1=\"" << y_right_top << "\" x2=\"" << x_right_bottom << "\" y2=\"" << y_right_bottom << "\"/>\n";
        
        // Horizontal boundaries at y=0 and y=4096
        auto [x_top_left, y_top_left] = transform(minX, 0, rightPanel);
        auto [x_top_right, y_top_right] = transform(maxX, 0, rightPanel);
        auto [x_bottom_left, y_bottom_left] = transform(minX, 4096, rightPanel);
        auto [x_bottom_right, y_bottom_right] = transform(maxX, 4096, rightPanel);
        
        file << "<line x1=\"" << x_top_left << "\" y1=\"" << y_top_left << "\" x2=\"" << x_top_right << "\" y2=\"" << y_top_right << "\"/>\n";
        file << "<line x1=\"" << x_bottom_left << "\" y1=\"" << y_bottom_left << "\" x2=\"" << x_bottom_right << "\" y2=\"" << y_bottom_right << "\"/>\n";
        
        file << "</g>\n";
        
        // Add tile boundary labels
        auto [label_x, label_y] = transform(0, maxY - 200, rightPanel);
        file << "<text x=\"" << label_x << "\" y=\"" << label_y << "\" font-family=\"Arial\" font-size=\"12\" fill=\"red\">Tile Edge (0)</text>\n";
        auto [label_x2, label_y2] = transform(4096, maxY - 200, rightPanel);
        file << "<text x=\"" << label_x2 << "\" y=\"" << label_y2 << "\" font-family=\"Arial\" font-size=\"12\" fill=\"red\">Tile Edge (4096)</text>\n";
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color, bool rightPanel = false) {
        if (polygon.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i + 1], rightPanel);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"0.2\" ";
        file << "stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
    }
    
    void drawTriangles(const std::vector<std::pair<double, double>>& vertices, 
                      const std::vector<uint32_t>& indices, 
                      const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.5\" opacity=\"0.7\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1]; 
            uint32_t i2 = indices[i + 2];
            
            if (i0 < vertices.size() && i1 < vertices.size() && i2 < vertices.size()) {
                auto [x0, y0] = transform(vertices[i0].first, vertices[i0].second, rightPanel);
                auto [x1, y1] = transform(vertices[i1].first, vertices[i1].second, rightPanel);
                auto [x2, y2] = transform(vertices[i2].first, vertices[i2].second, rightPanel);
                
                file << "<polygon points=\"" << x0 << "," << y0 << " " 
                     << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
            }
        }
        
        file << "</g>\n";
    }
    
    void drawEarcutTriangles(const std::vector<double>& polygon, 
                            const std::vector<uint32_t>& indices, 
                            const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.5\" opacity=\"0.7\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1]; 
            uint32_t i2 = indices[i + 2];
            
            if (i0*2+1 < polygon.size() && i1*2+1 < polygon.size() && i2*2+1 < polygon.size()) {
                auto [x0, y0] = transform(polygon[i0*2], polygon[i0*2+1], rightPanel);
                auto [x1, y1] = transform(polygon[i1*2], polygon[i1*2+1], rightPanel);
                auto [x2, y2] = transform(polygon[i2*2], polygon[i2*2+1], rightPanel);
                
                file << "<polygon points=\"" << x0 << "," << y0 << " " 
                     << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
            }
        }
        
        file << "</g>\n";
    }
    
    void addStats(const std::string& stats, bool rightPanel = false) {
        double statsX = rightPanel ? width * scale + 75 + 10 : 10;
        file << "<text x=\"" << (statsX + 50) << "\" y=\"" << (height * scale + 150) << "\" ";
        file << "font-family=\"Arial\" font-size=\"12\" fill=\"#666\">";
        file << stats << "</text>\n";
    }
};

// Create test polygons that simulate edge-crossing issues
std::vector<std::vector<double>> createEdgeCrossingTestCases() {
    std::vector<std::vector<double>> testCases;
    
    // Case 1: Polygon that crosses left edge (x=0)
    testCases.push_back({
        -200.0, 1500.0,    // Outside tile (left edge)
        800.0, 1400.0,     // Inside tile
        900.0, 2100.0,     // Inside tile  
        600.0, 2200.0,     // Inside tile
        -100.0, 1800.0     // Outside tile (left edge)
    });
    
    // Case 2: Polygon that crosses right edge (x=4096) 
    testCases.push_back({
        3500.0, 1000.0,    // Inside tile
        4200.0, 1100.0,    // Outside tile (right edge)
        4300.0, 1800.0,    // Outside tile (right edge)
        3600.0, 1900.0,    // Inside tile
        3400.0, 1400.0     // Inside tile
    });
    
    // Case 3: Large polygon spanning multiple grid cells and crossing edges
    testCases.push_back({
        -500.0, -300.0,    // Outside tile (bottom-left)
        4500.0, -200.0,    // Outside tile (bottom-right)
        4600.0, 4400.0,    // Outside tile (top-right)
        -400.0, 4500.0,    // Outside tile (top-left)
        -500.0, -300.0     // Close polygon
    });
    
    // Case 4: Complex polygon with edge touching
    testCases.push_back({
        2000.0, 2000.0,    // Center
        2500.0, 1500.0,    // Inside
        3000.0, 2000.0,    // Inside
        4096.0, 2200.0,    // Exactly on edge
        3500.0, 2800.0,    // Inside
        2500.0, 2500.0,    // Inside
        2000.0, 2000.0     // Close
    });
    
    return testCases;
}

int main() {
    std::cout << "=== MapLibre Edge-Crossing Artifact Analysis ===\n";
    std::cout << "Testing MeshCut vs Earcut with edge-crossing polygons\n\n";
    
    auto testCases = createEdgeCrossingTestCases();
    
    std::vector<std::string> caseNames = {
        "Left Edge Crossing",
        "Right Edge Crossing", 
        "Large Spanning Polygon",
        "Edge Touching Polygon"
    };
    
    for (size_t caseIdx = 0; caseIdx < testCases.size(); caseIdx++) {
        const auto& testPolygon = testCases[caseIdx];
        
        std::cout << "\\n=== Test Case " << (caseIdx + 1) << ": " << caseNames[caseIdx] << " ===\\n";
        std::cout << "Vertices: " << testPolygon.size() / 2 << "\\n";
        
        // Print coordinate range
        double minX = *std::min_element(testPolygon.begin(), testPolygon.end());
        double maxX = *std::max_element(testPolygon.begin(), testPolygon.end());
        double minY = testPolygon[1];
        double maxY = testPolygon[1];
        for (size_t i = 1; i < testPolygon.size(); i += 2) {
            minY = std::min(minY, testPolygon[i]);
            maxY = std::max(maxY, testPolygon[i]);
        }
        
        std::cout << "Coordinate range: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\\n";
        
        // Check if crosses tile boundaries
        bool crossesEdge = (minX < 0 || maxX > 4096 || minY < 0 || maxY > 4096);
        std::cout << "Crosses tile edge: " << (crossesEdge ? "YES" : "NO") << "\\n";
        
        // Test with both triangulation methods
        std::cout << "\\n";
        auto earcutResult = maplibreTriangulate(testPolygon, false);
        std::cout << "\\n";
        auto meshcutResult = maplibreTriangulate(testPolygon, true, 32);
        
        // Analyze differences
        double triangleDensity = (double)meshcutResult.indices.size() / (double)earcutResult.indices.size();
        std::cout << "\\nTriangle density increase: " << std::fixed << std::setprecision(2) << triangleDensity << "x\\n";
        
        // Create visualization
        std::string filename = "edge_artifact_case_" + std::to_string(caseIdx + 1) + ".svg";
        EdgeArtifactVisualizer viz(filename, testPolygon);
        
        // Left panel: Earcut
        viz.drawTitle("Earcut (MapLibre Default)", false);
        viz.drawTileBoundaries(false);
        viz.drawPolygon(testPolygon, "#4CAF50", false);
        viz.drawEarcutTriangles(testPolygon, earcutResult.indices, "#FF5722", false);
        viz.addStats("Earcut: " + std::to_string(earcutResult.indices.size() / 3) + " triangles", false);
        
        // Right panel: MeshCut
        viz.drawTitle("MeshCut (Terrain-Aware)", true);
        viz.drawTileBoundaries(true);
        viz.drawPolygon(testPolygon, "#4CAF50", true);
        viz.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", true);
        viz.addStats("MeshCut: " + std::to_string(meshcutResult.indices.size() / 3) + " triangles", true);
        
        std::cout << "Visualization saved: " << filename << "\\n";
    }
    
    std::cout << "\\n=== ANALYSIS COMPLETE ===\\n";
    std::cout << "Generated " << testCases.size() << " edge-crossing test visualizations\\n";
    std::cout << "\\nKey findings to look for:\\n";
    std::cout << "1. Coordinate overflow beyond int16_t range (-32768 to 32767)\\n";
    std::cout << "2. Different triangulation patterns near tile edges\\n";
    std::cout << "3. Grid artifacts when polygons cross boundaries\\n";
    std::cout << "4. Performance impact from increased triangle density\\n";
    std::cout << "\\nOpen the SVG files to see visual differences!\\n";
    
    return 0;
}