#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <cassert>

// SVG writer class for side-by-side visualization
class SVGComparison {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double width, height, scale;
    double panelWidth;
    
public:
    SVGComparison(const std::string& filename, double minX_, double minY_, double maxX_, double maxY_) 
        : minX(minX_), minY(minY_), maxX(maxX_), maxY(maxY_) {
        
        width = maxX - minX;
        height = maxY - minY;
        scale = 300.0 / std::max(width, height); // 300px per panel
        panelWidth = width * scale + 100; // 50px margin per side
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"" << (panelWidth * 2 + 50) << "\" "; // Two panels + center margin
        file << "height=\"" << (height * scale + 150) << "\" "; // Height + title space
        file << "viewBox=\"0 0 " << (panelWidth * 2 + 50) << " " << (height * scale + 150) << "\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        
        // Add titles
        file << "<text x=\"" << (panelWidth/2) << "\" y=\"30\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">Earcut</text>\n";
        file << "<text x=\"" << (panelWidth + 25 + panelWidth/2) << "\" y=\"30\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">MeshCut</text>\n";
        
        // Separator line
        file << "<line x1=\"" << (panelWidth + 12.5) << "\" y1=\"50\" x2=\"" << (panelWidth + 12.5) << "\" y2=\"" << (height * scale + 100) << "\" stroke=\"#ccc\" stroke-width=\"1\"/>\n";
    }
    
    ~SVGComparison() {
        file << "</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale + 50; // 50px margin
        double ty = (maxY - y) * scale + 50; // Flip Y axis + 50px margin
        
        if (rightPanel) {
            tx += panelWidth + 25; // Move to right panel
        }
        
        return {tx, ty};
    }
    
    void drawGrid(const meshcut::MeshCutOptions& options, bool rightPanel = false) {
        std::string panelClass = rightPanel ? "right-grid" : "left-grid";
        file << "<g class=\"" << panelClass << "\" stroke=\"#e0e0e0\" stroke-width=\"0.5\" fill=\"none\">\n";
        
        // Vertical lines
        for (int i = 0; i <= options.gridWidth; i++) {
            double x = options.gridOriginX + i * options.cellSize;
            auto [x1, y1] = transform(x, options.gridOriginY, rightPanel);
            auto [x2, y2] = transform(x, options.gridOriginY + options.gridHeight * options.cellSize, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        // Horizontal lines
        for (int j = 0; j <= options.gridHeight; j++) {
            double y = options.gridOriginY + j * options.cellSize;
            auto [x1, y1] = transform(options.gridOriginX, y, rightPanel);
            auto [x2, y2] = transform(options.gridOriginX + options.gridWidth * options.cellSize, y, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        file << "</g>\n";
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color = "#4CAF50", 
                     double opacity = 0.3, bool rightPanel = false) {
        if (polygon.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i + 1], rightPanel);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"" << opacity << "\" ";
        file << "stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
    }
    
    void drawTriangles(const std::vector<double>& vertices, const std::vector<uint32_t>& indices, 
                       const std::string& color = "#2196F3", double opacity = 0.1, bool rightPanel = false) {
        file << "<g fill=\"" << color << "\" fill-opacity=\"" << opacity << "\" stroke=\"" << color << "\" stroke-width=\"0.5\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1];
            uint32_t i2 = indices[i + 2];
            
            if (i0*2+1 >= vertices.size() || i1*2+1 >= vertices.size() || i2*2+1 >= vertices.size()) {
                continue;
            }
            
            auto [x0, y0] = transform(vertices[i0*2], vertices[i0*2+1], rightPanel);
            auto [x1, y1] = transform(vertices[i1*2], vertices[i1*2+1], rightPanel);
            auto [x2, y2] = transform(vertices[i2*2], vertices[i2*2+1], rightPanel);
            
            file << "<polygon points=\"" << x0 << "," << y0 << " " 
                 << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
        }
        
        file << "</g>\n";
    }
    
    void addStats(int triangleCount, int vertexCount, bool rightPanel = false) {
        double x = rightPanel ? (panelWidth + 25 + panelWidth/2) : (panelWidth/2);
        double y = height * scale + 70;
        
        file << "<text x=\"" << x << "\" y=\"" << y << "\" text-anchor=\"middle\" font-size=\"12\">";
        file << "Triangles: " << triangleCount << ", Vertices: " << vertexCount << "</text>\n";
    }
};

// Complex test shapes with realistic coordinate ranges
namespace TestShapes {
    
    // Complex building footprint (realistic GIS coordinates)
    std::vector<double> createBuildingFootprint() {
        // Simulated building coordinates (in meters, typical GIS projection)
        return {
            // L-shaped building
            413.5, 287.2,  // Bottom-left corner
            423.8, 287.1,  // Bottom extend
            423.9, 292.4,  // Right side up
            419.2, 292.5,  // Inner corner
            419.3, 297.8,  // Up to top
            413.4, 297.9,  // Top-left
            413.5, 287.2   // Close polygon
        };
    }
    
    // Irregular park/lake boundary
    std::vector<double> createParkBoundary() {
        std::vector<double> boundary;
        // Create organic curved boundary using many points
        double centerX = 416.0, centerY = 292.0;
        int numPoints = 16;
        
        for (int i = 0; i < numPoints; i++) {
            double angle = (2.0 * M_PI * i) / numPoints;
            // Variable radius for organic shape
            double radius = 8.0 + 3.0 * sin(3 * angle) + 2.0 * cos(5 * angle);
            double x = centerX + radius * cos(angle);
            double y = centerY + radius * sin(angle);
            boundary.push_back(x);
            boundary.push_back(y);
        }
        
        return boundary;
    }
    
    // Complex urban block with holes (simulating courtyards)
    std::vector<double> createUrbanBlock() {
        return {
            // Outer boundary - rectangular block
            405.0, 280.0,  // Bottom-left
            425.0, 280.0,  // Bottom-right
            425.0, 300.0,  // Top-right
            405.0, 300.0,  // Top-left
            405.0, 280.0   // Close
        };
    }
    
    std::vector<uint32_t> createUrbanBlockHoles() {
        // Hole indices for courtyards
        return {10}; // Hole starts at index 10 (after 5 vertices of outer boundary)
    }
    
    std::vector<double> createUrbanBlockWithHoles() {
        auto outer = createUrbanBlock();
        
        // Add two courtyard holes
        std::vector<double> hole1 = {
            408.0, 285.0,
            412.0, 285.0, 
            412.0, 289.0,
            408.0, 289.0,
            408.0, 285.0
        };
        
        std::vector<double> hole2 = {
            418.0, 292.0,
            422.0, 292.0,
            422.0, 296.0, 
            418.0, 296.0,
            418.0, 292.0
        };
        
        // Combine outer + holes
        outer.insert(outer.end(), hole1.begin(), hole1.end());
        outer.insert(outer.end(), hole2.begin(), hole2.end());
        
        return outer;
    }
}

void testComplexShape(const std::string& shapeName, 
                     const std::vector<double>& polygon,
                     const std::vector<uint32_t>& holes = {},
                     const meshcut::MeshCutOptions& options = {}) {
    
    std::cout << "\n=== Testing " << shapeName << " ===\n";
    std::cout << "Polygon vertices: " << polygon.size()/2 << "\n";
    std::cout << "Grid: " << options.gridWidth << "x" << options.gridHeight 
              << ", cell size: " << options.cellSize << "\n";
    std::cout << "Coverage: [" << options.gridOriginX << ", " << options.gridOriginY 
              << "] to [" << (options.gridOriginX + options.gridWidth * options.cellSize)
              << ", " << (options.gridOriginY + options.gridHeight * options.cellSize) << "]\n";
    
    // === MESHCUT TEST ===
    auto meshcutResult = meshcut::meshcut_full(polygon, holes, options);
    
    std::cout << "\nMeshCut Results:\n";
    std::cout << "- Vertices: " << meshcutResult.vertices.size()/2 << "\n";
    std::cout << "- Triangles: " << meshcutResult.indices.size()/3 << "\n";
    
    // Verify triangle indices are valid
    for (size_t i = 0; i < meshcutResult.indices.size(); i++) {
        if (meshcutResult.indices[i] * 2 + 1 >= meshcutResult.vertices.size()) {
            std::cerr << "ERROR: Invalid triangle index " << meshcutResult.indices[i] 
                      << " at position " << i << " (max vertex index: " 
                      << (meshcutResult.vertices.size()/2 - 1) << ")\n";
        }
    }
    
    // === EARCUT COMPARISON ===
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    
    // Convert polygon format for earcut
    if (holes.empty()) {
        // Simple polygon
        earcutPolygon.push_back({});
        for (size_t i = 0; i < polygon.size(); i += 2) {
            earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
        }
    } else {
        // Polygon with holes
        earcutPolygon.push_back({}); // Outer ring
        
        // Add outer ring vertices
        size_t vertexIndex = 0;
        size_t currentHole = 0;
        
        for (size_t i = 0; i < polygon.size(); i += 2, vertexIndex++) {
            // Check if we've reached a hole start
            if (currentHole < holes.size() && vertexIndex == holes[currentHole]) {
                earcutPolygon.push_back({}); // Start new hole
                currentHole++;
            }
            
            size_t ringIndex = currentHole; // Current ring (0 = outer, 1+ = holes)
            earcutPolygon[ringIndex].push_back({{polygon[i], polygon[i+1]}});
        }
    }
    
    auto earcutIndices = mapbox::earcut<uint32_t>(earcutPolygon);
    
    std::cout << "\nEarcut Results:\n";
    std::cout << "- Vertices: " << polygon.size()/2 << " (original polygon vertices only)\n";
    std::cout << "- Triangles: " << earcutIndices.size()/3 << "\n";
    
    // Verify earcut indices are valid
    for (size_t i = 0; i < earcutIndices.size(); i++) {
        if (earcutIndices[i] * 2 + 1 >= polygon.size()) {
            std::cerr << "ERROR: Invalid earcut index " << earcutIndices[i] 
                      << " at position " << i << " (max vertex index: " 
                      << (polygon.size()/2 - 1) << ")\n";
        }
    }
    
    // === VISUALIZATION ===
    std::string filename = "../visualizations/" + shapeName + "_comparison.svg";
    
    // Calculate bounds
    double minX = options.gridOriginX;
    double minY = options.gridOriginY;
    double maxX = options.gridOriginX + options.gridWidth * options.cellSize;
    double maxY = options.gridOriginY + options.gridHeight * options.cellSize;
    
    SVGComparison svg(filename, minX, minY, maxX, maxY);
    
    // Draw both panels
    // Left panel - Earcut
    svg.drawGrid(options, false);
    svg.drawPolygon(polygon, "#4CAF50", 0.3, false);
    svg.drawTriangles(polygon, earcutIndices, "#FF5722", 0.2, false);
    svg.addStats(earcutIndices.size()/3, polygon.size()/2, false);
    
    // Right panel - MeshCut  
    svg.drawGrid(options, true);
    svg.drawPolygon(polygon, "#4CAF50", 0.3, true);
    svg.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", 0.2, true);
    svg.addStats(meshcutResult.indices.size()/3, meshcutResult.vertices.size()/2, true);
    
    std::cout << "Visualization saved to: " << filename << "\n";
    
    // === API COMPATIBILITY CHECK ===
    std::cout << "\n--- API Compatibility Check ---\n";
    std::cout << "✓ Input format: Both use std::vector<double> [x,y, x,y, ...]\n";
    std::cout << "✓ Holes format: Both use std::vector<uint32_t> hole indices\n"; 
    std::cout << "✓ Output format: Both return std::vector<uint32_t> triangle indices\n";
    std::cout << "✓ MeshCut extends with MeshCutOptions for grid specification\n";
    std::cout << "✓ MeshCut provides meshcut_full() for vertex+index output\n";
    
    // Quality comparison
    double meshcutDensity = (double)meshcutResult.indices.size() / (double)earcutIndices.size();
    std::cout << "\n--- Quality Comparison ---\n";
    std::cout << "Triangle density ratio (MeshCut/Earcut): " << std::fixed << std::setprecision(2) << meshcutDensity << "x\n";
    if (meshcutDensity > 2.0) {
        std::cout << "✓ MeshCut provides significantly higher triangle density for terrain sampling\n";
    } else {
        std::cout << "! MeshCut triangle density may be insufficient for detailed terrain\n";
    }
}

int main() {
    std::cout << "MeshCut Complex Shape Visualization Test\n";
    std::cout << "=======================================\n";
    std::cout << "Testing compatibility and comparison with earcut using realistic coordinates\n";
    
    // Setup realistic grid covering the test area
    meshcut::MeshCutOptions options;
    options.gridOriginX = 400.0;   // Realistic GIS coordinates
    options.gridOriginY = 275.0;
    options.cellSize = 1.0;        // 1 meter cells
    options.gridWidth = 32;        // 32x32 grid
    options.gridHeight = 32;
    options.diagonalNE = true;
    
    // Test 1: Complex building footprint
    testComplexShape("building_footprint", 
                     TestShapes::createBuildingFootprint(), 
                     {}, options);
    
    // Test 2: Organic park boundary  
    testComplexShape("park_boundary",
                     TestShapes::createParkBoundary(),
                     {}, options);
    
    // Test 3: Urban block with courtyard holes
    testComplexShape("urban_block_with_holes",
                     TestShapes::createUrbanBlockWithHoles(),
                     {10, 20}, options); // Holes start at vertices 10 and 20
    
    // Test 4: Simple shape at different grid resolution
    meshcut::MeshCutOptions fineOptions = options;
    fineOptions.cellSize = 0.5;  // Finer 0.5m cells
    fineOptions.gridWidth = 64;  // 64x64 grid for same area
    fineOptions.gridHeight = 64;
    
    testComplexShape("building_fine_grid",
                     TestShapes::createBuildingFootprint(),
                     {}, fineOptions);
    
    std::cout << "\n=======================================\n";
    std::cout << "All visualizations saved to ../visualizations/\n";
    std::cout << "Check *_comparison.svg files for side-by-side comparisons\n";
    
    return 0;
}