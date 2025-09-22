#include "meshcut.hpp"
#include <iostream>
#include <iomanip>

int main() {
    std::cout << "🔍 MeshCut Algorithm Validation Test\n";
    std::cout << "====================================\n\n";
    
    // Use the same test polygon as the visualization test
    std::vector<double> polygon = {
        1.5, 1.0,  // Start point
        2.5, 0.5,  
        3.0, 1.8,
        2.2, 2.8,
        1.2, 2.5,
        0.8, 1.8,
        1.0, 1.2
    };
    
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0;
    options.cellSize = 0.5;
    options.gridWidth = 8;
    options.gridHeight = 8;
    options.diagonalNE = true;
    
    std::cout << "Test polygon: " << polygon.size()/2 << " vertices\n";
    std::cout << "Grid: " << options.gridWidth << "x" << options.gridHeight << " cells of size " << options.cellSize << "\n\n";
    
    // Test the optimized version
    auto result = meshcut::meshcut_full(polygon, {}, options);
    
    std::cout << "MeshCut Results:\n";
    std::cout << "- Vertices: " << result.vertices.size()/2 << "\n";
    std::cout << "- Triangles: " << result.indices.size()/3 << "\n\n";
    
    // Show some vertex coordinates to verify
    std::cout << "First 10 vertices:\n";
    for (size_t i = 0; i < std::min(result.vertices.size(), 20ul); i += 2) {
        std::cout << "  (" << std::fixed << std::setprecision(3) 
                  << result.vertices[i] << ", " << result.vertices[i+1] << ")\n";
    }
    
    std::cout << "\nFirst 10 triangles (indices):\n";
    for (size_t i = 0; i < std::min(result.indices.size(), 30ul); i += 3) {
        std::cout << "  [" << result.indices[i] << ", " << result.indices[i+1] << ", " << result.indices[i+2] << "]\n";
    }
    
    // Compare with earcut for reference
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({});
    for (size_t i = 0; i < polygon.size(); i += 2) {
        earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
    }
    
    auto earcutResult = mapbox::earcut<uint32_t>(earcutPolygon);
    std::cout << "\nEarcut comparison: " << earcutResult.size()/3 << " triangles\n";
    
    return 0;
}