#include "meshcut.hpp"
#include <iostream>
#include <vector>
#include <array>

int main() {
    std::cout << "MeshCut Basic Test\n";
    
    // Test 1: Simple square polygon
    std::vector<double> square = {
        0.0, 0.0,   // bottom-left
        2.0, 0.0,   // bottom-right  
        2.0, 2.0,   // top-right
        0.0, 2.0    // top-left
    };
    
    // Grid options - 4x4 grid covering the square
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0; 
    options.cellSize = 0.5;     // Each cell is 0.5x0.5
    options.gridWidth = 4;      
    options.gridHeight = 4;
    options.diagonalNE = true;  // Use / diagonal
    
    // Run meshcut
    std::vector<uint32_t> triangles = meshcut::meshcut(square, {}, options);
    
    std::cout << "Input polygon: 4 vertices (square)\n";
    std::cout << "Grid: " << options.gridWidth << "x" << options.gridHeight 
              << " cells of size " << options.cellSize << "\n";
    std::cout << "Output triangles: " << triangles.size() / 3 << " triangles\n";
    std::cout << "Triangle indices: ";
    for (size_t i = 0; i < triangles.size(); i++) {
        if (i % 3 == 0 && i > 0) std::cout << " | ";
        std::cout << triangles[i] << " ";
    }
    std::cout << "\n";
    
    // Test 2: Test with earcut-compatible polygon format
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({
        {{1.0, 1.0}},
        {{3.0, 1.0}}, 
        {{3.0, 3.0}},
        {{1.0, 3.0}}
    });
    
    auto triangles2 = meshcut::meshcut<uint32_t>(earcutPolygon, options);
    std::cout << "\nEarcut-format test: " << triangles2.size() / 3 << " triangles\n";
    
    std::cout << "\nBasic test completed successfully!\n";
    return 0;
}