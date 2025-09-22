// Test boundary detection fix
#include "../include/meshcut.hpp"
#include <iostream>
#include <vector>

int main() {
    // Test the boundary detection with the actual problematic polygon
    std::vector<double> polygon = {
        1.2, 1.1,   // Bottom-left-ish
        2.8, 1.0,   // Bottom-right-ish  
        3.0, 2.0,   // Right side
        2.5, 2.8,   // Top-right
        1.8, 3.1,   // Top
        1.0, 2.5,   // Top-left
        0.8, 1.8    // Left side
    };
    
    std::cout << "=== Testing boundary detection fix ===\n";
    
    // Create the polygon tester (internal class, need to access it)
    // Let me just call meshcut and see what happens
    MeshCutOptions options;
    options.gridSpacing = 0.5;
    options.enableSpatialIndex = true;
    
    auto result = meshcut(polygon, options);
    
    std::cout << "MeshCut result:\n";
    std::cout << "- Triangles: " << result.triangles.size() / 3 << "\n";
    std::cout << "- Vertices: " << result.vertices.size() / 2 << "\n";
    
    return 0;
}