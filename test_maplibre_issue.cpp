#include "include/meshcut.hpp"
#include <iostream>

int main() {
    // Simple triangle polygon - 3 vertices
    std::vector<double> triangle = {0.0, 0.0, 100.0, 0.0, 50.0, 100.0};
    
    meshcut::MeshCutOptions options;
    options.gridWidth = 4;
    options.gridHeight = 4;
    options.cellSize = 25.0;
    
    std::cout << "Original polygon vertices: " << triangle.size() / 2 << std::endl;
    
    // Test meshcut_full
    auto full_result = meshcut::meshcut_full(triangle, {}, options);
    std::cout << "MeshCut full result vertices: " << full_result.vertices.size() / 2 << std::endl;
    std::cout << "MeshCut indices count: " << full_result.indices.size() << std::endl;
    
    // Test meshcut (indices only)
    auto indices = meshcut::meshcut(triangle, {}, options);
    std::cout << "MeshCut indices-only result: " << indices.size() << " indices" << std::endl;
    
    // Check if any index is >= original vertex count
    size_t original_vertex_count = triangle.size() / 2;
    bool has_invalid_indices = false;
    for (auto idx : indices) {
        if (idx >= original_vertex_count) {
            has_invalid_indices = true;
            std::cout << "❌ INVALID INDEX: " << idx << " >= " << original_vertex_count << std::endl;
            break;
        }
    }
    
    if (!has_invalid_indices) {
        std::cout << "✅ All indices are valid for original polygon" << std::endl;
    }
    
    return 0;
}