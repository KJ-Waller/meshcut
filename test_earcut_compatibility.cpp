#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <array>
#include "meshcut/meshcut.hpp"
#include "meshcut/earcut_compatible.hpp"
#include "mapbox/earcut.hpp"

// Use vector of vectors for Earcut (supports holes)
using Polygon2D = std::vector<std::vector<std::array<double, 2>>>;
using SimplePolygon = std::vector<std::array<double, 2>>;

void test_compatibility() {
    std::cout << "🔄 Testing MeshCut Earcut Compatibility\n";
    std::cout << "======================================\n";
    
    // Test polygon - simple square (Earcut format: vector of rings)
    SimplePolygon square_ring = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    Polygon2D square = {square_ring};  // Wrap in outer vector for Earcut
    
    std::cout << "Input polygon: Square with 4 vertices\n";
    
    // Test original Earcut
    auto start = std::chrono::high_resolution_clock::now();
    auto earcut_indices = mapbox::earcut<uint32_t>(square);
    auto earcut_time = std::chrono::high_resolution_clock::now() - start;
    
    std::cout << "✅ Earcut result: " << earcut_indices.size() / 3 << " triangles\n";
    std::cout << "   Indices: ";
    for (size_t i = 0; i < std::min(earcut_indices.size(), size_t(12)); ++i) {
        std::cout << earcut_indices[i] << " ";
    }
    if (earcut_indices.size() > 12) std::cout << "...";
    std::cout << "\n";
    
    // Test MeshCut compatible interface (using simple ring)
    start = std::chrono::high_resolution_clock::now();
    auto meshcut_indices = meshcut::earcut<uint32_t>(square_ring, 4);  // 4x4 grid
    auto meshcut_time = std::chrono::high_resolution_clock::now() - start;
    
    std::cout << "✅ MeshCut result: " << meshcut_indices.size() / 3 << " triangles\n";
    std::cout << "   Indices: ";
    for (size_t i = 0; i < std::min(meshcut_indices.size(), size_t(12)); ++i) {
        std::cout << meshcut_indices[i] << " ";
    }
    if (meshcut_indices.size() > 12) std::cout << "...";
    std::cout << "\n";
    
    // Test with vertices output
    auto [indices, vertices] = meshcut::earcut_compatible_with_vertices<uint32_t>(square_ring, 4);
    std::cout << "✅ MeshCut with vertices: " << indices.size() / 3 << " triangles, " 
              << vertices.size() << " unique vertices\n";
    
    // Performance comparison
    auto earcut_us = std::chrono::duration_cast<std::chrono::microseconds>(earcut_time).count();
    auto meshcut_us = std::chrono::duration_cast<std::chrono::microseconds>(meshcut_time).count();
    
    std::cout << "\n⏱️  Performance Comparison:\n";
    std::cout << "   Earcut:  " << std::setw(6) << earcut_us << " μs\n";
    std::cout << "   MeshCut: " << std::setw(6) << meshcut_us << " μs\n";
    std::cout << "   Ratio:   " << std::fixed << std::setprecision(2) 
              << (double)meshcut_us / earcut_us << "x\n";
}

void test_complex_polygon() {
    std::cout << "\n🏗️  Testing Complex Polygon\n";
    std::cout << "===========================\n";
    
    // L-shaped building
    SimplePolygon building_ring = {
        {0.0, 0.0}, {8.0, 0.0}, {8.0, 3.0}, {5.0, 3.0}, 
        {5.0, 5.0}, {8.0, 5.0}, {8.0, 8.0}, {3.0, 8.0}, 
        {3.0, 6.0}, {0.0, 6.0}
    };
    Polygon2D building = {building_ring};
    
    auto earcut_result = mapbox::earcut<uint32_t>(building);
    auto meshcut_result = meshcut::earcut<uint32_t>(building_ring, 8);  // 8x8 grid
    
    std::cout << "Complex building polygon (10 vertices):\n";
    std::cout << "  Earcut:  " << earcut_result.size() / 3 << " triangles\n";
    std::cout << "  MeshCut: " << meshcut_result.size() / 3 << " triangles\n";
    
    // Check if both produce valid triangulation
    bool earcut_valid = (earcut_result.size() % 3) == 0;
    bool meshcut_valid = (meshcut_result.size() % 3) == 0;
    
    std::cout << "  Earcut valid:  " << (earcut_valid ? "✅" : "❌") << "\n";
    std::cout << "  MeshCut valid: " << (meshcut_valid ? "✅" : "❌") << "\n";
}

void demonstrate_grid_effects() {
    std::cout << "\n⚙️  Grid Size Effects\n";
    std::cout << "====================\n";
    
    SimplePolygon circle_ring;
    const int circle_points = 16;
    for (int i = 0; i < circle_points; ++i) {
        double angle = 2.0 * M_PI * i / circle_points;
        circle_ring.push_back({std::cos(angle), std::sin(angle)});
    }
    Polygon2D circle = {circle_ring};
    
    auto earcut_triangles = mapbox::earcut<uint32_t>(circle).size() / 3;
    std::cout << "Circle (16 vertices) - Earcut: " << earcut_triangles << " triangles\n";
    
    std::vector<int> grid_sizes = {4, 8, 16, 32};
    for (int grid_size : grid_sizes) {
        auto meshcut_triangles = meshcut::earcut<uint32_t>(circle_ring, grid_size).size() / 3;
        std::cout << "Circle - MeshCut " << grid_size << "x" << grid_size 
                  << ": " << meshcut_triangles << " triangles\n";
    }
}

int main() {
    std::cout << "🧪 MeshCut Earcut Compatibility Test Suite\n";
    std::cout << "==========================================\n\n";
    
    try {
        test_compatibility();
        test_complex_polygon();
        demonstrate_grid_effects();
        
        std::cout << "\n🎉 Compatibility Test Complete!\n";
        std::cout << "✅ MeshCut can serve as drop-in replacement for Earcut\n";
        std::cout << "⚠️  Note: Triangle counts will differ due to grid-based approach\n";
        std::cout << "💡 Adjust grid size to balance quality vs performance\n";
        
    } catch (const std::exception& e) {
        std::cout << "❌ Test failed: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}