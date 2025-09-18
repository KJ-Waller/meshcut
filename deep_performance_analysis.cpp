#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <array>
#include "meshcut/meshcut.hpp"
#include "mapbox/earcut.hpp"

class PerformanceProfiler {
    std::chrono::high_resolution_clock::time_point start_time;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double elapsed_us() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_time).count() / 1000.0;
    }
    
    void print_elapsed(const std::string& label) {
        std::cout << std::setw(30) << label << ": " << std::setw(10) << std::fixed 
                  << std::setprecision(2) << elapsed_us() << " μs\n";
    }
};

void profile_meshcut_components() {
    std::cout << "🔍 MeshCut Component-Level Performance Analysis\n";
    std::cout << "==============================================\n";
    
    // Test with a simple square that should be fast
    meshcut::Polygon square = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    int gridN = 4;  // 4x4 grid - should be fast
    
    PerformanceProfiler profiler;
    
    // Test individual components
    std::cout << "\n📊 Component Timing (4x4 grid, simple square):\n";
    
    // 1. Grid generation
    profiler.start();
    auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, gridN, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    profiler.print_elapsed("1. Build grid triangles");
    
    // 2. Classification 
    profiler.start();
    auto isInterior = meshcut::internal::classifyTriangles(gridTriangles, square);
    profiler.print_elapsed("2. Classify triangles");
    
    // 3. Find partial intersections
    profiler.start();
    auto isPartiallyIntersected = meshcut::internal::findPartiallyIntersectedTriangles(gridTriangles, square, isInterior);
    profiler.print_elapsed("3. Find partial intersections");
    
    // 4. Process partial triangles
    profiler.start();
    std::vector<meshcut::Triangle> partialTriangles;
    int earcut_calls = 0;
    for (size_t i = 0; i < gridTriangles.size(); ++i) {
        if (isPartiallyIntersected[i]) {
            auto intersection = meshcut::internal::computeTrianglePolygonIntersection(gridTriangles[i], square);
            if (intersection.size() >= 3) {
                auto triangulated = meshcut::internal::earcutTriangulate(intersection);
                partialTriangles.insert(partialTriangles.end(), triangulated.begin(), triangulated.end());
                earcut_calls++;
            }
        }
    }
    profiler.print_elapsed("4. Process partial triangles");
    
    // 5. Complete tessellation for comparison
    profiler.start();
    auto complete_result = meshcut::tessellate(square, bbox, gridN);
    profiler.print_elapsed("5. Complete tessellate() call");
    
    std::cout << "\n📈 Analysis:\n";
    std::cout << "  Grid triangles: " << gridTriangles.size() << "\n";
    std::cout << "  Interior triangles: " << std::count(isInterior.begin(), isInterior.end(), true) << "\n";
    std::cout << "  Partial triangles: " << std::count(isPartiallyIntersected.begin(), isPartiallyIntersected.end(), true) << "\n";
    std::cout << "  Earcut calls made: " << earcut_calls << "\n";
    std::cout << "  Final triangles: " << complete_result.size() << "\n";
}

void test_earcut_overhead() {
    std::cout << "\n🔬 Earcut Call Overhead Analysis\n";
    std::cout << "================================\n";
    
    // Test single large vs multiple small Earcut calls
    meshcut::Polygon large_polygon = {{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}};
    
    PerformanceProfiler profiler;
    
    // Single large Earcut call
    profiler.start();
    std::vector<std::vector<std::array<double, 2>>> earcut_polygon = {
        {{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}}
    };
    auto single_result = mapbox::earcut<uint32_t>(earcut_polygon);
    profiler.print_elapsed("Single large Earcut call");
    
    // Multiple small Earcut calls (simulate MeshCut's approach)
    profiler.start();
    int num_calls = 16;  // Simulate processing 16 partial triangles
    std::vector<uint32_t> multi_result;
    for (int i = 0; i < num_calls; ++i) {
        std::vector<std::vector<std::array<double, 2>>> small_polygon = {
            {{0.0, 0.0}, {1.0, 0.0}, {0.5, 1.0}}  // Small triangle
        };
        auto result = mapbox::earcut<uint32_t>(small_polygon);
        multi_result.insert(multi_result.end(), result.begin(), result.end());
    }
    profiler.print_elapsed("Multiple small Earcut calls");
    
    std::cout << "  Single call result: " << single_result.size() / 3 << " triangles\n";
    std::cout << "  Multi call result: " << multi_result.size() / 3 << " triangles\n";
}

void test_realistic_tile_scenario() {
    std::cout << "\n🗺️  Realistic MapLibre Tile Scenario\n";
    std::cout << "===================================\n";
    
    // Create a polygon that spans part of a tile (realistic building footprint)
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.4, 0.1}, {0.4, 0.3}, {0.3, 0.3}, 
        {0.3, 0.4}, {0.4, 0.4}, {0.4, 0.6}, {0.1, 0.6}
    };
    
    // Tile bounds (0,0) to (1,1) with 32x32 grid
    meshcut::Polygon tile_bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    PerformanceProfiler profiler;
    
    std::cout << "Testing building footprint on 32x32 tile grid:\n";
    
    // Test with different grid sizes
    std::vector<int> grid_sizes = {4, 8, 16, 32};
    
    for (int grid_size : grid_sizes) {
        profiler.start();
        auto result = meshcut::tessellate(building, tile_bbox, grid_size);
        double time_us = profiler.elapsed_us();
        
        std::cout << "  " << grid_size << "x" << grid_size << " grid: " 
                  << std::setw(8) << std::fixed << std::setprecision(2) << time_us 
                  << " μs, " << result.size() << " triangles\n";
    }
    
    // Compare with single Earcut call
    profiler.start();
    std::vector<std::vector<std::array<double, 2>>> earcut_building;
    earcut_building.push_back({});
    for (const auto& pt : building) {
        earcut_building[0].push_back({pt.x, pt.y});
    }
    auto earcut_result = mapbox::earcut<uint32_t>(earcut_building);
    profiler.print_elapsed("Earcut (same polygon)");
    std::cout << "                     " << earcut_result.size() / 3 << " triangles\n";
}

void analyze_algorithm_complexity() {
    std::cout << "\n📐 Algorithm Complexity Analysis\n";
    std::cout << "===============================\n";
    
    PerformanceProfiler profiler;
    
    // Test how performance scales with polygon complexity
    std::vector<int> polygon_sizes = {4, 8, 16, 32, 64};
    
    std::cout << "Performance vs Polygon Complexity (8x8 grid):\n";
    std::cout << "Vertices | MeshCut (μs) | Earcut (μs) | Ratio\n";
    std::cout << "---------|--------------|-------------|------\n";
    
    for (int n : polygon_sizes) {
        // Create regular n-gon
        meshcut::Polygon polygon;
        std::vector<std::vector<std::array<double, 2>>> earcut_polygon = {{}};
        
        for (int i = 0; i < n; ++i) {
            double angle = 2.0 * M_PI * i / n;
            double x = 0.5 + 0.4 * std::cos(angle);
            double y = 0.5 + 0.4 * std::sin(angle);
            polygon.push_back({x, y});
            earcut_polygon[0].push_back({x, y});
        }
        
        meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
        
        // Test MeshCut
        profiler.start();
        auto meshcut_result = meshcut::tessellate(polygon, bbox, 8);
        double meshcut_time = profiler.elapsed_us();
        
        // Test Earcut
        profiler.start();
        auto earcut_result = mapbox::earcut<uint32_t>(earcut_polygon);
        double earcut_time = profiler.elapsed_us();
        
        double ratio = meshcut_time / earcut_time;
        
        std::cout << std::setw(8) << n << " | " 
                  << std::setw(12) << std::fixed << std::setprecision(1) << meshcut_time << " | "
                  << std::setw(11) << std::fixed << std::setprecision(1) << earcut_time << " | "
                  << std::setw(5) << std::fixed << std::setprecision(1) << ratio << "x\n";
    }
}

int main() {
    std::cout << "🚀 Deep MeshCut Performance Investigation\n";
    std::cout << "========================================\n";
    std::cout << "Investigating why MeshCut appears slow for MapLibre terrain use case\n\n";
    
    profile_meshcut_components();
    test_earcut_overhead();
    test_realistic_tile_scenario();
    analyze_algorithm_complexity();
    
    std::cout << "\n🎯 Summary:\n";
    std::cout << "This analysis should reveal where the performance bottleneck is:\n";
    std::cout << "1. Grid generation (should be fast)\n";
    std::cout << "2. Triangle classification (should be fast)\n"; 
    std::cout << "3. Partial intersection processing (potential bottleneck)\n";
    std::cout << "4. Multiple Earcut calls (potential bottleneck)\n";
    std::cout << "5. Overall algorithm scaling\n";
    
    return 0;
}