#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include "meshcut/meshcut.hpp"
#include "mapbox/earcut.hpp"

class Timer {
    std::chrono::high_resolution_clock::time_point start;
public:
    void reset() { start = std::chrono::high_resolution_clock::now(); }
    double elapsed_us() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000.0;
    }
};

void test_earcut_overhead() {
    std::cout << "🔬 Testing Earcut Function Call Overhead\n";
    std::cout << "========================================\n";
    
    Timer timer;
    
    // Create one moderately complex polygon
    std::vector<std::vector<std::array<double, 2>>> large_poly = {
        {{0.0, 0.0}, {1.0, 0.0}, {1.0, 0.5}, {0.8, 0.5}, {0.8, 0.8}, {1.0, 0.8}, {1.0, 1.0}, {0.0, 1.0}}
    };
    
    // Create 20 tiny triangles (simulating our partial triangles)
    std::vector<std::vector<std::vector<std::array<double, 2>>>> tiny_polys;
    for (int i = 0; i < 20; ++i) {
        double offset = i * 0.05;
        tiny_polys.push_back({{{offset, 0.0}, {offset + 0.04, 0.0}, {offset + 0.02, 0.04}}});
    }
    
    // Test single large call
    timer.reset();
    auto large_result = mapbox::earcut<uint32_t>(large_poly);
    double large_time = timer.elapsed_us();
    
    // Test many small calls
    timer.reset();
    std::vector<uint32_t> small_results;
    for (auto& poly : tiny_polys) {
        auto result = mapbox::earcut<uint32_t>(poly);
        small_results.insert(small_results.end(), result.begin(), result.end());
    }
    double small_time = timer.elapsed_us();
    
    std::cout << "One large polygon (8 vertices): " << large_time << " μs → " << large_result.size()/3 << " triangles\n";
    std::cout << "20 tiny polygons (3 vertices each): " << small_time << " μs → " << small_results.size()/3 << " triangles\n";
    std::cout << "Overhead ratio: " << small_time / large_time << "x\n\n";
    
    // Test pure function call overhead
    timer.reset();
    for (int i = 0; i < 1000; ++i) {
        mapbox::earcut<uint32_t>(tiny_polys[0]);
    }
    double overhead_time = timer.elapsed_us();
    std::cout << "1000 identical tiny Earcut calls: " << overhead_time << " μs\n";
    std::cout << "Average per call: " << overhead_time / 1000.0 << " μs overhead\n\n";
}

void test_clipper_overhead() {
    std::cout << "⚙️ Testing Clipper2 Intersection Overhead\n";
    std::cout << "=========================================\n";
    
    Timer timer;
    meshcut::Polygon square = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    // Create a single triangle for intersection
    meshcut::Triangle test_triangle{{0.1, 0.1}, {0.2, 0.1}, {0.15, 0.2}};
    
    // Test single intersection
    timer.reset();
    auto intersection = meshcut::internal::computeTrianglePolygonIntersection(test_triangle, square);
    double single_time = timer.elapsed_us();
    
    // Test 100 identical intersections 
    timer.reset();
    for (int i = 0; i < 100; ++i) {
        auto result = meshcut::internal::computeTrianglePolygonIntersection(test_triangle, square);
    }
    double multiple_time = timer.elapsed_us();
    
    std::cout << "Single Clipper intersection: " << single_time << " μs\n";
    std::cout << "100 identical intersections: " << multiple_time << " μs\n";
    std::cout << "Average per intersection: " << multiple_time / 100.0 << " μs\n";
    std::cout << "Overhead per call: " << (multiple_time / 100.0) - (single_time / 1.0) << " μs\n\n";
}

void compare_algorithmic_approaches() {
    std::cout << "🏎️ Algorithmic Approach Comparison\n";
    std::cout << "==================================\n";
    
    Timer timer;
    
    // Test polygon - building footprint
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    // MeshCut approach (current implementation)
    timer.reset();
    auto meshcut_result = meshcut::tessellate(building, bbox, 8);
    double meshcut_time = timer.elapsed_us();
    
    // Pure Earcut approach
    timer.reset();
    std::vector<std::vector<std::array<double, 2>>> earcut_poly;
    earcut_poly.push_back({});
    for (auto& pt : building) {
        earcut_poly[0].push_back({pt.x, pt.y});
    }
    auto earcut_result = mapbox::earcut<uint32_t>(earcut_poly);
    double earcut_time = timer.elapsed_us();
    
    std::cout << "Building footprint (8 vertices):\n";
    std::cout << "  MeshCut (8x8 grid): " << meshcut_time << " μs → " << meshcut_result.size() << " triangles\n";
    std::cout << "  Pure Earcut:        " << earcut_time << " μs → " << earcut_result.size()/3 << " triangles\n";
    std::cout << "  Time ratio:          " << meshcut_time / earcut_time << "x slower\n";
    std::cout << "  Triangle ratio:      " << (double)meshcut_result.size() / (earcut_result.size()/3) << "x more triangles\n\n";
    
    // Break down where MeshCut time goes
    timer.reset();
    auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, 8, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    double grid_time = timer.elapsed_us();
    
    timer.reset();
    auto isInterior = meshcut::internal::classifyTriangles(gridTriangles, building);
    double classify_time = timer.elapsed_us();
    
    timer.reset();
    auto isPartial = meshcut::internal::findPartiallyIntersectedTriangles(gridTriangles, building, isInterior);
    double partial_time = timer.elapsed_us();
    
    int partial_count = std::count(isPartial.begin(), isPartial.end(), true);
    int interior_count = std::count(isInterior.begin(), isInterior.end(), true);
    
    std::cout << "MeshCut breakdown:\n";
    std::cout << "  Grid generation:     " << grid_time << " μs\n";
    std::cout << "  Triangle classify:   " << classify_time << " μs\n";
    std::cout << "  Find partials:       " << partial_time << " μs\n";
    std::cout << "  Remaining time:      " << meshcut_time - grid_time - classify_time - partial_time << " μs (intersection + earcut)\n";
    std::cout << "  Interior triangles:  " << interior_count << "\n";
    std::cout << "  Partial triangles:   " << partial_count << " (each needs Clipper + Earcut)\n";
}

int main() {
    std::cout << "🧪 MeshCut Performance Deep Dive\n";
    std::cout << "=================================\n\n";
    
    test_earcut_overhead();
    test_clipper_overhead(); 
    compare_algorithmic_approaches();
    
    std::cout << "💡 Analysis:\n";
    std::cout << "The performance issue likely comes from:\n";
    std::cout << "1. Function call overhead (setup/teardown)\n";
    std::cout << "2. Memory allocation overhead per call\n";
    std::cout << "3. Clipper2 intersection computation costs\n";
    std::cout << "4. Many small operations vs few large ones\n\n";
    
    std::cout << "🎯 This is a classic 'death by a thousand cuts' scenario\n";
    std::cout << "The algorithm is correct, but the implementation has overhead\n";
    
    return 0;
}