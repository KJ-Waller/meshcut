#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "meshcut/meshcut.hpp"

class PreciseProfiler {
    std::chrono::high_resolution_clock::time_point start_time;
    
public:
    void start() {
        start_time = std::chrono::high_resolution_clock::now();
    }
    
    double elapsed_us() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_time).count() / 1000.0;
    }
};

void profile_detailed_components() {
    std::cout << "🔬 Detailed MeshCut Component Profiling\n";
    std::cout << "======================================\n";
    
    // Realistic building footprint
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    // Test with 16x16 grid (your minimum requirement)
    int gridN = 16;
    
    PreciseProfiler profiler;
    
    std::cout << "Testing with 16x16 grid on building footprint:\n";
    std::cout << "Polygon vertices: " << building.size() << "\n\n";
    
    // 1. Grid Generation
    profiler.start();
    auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, gridN, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    double grid_time = profiler.elapsed_us();
    
    // 2. Classification
    profiler.start();
    auto isInterior = meshcut::internal::classifyTriangles(gridTriangles, building);
    double classify_time = profiler.elapsed_us();
    
    int interior_count = std::count(isInterior.begin(), isInterior.end(), true);
    
    // 3. Find partials
    profiler.start();
    auto isPartiallyIntersected = meshcut::internal::findPartiallyIntersectedTriangles(gridTriangles, building, isInterior);
    double partial_find_time = profiler.elapsed_us();
    
    int partial_count = std::count(isPartiallyIntersected.begin(), isPartiallyIntersected.end(), true);
    
    // 4. Process partials (the suspected bottleneck)
    profiler.start();
    std::vector<meshcut::Triangle> partialTriangles;
    int successful_intersections = 0;
    int earcut_calls = 0;
    double total_intersection_time = 0;
    double total_earcut_time = 0;
    
    for (size_t i = 0; i < gridTriangles.size(); ++i) {
        if (isPartiallyIntersected[i]) {
            // Time the intersection computation
            PreciseProfiler intersection_profiler;
            intersection_profiler.start();
            auto intersection = meshcut::internal::computeTrianglePolygonIntersection(gridTriangles[i], building);
            double intersection_time = intersection_profiler.elapsed_us();
            total_intersection_time += intersection_time;
            
            if (intersection.size() >= 3) {
                successful_intersections++;
                
                // Time the Earcut call
                PreciseProfiler earcut_profiler;
                earcut_profiler.start();
                auto triangulated = meshcut::internal::earcutTriangulate(intersection);
                double earcut_time = earcut_profiler.elapsed_us();
                total_earcut_time += earcut_time;
                earcut_calls++;
                
                partialTriangles.insert(partialTriangles.end(), triangulated.begin(), triangulated.end());
            }
        }
    }
    
    double partial_process_time = profiler.elapsed_us();
    
    // 5. Complete call for comparison
    profiler.start();
    auto complete_result = meshcut::tessellate(building, bbox, gridN);
    double complete_time = profiler.elapsed_us();
    
    // Results
    std::cout << "📊 Detailed Timing Breakdown:\n";
    std::cout << "────────────────────────────────\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Grid generation:        " << std::setw(8) << grid_time << " μs (" << std::setw(4) << (grid_time/complete_time)*100 << "%)\n";
    std::cout << "Triangle classification:" << std::setw(8) << classify_time << " μs (" << std::setw(4) << (classify_time/complete_time)*100 << "%)\n";
    std::cout << "Find partial triangles: " << std::setw(8) << partial_find_time << " μs (" << std::setw(4) << (partial_find_time/complete_time)*100 << "%)\n";
    std::cout << "Process partial triangles:" << std::setw(6) << partial_process_time << " μs (" << std::setw(4) << (partial_process_time/complete_time)*100 << "%)\n";
    std::cout << "  └─ Total intersections: " << std::setw(8) << total_intersection_time << " μs (" << std::setw(4) << (total_intersection_time/complete_time)*100 << "%)\n";
    std::cout << "  └─ Total Earcut calls:  " << std::setw(8) << total_earcut_time << " μs (" << std::setw(4) << (total_earcut_time/complete_time)*100 << "%)\n";
    std::cout << "TOTAL:                  " << std::setw(8) << complete_time << " μs\n";
    
    std::cout << "\n📈 Algorithm Analysis:\n";
    std::cout << "Grid triangles:         " << gridTriangles.size() << "\n";
    std::cout << "Interior triangles:     " << interior_count << "\n";
    std::cout << "Partial triangles:      " << partial_count << "\n";
    std::cout << "Successful intersections:" << successful_intersections << "\n";
    std::cout << "Earcut calls made:      " << earcut_calls << "\n";
    std::cout << "Final triangles:        " << complete_result.size() << "\n";
    
    std::cout << "\n⚡ Per-Operation Costs:\n";
    if (successful_intersections > 0) {
        std::cout << "Avg intersection time:  " << std::setw(8) << total_intersection_time/successful_intersections << " μs\n";
    }
    if (earcut_calls > 0) {
        std::cout << "Avg Earcut time:        " << std::setw(8) << total_earcut_time/earcut_calls << " μs\n";
    }
    std::cout << "Grid triangle classify: " << std::setw(8) << classify_time/gridTriangles.size() << " μs/triangle\n";
}

void analyze_classification_bottleneck() {
    std::cout << "\n🎯 Classification Performance Analysis\n";
    std::cout << "=====================================\n";
    
    meshcut::Polygon simple_square = {{0.2, 0.2}, {0.8, 0.2}, {0.8, 0.8}, {0.2, 0.8}};
    meshcut::Polygon complex_star;
    
    // Create complex 20-point star
    for (int i = 0; i < 20; ++i) {
        double angle = i * M_PI / 10;
        double radius = (i % 2 == 0) ? 0.4 : 0.2;
        complex_star.push_back(meshcut::Point{0.5 + radius * cos(angle), 0.5 + radius * sin(angle)});
    }
    
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, 16, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    
    PreciseProfiler profiler;
    
    // Test classification with simple vs complex polygons
    profiler.start();
    auto simple_classification = meshcut::internal::classifyTriangles(gridTriangles, simple_square);
    double simple_time = profiler.elapsed_us();
    
    profiler.start();
    auto complex_classification = meshcut::internal::classifyTriangles(gridTriangles, complex_star);
    double complex_time = profiler.elapsed_us();
    
    std::cout << "Classification timing:\n";
    std::cout << "Simple square (4 vertices):  " << simple_time << " μs\n";
    std::cout << "Complex star (20 vertices):  " << complex_time << " μs\n";
    std::cout << "Time scales with polygon complexity: " << complex_time/simple_time << "x\n";
    
    // Analyze what happens inside classification
    int grid_size = gridTriangles.size();
    std::cout << "\nPer-triangle classification cost:\n";
    std::cout << "Simple: " << simple_time/grid_size << " μs/triangle\n";
    std::cout << "Complex: " << complex_time/grid_size << " μs/triangle\n";
}

void test_earcut_batching_potential() {
    std::cout << "\n🚀 Earcut Batching Analysis\n";
    std::cout << "===========================\n";
    
    // Create multiple small polygons (simulating partial triangles)
    std::vector<meshcut::Polygon> small_polygons;
    for (int i = 0; i < 20; ++i) {
        meshcut::Polygon tri;
        double base_x = i * 0.1;
        tri.push_back(meshcut::Point{base_x, 0.0});
        tri.push_back(meshcut::Point{base_x + 0.05, 0.0});
        tri.push_back(meshcut::Point{base_x + 0.025, 0.05});
        small_polygons.push_back(tri);
    }
    
    // Create one large polygon with equivalent area
    meshcut::Polygon large_polygon;
    for (int i = 0; i <= 20; ++i) {
        large_polygon.push_back(meshcut::Point{i * 0.1, 0.0});
    }
    for (int i = 20; i >= 0; --i) {
        large_polygon.push_back(meshcut::Point{i * 0.1, 0.05});
    }
    
    PreciseProfiler profiler;
    
    // Test many small Earcut calls
    profiler.start();
    std::vector<meshcut::Triangle> small_results;
    for (const auto& poly : small_polygons) {
        auto triangles = meshcut::internal::earcutTriangulate(poly);
        small_results.insert(small_results.end(), triangles.begin(), triangles.end());
    }
    double small_time = profiler.elapsed_us();
    
    // Test one large Earcut call
    profiler.start();
    auto large_result = meshcut::internal::earcutTriangulate(large_polygon);
    double large_time = profiler.elapsed_us();
    
    std::cout << "Earcut call comparison:\n";
    std::cout << "20 small calls:  " << small_time << " μs → " << small_results.size() << " triangles\n";
    std::cout << "1 large call:    " << large_time << " μs → " << large_result.size() << " triangles\n";
    std::cout << "Small/Large ratio: " << small_time/large_time << "x (overhead factor)\n";
    
    if (small_time > large_time) {
        std::cout << "✅ Batching could help! Overhead factor: " << small_time/large_time << "x\n";
    } else {
        std::cout << "❌ Batching unlikely to help much\n";
    }
}

int main() {
    std::cout << "🕵️ Advanced MeshCut Bottleneck Investigation\n";
    std::cout << "===========================================\n";
    std::cout << "Target: 16x16 grid performance for MapLibre terrain\n\n";
    
    profile_detailed_components();
    analyze_classification_bottleneck();
    test_earcut_batching_potential();
    
    std::cout << "\n🎯 Optimization Recommendations:\n";
    std::cout << "1. If classification is >30%: Optimize point-in-polygon tests\n";
    std::cout << "2. If intersections are >40%: Further optimize geometric intersection\n";
    std::cout << "3. If Earcut calls are >30%: Consider batching or caching\n";
    std::cout << "4. If overall distribution is even: Need algorithmic changes\n";
    
    return 0;
}