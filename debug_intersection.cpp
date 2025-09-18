#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "meshcut/meshcut.hpp"

void debug_intersection_issue() {
    std::cout << "🔍 Debugging Intersection Detection Issue\n";
    std::cout << "=========================================\n";
    
    // Building footprint that should definitely have partial intersections
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, 8, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    auto isInterior = meshcut::internal::classifyTriangles(gridTriangles, building);
    auto isPartiallyIntersected = meshcut::internal::findPartiallyIntersectedTriangles(gridTriangles, building, isInterior);
    
    int interior_count = std::count(isInterior.begin(), isInterior.end(), true);
    int partial_count = std::count(isPartiallyIntersected.begin(), isPartiallyIntersected.end(), true);
    
    std::cout << "8x8 Grid Results:\n";
    std::cout << "Total triangles: " << gridTriangles.size() << "\n";
    std::cout << "Interior triangles: " << interior_count << "\n";
    std::cout << "Partial triangles: " << partial_count << "\n";
    
    // Test some specific partial triangles
    int successful_intersections = 0;
    int failed_intersections = 0;
    
    for (size_t i = 0; i < gridTriangles.size(); ++i) {
        if (isPartiallyIntersected[i]) {
            auto intersection = meshcut::internal::computeTrianglePolygonIntersection(gridTriangles[i], building);
            
            if (intersection.size() >= 3) {
                successful_intersections++;
                if (successful_intersections <= 3) {  // Show first 3
                    std::cout << "\n✅ Triangle " << i << " intersection (size: " << intersection.size() << "):\n";
                    std::cout << "Triangle: [" << gridTriangles[i].a.x << "," << gridTriangles[i].a.y << "] [" 
                              << gridTriangles[i].b.x << "," << gridTriangles[i].b.y << "] ["
                              << gridTriangles[i].c.x << "," << gridTriangles[i].c.y << "]\n";
                    std::cout << "Intersection: ";
                    for (const auto& p : intersection) {
                        std::cout << "[" << p.x << "," << p.y << "] ";
                    }
                    std::cout << "\n";
                }
            } else {
                failed_intersections++;
                if (failed_intersections <= 3) {  // Show first 3 failures
                    std::cout << "\n❌ Triangle " << i << " failed intersection (size: " << intersection.size() << "):\n";
                    std::cout << "Triangle: [" << gridTriangles[i].a.x << "," << gridTriangles[i].a.y << "] [" 
                              << gridTriangles[i].b.x << "," << gridTriangles[i].b.y << "] ["
                              << gridTriangles[i].c.x << "," << gridTriangles[i].c.y << "]\n";
                    std::cout << "Result: ";
                    for (const auto& p : intersection) {
                        std::cout << "[" << p.x << "," << p.y << "] ";
                    }
                    std::cout << "\n";
                }
            }
        }
    }
    
    std::cout << "\n📊 Intersection Results:\n";
    std::cout << "Successful intersections: " << successful_intersections << "\n";
    std::cout << "Failed intersections: " << failed_intersections << "\n";
    std::cout << "Success rate: " << (successful_intersections * 100.0) / (successful_intersections + failed_intersections) << "%\n";
    
    // Test full tessellation to see expected result
    auto complete_result = meshcut::tessellate(building, bbox, 8);
    std::cout << "\nComplete tessellation: " << complete_result.size() << " triangles\n";
    std::cout << "Expected: ~" << interior_count + (successful_intersections * 2) << " triangles\n";
}

void optimize_classification_speed() {
    std::cout << "\n⚡ Classification Speed Optimization Analysis\n";
    std::cout << "=============================================\n";
    
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    std::vector<int> grid_sizes = {4, 8, 16, 32};
    
    for (int gridN : grid_sizes) {
        auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, gridN, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
        
        auto start = std::chrono::high_resolution_clock::now();
        auto isInterior = meshcut::internal::classifyTriangles(gridTriangles, building);
        auto end = std::chrono::high_resolution_clock::now();
        
        double classify_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() / 1000.0;
        int interior_count = std::count(isInterior.begin(), isInterior.end(), true);
        
        std::cout << gridN << "x" << gridN << " grid: " << classify_time << " μs, " 
                  << classify_time / gridTriangles.size() << " μs/triangle, "
                  << interior_count << " interior\n";
    }
    
    std::cout << "\n💡 Classification scales as O(n * m) where n=grid_triangles, m=polygon_vertices\n";
    std::cout << "For 32x32 grid (2048 triangles) × 8 vertices = 16,384 point-in-polygon tests\n";
    std::cout << "At 0.34 μs/triangle × 2048 triangles = ~696 μs just for classification\n";
}

int main() {
    debug_intersection_issue();
    optimize_classification_speed();
    
    std::cout << "\n🚀 Key Findings & Next Steps:\n";
    std::cout << "1. Check if intersection algorithm is working correctly\n";
    std::cout << "2. Classification is O(n*m) - needs spatial optimization\n";
    std::cout << "3. Consider grid-aware point-in-polygon algorithms\n";
    std::cout << "4. Batch processing or early termination strategies\n";
    
    return 0;
}