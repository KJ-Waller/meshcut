#include <iostream>
#include <vector>
#include <iomanip>
#include "meshcut/meshcut.hpp"

int main() {
    std::cout << "🔧 Testing Complete Intersection + Earcut Pipeline\n";
    std::cout << "=================================================\n";
    
    // Use the same test case
    meshcut::Triangle triangle = {
        {0.05, 0.05},   // inside
        {0.15, 0.05},   // outside  
        {0.10, 0.15}    // outside
    };
    
    meshcut::Polygon square = {
        {0.1, 0.1},     // bottom-left
        {0.2, 0.1},     // bottom-right
        {0.2, 0.2},     // top-right
        {0.1, 0.2}      // top-left
    };
    
    std::cout << "Testing intersection + Earcut pipeline:\n";
    
    auto intersection = meshcut::internal::fastTrianglePolygonIntersection(triangle, square);
    
    std::cout << "Intersection (" << intersection.size() << " points): ";
    for (const auto& p : intersection) {
        std::cout << "[" << std::fixed << std::setprecision(3) << p.x << "," << p.y << "] ";
    }
    std::cout << "\n";
    
    if (intersection.size() >= 3) {
        std::cout << "✅ Size check passed, attempting Earcut triangulation...\n";
        
        try {
            auto triangulated = meshcut::internal::earcutTriangulate(intersection);
            std::cout << "✅ Earcut succeeded! Generated " << triangulated.size() << " triangles\n";
            
            for (size_t i = 0; i < triangulated.size(); ++i) {
                const auto& tri = triangulated[i];
                std::cout << "  Triangle " << i << ": [" << tri.a.x << "," << tri.a.y << "] [" 
                          << tri.b.x << "," << tri.b.y << "] [" << tri.c.x << "," << tri.c.y << "]\n";
            }
        } catch (const std::exception& e) {
            std::cout << "❌ Earcut failed: " << e.what() << "\n";
        }
    } else {
        std::cout << "❌ Size check failed (" << intersection.size() << " < 3)\n";
    }
    
    // Test with a building footprint case that was failing
    std::cout << "\n🏢 Testing with building footprint edge case:\n";
    
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, 8, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    auto isInterior = meshcut::internal::classifyTriangles(gridTriangles, building);
    auto isPartiallyIntersected = meshcut::internal::findPartiallyIntersectedTriangles(gridTriangles, building, isInterior);
    
    // Find first partial triangle and test it
    for (size_t i = 0; i < gridTriangles.size() && i < 10; ++i) {
        if (isPartiallyIntersected[i]) {
            std::cout << "\nTesting partial triangle " << i << ":\n";
            const auto& tri = gridTriangles[i];
            std::cout << "Triangle: [" << tri.a.x << "," << tri.a.y << "] [" 
                      << tri.b.x << "," << tri.b.y << "] [" << tri.c.x << "," << tri.c.y << "]\n";
            
            auto result = meshcut::internal::fastTrianglePolygonIntersection(tri, building);
            std::cout << "Intersection (" << result.size() << " points): ";
            for (const auto& p : result) {
                std::cout << "[" << p.x << "," << p.y << "] ";
            }
            std::cout << "\n";
            
            if (result.size() >= 3) {
                auto earcut_result = meshcut::internal::earcutTriangulate(result);
                std::cout << "✅ Earcut: " << earcut_result.size() << " triangles\n";
                break;
            } else {
                std::cout << "❌ Not enough points for triangulation\n";
            }
        }
    }
    
    return 0;
}