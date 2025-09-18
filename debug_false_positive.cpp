#include <iostream>
#include <vector>
#include <iomanip>
#include "meshcut/meshcut.hpp"

void test_specific_triangle() {
    std::cout << "🔍 Testing Specific Triangle That Should NOT Intersect\n";
    std::cout << "=====================================================\n";
    
    // This is one of the failing triangles - should NOT intersect building
    meshcut::Triangle triangle = {
        {0.125, 0.000},  // bottom-right of grid cell
        {0.125, 0.125},  // top-right of grid cell  
        {0.000, 0.125}   // top-left of grid cell
    };
    
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    
    std::cout << "Triangle: [" << triangle.a.x << "," << triangle.a.y << "] [" 
              << triangle.b.x << "," << triangle.b.y << "] [" << triangle.c.x << "," << triangle.c.y << "]\n";
    
    // Test individual conditions
    bool hasVertexInside = meshcut::internal::pointInPolygon(triangle.a, building) ||
                          meshcut::internal::pointInPolygon(triangle.b, building) ||
                          meshcut::internal::pointInPolygon(triangle.c, building);
    
    bool intersectsEdges = meshcut::internal::triangleIntersectsPolygonEdge(triangle, building);
    
    bool hasPolygonVertexInside = false;
    for (const auto& polyPoint : building) {
        if (meshcut::internal::pointInTriangle(polyPoint, triangle)) {
            hasPolygonVertexInside = true;
            std::cout << "Found polygon vertex inside triangle: [" << polyPoint.x << "," << polyPoint.y << "]\n";
            break;
        }
    }
    
    std::cout << "hasVertexInside: " << hasVertexInside << "\n";
    std::cout << "intersectsEdges: " << intersectsEdges << "\n"; 
    std::cout << "hasPolygonVertexInside: " << hasPolygonVertexInside << "\n";
    std::cout << "Overall partial: " << (hasVertexInside || intersectsEdges || hasPolygonVertexInside) << "\n";
    
    if (intersectsEdges) {
        std::cout << "\n🔍 Debugging edge intersections (this should be FALSE):\n";
        
        std::vector<std::pair<meshcut::Point, meshcut::Point>> triangleEdges = {
            {triangle.a, triangle.b},
            {triangle.b, triangle.c},
            {triangle.c, triangle.a}
        };
        
        for (size_t i = 0; i < building.size(); ++i) {
            meshcut::Point polyStart = building[i];
            meshcut::Point polyEnd = building[(i + 1) % building.size()];
            
            std::cout << "  Building edge " << i << ": [" << polyStart.x << "," << polyStart.y << "] -> [" << polyEnd.x << "," << polyEnd.y << "]\n";
            
            for (size_t j = 0; j < triangleEdges.size(); ++j) {
                bool intersects = meshcut::internal::lineSegmentsIntersect(
                    triangleEdges[j].first, triangleEdges[j].second, polyStart, polyEnd);
                
                if (intersects) {
                    std::cout << "    ❌ FALSE POSITIVE: Triangle edge " << j << " intersects!\n";
                    std::cout << "       Triangle edge: [" << triangleEdges[j].first.x << "," << triangleEdges[j].first.y 
                              << "] -> [" << triangleEdges[j].second.x << "," << triangleEdges[j].second.y << "]\n";
                }
            }
        }
    }
}

int main() {
    test_specific_triangle();
    return 0;
}