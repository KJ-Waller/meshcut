#include <iostream>
#include <vector>
#include "meshcut/meshcut.hpp"

int main() {
    std::cout << "🔧 Testing Intersection Algorithm Step by Step\n";
    std::cout << "==============================================\n";
    
    // Simple test case: triangle intersecting simple square
    meshcut::Triangle triangle = {
        {0.05, 0.05},   // inside
        {0.15, 0.05},   // outside  
        {0.10, 0.15}    // outside
    };
    
    // Simple square polygon
    meshcut::Polygon square = {
        {0.1, 0.1},     // bottom-left
        {0.2, 0.1},     // bottom-right
        {0.2, 0.2},     // top-right
        {0.1, 0.2}      // top-left
    };
    
    std::cout << "Triangle: [" << triangle.a.x << "," << triangle.a.y << "] [" 
              << triangle.b.x << "," << triangle.b.y << "] [" << triangle.c.x << "," << triangle.c.y << "]\n";
    
    std::cout << "Square: ";
    for (const auto& p : square) {
        std::cout << "[" << p.x << "," << p.y << "] ";
    }
    std::cout << "\n\n";
    
    // Test the intersection
    auto result = meshcut::internal::fastTrianglePolygonIntersection(triangle, square);
    
    std::cout << "Intersection result (" << result.size() << " points): ";
    for (const auto& p : result) {
        std::cout << "[" << p.x << "," << p.y << "] ";
    }
    std::cout << "\n\n";
    
    // Test individual edge checks
    std::cout << "Testing isInsideEdge function:\n";
    for (size_t i = 0; i < square.size(); ++i) {
        size_t j = (i + 1) % square.size();
        meshcut::Point edgeStart = square[i];
        meshcut::Point edgeEnd = square[j];
        
        std::cout << "Edge " << i << ": [" << edgeStart.x << "," << edgeStart.y << "] -> [" << edgeEnd.x << "," << edgeEnd.y << "]\n";
        
        // Test triangle vertices against this edge
        bool a_inside = meshcut::internal::isInsideEdge(triangle.a, edgeStart, edgeEnd);
        bool b_inside = meshcut::internal::isInsideEdge(triangle.b, edgeStart, edgeEnd);
        bool c_inside = meshcut::internal::isInsideEdge(triangle.c, edgeStart, edgeEnd);
        
        std::cout << "  Triangle.a inside: " << a_inside << "\n";
        std::cout << "  Triangle.b inside: " << b_inside << "\n";
        std::cout << "  Triangle.c inside: " << c_inside << "\n";
    }
    
    // Test a point we know should be inside
    meshcut::Point inside_point = {0.15, 0.15};
    meshcut::Point outside_point = {0.05, 0.05};
    
    std::cout << "\nTesting known points:\n";
    std::cout << "Point [0.15,0.15] (should be inside square):\n";
    for (size_t i = 0; i < square.size(); ++i) {
        size_t j = (i + 1) % square.size();
        bool inside = meshcut::internal::isInsideEdge(inside_point, square[i], square[j]);
        std::cout << "  Edge " << i << ": " << inside << "\n";
    }
    
    std::cout << "\nPoint [0.05,0.05] (should be outside square):\n";
    for (size_t i = 0; i < square.size(); ++i) {
        size_t j = (i + 1) % square.size();
        bool inside = meshcut::internal::isInsideEdge(outside_point, square[i], square[j]);
        std::cout << "  Edge " << i << ": " << inside << "\n";
    }
    
    return 0;
}