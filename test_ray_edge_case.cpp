// Test ray casting edge case for vertex at (6,4)
#include <iostream>
#include <vector>
#include <cmath>

struct Point {
    double x, y;
    Point(double x, double y) : x(x), y(y) {}
};

class Edge {
public:
    double x1, y1, x2, y2;
    double minY, maxY;
    
    Edge(double x1, double y1, double x2, double y2) 
        : x1(x1), y1(y1), x2(x2), y2(y2), minY(std::min(y1, y2)), maxY(std::max(y1, y2)) {}
    
    double getXAtY(double y) const {
        if (std::abs(y2 - y1) < 1e-10) return x1; // Horizontal edge
        return x1 + (x2 - x1) * (y - y1) / (y2 - y1);
    }
    
    bool intersectsRay(double x, double y) const {
        // Quick Y bounds check
        if (y < minY || y > maxY) return false;
        
        // Handle edge cases - endpoint exactly on ray
        if (std::abs(y - y1) < 1e-10) {
            // Only count if other endpoint is above ray
            return y2 > y;
        }
        if (std::abs(y - y2) < 1e-10) {
            // Only count if other endpoint is above ray  
            return y1 > y;
        }
        
        // Compute intersection X
        double intersectX = getXAtY(y);
        return intersectX > x;
    }
};

bool isInside(const std::vector<Edge>& edges, double x, double y) {
    int intersections = 0;
    
    std::cout << "Testing point (" << x << "," << y << "):\n";
    for (size_t i = 0; i < edges.size(); i++) {
        const auto& edge = edges[i];
        if (edge.intersectsRay(x, y)) {
            intersections++;
            std::cout << "  Edge " << i << " [(" << edge.x1 << "," << edge.y1 << ") -> (" 
                     << edge.x2 << "," << edge.y2 << ")] intersects ray\n";
        }
    }
    
    std::cout << "  Total intersections: " << intersections << " -> " 
              << (intersections % 2 == 1 ? "INSIDE" : "OUTSIDE") << "\n\n";
    return (intersections % 2) == 1;
}

int main() {
    // Create a simple square polygon with vertex at (6,4)
    // Square from (4,2) to (8,6) - vertex at (6,4) is on the bottom edge
    std::vector<Edge> edges = {
        Edge(4, 2, 8, 2),  // bottom edge (6,4) is NOT a vertex here
        Edge(8, 2, 8, 6),  // right edge
        Edge(8, 6, 4, 6),  // top edge
        Edge(4, 6, 4, 2)   // left edge
    };
    
    std::cout << "=== Testing square polygon (4,2) to (8,6) ===\n";
    
    // Test the problematic point (6,4)
    bool result1 = isInside(edges, 6, 4);
    
    // Now test with a polygon that actually has vertex at (6,4)
    std::vector<Edge> edges2 = {
        Edge(4, 2, 6, 4),  // edge to vertex (6,4)
        Edge(6, 4, 8, 2),  // edge from vertex (6,4) 
        Edge(8, 2, 8, 6),  // right edge
        Edge(8, 6, 4, 6),  // top edge
        Edge(4, 6, 4, 2)   // left edge
    };
    
    std::cout << "=== Testing polygon WITH vertex at (6,4) ===\n";
    bool result2 = isInside(edges2, 6, 4);
    
    return 0;
}