// Test the specific case of ray starting from polygon vertex
#include <iostream>
#include <vector>
#include <cmath>

struct Edge {
    double x1, y1, x2, y2;
    double minY, maxY;
    
    Edge(double x1, double y1, double x2, double y2) 
        : x1(x1), y1(y1), x2(x2), y2(y2), minY(std::min(y1, y2)), maxY(std::max(y1, y2)) {}
    
    double getXAtY(double y) const {
        if (std::abs(y2 - y1) < 1e-10) return x1;
        return x1 + (x2 - x1) * (y - y1) / (y2 - y1);
    }
    
    bool intersectsRay(double x, double y) const {
        if (y < minY || y > maxY) return false;
        
        // Handle edge cases - endpoint exactly on ray
        if (std::abs(y - y1) < 1e-10) {
            return y2 > y;
        }
        if (std::abs(y - y2) < 1e-10) {
            return y1 > y;
        }
        
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
                     << edge.x2 << "," << edge.y2 << ")] intersects\n";
        } else {
            std::cout << "  Edge " << i << " [(" << edge.x1 << "," << edge.y1 << ") -> (" 
                     << edge.x2 << "," << edge.y2 << ")] NO intersection\n";
        }
    }
    
    std::cout << "  Total intersections: " << intersections << " -> " 
              << (intersections % 2 == 1 ? "INSIDE" : "OUTSIDE") << "\n\n";
    return (intersections % 2) == 1;
}

int main() {
    // Create the actual problematic polygon from SVG (scaled down to match grid coords)
    // SVG polygon: "120,290 280,300 300,200 250,120 180,90 100,150 80,220"
    // Grid coordinate (6,4) corresponds to SVG (300,200)
    std::vector<Edge> edges = {
        Edge(2.4, 1.0, 5.6, 1.2),  // 120,290 -> 280,300 (scaled by 1/50)
        Edge(5.6, 1.2, 6.0, 4.0),  // 280,300 -> 300,200 *** VERTEX AT (6,4) ***
        Edge(6.0, 4.0, 5.0, 6.4),  // 300,200 -> 250,120 *** FROM VERTEX AT (6,4) ***
        Edge(5.0, 6.4, 3.6, 7.2),  // 250,120 -> 180,90
        Edge(3.6, 7.2, 2.0, 5.0),  // 180,90 -> 100,150
        Edge(2.0, 5.0, 1.6, 2.4),  // 100,150 -> 80,220
        Edge(1.6, 2.4, 2.4, 1.0)   // 80,220 -> 120,290 (close polygon)
    };
    
    std::cout << "=== Testing ray from polygon vertex (6,4) ===\n";
    
    // Test the exact vertex point
    bool result = isInside(edges, 6.0, 4.0);
    
    // Test nearby points
    std::cout << "\n=== Testing nearby points ===\n";
    isInside(edges, 5.99, 4.0);
    isInside(edges, 6.01, 4.0);
    isInside(edges, 6.0, 3.99);
    isInside(edges, 6.0, 4.01);
    
    return 0;
}