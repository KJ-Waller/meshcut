#include <iostream>
#include <vector>
#include <fstream>
#include "../test/test_comprehensive_tile_suite.cpp"

int main() {
    TileReader reader;
    std::cout << "Downloading Innsbruck tile data...\n";
    
    bool success = reader.downloadTile(14, 8712, 5741);
    if (!success) {
        std::cerr << "Failed to download tile\n";
        return 1;
    }
    
    auto polygons = reader.extractTestPolygons(250);
    
    // Find polygon 202 (index 201 since 0-based)
    if (polygons.size() > 201) {
        const auto& poly202 = polygons[201]; // 0-based index
        
        std::cout << "\n=== POLYGON 202 DETAILED ANALYSIS ===\n";
        std::cout << "Name: " << poly202.name << "\n";
        std::cout << "ID: " << poly202.id << "\n";
        std::cout << "Rings: " << poly202.rings.size() << "\n\n";
        
        for (size_t i = 0; i < poly202.rings.size(); i++) {
            const auto& ring = poly202.rings[i];
            std::cout << "Ring " << i << " (" << (i == 0 ? "outer" : "hole") << "):\n";
            std::cout << "  Vertices: " << ring.size() << "\n";
            
            if (!ring.empty()) {
                std::cout << "  First point: (" << ring[0].x << ", " << ring[0].y << ")\n";
                std::cout << "  Last point: (" << ring.back().x << ", " << ring.back().y << ")\n";
                std::cout << "  Closed: " << (ring.size() > 2 && 
                    std::abs(ring[0].x - ring.back().x) < 0.01 && 
                    std::abs(ring[0].y - ring.back().y) < 0.01 ? "YES" : "NO") << "\n";
                
                // Print first few points to see pattern
                std::cout << "  First 5 points: ";
                for (size_t j = 0; j < std::min(5UL, ring.size()); j++) {
                    std::cout << "(" << ring[j].x << ", " << ring[j].y << ")";
                    if (j < std::min(5UL, ring.size()) - 1) std::cout << " -> ";
                }
                std::cout << "\n";
                
                if (ring.size() > 5) {
                    std::cout << "  Last 5 points: ";
                    for (size_t j = std::max(0UL, ring.size() - 5); j < ring.size(); j++) {
                        std::cout << "(" << ring[j].x << ", " << ring[j].y << ")";
                        if (j < ring.size() - 1) std::cout << " -> ";
                    }
                    std::cout << "\n";
                }
            }
            std::cout << "\n";
        }
        
        // Test Earcut format
        auto earcutFormat = poly202.toEarcutFormat();
        std::cout << "Earcut format rings: " << earcutFormat.size() << "\n";
        for (size_t i = 0; i < earcutFormat.size(); i++) {
            std::cout << "  Earcut ring " << i << ": " << earcutFormat[i].size() << " vertices\n";
        }
    }
    
    return 0;
}