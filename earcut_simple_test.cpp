#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include "../third_party/earcut.hpp"

int main() {
    // Simple test: outer square with inner square hole
    std::vector<std::vector<std::array<double, 2>>> polygon = {
        // Outer ring (counter-clockwise)
        {{
            {{0, 0}}, {{100, 0}}, {{100, 100}}, {{0, 100}}
        }},
        // Inner hole (clockwise)  
        {{
            {{20, 20}}, {{20, 80}}, {{80, 80}}, {{80, 20}}
        }}
    };
    
    auto result = mapbox::earcut<uint32_t>(polygon);
    
    std::cout << "Simple hole test:\n";
    std::cout << "Input: outer ring (4 pts) + hole (4 pts) = 8 total vertices\n";
    std::cout << "Output: " << result.size()/3 << " triangles\n";
    std::cout << "Triangles: ";
    for (size_t i = 0; i < result.size(); i += 3) {
        std::cout << "(" << result[i] << "," << result[i+1] << "," << result[i+2] << ") ";
    }
    std::cout << "\n\n";
    
    // Generate SVG to visualize
    std::ofstream svg("earcut_test.svg");
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg width=\"200\" height=\"200\" xmlns=\"http://www.w3.org/2000/svg\">\n";
    svg << "  <rect width=\"200\" height=\"200\" fill=\"white\" stroke=\"black\"/>\n";
    
    // Flatten coordinates for visualization
    std::vector<double> coords;
    for (const auto& ring : polygon) {
        for (const auto& pt : ring) {
            coords.push_back(pt[0] + 50); // offset for viewing
            coords.push_back(pt[1] + 50);
        }
    }
    
    // Draw triangulation
    svg << "  <g fill=\"none\" stroke=\"red\" stroke-width=\"1\">\n";
    for (size_t i = 0; i < result.size(); i += 3) {
        uint32_t i0 = result[i];
        uint32_t i1 = result[i + 1]; 
        uint32_t i2 = result[i + 2];
        
        if (i0*2+1 < coords.size() && i1*2+1 < coords.size() && i2*2+1 < coords.size()) {
            double x0 = coords[i0*2], y0 = coords[i0*2+1];
            double x1 = coords[i1*2], y1 = coords[i1*2+1];
            double x2 = coords[i2*2], y2 = coords[i2*2+1];
            
            svg << "    <polygon points=\"" << x0 << "," << y0 << " " 
                << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
        }
    }
    svg << "  </g>\n";
    
    // Draw outline
    svg << "  <polygon points=\"50,50 150,50 150,150 50,150\" fill=\"none\" stroke=\"blue\" stroke-width=\"2\"/>\n";
    svg << "  <polygon points=\"70,70 70,130 130,130 130,70\" fill=\"none\" stroke=\"green\" stroke-width=\"2\"/>\n";
    
    svg << "</svg>\n";
    svg.close();
    
    std::cout << "Visualization saved to earcut_test.svg\n";
    return 0;
}