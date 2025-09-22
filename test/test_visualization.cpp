#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>

// Helper function to create a simple but interesting polygon
std::vector<double> createTestPolygon() {
    // Create a polygon that will intersect grid boundaries nicely
    return {
        1.2, 1.1,   // Bottom-left-ish
        2.8, 1.0,   // Bottom-right-ish  
        3.0, 2.0,   // Right side
        2.5, 2.8,   // Top-right
        1.8, 3.1,   // Top
        1.0, 2.5,   // Top-left
        0.8, 1.8    // Left side
    };
}

// SVG writer class for visualization
class SVGWriter {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double width, height, scale;
    
public:
    SVGWriter(const std::string& filename, double minX_, double minY_, double maxX_, double maxY_) 
        : minX(minX_), minY(minY_), maxX(maxX_), maxY(maxY_) {
        
        width = maxX - minX;
        height = maxY - minY;
        scale = 400.0 / std::max(width, height); // 400px canvas
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"" << (width * scale + 100) << "\" ";
        file << "height=\"" << (height * scale + 100) << "\" ";
        file << "viewBox=\"0 0 " << (width * scale + 100) << " " << (height * scale + 100) << "\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        file << "<g transform=\"translate(50,50)\">\n"; // 50px margin
    }
    
    ~SVGWriter() {
        file << "</g>\n</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y) {
        double tx = (x - minX) * scale;
        double ty = (maxY - y) * scale; // Flip Y axis
        return {tx, ty};
    }
    
    void drawGrid(const meshcut::MeshCutOptions& options) {
        file << "<g stroke=\"#e0e0e0\" stroke-width=\"0.5\" fill=\"none\">\n";
        
        // Vertical lines
        for (int i = 0; i <= options.gridWidth; i++) {
            double x = options.gridOriginX + i * options.cellSize;
            auto [x1, y1] = transform(x, options.gridOriginY);
            auto [x2, y2] = transform(x, options.gridOriginY + options.gridHeight * options.cellSize);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        // Horizontal lines
        for (int j = 0; j <= options.gridHeight; j++) {
            double y = options.gridOriginY + j * options.cellSize;
            auto [x1, y1] = transform(options.gridOriginX, y);
            auto [x2, y2] = transform(options.gridOriginX + options.gridWidth * options.cellSize, y);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        file << "</g>\n";
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color = "#4CAF50", double opacity = 0.3) {
        if (polygon.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i + 1]);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"" << opacity << "\" ";
        file << "stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
    }
    
    void drawTriangles(const std::vector<double>& vertices, const std::vector<uint32_t>& indices, 
                       const std::string& color = "#2196F3", double opacity = 0.1) {
        file << "<g fill=\"" << color << "\" fill-opacity=\"" << opacity << "\" stroke=\"" << color << "\" stroke-width=\"0.5\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1]; 
            uint32_t i2 = indices[i + 2];
            
            if (i0*2+1 < vertices.size() && i1*2+1 < vertices.size() && i2*2+1 < vertices.size()) {
                auto [x0, y0] = transform(vertices[i0*2], vertices[i0*2+1]);
                auto [x1, y1] = transform(vertices[i1*2], vertices[i1*2+1]);
                auto [x2, y2] = transform(vertices[i2*2], vertices[i2*2+1]);
                
                file << "<polygon points=\"" << x0 << "," << y0 << " " 
                     << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
            }
        }
        
        file << "</g>\n";
    }
    
    void drawVertices(const std::vector<double>& vertices, const std::string& color = "#F44336") {
        file << "<g fill=\"" << color << "\" stroke=\"white\" stroke-width=\"1\">\n";
        
        for (size_t i = 0; i < vertices.size(); i += 2) {
            auto [x, y] = transform(vertices[i], vertices[i + 1]);
            file << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"2\"/>\n";
        }
        
        file << "</g>\n";
    }
    
    void addTitle(const std::string& title) {
        file << "<text x=\"" << (width * scale / 2) << "\" y=\"-20\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"16\" font-weight=\"bold\">" << title << "</text>\n";
    }
    
    void addStats(int triangles, int vertices, int cells) {
        file << "<text x=\"0\" y=\"" << (height * scale + 30) << "\" ";
        file << "font-family=\"Arial\" font-size=\"12\" fill=\"#666\">";
        file << "Triangles: " << triangles << " | Vertices: " << vertices << " | Cells processed: " << cells;
        file << "</text>\n";
    }
};

int main() {
    std::cout << "MeshCut Visualization Test\n";
    std::cout << "==========================\n";
    
    // Create test polygon
    auto polygon = createTestPolygon();
    std::cout << "Test polygon: " << polygon.size()/2 << " vertices\n";
    
    // MeshCut options for good visualization
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0;
    options.cellSize = 0.5;
    options.gridWidth = 8;
    options.gridHeight = 8;
    options.diagonalNE = true;
    
    std::cout << "Grid: " << options.gridWidth << "x" << options.gridHeight 
              << " cells of size " << options.cellSize << "\n";
    
    // Run MeshCut
    auto meshcutResult = meshcut::meshcut_full(polygon, {}, options);
    
    std::cout << "MeshCut generated:\n";
    std::cout << "- " << meshcutResult.indices.size()/3 << " triangles\n";
    std::cout << "- " << meshcutResult.vertices.size()/2 << " vertices\n";
    
    // Also run regular earcut for comparison
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({});
    for (size_t i = 0; i < polygon.size(); i += 2) {
        earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
    }
    auto earcutIndices = mapbox::earcut<uint32_t>(earcutPolygon);
    
    std::cout << "Earcut generated: " << earcutIndices.size()/3 << " triangles\n";
    
    // Calculate bounds for SVG
    double minX = 0.0, maxX = 4.0;
    double minY = 0.0, maxY = 4.0;
    
    // Create SVG visualization 
    {
        SVGWriter svg("meshcut_visualization.svg", minX, minY, maxX, maxY);
        svg.addTitle("MeshCut Grid-Constrained Triangulation");
        
        // Draw grid first (background)
        svg.drawGrid(options);
        
        // Draw original polygon
        svg.drawPolygon(polygon, "#4CAF50", 0.2);
        
        // Draw triangles
        svg.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", 0.15);
        
        // Draw vertices
        svg.drawVertices(meshcutResult.vertices, "#F44336");
        
        // Add statistics
        svg.addStats(meshcutResult.indices.size()/3, meshcutResult.vertices.size()/2, 
                    options.gridWidth * options.gridHeight);
    }
    
    // Create comparison SVG with earcut
    {
        SVGWriter svg("earcut_comparison.svg", minX, minY, maxX, maxY);
        svg.addTitle("Earcut Standard Triangulation");
        
        // Draw original polygon
        svg.drawPolygon(polygon, "#4CAF50", 0.2);
        
        // Draw earcut triangles (using original polygon vertices)
        svg.drawTriangles(polygon, earcutIndices, "#FF9800", 0.15);
        
        // Draw original vertices
        std::vector<double> originalVertices = polygon;
        svg.drawVertices(originalVertices, "#F44336");
        
        svg.addStats(earcutIndices.size()/3, polygon.size()/2, 0);
    }
    
    std::cout << "\nVisualization files created:\n";
    std::cout << "- meshcut_visualization.svg (grid-constrained)\n";
    std::cout << "- earcut_comparison.svg (standard triangulation)\n";
    std::cout << "\nOpen these SVG files in a web browser to see the results!\n";
    
    return 0;
}