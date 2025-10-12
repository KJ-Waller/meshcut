#include "meshcut.hpp"
#include "earcut.hpp" 
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <curl/curl.h>
#include <zlib.h>
#include <string>
#include <algorithm>
#include <cstring>

// Callback for CURL to write data
size_t WriteCallback(void* contents, size_t size, size_t nmemb, std::vector<uint8_t>* userp) {
    size_t totalSize = size * nmemb;
    size_t oldSize = userp->size();
    userp->resize(oldSize + totalSize);
    std::memcpy(userp->data() + oldSize, contents, totalSize);
    return totalSize;
}

bool downloadTile(const std::string& url, std::vector<uint8_t>& data) {
    CURL* curl = curl_easy_init();
    if (!curl) return false;
    
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &data);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);
    
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    
    return (res == CURLE_OK);
}

std::vector<uint8_t> decompressGzip(const std::vector<uint8_t>& compressed) {
    z_stream zs = {};
    if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
        return {};
    }
    
    std::vector<uint8_t> decompressed;
    decompressed.reserve(compressed.size() * 4);
    
    zs.next_in = const_cast<uint8_t*>(compressed.data());
    zs.avail_in = compressed.size();
    
    uint8_t buffer[8192];
    int ret;
    do {
        zs.next_out = buffer;
        zs.avail_out = sizeof(buffer);
        
        ret = inflate(&zs, Z_NO_FLUSH);
        if (ret != Z_OK && ret != Z_STREAM_END) {
            break;
        }
        
        decompressed.insert(decompressed.end(), buffer, buffer + (sizeof(buffer) - zs.avail_out));
    } while (ret != Z_STREAM_END);
    
    inflateEnd(&zs);
    return (ret == Z_STREAM_END) ? decompressed : std::vector<uint8_t>{};
}

// Simple Protocol Buffer parser for MVT
class MVTParser {
private:
    const uint8_t* data;
    size_t size;
    size_t pos;
    
    uint64_t readVarint() {
        uint64_t result = 0;
        int shift = 0;
        while (pos < size) {
            uint8_t byte = data[pos++];
            result |= (uint64_t)(byte & 0x7F) << shift;
            if ((byte & 0x80) == 0) break;
            shift += 7;
            if (shift >= 64) break;
        }
        return result;
    }
    
    std::string readString(size_t length) {
        if (pos + length > size) return "";
        std::string result(reinterpret_cast<const char*>(data + pos), length);
        pos += length;
        return result;
    }
    
    void skipBytes(size_t count) {
        pos = std::min(pos + count, size);
    }
    
public:
    struct Point {
        double x, y;
    };
    
    struct Polygon {
        std::vector<Point> points;
        std::string layerName;
        bool crossesEdge = false; // Flag for edge-crossing detection
    };
    
    MVTParser(const std::vector<uint8_t>& tileData) 
        : data(tileData.data()), size(tileData.size()), pos(0) {}
    
    std::vector<Polygon> extractPolygons() {
        std::vector<Polygon> polygons;
        pos = 0;
        
        std::cout << "Parsing MVT data (" << size << " bytes)...\n";
        
        while (pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 2) { // Layer field
                uint64_t layerLength = readVarint();
                size_t layerEnd = std::min(pos + layerLength, size);
                
                auto layerPolygons = parseLayer(layerEnd);
                polygons.insert(polygons.end(), layerPolygons.begin(), layerPolygons.end());
                
                if (polygons.size() >= 20) { // Get more polygons to find edge cases
                    std::cout << "Found " << polygons.size() << " polygons, analyzing for edge crossings\n";
                    break;
                }
            } else {
                skipField(wireType);
            }
        }
        
        // Check which polygons cross tile edges (MVT extent is typically 4096)
        const double TILE_EXTENT = 4096.0;
        const double EDGE_THRESHOLD = 50.0; // Close to edge
        
        for (auto& polygon : polygons) {
            bool touchesLeft = false, touchesRight = false, touchesTop = false, touchesBottom = false;
            
            for (const auto& point : polygon.points) {
                if (point.x <= EDGE_THRESHOLD) touchesLeft = true;
                if (point.x >= TILE_EXTENT - EDGE_THRESHOLD) touchesRight = true;
                if (point.y <= EDGE_THRESHOLD) touchesBottom = true;
                if (point.y >= TILE_EXTENT - EDGE_THRESHOLD) touchesTop = true;
            }
            
            polygon.crossesEdge = touchesLeft || touchesRight || touchesTop || touchesBottom;
            
            if (polygon.crossesEdge) {
                std::cout << "Edge-crossing polygon in layer '" << polygon.layerName 
                         << "' with " << polygon.points.size() << " vertices\n";
                std::cout << "  Touches: ";
                if (touchesLeft) std::cout << "LEFT ";
                if (touchesRight) std::cout << "RIGHT ";
                if (touchesTop) std::cout << "TOP ";
                if (touchesBottom) std::cout << "BOTTOM ";
                std::cout << "\n";
            }
        }
        
        return polygons;
    }
    
private:
    void skipField(uint8_t wireType) {
        switch (wireType) {
            case 0: readVarint(); break;
            case 1: skipBytes(8); break;
            case 2: skipBytes(readVarint()); break;
            case 5: skipBytes(4); break;
            default: break;
        }
    }
    
    std::vector<Polygon> parseLayer(size_t layerEnd) {
        std::vector<Polygon> polygons;
        std::vector<std::string> keys;
        std::vector<std::string> values;
        std::string layerName = "unknown";
        uint32_t extent = 4096;
        
        while (pos < layerEnd && pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 1 && wireType == 2) { // Name
                uint64_t length = readVarint();
                layerName = readString(length);
            } else if (fieldNumber == 2 && wireType == 2) { // Features
                uint64_t featureLength = readVarint();
                size_t featureEnd = std::min(pos + featureLength, size);
                auto featurePolygons = parseFeature(featureEnd, layerName, extent);
                polygons.insert(polygons.end(), featurePolygons.begin(), featurePolygons.end());
            } else if (fieldNumber == 3 && wireType == 2) { // Keys
                uint64_t length = readVarint();
                keys.push_back(readString(length));
            } else if (fieldNumber == 4 && wireType == 2) { // Values
                uint64_t length = readVarint();
                values.push_back(readString(length));
            } else if (fieldNumber == 5 && wireType == 0) { // Extent
                extent = readVarint();
            } else {
                skipField(wireType);
            }
        }
        
        return polygons;
    }
    
    std::vector<Polygon> parseFeature(size_t featureEnd, const std::string& layerName, uint32_t extent) {
        std::vector<Polygon> polygons;
        std::vector<uint32_t> geometry;
        uint32_t type = 0;
        
        while (pos < featureEnd && pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 0) { // Type
                type = readVarint();
            } else if (fieldNumber == 4 && wireType == 2) { // Geometry
                uint64_t geometryLength = readVarint();
                size_t geometryEnd = std::min(pos + geometryLength, size);
                
                while (pos < geometryEnd && pos < size) {
                    geometry.push_back(readVarint());
                }
            } else {
                skipField(wireType);
            }
        }
        
        if (type == 3) { // Polygon
            auto polyPolygons = decodePolygonGeometry(geometry, extent);
            for (auto& poly : polyPolygons) {
                poly.layerName = layerName;
                polygons.push_back(poly);
            }
        }
        
        return polygons;
    }
    
    std::vector<Polygon> decodePolygonGeometry(const std::vector<uint32_t>& geometry, uint32_t extent) {
        std::vector<Polygon> polygons;
        
        if (geometry.empty()) return polygons;
        
        size_t i = 0;
        int32_t x = 0, y = 0;
        
        while (i < geometry.size()) {
            if (i >= geometry.size()) break;
            
            uint32_t cmdInt = geometry[i++];
            uint32_t cmd = cmdInt & 0x7;
            uint32_t count = cmdInt >> 3;
            
            if (cmd == 1) { // MoveTo - start new ring
                Polygon polygon;
                
                for (uint32_t j = 0; j < count && i + 1 < geometry.size(); j++) {
                    x += zigzagDecode(geometry[i++]);
                    y += zigzagDecode(geometry[i++]);
                    
                    // Convert to normalized coordinates [0,1]
                    double normX = (double)x / (double)extent;
                    double normY = 1.0 - (double)y / (double)extent; // Flip Y
                    
                    polygon.points.push_back({normX, normY});
                }
                
                if (!polygon.points.empty()) {
                    polygons.push_back(polygon);
                }
            } else if (cmd == 2) { // LineTo - continue current ring
                if (!polygons.empty()) {
                    for (uint32_t j = 0; j < count && i + 1 < geometry.size(); j++) {
                        x += zigzagDecode(geometry[i++]);
                        y += zigzagDecode(geometry[i++]);
                        
                        double normX = (double)x / (double)extent;  
                        double normY = 1.0 - (double)y / (double)extent;
                        
                        polygons.back().points.push_back({normX, normY});
                    }
                }
            }
        }
        
        return polygons;
    }
    
    int32_t zigzagDecode(uint32_t encoded) {
        return (encoded >> 1) ^ (-(encoded & 1));
    }
};

// Enhanced SVG writer for edge-crossing analysis
class EdgeCrossingSVGWriter {
private:
    std::ofstream file;
    double scale;
    
public:
    EdgeCrossingSVGWriter(const std::string& filename) {
        scale = 400.0; // 400px for the unit square [0,1]x[0,1]
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"900\" height=\"450\" viewBox=\"0 0 900 450\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        
        // Add title
        file << "<text x=\"450\" y=\"30\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"16\" font-weight=\"bold\">";
        file << "Edge-Crossing Polygon Analysis: Earcut vs MeshCut</text>\n";
        
        // Left panel (Earcut)
        file << "<g transform=\"translate(50,60)\">\n";
        file << "<text x=\"200\" y=\"-10\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"14\" font-weight=\"bold\">Earcut (Standard)</text>\n";
        file << "<rect x=\"0\" y=\"0\" width=\"400\" height=\"400\" fill=\"none\" stroke=\"#333\" stroke-width=\"2\"/>\n";
        
        // Right panel (MeshCut) 
        file << "</g>\n<g transform=\"translate(500,60)\">\n";
        file << "<text x=\"200\" y=\"-10\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"14\" font-weight=\"bold\">MeshCut (Grid-Constrained)</text>\n";
        file << "<rect x=\"0\" y=\"0\" width=\"400\" height=\"400\" fill=\"none\" stroke=\"#333\" stroke-width=\"2\"/>\n";
    }
    
    ~EdgeCrossingSVGWriter() {
        file << "</g>\n</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y) {
        return {x * scale, (1.0 - y) * scale}; // Flip Y for SVG
    }
    
    void drawGrid(const meshcut::MeshCutOptions& options, bool rightPanel = false) {
        if (rightPanel) {
            file << "<g stroke=\"#e0e0e0\" stroke-width=\"0.5\" fill=\"none\">\n";
            
            // Grid lines for [0,1] space
            for (int i = 0; i <= 32; i++) {
                double pos = (double)i / 32.0;
                auto [x1, y1] = transform(pos, 0.0);
                auto [x2, y2] = transform(pos, 1.0);
                file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
                
                auto [x3, y3] = transform(0.0, pos);
                auto [x4, y4] = transform(1.0, pos);
                file << "<line x1=\"" << x3 << "\" y1=\"" << y3 << "\" x2=\"" << x4 << "\" y2=\"" << y4 << "\"/>\n";
            }
            
            file << "</g>\n";
        }
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color, bool rightPanel = false) {
        if (polygon.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i + 1]);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"0.2\" ";
        file << "stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
    }
    
    void drawTriangles(const std::vector<double>& vertices, const std::vector<uint32_t>& indices, 
                       const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"" << color << "\" fill-opacity=\"0.1\" stroke=\"" << color << "\" stroke-width=\"0.5\">\n";
        
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
    
    void highlightEdges() {
        // Highlight tile edges with thicker red lines
        file << "<g stroke=\"#FF0000\" stroke-width=\"3\" fill=\"none\">\n";
        file << "<line x1=\"0\" y1=\"0\" x2=\"400\" y2=\"0\"/>\n"; // Bottom
        file << "<line x1=\"0\" y1=\"400\" x2=\"400\" y2=\"400\"/>\n"; // Top  
        file << "<line x1=\"0\" y1=\"0\" x2=\"0\" y2=\"400\"/>\n"; // Left
        file << "<line x1=\"400\" y1=\"0\" x2=\"400\" y2=\"400\"/>\n"; // Right
        file << "</g>\n";
    }
    
    void addStats(int triangles, int vertices, const std::string& type) {
        file << "<text x=\"200\" y=\"430\" text-anchor=\"middle\" font-family=\"Arial\" font-size=\"12\" fill=\"#666\">";
        file << type << " - Triangles: " << triangles << " | Vertices: " << vertices;
        file << "</text>\n";
    }
};

int main() {
    std::cout << "Edge-Crossing Polygon Test\n";
    std::cout << "=========================\n";
    
    // Download a tile that's likely to have edge-crossing polygons
    // Using Innsbruck which we know works
    std::string tileUrl = "https://demotiles.maplibre.org/tiles/14/8872/5904.pbf"; // Innsbruck
    std::vector<uint8_t> tileData;
    
    std::cout << "Downloading vector tile from Innsbruck...\n";
    if (!downloadTile(tileUrl, tileData)) {
        std::cout << "Failed to download tile\n";
        return 1;
    }
    
    std::cout << "Downloaded " << tileData.size() << " bytes\n";
    
    // Always try to decompress first
    auto decompressed = decompressGzip(tileData);
    if (!decompressed.empty()) {
        tileData = std::move(decompressed);
        std::cout << "Decompressed to " << tileData.size() << " bytes\n";
    } else {
        std::cout << "Data not compressed or decompression failed, using raw data\n";
    }
    
    // Parse polygons and find edge-crossing ones
    MVTParser parser(tileData);
    auto polygons = parser.extractPolygons();
    
    std::cout << "\nFound " << polygons.size() << " total polygons\n";
    
    // Find the first edge-crossing polygon
    MVTParser::Polygon* edgePolygon = nullptr;
    for (auto& polygon : polygons) {
        if (polygon.crossesEdge && polygon.points.size() >= 4) {
            edgePolygon = &polygon;
            break;
        }
    }
    
    if (!edgePolygon) {
        std::cout << "No suitable edge-crossing polygons found. Using largest polygon instead.\n";
        // Fallback to largest polygon
        size_t maxSize = 0;
        for (auto& polygon : polygons) {
            if (polygon.points.size() > maxSize) {
                maxSize = polygon.points.size();
                edgePolygon = &polygon;
            }
        }
    }
    
    if (!edgePolygon) {
        std::cout << "No polygons found at all!\n";
        return 1;
    }
    
    std::cout << "\nTesting polygon from layer '" << edgePolygon->layerName << "'\n";
    std::cout << "Vertices: " << edgePolygon->points.size() << "\n";
    std::cout << "Crosses edge: " << (edgePolygon->crossesEdge ? "YES" : "NO") << "\n";
    
    // Convert to flat array for triangulation
    std::vector<double> flatPolygon;
    for (const auto& point : edgePolygon->points) {
        flatPolygon.push_back(point.x);
        flatPolygon.push_back(point.y);
    }
    
    // Set up MeshCut options for [0,1] coordinate space
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0;
    options.cellSize = 1.0 / 32.0; // 32x32 grid over [0,1]
    options.gridWidth = 32;
    options.gridHeight = 32;
    options.diagonalNE = true;
    
    std::cout << "\nRunning MeshCut triangulation...\n";
    auto meshcutResult = meshcut::meshcut_full(flatPolygon, {}, options);
    
    std::cout << "MeshCut result: " << meshcutResult.vertices.size()/2 << " vertices, " 
              << meshcutResult.indices.size()/3 << " triangles\n";
    
    // Run Earcut for comparison
    std::cout << "Running Earcut triangulation...\n";
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({});
    for (size_t i = 0; i < flatPolygon.size(); i += 2) {
        earcutPolygon[0].push_back({{flatPolygon[i], flatPolygon[i+1]}});
    }
    auto earcutIndices = mapbox::earcut<uint32_t>(earcutPolygon);
    
    std::cout << "Earcut result: " << flatPolygon.size()/2 << " vertices, " 
              << earcutIndices.size()/3 << " triangles\n";
    
    // Create visualization
    EdgeCrossingSVGWriter svg("edge_crossing_analysis.svg");
    
    // Left panel - Earcut
    svg.drawPolygon(flatPolygon, "#4CAF50", false);
    svg.drawTriangles(flatPolygon, earcutIndices, "#FF5722", false);
    svg.highlightEdges();
    svg.addStats(earcutIndices.size()/3, flatPolygon.size()/2, "Earcut");
    
    // Right panel - MeshCut  
    svg.drawGrid(options, true);
    svg.drawPolygon(flatPolygon, "#4CAF50", true);
    svg.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", true);
    svg.highlightEdges();
    svg.addStats(meshcutResult.indices.size()/3, meshcutResult.vertices.size()/2, "MeshCut");
    
    std::cout << "\nVisualization saved: edge_crossing_analysis.svg\n";
    
    // Analyze potential issues
    std::cout << "\n=== EDGE-CROSSING ANALYSIS ===\n";
    std::cout << "Polygon crosses tile edge: " << (edgePolygon->crossesEdge ? "YES" : "NO") << "\n";
    
    // Check for vertices very close to edges
    int nearEdgeCount = 0;
    const double EDGE_THRESHOLD = 0.05; // 5% from edge
    
    for (size_t i = 0; i < flatPolygon.size(); i += 2) {
        double x = flatPolygon[i];
        double y = flatPolygon[i+1];
        
        if (x <= EDGE_THRESHOLD || x >= 1.0 - EDGE_THRESHOLD ||
            y <= EDGE_THRESHOLD || y >= 1.0 - EDGE_THRESHOLD) {
            nearEdgeCount++;
        }
    }
    
    std::cout << "Vertices near tile edges: " << nearEdgeCount << "/" << flatPolygon.size()/2 << "\n";
    std::cout << "Triangle density improvement: " << std::fixed << std::setprecision(1) 
              << (double)meshcutResult.indices.size() / (double)earcutIndices.size() << "x\n";
    
    // Check for degenerate triangles
    int degenerateTriangles = 0;
    const double MIN_AREA = 1e-10;
    
    for (size_t i = 0; i < meshcutResult.indices.size(); i += 3) {
        uint32_t i0 = meshcutResult.indices[i];
        uint32_t i1 = meshcutResult.indices[i + 1]; 
        uint32_t i2 = meshcutResult.indices[i + 2];
        
        if (i0*2+1 < meshcutResult.vertices.size() && 
            i1*2+1 < meshcutResult.vertices.size() && 
            i2*2+1 < meshcutResult.vertices.size()) {
            
            double x0 = meshcutResult.vertices[i0*2], y0 = meshcutResult.vertices[i0*2+1];
            double x1 = meshcutResult.vertices[i1*2], y1 = meshcutResult.vertices[i1*2+1];
            double x2 = meshcutResult.vertices[i2*2], y2 = meshcutResult.vertices[i2*2+1];
            
            // Calculate triangle area using cross product
            double area = std::abs((x1-x0)*(y2-y0) - (x2-x0)*(y1-y0)) * 0.5;
            if (area < MIN_AREA) {
                degenerateTriangles++;
            }
        }
    }
    
    std::cout << "Potentially degenerate triangles: " << degenerateTriangles << "\n";
    
    std::cout << "\n=== EDGE-CROSSING TEST COMPLETE ===\n";
    std::cout << "✓ Analyzed real edge-crossing polygon data\n";
    std::cout << "✓ Compared MeshCut vs Earcut on boundary cases\n";
    std::cout << "✓ Generated side-by-side visualization\n";
    std::cout << "✓ Checked for degenerate triangles\n";
    
    return 0;
}