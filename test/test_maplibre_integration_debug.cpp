#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include <limits>
#include <curl/curl.h>
#include <zlib.h>

// MapLibre coordinate system constants (from your integration)
constexpr int16_t EXTENT = 8192;  // From util::EXTENT in MapLibre

// Simulate MapLibre's GeometryCoordinate type
struct GeometryCoordinate {
    int16_t x, y;
    GeometryCoordinate(int16_t x_, int16_t y_) : x(x_), y(y_) {}
};

// Simulate MapLibre's coordinate conversion
GeometryCoordinate doubleToGeometry(double x, double y) {
    return GeometryCoordinate(
        static_cast<int16_t>(std::round(x)),
        static_cast<int16_t>(std::round(y))
    );
}

// CURL data structure
struct HTTPResponse {
    std::vector<uint8_t> data;
};

size_t writeCallback(void* contents, size_t size, size_t nmemb, HTTPResponse* response) {
    size_t totalSize = size * nmemb;
    uint8_t* data = static_cast<uint8_t*>(contents);
    response->data.insert(response->data.end(), data, data + totalSize);
    return totalSize;
}

bool downloadTile(const std::string& url, std::vector<uint8_t>& tileData) {
    CURL* curl = curl_easy_init();
    if (!curl) return false;
    
    HTTPResponse response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);
    
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    
    if (res == CURLE_OK) {
        tileData = std::move(response.data);
        return true;
    }
    return false;
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

// Enhanced MVT parser that finds edge-crossing polygons
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
        bool crossesEdge;  // NEW: tracks if polygon crosses tile boundary
        
        // Helper to check if polygon crosses tile edges
        void analyzeBoundaries() {
            crossesEdge = false;
            for (const auto& pt : points) {
                // Check if any vertex is on or very close to tile boundaries
                if (pt.x <= 0.0 || pt.x >= 4096.0 || pt.y <= 0.0 || pt.y >= 4096.0) {
                    crossesEdge = true;
                    break;
                }
                // Also check for vertices very close to edges (tolerance for edge detection)
                if (pt.x <= 10.0 || pt.x >= 4086.0 || pt.y <= 10.0 || pt.y >= 4086.0) {
                    crossesEdge = true;
                    break;
                }
            }
        }
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
            
            if (fieldNumber == 3 && wireType == 2) {
                uint64_t layerLength = readVarint();
                size_t layerEnd = std::min(pos + layerLength, size);
                
                auto layerPolygons = parseLayer(layerEnd);
                for (auto& poly : layerPolygons) {
                    poly.analyzeBoundaries();  // Analyze edge crossing
                    polygons.push_back(std::move(poly));
                }
                
                if (polygons.size() >= 20) {
                    std::cout << "Found " << polygons.size() << " polygons, stopping for analysis\n";
                    break;
                }
            } else {
                skipField(wireType);
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
            
            if (fieldNumber == 1 && wireType == 2) { // Layer name
                uint64_t length = readVarint();
                layerName = readString(length);
            } else if (fieldNumber == 2 && wireType == 2) { // Feature
                uint64_t featureLength = readVarint();
                size_t featureEnd = std::min(pos + featureLength, size);
                auto polygon = parseFeature(featureEnd, keys, values, layerName, extent);
                if (!polygon.points.empty()) {
                    polygons.push_back(std::move(polygon));
                }
            } else if (fieldNumber == 3 && wireType == 2) { // Keys
                uint64_t length = readVarint();
                keys.push_back(readString(length));
            } else if (fieldNumber == 4 && wireType == 2) { // Values
                uint64_t length = readVarint();
                values.push_back(readString(length));
            } else if (fieldNumber == 5 && wireType == 0) { // Extent
                extent = static_cast<uint32_t>(readVarint());
            } else {
                skipField(wireType);
            }
        }
        
        return polygons;
    }
    
    Polygon parseFeature(size_t featureEnd, const std::vector<std::string>& keys, 
                        const std::vector<std::string>& values, 
                        const std::string& layerName, uint32_t extent) {
        Polygon polygon;
        polygon.layerName = layerName;
        
        while (pos < featureEnd && pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 4 && wireType == 2) { // Geometry
                uint64_t geomLength = readVarint();
                size_t geomEnd = std::min(pos + geomLength, size);
                polygon.points = parseGeometry(geomEnd, extent);
                break;
            } else {
                skipField(wireType);
            }
        }
        
        return polygon;
    }
    
    std::vector<Point> parseGeometry(size_t geomEnd, uint32_t extent) {
        std::vector<Point> points;
        
        int32_t x = 0, y = 0;
        bool inRing = false;
        
        while (pos < geomEnd && pos < size) {
            uint32_t command = static_cast<uint32_t>(readVarint());
            uint32_t id = command & 0x7;
            uint32_t count = command >> 3;
            
            if (id == 1) { // MoveTo
                inRing = true;
                for (uint32_t i = 0; i < count && pos < geomEnd; i++) {
                    int32_t dx = static_cast<int32_t>(readVarint());
                    int32_t dy = static_cast<int32_t>(readVarint());
                    dx = (dx >> 1) ^ (-(dx & 1));
                    dy = (dy >> 1) ^ (-(dy & 1));
                    x += dx;
                    y += dy;
                    
                    points.push_back({static_cast<double>(x), static_cast<double>(y)});
                }
            } else if (id == 2 && inRing) { // LineTo
                for (uint32_t i = 0; i < count && pos < geomEnd; i++) {
                    int32_t dx = static_cast<int32_t>(readVarint());
                    int32_t dy = static_cast<int32_t>(readVarint());
                    dx = (dx >> 1) ^ (-(dx & 1));
                    dy = (dy >> 1) ^ (-(dy & 1));
                    x += dx;
                    y += dy;
                    
                    points.push_back({static_cast<double>(x), static_cast<double>(y)});
                }
            } else if (id == 7) { // ClosePath
                if (inRing && !points.empty()) {
                    // Close the ring if not already closed
                    if (points.back().x != points[0].x || points.back().y != points[0].y) {
                        points.push_back(points[0]);
                    }
                }
                inRing = false;
            }
        }
        
        return points;
    }
};

// MapLibre-style triangulation function that mimics your integration
struct TriangulationResult {
    std::vector<uint32_t> indices;
    std::vector<std::pair<double, double>> vertices;
    bool usedMeshcut;
    
    TriangulationResult(std::vector<uint32_t> idx, bool meshcut = false) 
        : indices(std::move(idx)), usedMeshcut(meshcut) {}
    TriangulationResult(std::vector<uint32_t> idx, std::vector<std::pair<double, double>> verts, bool meshcut = true) 
        : indices(std::move(idx)), vertices(std::move(verts)), usedMeshcut(meshcut) {}
};

TriangulationResult maplibreTriangulate(const std::vector<double>& polygon, bool useMeshcut = true, uint32_t gridSize = 32) {
    if (useMeshcut) {
        // Replicate your MapLibre meshcut integration
        meshcut::MeshCutOptions options;
        
        // Configure grid exactly like in your integration
        const double extentDouble = static_cast<double>(EXTENT);
        options.gridWidth = gridSize;
        options.gridHeight = gridSize;
        options.cellSize = (2.0 * extentDouble) / gridSize;
        options.gridOriginX = -extentDouble;
        options.gridOriginY = -extentDouble;
        options.diagonalNE = true;
        
        try {
            auto meshResult = meshcut::meshcut_full(polygon, {}, options);
            
            if (!meshResult.indices.empty() && !meshResult.vertices.empty()) {
                // Convert to coordinate pairs
                std::vector<std::pair<double, double>> allVertices;
                allVertices.reserve(meshResult.vertices.size() / 2);
                
                for (size_t i = 0; i < meshResult.vertices.size(); i += 2) {
                    allVertices.emplace_back(meshResult.vertices[i], meshResult.vertices[i + 1]);
                }
                
                return TriangulationResult(std::move(meshResult.indices), std::move(allVertices), true);
            }
        } catch (const std::exception& e) {
            std::cout << "MeshCut failed: " << e.what() << ", falling back to earcut\n";
        }
    }
    
    // Earcut fallback
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({});
    for (size_t i = 0; i < polygon.size(); i += 2) {
        earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
    }
    auto earcutIndices = mapbox::earcut<uint32_t>(earcutPolygon);
    
    return TriangulationResult(std::move(earcutIndices), false);
}

// Analyze coordinate overflow issues (like in your integration)
struct CoordinateAnalysis {
    size_t overflowCount = 0;
    size_t totalVertices = 0;
    double minX = 999999, maxX = -999999, minY = 999999, maxY = -999999;
    std::vector<std::pair<double, double>> overflowVertices;
    
    void analyze(const std::vector<std::pair<double, double>>& vertices) {
        totalVertices = vertices.size();
        
        for (const auto& vertex : vertices) {
            double x = vertex.first;
            double y = vertex.second;
            
            minX = std::min(minX, x);
            maxX = std::max(maxX, x);
            minY = std::min(minY, y);
            maxY = std::max(maxY, y);
            
            // Check for int16_t overflow (like in your MapLibre integration)
            if (x < std::numeric_limits<int16_t>::min() || x > std::numeric_limits<int16_t>::max() ||
                y < std::numeric_limits<int16_t>::min() || y > std::numeric_limits<int16_t>::max()) {
                overflowCount++;
                overflowVertices.emplace_back(x, y);
            }
        }
    }
    
    void print() {
        std::cout << "COORDINATE ANALYSIS:\n";
        std::cout << "- Total vertices: " << totalVertices << "\n";
        std::cout << "- Range: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\n";
        std::cout << "- Overflow vertices: " << overflowCount << " (out of int16_t range)\n";
        
        if (overflowCount > 0 && overflowCount <= 5) {
            std::cout << "- Overflow examples: ";
            for (size_t i = 0; i < std::min(overflowVertices.size(), size_t(3)); i++) {
                std::cout << "(" << overflowVertices[i].first << "," << overflowVertices[i].second << ") ";
            }
            std::cout << "\n";
        }
    }
};

// SVG visualization for edge artifacts
class EdgeArtifactVisualizer {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double width, height, scale;
    
public:
    EdgeArtifactVisualizer(const std::string& filename, const std::vector<double>& polygon) {
        // Calculate bounds
        minX = maxX = polygon[0];
        minY = maxY = polygon[1];
        for (size_t i = 0; i < polygon.size(); i += 2) {
            minX = std::min(minX, polygon[i]);
            maxX = std::max(maxX, polygon[i]);
            minY = std::min(minY, polygon[i+1]);
            maxY = std::max(maxY, polygon[i+1]);
        }
        
        // Add some padding
        double padX = (maxX - minX) * 0.1;
        double padY = (maxY - minY) * 0.1;
        minX -= padX; maxX += padX;
        minY -= padY; maxY += padY;
        
        width = maxX - minX;
        height = maxY - minY;
        scale = 600.0 / std::max(width, height);
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"" << (width * scale * 2 + 150) << "\" ";  // Side by side
        file << "height=\"" << (height * scale + 200) << "\" ";
        file << "viewBox=\"0 0 " << (width * scale * 2 + 150) << " " << (height * scale + 200) << "\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    }
    
    ~EdgeArtifactVisualizer() {
        file << "</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale;
        double ty = (maxY - y) * scale;
        if (rightPanel) tx += width * scale + 75;  // Move to right panel
        return {tx + 50, ty + 100};  // Add margins
    }
    
    void drawTitle(const std::string& title, bool rightPanel = false) {
        double titleX = rightPanel ? width * scale + 75 + (width * scale / 2) : (width * scale / 2);
        file << "<text x=\"" << (titleX + 50) << "\" y=\"50\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"18\" font-weight=\"bold\">" << title << "</text>\n";
    }
    
    void drawTileBoundaries(bool rightPanel = false) {
        // Draw tile boundaries at 0, 4096
        file << "<g stroke=\"#ff0000\" stroke-width=\"3\" stroke-dasharray=\"5,5\" fill=\"none\" opacity=\"0.7\">\n";
        
        // Vertical boundaries
        auto [x_left_top, y_left_top] = transform(0, maxY, rightPanel);
        auto [x_left_bottom, y_left_bottom] = transform(0, minY, rightPanel);
        auto [x_right_top, y_right_top] = transform(4096, maxY, rightPanel);
        auto [x_right_bottom, y_right_bottom] = transform(4096, minY, rightPanel);
        
        file << "<line x1=\"" << x_left_top << "\" y1=\"" << y_left_top << "\" x2=\"" << x_left_bottom << "\" y2=\"" << y_left_bottom << "\"/>\n";
        file << "<line x1=\"" << x_right_top << "\" y1=\"" << y_right_top << "\" x2=\"" << x_right_bottom << "\" y2=\"" << y_right_bottom << "\"/>\n";
        
        // Horizontal boundaries  
        auto [x_top_left, y_top_left] = transform(minX, 0, rightPanel);
        auto [x_top_right, y_top_right] = transform(maxX, 0, rightPanel);
        auto [x_bottom_left, y_bottom_left] = transform(minX, 4096, rightPanel);
        auto [x_bottom_right, y_bottom_right] = transform(maxX, 4096, rightPanel);
        
        file << "<line x1=\"" << x_top_left << "\" y1=\"" << y_top_left << "\" x2=\"" << x_top_right << "\" y2=\"" << y_top_right << "\"/>\n";
        file << "<line x1=\"" << x_bottom_left << "\" y1=\"" << y_bottom_left << "\" x2=\"" << x_bottom_right << "\" y2=\"" << y_bottom_right << "\"/>\n";
        
        file << "</g>\n";
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color, bool rightPanel = false) {
        if (polygon.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i + 1], rightPanel);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"0.2\" ";
        file << "stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
    }
    
    void drawTriangles(const std::vector<std::pair<double, double>>& vertices, 
                      const std::vector<uint32_t>& indices, 
                      const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.8\" opacity=\"0.6\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1]; 
            uint32_t i2 = indices[i + 2];
            
            if (i0 < vertices.size() && i1 < vertices.size() && i2 < vertices.size()) {
                auto [x0, y0] = transform(vertices[i0].first, vertices[i0].second, rightPanel);
                auto [x1, y1] = transform(vertices[i1].first, vertices[i1].second, rightPanel);
                auto [x2, y2] = transform(vertices[i2].first, vertices[i2].second, rightPanel);
                
                file << "<polygon points=\"" << x0 << "," << y0 << " " 
                     << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
            }
        }
        
        file << "</g>\n";
    }
    
    void drawEarcutTriangles(const std::vector<double>& polygon, 
                            const std::vector<uint32_t>& indices, 
                            const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.8\" opacity=\"0.6\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1]; 
            uint32_t i2 = indices[i + 2];
            
            if (i0*2+1 < polygon.size() && i1*2+1 < polygon.size() && i2*2+1 < polygon.size()) {
                auto [x0, y0] = transform(polygon[i0*2], polygon[i0*2+1], rightPanel);
                auto [x1, y1] = transform(polygon[i1*2], polygon[i1*2+1], rightPanel);
                auto [x2, y2] = transform(polygon[i2*2], polygon[i2*2+1], rightPanel);
                
                file << "<polygon points=\"" << x0 << "," << y0 << " " 
                     << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
            }
        }
        
        file << "</g>\n";
    }
    
    void addStats(const std::string& stats, bool rightPanel = false) {
        double statsX = rightPanel ? width * scale + 75 + 10 : 10;
        file << "<text x=\"" << (statsX + 50) << "\" y=\"" << (height * scale + 150) << "\" ";
        file << "font-family=\"Arial\" font-size=\"12\" fill=\"#666\">";
        file << stats << "</text>\n";
    }
};

int main() {
    std::cout << "=== MapLibre Integration Debug Test ===\n";
    std::cout << "Analyzing edge-crossing artifacts in MeshCut vs Earcut\n\n";
    
    // Try multiple tiles to find one with polygons
    std::vector<std::string> tileUrls = {
        "https://demotiles.maplibre.org/data/v3/11/1083/724.pbf",  // Innsbruck (known to work)
        "https://demotiles.maplibre.org/data/v3/8/136/91.pbf",     // Zurich area  
        "https://demotiles.maplibre.org/data/v3/10/544/362.pbf"    // Another European tile
    };
    
    std::vector<uint8_t> tileData;
    std::string workingUrl;
    
    for (const auto& tileUrl : tileUrls) {
        std::cout << "Trying tile: " << tileUrl << "\n";
        if (downloadTile(tileUrl, tileData)) {
            std::cout << "Downloaded " << tileData.size() << " bytes\n";
            
            // Try decompression if needed
            std::vector<uint8_t> decompressed = decompressGzip(tileData);
            if (!decompressed.empty()) {
                std::cout << "Decompressed to " << decompressed.size() << " bytes\n";
                tileData = std::move(decompressed);
            }
            
            // Quick parse to check if we have polygons
            MVTParser parser(tileData);
            auto polygons = parser.extractPolygons();
            std::cout << "Found " << polygons.size() << " polygons\n";
            
            if (!polygons.empty()) {
                workingUrl = tileUrl;
                std::cout << "✓ Using this tile for testing\n\n";
                break;
            }
        }
        std::cout << "No polygons in this tile, trying next...\n\n";
    }
    
    if (tileData.empty()) {
        std::cerr << "Failed to download any usable tile\n";
        return 1;
    }
    
    // Parse MVT data
    MVTParser parser(tileData);
    auto polygons = parser.extractPolygons();
    
    std::cout << "Found " << polygons.size() << " polygons\n";
    
    // Find edge-crossing polygons
    std::vector<const MVTParser::Polygon*> edgePolygons;
    std::vector<const MVTParser::Polygon*> interiorPolygons;
    
    for (const auto& poly : polygons) {
        if (poly.crossesEdge) {
            edgePolygons.push_back(&poly);
        } else {
            interiorPolygons.push_back(&poly);
        }
    }
    
    std::cout << "Edge-crossing polygons: " << edgePolygons.size() << "\n";
    std::cout << "Interior polygons: " << interiorPolygons.size() << "\n\n";
    
    if (edgePolygons.empty()) {
        std::cout << "No edge-crossing polygons found in this tile. Trying with first available polygon...\n";
        if (!polygons.empty()) {
            edgePolygons.push_back(&polygons[0]);
        }
    }
    
    if (edgePolygons.empty()) {
        std::cerr << "No polygons found to test\n";
        return 1;
    }
    
    // Test each edge-crossing polygon
    for (size_t polyIdx = 0; polyIdx < std::min(edgePolygons.size(), size_t(3)); polyIdx++) {
        const auto& testPoly = *edgePolygons[polyIdx];
        
        std::cout << "=== Testing Edge Polygon " << (polyIdx + 1) << " ===\n";
        std::cout << "Layer: " << testPoly.layerName << "\n";
        std::cout << "Vertices: " << testPoly.points.size() << "\n";
        std::cout << "Crosses edge: " << (testPoly.crossesEdge ? "YES" : "NO") << "\n";
        
        // Convert to flat coordinate array
        std::vector<double> flatPolygon;
        for (const auto& pt : testPoly.points) {
            flatPolygon.push_back(pt.x);
            flatPolygon.push_back(pt.y);
        }
        
        // Print coordinate range
        double minX = *std::min_element(flatPolygon.begin(), flatPolygon.begin() + flatPolygon.size());
        double maxX = *std::max_element(flatPolygon.begin(), flatPolygon.begin() + flatPolygon.size());
        auto minYit = std::min_element(flatPolygon.begin() + 1, flatPolygon.end(), [](double a, double b) {
            return a < b;
        });
        auto maxYit = std::max_element(flatPolygon.begin() + 1, flatPolygon.end(), [](double a, double b) {
            return a < b;
        });
        double minY = (minYit != flatPolygon.end()) ? *minYit : 0;
        double maxY = (maxYit != flatPolygon.end()) ? *maxYit : 0;
        
        std::cout << "Coordinate range: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\n";
        
        // Test with both triangulation methods
        std::cout << "\n--- Earcut Triangulation ---\n";
        auto earcutResult = maplibreTriangulate(flatPolygon, false);
        std::cout << "Triangles: " << earcutResult.indices.size() / 3 << "\n";
        std::cout << "Uses original vertices: " << flatPolygon.size() / 2 << "\n";
        
        std::cout << "\n--- MeshCut Triangulation ---\n";
        auto meshcutResult = maplibreTriangulate(flatPolygon, true, 32);
        std::cout << "Triangles: " << meshcutResult.indices.size() / 3 << "\n";
        std::cout << "Total vertices: " << meshcutResult.vertices.size() << "\n";
        
        // Analyze coordinate overflow (critical for MapLibre)
        CoordinateAnalysis analysis;
        analysis.analyze(meshcutResult.vertices);
        analysis.print();
        
        // Create visualization
        std::string filename = "edge_polygon_debug_" + std::to_string(polyIdx + 1) + ".svg";
        EdgeArtifactVisualizer viz(filename, flatPolygon);
        
        // Left panel: Earcut
        viz.drawTitle("Earcut (Original)", false);
        viz.drawTileBoundaries(false);
        viz.drawPolygon(flatPolygon, "#4CAF50", false);
        viz.drawEarcutTriangles(flatPolygon, earcutResult.indices, "#FF5722", false);
        viz.addStats("Earcut: " + std::to_string(earcutResult.indices.size() / 3) + " triangles", false);
        
        // Right panel: MeshCut
        viz.drawTitle("MeshCut (Grid-Enhanced)", true);
        viz.drawTileBoundaries(true);
        viz.drawPolygon(flatPolygon, "#4CAF50", true);
        viz.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", true);
        viz.addStats("MeshCut: " + std::to_string(meshcutResult.indices.size() / 3) + " triangles, " + 
                    std::to_string(analysis.overflowCount) + " overflow", true);
        
        std::cout << "Visualization saved: " << filename << "\n";
        
        // Print detailed analysis
        std::cout << "\n--- Potential Issues ---\n";
        if (analysis.overflowCount > 0) {
            std::cout << "⚠️  COORDINATE OVERFLOW: " << analysis.overflowCount << " vertices exceed int16_t range\n";
            std::cout << "   This will cause truncation when converting to GeometryCoordinate\n";
        }
        
        if (testPoly.crossesEdge) {
            std::cout << "⚠️  EDGE CROSSING: Polygon touches tile boundaries\n";
            std::cout << "   May cause artifacts when tiles are combined\n";
        }
        
        double triangleDensity = (double)meshcutResult.indices.size() / (double)earcutResult.indices.size();
        if (triangleDensity > 5.0) {
            std::cout << "⚠️  HIGH TRIANGLE DENSITY: " << std::fixed << std::setprecision(1) << triangleDensity << "x increase\n";
            std::cout << "   May impact rendering performance\n";
        }
        
        std::cout << "\n";
    }
    
    std::cout << "=== Summary ===\n";
    std::cout << "Generated edge polygon debug visualizations\n";
    std::cout << "Key findings:\n";
    std::cout << "1. Check for coordinate overflow (int16_t truncation)\n";
    std::cout << "2. Analyze edge-crossing behavior differences\n";
    std::cout << "3. Compare triangle density and quality\n";
    std::cout << "\nOpen the generated SVG files to see visual differences!\n";
    
    return 0;
}