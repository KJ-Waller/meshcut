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

class TileTestSuite {
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
    
    struct TileInfo {
        int z, x, y;
        uint32_t extent;
        std::string source;
    };
    
    struct TestPolygon {
        std::vector<std::vector<Point>> rings;  // First ring = outer boundary, rest = holes
        std::string layerName;
        TileInfo tile;
        int id;
        
        // Get all points from all rings as a flat coordinate array (for MeshCut compatibility)
        std::vector<double> toFlatCoords() const {
            std::vector<double> coords;
            for (const auto& ring : rings) {
                for (const auto& pt : ring) {
                    coords.push_back(pt.x);
                    coords.push_back(pt.y);
                }
            }
            return coords;
        }
        
        // Get rings in Earcut-compatible format
        std::vector<std::vector<std::array<double, 2>>> toEarcutFormat() const {
            std::vector<std::vector<std::array<double, 2>>> earcutPolygon;
            for (const auto& ring : rings) {
                earcutPolygon.push_back({});
                for (const auto& pt : ring) {
                    earcutPolygon.back().push_back({{pt.x, pt.y}});
                }
            }
            return earcutPolygon;
        }
        
        // Get flattened coordinates in Earcut order (for visualization of Earcut results)
        std::vector<double> toEarcutCoords() const {
            std::vector<double> coords;
            for (const auto& ring : rings) {
                for (const auto& pt : ring) {
                    coords.push_back(pt.x);
                    coords.push_back(pt.y);
                }
            }
            return coords;
        }
        
        // Get coordinates and hole indices for MeshCut (earcut-compatible format)
        std::pair<std::vector<double>, std::vector<uint32_t>> toMeshCutFormat() const {
            std::vector<double> coords;
            std::vector<uint32_t> holes;
            
            uint32_t vertexIndex = 0;
            
            // Add outer ring
            if (!rings.empty()) {
                for (const auto& pt : rings[0]) {
                    coords.push_back(pt.x);
                    coords.push_back(pt.y);
                    vertexIndex++;
                }
                
                // Add holes
                for (size_t i = 1; i < rings.size(); i++) {
                    holes.push_back(vertexIndex); // Start index of this hole
                    for (const auto& pt : rings[i]) {
                        coords.push_back(pt.x);
                        coords.push_back(pt.y);
                        vertexIndex++;
                    }
                }
            }
            
            return {std::move(coords), std::move(holes)};
        }
        
        void printInfo() const {
            if (rings.empty() || rings[0].empty()) return;
            
            double minX = rings[0][0].x, maxX = rings[0][0].x;
            double minY = rings[0][0].y, maxY = rings[0][0].y;
            
            for (const auto& ring : rings) {
                for (const auto& pt : ring) {
                    minX = std::min(minX, pt.x);
                    maxX = std::max(maxX, pt.x);
                    minY = std::min(minY, pt.y);
                    maxY = std::max(maxY, pt.y);
                }
            }
            
            int totalVertices = 0;
            for (const auto& ring : rings) {
                totalVertices += ring.size();
            }
            
            std::cout << "Polygon " << id << " (" << layerName << "):\n";
            std::cout << "  Rings: " << rings.size() << " (1 outer";
            if (rings.size() > 1) std::cout << ", " << (rings.size() - 1) << " holes";
            std::cout << ")\n";
            std::cout << "  Total vertices: " << totalVertices << "\n";
            std::cout << "  Bounds: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\n";
            
            bool extendsOutside = (minX < 0 || maxX > tile.extent || minY < 0 || maxY > tile.extent);
            std::cout << "  Extends beyond tile: " << (extendsOutside ? "YES" : "NO") << "\n";
            
            if (extendsOutside) {
                std::cout << "  Buffer zones: ";
                if (minX < 0) std::cout << "Left(" << (-minX) << ") ";
                if (maxX > tile.extent) std::cout << "Right(" << (maxX - tile.extent) << ") ";
                if (minY < 0) std::cout << "Bottom(" << (-minY) << ") ";
                if (maxY > tile.extent) std::cout << "Top(" << (maxY - tile.extent) << ") ";
                std::cout << "\n";
            }
            
            double coverage = ((maxX - minX) * (maxY - minY)) / (tile.extent * tile.extent) * 100.0;
            std::cout << "  Tile coverage: " << std::fixed << std::setprecision(1) << coverage << "%\n";
            
            // Debug output for polygon 202
            if (id == 202) {
                std::cout << "\n=== DETAILED POLYGON 202 ANALYSIS ===\n";
                for (size_t i = 0; i < rings.size(); i++) {
                    const auto& ring = rings[i];
                    std::cout << "Ring " << i << " (" << (i == 0 ? "outer" : "hole") << "):\n";
                    std::cout << "  Vertices: " << ring.size() << "\n";
                    
                    if (!ring.empty()) {
                        std::cout << "  First: (" << ring[0].x << ", " << ring[0].y << ")\n";
                        std::cout << "  Last: (" << ring.back().x << ", " << ring.back().y << ")\n";
                        
                        // Calculate winding order (shoelace formula)
                        double area = 0.0;
                        for (size_t j = 0; j < ring.size(); j++) {
                            size_t next = (j + 1) % ring.size();
                            area += (ring[next].x - ring[j].x) * (ring[next].y + ring[j].y);
                        }
                        bool clockwise = area > 0;
                        std::cout << "  Winding: " << (clockwise ? "CLOCKWISE" : "COUNTER-CLOCKWISE") << " (area: " << area << ")\n";
                        
                        // For holes, check if they're actually inside the outer ring
                        if (i > 0 && !rings[0].empty()) {
                            // Simple point-in-polygon test for first point of hole
                            auto& pt = ring[0];
                            int crossings = 0;
                            const auto& outer = rings[0];
                            for (size_t j = 0; j < outer.size() - 1; j++) {
                                if (((outer[j].y <= pt.y) && (pt.y < outer[j+1].y)) ||
                                    ((outer[j+1].y <= pt.y) && (pt.y < outer[j].y))) {
                                    double intersectX = outer[j].x + (pt.y - outer[j].y) * (outer[j+1].x - outer[j].x) / (outer[j+1].y - outer[j].y);
                                    if (pt.x < intersectX) crossings++;
                                }
                            }
                            bool inside = (crossings % 2) == 1;
                            std::cout << "  Inside outer ring: " << (inside ? "YES" : "NO") << " (crossings: " << crossings << ")\n";
                        }
                        
                        // Print first few and last few points to see if there are connections
                        std::cout << "  All points: ";
                        for (size_t j = 0; j < std::min(ring.size(), 20UL); j++) {
                            std::cout << "(" << ring[j].x << ", " << ring[j].y << ")";
                            if (j < std::min(ring.size(), 20UL) - 1) std::cout << " -> ";
                        }
                        if (ring.size() > 20) std::cout << " ... [" << (ring.size() - 20) << " more]";
                        std::cout << "\n";
                    }
                    std::cout << "\n";
                }
                std::cout << "=== END POLYGON 202 ANALYSIS ===\n\n";
            }
            
            std::cout << "\n";
        }
    };
    
    TileTestSuite(const std::vector<uint8_t>& tileData, int z, int x, int y, const std::string& source) 
        : data(tileData.data()), size(tileData.size()), pos(0) {
        tileInfo.z = z;
        tileInfo.x = x;
        tileInfo.y = y;
        tileInfo.extent = 4096;
        tileInfo.source = source;
    }
    
    TileInfo tileInfo;
    
    std::vector<TestPolygon> extractTestPolygons(int maxPolygons = 1000) {
        std::vector<TestPolygon> polygons;
        pos = 0;
        
        std::cout << "=== TILE TEST SUITE ===\n";
        std::cout << "Source: " << tileInfo.source << "\n";
        std::cout << "Tile: z=" << tileInfo.z << ", x=" << tileInfo.x << ", y=" << tileInfo.y << "\n";
        std::cout << "Data size: " << size << " bytes\n\n";
        
        int polygonId = 1;
        
        while (pos < size && polygons.size() < maxPolygons) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 2) {
                uint64_t layerLength = readVarint();
                size_t layerEnd = std::min(pos + layerLength, size);
                
                auto layerPolygons = parseLayer(layerEnd, polygonId);
                for (auto& poly : layerPolygons) {
                    poly.tile = tileInfo;
                    polygons.push_back(std::move(poly));
                    if (polygons.size() >= maxPolygons) break;
                }
            } else {
                skipField(wireType);
            }
        }
        
        std::cout << "Extracted " << polygons.size() << " test polygons\n\n";
        
        for (auto& polygon : polygons) {
            polygon.printInfo();
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
    
    std::vector<TestPolygon> parseLayer(size_t layerEnd, int& polygonId) {
        std::vector<TestPolygon> polygons;
        std::vector<std::string> keys;
        std::vector<std::string> values;
        std::string layerName = "unknown";
        
        int pointCount = 0, lineCount = 0, polygonCount = 0, unknownCount = 0;
        
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
                auto [polygon, geomType] = parseFeature(featureEnd, keys, values, layerName, polygonId);
                
                // Count geometry types
                if (geomType == 1) pointCount++;
                else if (geomType == 2) lineCount++;
                else if (geomType == 3) polygonCount++;
                else unknownCount++;
                
                if (!polygon.rings.empty() && !polygon.rings[0].empty() && polygon.rings[0].size() >= 4) {
                    polygons.push_back(std::move(polygon));
                    polygonId++;
                }
            } else if (fieldNumber == 3 && wireType == 2) { // Keys
                uint64_t length = readVarint();
                keys.push_back(readString(length));
            } else if (fieldNumber == 4 && wireType == 2) { // Values
                uint64_t length = readVarint();
                values.push_back(readString(length));
            } else if (fieldNumber == 5 && wireType == 0) { // Extent
                tileInfo.extent = static_cast<uint32_t>(readVarint());
            } else {
                skipField(wireType);
            }
        }
        
        if (pointCount > 0 || lineCount > 0 || polygonCount > 0 || unknownCount > 0) {
            std::cout << "Layer '" << layerName << "': ";
            if (pointCount > 0) std::cout << pointCount << " points, ";
            if (lineCount > 0) std::cout << lineCount << " lines, ";
            if (polygonCount > 0) std::cout << polygonCount << " polygons, ";
            if (unknownCount > 0) std::cout << unknownCount << " unknown";
            std::cout << " -> " << polygons.size() << " extracted\n";
        }
        
        return polygons;
    }
    
    std::pair<TestPolygon, uint32_t> parseFeature(size_t featureEnd, const std::vector<std::string>& keys, 
                                                 const std::vector<std::string>& values, 
                                                 const std::string& layerName, int id) {
        TestPolygon polygon;
        polygon.layerName = layerName;
        polygon.id = id;
        
        uint32_t geometryType = 0;  // 0=Unknown, 1=Point, 2=LineString, 3=Polygon
        
        while (pos < featureEnd && pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 0) { // Geometry type
                geometryType = static_cast<uint32_t>(readVarint());
            } else if (fieldNumber == 4 && wireType == 2) { // Geometry
                uint64_t geomLength = readVarint();
                size_t geomEnd = std::min(pos + geomLength, size);
                
                // Only process Polygon geometries (type 3), skip Points (1) and LineStrings (2)
                if (geometryType == 3) {
                    polygon.rings = parseGeometry(geomEnd);
                } else {
                    // Skip geometry data for non-polygon types
                    pos = geomEnd;
                }
                break;
            } else {
                skipField(wireType);
            }
        }
        
        return std::make_pair(std::move(polygon), geometryType);
    }
    
    std::vector<std::vector<Point>> parseGeometry(size_t geomEnd) {
        std::vector<std::vector<Point>> rings;
        std::vector<Point> currentRing;
        
        int32_t x = 0, y = 0;
        bool inRing = false;
        
        while (pos < geomEnd && pos < size) {
            uint32_t command = static_cast<uint32_t>(readVarint());
            uint32_t id = command & 0x7;
            uint32_t count = command >> 3;
            
            if (id == 1) { // MoveTo - starts a new ring
                // If we were in a ring, close it and add to rings
                if (inRing && !currentRing.empty()) {
                    rings.push_back(std::move(currentRing));
                    currentRing.clear();
                }
                
                inRing = true;
                for (uint32_t i = 0; i < count && pos < geomEnd; i++) {
                    int32_t dx = static_cast<int32_t>(readVarint());
                    int32_t dy = static_cast<int32_t>(readVarint());
                    dx = (dx >> 1) ^ (-(dx & 1));
                    dy = (dy >> 1) ^ (-(dy & 1));
                    x += dx;
                    y += dy;
                    
                    currentRing.push_back({static_cast<double>(x), static_cast<double>(y)});
                }
            } else if (id == 2 && inRing) { // LineTo
                for (uint32_t i = 0; i < count && pos < geomEnd; i++) {
                    int32_t dx = static_cast<int32_t>(readVarint());
                    int32_t dy = static_cast<int32_t>(readVarint());
                    dx = (dx >> 1) ^ (-(dx & 1));
                    dy = (dy >> 1) ^ (-(dy & 1));
                    x += dx;
                    y += dy;
                    
                    currentRing.push_back({static_cast<double>(x), static_cast<double>(y)});
                }
            } else if (id == 7) { // ClosePath
                if (inRing && !currentRing.empty()) {
                    // Ensure ring is closed
                    if (currentRing.back().x != currentRing[0].x || currentRing.back().y != currentRing[0].y) {
                        currentRing.push_back(currentRing[0]);
                    }
                    rings.push_back(std::move(currentRing));
                    currentRing.clear();
                }
                inRing = false;
            }
        }
        
        // Handle case where final ring wasn't closed with ClosePath
        if (inRing && !currentRing.empty()) {
            if (currentRing.back().x != currentRing[0].x || currentRing.back().y != currentRing[0].y) {
                currentRing.push_back(currentRing[0]);
            }
            rings.push_back(std::move(currentRing));
        }
        
        return rings;
    }
};

class VisualTestSuite {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double width, height, scale;
    uint32_t tileExtent;
    
public:
    VisualTestSuite(const std::string& filename, const std::vector<TileTestSuite::TestPolygon>& polygons) {
        if (polygons.empty()) return;
        
        tileExtent = polygons[0].tile.extent;
        
        // Calculate bounds from all polygons, expand for better visualization
        if (!polygons[0].rings.empty() && !polygons[0].rings[0].empty()) {
            minX = maxX = polygons[0].rings[0][0].x;
            minY = maxY = polygons[0].rings[0][0].y;
        }
        
        for (const auto& poly : polygons) {
            for (const auto& ring : poly.rings) {
                for (const auto& pt : ring) {
                    minX = std::min(minX, pt.x);
                    maxX = std::max(maxX, pt.x);
                    minY = std::min(minY, pt.y);
                    maxY = std::max(maxY, pt.y);
                }
            }
        }
        
        // Expand bounds to show tile context
        minX = std::min(minX, -200.0);
        maxX = std::max(maxX, tileExtent + 200.0);
        minY = std::min(minY, -200.0);
        maxY = std::max(maxY, tileExtent + 200.0);
        
        width = maxX - minX;
        height = maxY - minY;
        scale = 800.0 / std::max(width, height);
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"" << (width * scale * 2 + 200) << "\" ";
        file << "height=\"" << (height * scale + 200) << "\" ";
        file << "viewBox=\"0 0 " << (width * scale * 2 + 200) << " " << (height * scale + 200) << "\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        
        // Add title with tile info
        file << "<text x=\"" << (width * scale + 100) << "\" y=\"30\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"18\" font-weight=\"bold\">";
        file << "Tile Test Suite: " << polygons[0].tile.source << "</text>\n";
        
        file << "<text x=\"" << (width * scale + 100) << "\" y=\"50\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"14\" fill=\"#666\">";
        file << "z=" << polygons[0].tile.z << ", x=" << polygons[0].tile.x << ", y=" << polygons[0].tile.y << "</text>\n";
    }
    
    ~VisualTestSuite() {
        file << "</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale;
        double ty = (y - minY) * scale;  // Fixed: don't flip Y coordinate
        if (rightPanel) tx += width * scale + 100;
        return {tx + 50, ty + 100};
    }
    
    void drawTileBoundary(bool rightPanel = false) {
        file << "<g stroke=\"#ff0000\" stroke-width=\"2\" fill=\"none\" opacity=\"0.8\">\n";
        
        auto [x1, y1] = transform(0, 0, rightPanel);
        auto [x2, y2] = transform(tileExtent, 0, rightPanel);
        auto [x3, y3] = transform(tileExtent, tileExtent, rightPanel);
        auto [x4, y4] = transform(0, tileExtent, rightPanel);
        
        file << "<rect x=\"" << x1 << "\" y=\"" << y4 << "\" width=\"" << (x2-x1) << "\" height=\"" << (y1-y4) << "\"/>\n";
        file << "</g>\n";
        
        // Add boundary label
        auto [labelX, labelY] = transform(tileExtent/2, -50, rightPanel);
        file << "<text x=\"" << labelX << "\" y=\"" << labelY << "\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"11\" fill=\"red\">Tile Boundary (0-" << tileExtent << ")</text>\n";
    }
    
    void drawGrid(int gridSize, bool rightPanel = false) {
        double cellSize = (tileExtent + 256.0) / gridSize;
        double gridOriginX = -128.0;
        double gridOriginY = -128.0;
        
        file << "<g stroke=\"#cccccc\" stroke-width=\"0.5\" fill=\"none\" opacity=\"0.6\">\n";
        
        // Vertical grid lines
        for (int i = 0; i <= gridSize; i++) {
            double x = gridOriginX + i * cellSize;
            auto [x1, y1] = transform(x, gridOriginY, rightPanel);
            auto [x2, y2] = transform(x, gridOriginY + gridSize * cellSize, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        // Horizontal grid lines
        for (int j = 0; j <= gridSize; j++) {
            double y = gridOriginY + j * cellSize;
            auto [x1, y1] = transform(gridOriginX, y, rightPanel);
            auto [x2, y2] = transform(gridOriginX + gridSize * cellSize, y, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        file << "</g>\n";
    }
    
    void drawPolygon(const std::vector<double>& coords, const std::string& color, bool rightPanel = false, double opacity = 0.3) {
        if (coords.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < coords.size(); i += 2) {
            auto [x, y] = transform(coords[i], coords[i + 1], rightPanel);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"" << opacity << "\" ";
        file << "stroke=\"" << color << "\" stroke-width=\"1.5\"/>\n";
    }
    
    // Draw polygon with proper hole support using SVG path fill-rule
    void drawPolygonWithHoles(const std::vector<std::vector<TileTestSuite::Point>>& rings, const std::string& color, bool rightPanel = false, double opacity = 0.3) {
        if (rings.empty() || rings[0].empty()) return;
        
        file << "<path d=\"";
        
        // Draw outer ring
        const auto& outer = rings[0];
        auto [x0, y0] = transform(outer[0].x, outer[0].y, rightPanel);
        file << "M " << x0 << " " << y0 << " ";
        for (size_t i = 1; i < outer.size(); i++) {
            auto [x, y] = transform(outer[i].x, outer[i].y, rightPanel);
            file << "L " << x << " " << y << " ";
        }
        file << "Z ";
        
        // Draw holes
        for (size_t r = 1; r < rings.size(); r++) {
            const auto& hole = rings[r];
            if (hole.empty()) continue;
            
            auto [hx0, hy0] = transform(hole[0].x, hole[0].y, rightPanel);
            file << "M " << hx0 << " " << hy0 << " ";
            for (size_t i = 1; i < hole.size(); i++) {
                auto [hx, hy] = transform(hole[i].x, hole[i].y, rightPanel);
                file << "L " << hx << " " << hy << " ";
            }
            file << "Z ";
        }
        
        file << "\" fill=\"" << color << "\" fill-opacity=\"" << opacity << "\" ";
        file << "fill-rule=\"evenodd\" stroke=\"" << color << "\" stroke-width=\"1.5\"/>\n";
    }
    
    void drawTriangulation(const std::vector<double>& coords, const std::vector<uint32_t>& indices, 
                          const std::string& color, bool rightPanel = false, const std::string& method = "") {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.4\" opacity=\"0.7\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i];
            uint32_t i1 = indices[i + 1]; 
            uint32_t i2 = indices[i + 2];
            
            if (i0*2+1 < coords.size() && i1*2+1 < coords.size() && i2*2+1 < coords.size()) {
                auto [x0, y0] = transform(coords[i0*2], coords[i0*2+1], rightPanel);
                auto [x1, y1] = transform(coords[i1*2], coords[i1*2+1], rightPanel);
                auto [x2, y2] = transform(coords[i2*2], coords[i2*2+1], rightPanel);
                
                file << "<polygon points=\"" << x0 << "," << y0 << " " 
                     << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
            }
        }
        
        file << "</g>\n";
    }
    
    void drawMeshCutTriangulation(const std::vector<std::pair<double, double>>& vertices, 
                                 const std::vector<uint32_t>& indices, 
                                 const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.4\" opacity=\"0.7\">\n";
        
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
    
    void addPanelTitle(const std::string& title, bool rightPanel = false) {
        double titleX = rightPanel ? width * scale + 100 + (width * scale / 2) : (width * scale / 2);
        file << "<text x=\"" << (titleX + 50) << "\" y=\"90\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"16\" font-weight=\"bold\">" << title << "</text>\n";
    }
    
    void addStats(const std::string& stats, bool rightPanel = false) {
        double statsX = rightPanel ? width * scale + 100 + 10 : 10;
        file << "<text x=\"" << (statsX + 50) << "\" y=\"" << (height * scale + 150) << "\" ";
        file << "font-family=\"Arial\" font-size=\"12\" fill=\"#666\">" << stats << "</text>\n";
    }
};

meshcut::MeshCutOptions createMeshCutOptions(uint32_t tileExtent, uint32_t gridSize = 32) {
    meshcut::MeshCutOptions options;
    options.gridOriginX = -128.0;  
    options.gridOriginY = -128.0;  
    options.gridWidth = gridSize;
    options.gridHeight = gridSize;
    options.cellSize = (tileExtent + 256.0) / gridSize;
    options.diagonalNE = true;
    return options;
}

int main() {
    std::cout << "=== COMPREHENSIVE TILE TEST SUITE ===\n";
    std::cout << "Testing the exact same tile you're using in MapLibre Native\n\n";
    
    // Use the exact same tile source and coordinates you're testing
    int z = 14, x = 8712, y = 5741;  // Innsbruck coordinates from your testing
    std::string tileUrl = "https://demotiles.maplibre.org/tiles-omt/" + std::to_string(z) + "/" + 
                         std::to_string(x) + "/" + std::to_string(y) + ".pbf";
    
    std::cout << "Loading tile from your MapLibre source...\n";
    std::cout << "URL: " << tileUrl << "\n\n";
    
    std::vector<uint8_t> tileData;
    if (!downloadTile(tileUrl, tileData)) {
        std::cerr << "Failed to download tile\n";
        return 1;
    }
    
    std::cout << "Downloaded " << tileData.size() << " bytes\n";
    
    // Decompress if needed
    std::vector<uint8_t> decompressed = decompressGzip(tileData);
    if (!decompressed.empty()) {
        std::cout << "Decompressed to " << decompressed.size() << " bytes\n";
        tileData = std::move(decompressed);
    }
    
    // Extract polygons for testing
    TileTestSuite testSuite(tileData, z, x, y, "MapLibre Demo Tiles");
    auto testPolygons = testSuite.extractTestPolygons(1000);  // Get ALL polygons (up to 1000)
    
    if (testPolygons.empty()) {
        std::cerr << "No test polygons found\n";
        return 1;
    }
    
    // Set up MeshCut options
    auto meshcutOptions = createMeshCutOptions(testSuite.tileInfo.extent, 32);
    
    std::cout << "=== GRID CONFIGURATION ===\n";
    std::cout << "Grid size: 32x32 cells\n";
    std::cout << "Cell size: " << meshcutOptions.cellSize << " units\n";
    std::cout << "Grid covers: " << meshcutOptions.gridOriginX << " to " 
              << (meshcutOptions.gridOriginX + meshcutOptions.gridWidth * meshcutOptions.cellSize) << "\n";
    std::cout << "Tile extent: 0 to " << testSuite.tileInfo.extent << "\n\n";
    
    // Create combined visualization showing all polygons FIRST (before individual processing)
    std::cout << "=== COMBINED TILE VISUALIZATION ===\n";
    std::string combinedFilename = "test_tile_combined_analysis.svg";
    VisualTestSuite combinedViz(combinedFilename, testPolygons);
    
    // Left panel: All polygons with Earcut
    combinedViz.addPanelTitle("Complete Tile - Earcut", false);
    combinedViz.drawTileBoundary(false);
    
    // Right panel: All polygons with MeshCut
    combinedViz.addPanelTitle("Complete Tile - MeshCut", true);
    combinedViz.drawTileBoundary(true);
    combinedViz.drawGrid(32, true);
    
    // Draw all polygons in different colors (cycle through colors for all polygons)
    std::vector<std::string> colors = {"#4CAF50", "#2196F3", "#FF9800", "#9C27B0", "#F44336", 
                                      "#00BCD4", "#8BC34A", "#FFC107", "#795548", "#607D8B",
                                      "#E91E63", "#3F51B5", "#FF5722", "#673AB7", "#009688",
                                      "#CDDC39", "#FF6F00", "#E91E63", "#3F51B5", "#795548"};
    
    for (size_t i = 0; i < testPolygons.size(); i++) {
        // Cycle through colors if we have more polygons than colors
        std::string color = colors[i % colors.size()];
        combinedViz.drawPolygonWithHoles(testPolygons[i].rings, color, false, 0.15);  // Left panel (Earcut)
        combinedViz.drawPolygonWithHoles(testPolygons[i].rings, color, true, 0.15);   // Right panel (MeshCut now supports holes!)
    }
    
    combinedViz.addStats("All polygons: " + std::to_string(testPolygons.size()) + " from tile", false);
    combinedViz.addStats("32x32 terrain grid overlay", true);
    
    std::cout << "Combined tile visualization saved: " << combinedFilename << "\n";
    std::cout << "SUCCESS: Generated combined view with ALL " << testPolygons.size() << " polygons!\n\n";
    
    // Test each polygon individually (with error handling to prevent crashes)
    std::cout << "=== INDIVIDUAL POLYGON TESTS ===\n";
    int successCount = 0;
    int errorCount = 0;
    
    for (size_t i = 0; i < testPolygons.size(); i++) {
        const auto& poly = testPolygons[i];
        
        try {
            auto earcutCoords = poly.toEarcutCoords(); // For Earcut visualization (proper ring order)
            auto [meshcutCoords, meshcutHoles] = poly.toMeshCutFormat(); // For MeshCut (proper hole support)
            
            std::cout << "Testing Polygon " << poly.id << " (" << poly.layerName << ")...\n";
            
            // Test with Earcut (properly handling holes)
            auto earcutPolygon = poly.toEarcutFormat();
            
            // Debug output for polygon 202 inputs
            if (poly.id == 202) {
                std::cout << "\n=== POLYGON 202 INPUT DEBUG ===\n";
                std::cout << "Earcut polygon rings: " << earcutPolygon.size() << "\n";
                std::cout << "MeshCut coords: " << meshcutCoords.size()/2 << " vertices, holes at indices: [";
                for (size_t i = 0; i < meshcutHoles.size(); i++) {
                    std::cout << meshcutHoles[i];
                    if (i < meshcutHoles.size() - 1) std::cout << ", ";
                }
                std::cout << "]\n";
                std::cout << "=== END INPUT DEBUG ===\n\n";
            }
            
            auto earcutResult = mapbox::earcut<uint32_t>(earcutPolygon);
            
            // Test with MeshCut (now with proper hole support)
            auto meshcutResult = meshcut::meshcut_full(meshcutCoords, meshcutHoles, meshcutOptions);
            
            std::cout << "  Earcut: " << earcutCoords.size()/2 << " vertices → " << earcutResult.size()/3 << " triangles\n";
            std::cout << "  MeshCut: " << meshcutCoords.size()/2 << " input vertices → " 
                      << meshcutResult.vertices.size()/2 << " output vertices → " 
                      << meshcutResult.indices.size()/3 << " triangles\n";
            
            double improvement = (double)meshcutResult.indices.size() / (double)earcutResult.size();
            std::cout << "  Triangle ratio: " << std::fixed << std::setprecision(1) << improvement << "x\n";
            
            // Create individual comparison visualization
            std::string filename = "test_polygon_" + std::to_string(poly.id) + "_comparison.svg";
            VisualTestSuite viz(filename, {poly});
            
            // Left panel: Earcut
            viz.addPanelTitle("Earcut (Current MapLibre)", false);
            viz.drawTileBoundary(false);
            viz.drawPolygonWithHoles(poly.rings, "#4CAF50", false, 0.3);
            viz.drawTriangulation(earcutCoords, earcutResult, "#FF5722", false);
            viz.addStats("Earcut: " + std::to_string(earcutResult.size()/3) + " triangles", false);
            
            // Right panel: MeshCut
            viz.addPanelTitle("MeshCut (Terrain-Aware)", true);
            viz.drawTileBoundary(true);
            viz.drawGrid(32, true);
            viz.drawPolygonWithHoles(poly.rings, "#4CAF50", true, 0.3);  // MeshCut now supports holes!
            
            // Convert MeshCut result for visualization
            std::vector<std::pair<double, double>> meshVertices;
            for (size_t j = 0; j < meshcutResult.vertices.size(); j += 2) {
                meshVertices.emplace_back(meshcutResult.vertices[j], meshcutResult.vertices[j+1]);
            }
            viz.drawMeshCutTriangulation(meshVertices, meshcutResult.indices, "#2196F3", true);
            viz.addStats("MeshCut: " + std::to_string(meshcutResult.indices.size()/3) + " triangles, " + 
                        std::to_string(meshVertices.size()) + " vertices", true);
            
            std::cout << "  Individual comparison saved: " << filename << "\n\n";
            successCount++;
            
        } catch (const std::exception& e) {
            std::cout << "  ERROR processing polygon " << poly.id << ": " << e.what() << "\n";
            std::cout << "  Skipping this polygon and continuing...\n\n";
            errorCount++;
        }
    }
    
    std::cout << "\n=== TEST SUITE COMPLETE ===\n";
    std::cout << "Generated visualizations for visual evaluation:\n";
    std::cout << "- Combined tile analysis (1 file) - showing ALL " << testPolygons.size() << " polygons\n";
    std::cout << "- Individual polygon comparisons (" << successCount << " successful, " << errorCount << " errors)\n";
    std::cout << "- All using the exact same tile source as your MapLibre testing\n";
    std::cout << "\nNow you can visually evaluate the differences without building for Android!\n";
    std::cout << "The combined visualization shows all fill areas as requested.\n";
    
    return 0;
}