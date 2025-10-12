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

// Enhanced MVT parser that analyzes coordinate ranges
class MVTAnalyzer {
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
        uint32_t extent;
        
        void printCoordinateAnalysis() const {
            if (points.empty()) return;
            
            double minX = points[0].x, maxX = points[0].x;
            double minY = points[0].y, maxY = points[0].y;
            
            for (const auto& pt : points) {
                minX = std::min(minX, pt.x);
                maxX = std::max(maxX, pt.x);
                minY = std::min(minY, pt.y);
                maxY = std::max(maxY, pt.y);
            }
            
            std::cout << "  Layer: " << layerName << std::endl;
            std::cout << "  Vertices: " << points.size() << std::endl;
            std::cout << "  Tile extent: " << extent << std::endl;
            std::cout << "  Raw coordinate range: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")" << std::endl;
            
            // Check if coordinates are outside 0-extent range
            bool hasOutOfBounds = minX < 0 || maxX > extent || minY < 0 || maxY > extent;
            std::cout << "  ⚠️  Contains out-of-bounds coordinates: " << (hasOutOfBounds ? "YES" : "NO") << std::endl;
            
            if (hasOutOfBounds) {
                std::cout << "  📍 Out-of-bounds details:" << std::endl;
                if (minX < 0) std::cout << "    - Min X below 0: " << minX << std::endl;
                if (maxX > extent) std::cout << "    - Max X above extent: " << maxX << " (extent=" << extent << ")" << std::endl;
                if (minY < 0) std::cout << "    - Min Y below 0: " << minY << std::endl;
                if (maxY > extent) std::cout << "    - Max Y above extent: " << maxY << " (extent=" << extent << ")" << std::endl;
            }
            
            std::cout << "  📊 Sample coordinates:" << std::endl;
            size_t sampleCount = std::min(points.size(), size_t(5));
            for (size_t i = 0; i < sampleCount; i++) {
                std::cout << "    [" << i << "] (" << points[i].x << "," << points[i].y << ")" << std::endl;
            }
            std::cout << std::endl;
        }
    };
    
    MVTAnalyzer(const std::vector<uint8_t>& tileData) 
        : data(tileData.data()), size(tileData.size()), pos(0) {}
    
    std::vector<Polygon> analyzePolygons() {
        std::vector<Polygon> polygons;
        pos = 0;
        
        std::cout << "Analyzing MVT coordinate ranges (" << size << " bytes)...\n\n";
        
        while (pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 2) {
                uint64_t layerLength = readVarint();
                size_t layerEnd = std::min(pos + layerLength, size);
                
                auto layerPolygons = parseLayer(layerEnd);
                for (auto& poly : layerPolygons) {
                    polygons.push_back(std::move(poly));
                }
                
                if (polygons.size() >= 10) {
                    std::cout << "Found " << polygons.size() << " polygons, stopping for analysis\n\n";
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
        uint32_t extent = 4096;  // Default MVT extent
        
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
        polygon.extent = extent;
        
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
                    
                    // Store RAW coordinates (not normalized)
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
                    
                    // Store RAW coordinates (not normalized)
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

int main() {
    std::cout << "=== REAL TILE COORDINATE ANALYSIS ===\n";
    std::cout << "Analyzing actual vector tile data to understand coordinate ranges\n\n";
    
    // Test multiple real tiles to see coordinate patterns
    std::vector<std::pair<std::string, std::string>> testTiles = {
        {"https://demotiles.maplibre.org/tiles-omt/14/8710/5744.pbf", "Innsbruck (z14)"},
        {"https://demotiles.maplibre.org/tiles-omt/10/544/362.pbf", "Europe (z10)"},
        {"https://demotiles.maplibre.org/tiles-omt/12/2044/1361.pbf", "Germany (z12)"}
    };
    
    for (const auto& [url, description] : testTiles) {
        std::cout << "=== " << description << " ===\n";
        std::cout << "URL: " << url << "\n\n";
        
        std::vector<uint8_t> tileData;
        if (!downloadTile(url, tileData)) {
            std::cout << "❌ Failed to download tile\n\n";
            continue;
        }
        
        std::cout << "Downloaded " << tileData.size() << " bytes\n";
        
        // Try decompression if needed
        std::vector<uint8_t> decompressed = decompressGzip(tileData);
        if (!decompressed.empty()) {
            std::cout << "Decompressed to " << decompressed.size() << " bytes\n";
            tileData = std::move(decompressed);
        }
        
        // Analyze coordinate ranges
        MVTAnalyzer analyzer(tileData);
        auto polygons = analyzer.analyzePolygons();
        
        std::cout << "=== COORDINATE ANALYSIS RESULTS ===\n";
        for (const auto& poly : polygons) {
            poly.printCoordinateAnalysis();
        }
        
        std::cout << "================================\n\n";
    }
    
    std::cout << "=== CONCLUSION ===\n";
    std::cout << "This analysis shows the actual coordinate ranges in real MVT tiles.\n";
    std::cout << "Key questions answered:\n";
    std::cout << "1. Do polygons actually cross tile boundaries?\n";
    std::cout << "2. What are the real coordinate ranges?\n";
    std::cout << "3. Are coordinates normalized or in tile units?\n";
    std::cout << "4. Is the 0-4096 assumption correct?\n";
    
    return 0;
}