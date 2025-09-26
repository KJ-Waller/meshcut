#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>
#include <cmath>
#include <cassert>
#include <string>
#include <map>
#include <curl/curl.h>
#include <zlib.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Callback function for curl to write data
struct CurlResponse {
    std::vector<uint8_t> data;
};

size_t WriteCallback(void* contents, size_t size, size_t nmemb, CurlResponse* response) {
    size_t totalSize = size * nmemb;
    uint8_t* bytes = (uint8_t*)contents;
    response->data.insert(response->data.end(), bytes, bytes + totalSize);
    return totalSize;
}

// Simple Protocol Buffer decoder for MVT (Mapbox Vector Tiles)
class MVTDecoder {
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
        }
        return result;
    }
    
    std::string readString(size_t length) {
        if (pos + length > size) return "";
        std::string result((char*)(data + pos), length);
        pos += length;
        return result;
    }
    
    std::vector<uint8_t> readBytes(size_t length) {
        if (pos + length > size) return {};
        std::vector<uint8_t> result(data + pos, data + pos + length);
        pos += length;
        return result;
    }
    
public:
    struct Feature {
        std::string type;
        std::vector<std::vector<std::pair<double, double>>> polygons;
        std::map<std::string, std::string> properties;
    };
    
    struct Layer {
        std::string name;
        std::vector<Feature> features;
        uint32_t extent = 4096;
    };
    
    MVTDecoder(const std::vector<uint8_t>& tileData) {
        data = tileData.data();
        size = tileData.size();
        pos = 0;
    }
    
    std::vector<Layer> decode() {
        std::vector<Layer> layers;
        int iterations = 0;
        const int MAX_ITERATIONS = 10000; // Safety limit
        
        while (pos < size && iterations < MAX_ITERATIONS) {
            iterations++;
            
            if (iterations % 1000 == 0) {
                std::cout << "  Decode iteration " << iterations << ", pos=" << pos << "/" << size << "\n";
            }
            
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 2) { // Layer field
                uint64_t length = readVarint();
                size_t layerEnd = pos + length;
                
                std::cout << "  Decoding layer at pos " << pos << ", length=" << length << "\n";
                Layer layer = decodeLayer(layerEnd);
                if (!layer.name.empty()) {
                    layers.push_back(layer);
                    std::cout << "  Found layer: " << layer.name << " with " << layer.features.size() << " features\n";
                }
            } else {
                skipField(wireType);
            }
        }
        
        if (iterations >= MAX_ITERATIONS) {
            std::cout << "WARNING: Hit iteration limit during decode\n";
        }
        
        return layers;
    }
    
private:
    void skipField(uint8_t wireType) {
        switch (wireType) {
            case 0: readVarint(); break; // Varint
            case 1: pos += 8; break; // 64-bit
            case 2: pos += readVarint(); break; // Length-delimited
            case 5: pos += 4; break; // 32-bit
        }
    }
    
    Layer decodeLayer(size_t layerEnd) {
        Layer layer;
        std::vector<std::string> keys;
        std::vector<std::string> values;
        
        int layerIterations = 0;
        const int MAX_LAYER_ITERATIONS = 5000;
        
        while (pos < layerEnd && layerIterations < MAX_LAYER_ITERATIONS) {
            layerIterations++;
            
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            switch (fieldNumber) {
                case 1: { // version
                    readVarint();
                    break;
                }
                case 2: { // name
                    uint64_t length = readVarint();
                    layer.name = readString(length);
                    break;
                }
                case 3: { // features
                    uint64_t length = readVarint();
                    size_t featureEnd = pos + length;
                    Feature feature = decodeFeature(featureEnd, keys, values);
                    if (!feature.polygons.empty()) {
                        layer.features.push_back(feature);
                    }
                    break;
                }
                case 4: { // keys
                    uint64_t length = readVarint();
                    keys.push_back(readString(length));
                    break;
                }
                case 5: { // values
                    uint64_t length = readVarint();
                    values.push_back(decodeValue(pos + length));
                    break;
                }
                case 15: { // extent
                    layer.extent = readVarint();
                    break;
                }
                default:
                    skipField(wireType);
                    break;
            }
        }
        
        return layer;
    }
    
    std::string decodeValue(size_t valueEnd) {
        while (pos < valueEnd) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 1 && wireType == 2) { // string value
                uint64_t length = readVarint();
                return readString(length);
            } else {
                skipField(wireType);
            }
        }
        return "";
    }
    
    Feature decodeFeature(size_t featureEnd, const std::vector<std::string>& keys, const std::vector<std::string>& values) {
        Feature feature;
        std::vector<uint32_t> tags;
        
        while (pos < featureEnd) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            switch (fieldNumber) {
                case 1: { // id
                    readVarint();
                    break;
                }
                case 2: { // tags
                    uint64_t length = readVarint();
                    size_t tagsEnd = pos + length;
                    while (pos < tagsEnd) {
                        tags.push_back(readVarint());
                    }
                    break;
                }
                case 3: { // type (1=Point, 2=LineString, 3=Polygon)
                    uint32_t geomType = readVarint();
                    if (geomType == 3) {
                        feature.type = "Polygon";
                    }
                    break;
                }
                case 4: { // geometry
                    uint64_t length = readVarint();
                    size_t geomEnd = pos + length;
                    if (feature.type == "Polygon") {
                        feature.polygons = decodePolygonGeometry(geomEnd);
                    } else {
                        pos = geomEnd; // Skip non-polygon geometry
                    }
                    break;
                }
                default:
                    skipField(wireType);
                    break;
            }
        }
        
        // Decode tags into properties
        for (size_t i = 0; i < tags.size(); i += 2) {
            if (i + 1 < tags.size() && tags[i] < keys.size() && tags[i + 1] < values.size()) {
                feature.properties[keys[tags[i]]] = values[tags[i + 1]];
            }
        }
        
        return feature;
    }
    
    std::vector<std::vector<std::pair<double, double>>> decodePolygonGeometry(size_t geomEnd) {
        std::vector<std::vector<std::pair<double, double>>> polygons;
        std::vector<std::pair<double, double>> currentRing;
        
        int32_t x = 0, y = 0;
        
        while (pos < geomEnd) {
            uint32_t command = readVarint();
            uint32_t cmdId = command & 0x7;
            uint32_t cmdCount = command >> 3;
            
            if (cmdId == 1) { // MoveTo
                if (!currentRing.empty()) {
                    polygons.push_back(currentRing);
                    currentRing.clear();
                }
                for (uint32_t i = 0; i < cmdCount; i++) {
                    int32_t dx = zigzagDecode(readVarint());
                    int32_t dy = zigzagDecode(readVarint());
                    x += dx;
                    y += dy;
                    currentRing.push_back({(double)x, (double)y});
                }
            } else if (cmdId == 2) { // LineTo
                for (uint32_t i = 0; i < cmdCount; i++) {
                    int32_t dx = zigzagDecode(readVarint());
                    int32_t dy = zigzagDecode(readVarint());
                    x += dx;
                    y += dy;
                    currentRing.push_back({(double)x, (double)y});
                }
            } else if (cmdId == 7) { // ClosePath
                if (!currentRing.empty() && currentRing.size() > 2) {
                    // Close the ring by adding the first point again if needed
                    if (currentRing.front() != currentRing.back()) {
                        currentRing.push_back(currentRing.front());
                    }
                }
            }
        }
        
        if (!currentRing.empty()) {
            polygons.push_back(currentRing);
        }
        
        return polygons;
    }
    
    int32_t zigzagDecode(uint32_t n) {
        return (n >> 1) ^ (-(n & 1));
    }
};

// Coordinate conversion utilities
namespace CoordUtils {
    // Convert tile coordinates to Web Mercator
    std::pair<double, double> tileToWebMercator(double tileX, double tileY, int z, int x, int y, uint32_t extent) {
        double n = pow(2.0, z);
        double lon_min = x / n * 360.0 - 180.0;
        double lon_max = (x + 1) / n * 360.0 - 180.0;
        double lat_min = atan(sinh(M_PI * (1 - 2 * (y + 1) / n))) * 180.0 / M_PI;
        double lat_max = atan(sinh(M_PI * (1 - 2 * y / n))) * 180.0 / M_PI;
        
        // Convert lat/lon to Web Mercator
        auto latLonToWebMercator = [](double lat, double lon) -> std::pair<double, double> {
            double x = lon * 20037508.34 / 180.0;
            double y = log(tan((90.0 + lat) * M_PI / 360.0)) / (M_PI / 180.0);
            y = y * 20037508.34 / 180.0;
            return {x, y};
        };
        
        auto [minX_merc, minY_merc] = latLonToWebMercator(lat_min, lon_min);
        auto [maxX_merc, maxY_merc] = latLonToWebMercator(lat_max, lon_max);
        
        // Scale tile coordinates to Web Mercator
        double webMercX = minX_merc + (tileX / extent) * (maxX_merc - minX_merc);
        double webMercY = maxY_merc - (tileY / extent) * (maxY_merc - minY_merc); // Flip Y
        
        return {webMercX, webMercY};
    }
    
    // Convert Web Mercator back to lat/lon for verification
    std::pair<double, double> webMercatorToLatLon(double x, double y) {
        double lon = x / 20037508.34 * 180.0;
        double lat = (2.0 * atan(exp(y / 20037508.34 * M_PI)) - M_PI / 2.0) * 180.0 / M_PI;
        return {lat, lon};
    }
}

bool downloadTile(int z, int x, int y, std::vector<uint8_t>& tileData) {
    CURL* curl;
    CURLcode res;
    CurlResponse response;
    
    curl = curl_easy_init();
    if (!curl) {
        std::cerr << "Failed to initialize curl\n";
        return false;
    }
    
    std::string url = "https://demotiles.maplibre.org/tiles-omt/" + 
                      std::to_string(z) + "/" + 
                      std::to_string(x) + "/" + 
                      std::to_string(y) + ".pbf";
    
    std::cout << "Downloading tile from: " << url << std::endl;
    
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 30L);
    
    res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    
    if (res != CURLE_OK) {
        std::cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << "\n";
        return false;
    }
    
    if (response.data.empty()) {
        std::cerr << "Downloaded tile is empty\n";
        return false;
    }
    
    std::cout << "Downloaded " << response.data.size() << " bytes\n";
    
    // Copy data to output vector
    tileData.assign(response.data.begin(), response.data.end());
    
    return true;
}

// SVG visualization class (simplified for real data)
class TileVisualization {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double scale;
    double panelWidth;
    
public:
    TileVisualization(const std::string& filename, double minX_, double minY_, double maxX_, double maxY_) 
        : minX(minX_), minY(minY_), maxX(maxX_), maxY(maxY_) {
        
        double width = maxX - minX;
        double height = maxY - minY;
        scale = 400.0 / std::max(width, height);
        panelWidth = width * scale + 100;
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" ";
        file << "width=\"" << (panelWidth * 2 + 50) << "\" ";
        file << "height=\"" << (height * scale + 200) << "\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        
        file << "<text x=\"" << (panelWidth/2) << "\" y=\"25\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">Earcut</text>\n";
        file << "<text x=\"" << (panelWidth + 25 + panelWidth/2) << "\" y=\"25\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">MeshCut</text>\n";
    }
    
    ~TileVisualization() {
        file << "</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale + 50;
        double ty = (maxY - y) * scale + 50;
        if (rightPanel) tx += panelWidth + 25;
        return {tx, ty};
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color, bool rightPanel = false) {
        if (polygon.size() < 6) return;
        
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i + 1], rightPanel);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"0.4\" stroke=\"" << color << "\"/>\n";
    }
    
    void drawTriangles(const std::vector<double>& vertices, const std::vector<uint32_t>& indices, 
                       const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.3\">\n";
        
        for (size_t i = 0; i < indices.size(); i += 3) {
            uint32_t i0 = indices[i] * 2;
            uint32_t i1 = indices[i + 1] * 2;
            uint32_t i2 = indices[i + 2] * 2;
            
            if (i0 + 1 >= vertices.size() || i1 + 1 >= vertices.size() || i2 + 1 >= vertices.size()) continue;
            
            auto [x0, y0] = transform(vertices[i0], vertices[i0 + 1], rightPanel);
            auto [x1, y1] = transform(vertices[i1], vertices[i1 + 1], rightPanel);
            auto [x2, y2] = transform(vertices[i2], vertices[i2 + 1], rightPanel);
            
            file << "<polygon points=\"" << x0 << "," << y0 << " " 
                 << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
        }
        
        file << "</g>\n";
    }
};

int main() {
    std::cout << "MeshCut Real Vector Tile Test - ACTUAL DATA\n";
    std::cout << "==========================================\n";
    
    // Download the actual Innsbruck tile
    int z = 14, x = 8710, y = 5744;
    std::vector<uint8_t> tileData;
    
    if (!downloadTile(z, x, y, tileData)) {
        std::cerr << "Failed to download tile\n";
        return 1;
    }
    
    // Decode the MVT data
    std::cout << "Starting MVT decoding...\n";
    MVTDecoder decoder(tileData);
    auto layers = decoder.decode();
    std::cout << "MVT decoding completed!\n";
    
    std::cout << "\nFound " << layers.size() << " layers in tile:\n";
    for (const auto& layer : layers) {
        std::cout << "- Layer '" << layer.name << "': " << layer.features.size() << " polygon features\n";
    }
    
    // Find fill layers with polygons
    std::vector<MVTDecoder::Feature> fillPolygons;
    for (const auto& layer : layers) {
        if (layer.name.find("fill") != std::string::npos || 
            layer.name == "building" || layer.name == "water" || 
            layer.name == "landuse" || layer.name == "landcover") {
            
            for (const auto& feature : layer.features) {
                if (feature.type == "Polygon" && !feature.polygons.empty()) {
                    fillPolygons.push_back(feature);
                }
            }
        }
    }
    
    std::cout << "\nFound " << fillPolygons.size() << " fill polygons total\n";
    
    if (fillPolygons.empty()) {
        std::cout << "No fill polygons found in tile\n";
        return 1;
    }
    
    // Just test the first polygon for demo
    if (fillPolygons.empty()) {
        std::cout << "No fill polygons found, checking first layer with features...\n";
        for (const auto& layer : layers) {
            if (!layer.features.empty()) {
                fillPolygons.push_back(layer.features[0]);
                std::cout << "Using first polygon from layer: " << layer.name << "\n";
                break;
            }
        }
    }
    
    if (fillPolygons.empty()) {
        std::cout << "No polygon features found in any layer\n";
        return 1;
    }
    
    // Test just the first polygon
    int tested = 0;
    double tileMinX = 1267020, tileMinY = 5985325, tileMaxX = 1269466, tileMaxY = 5987771;
    
    for (const auto& feature : fillPolygons) {
        if (tested >= 5) break; // Limit to first 5 polygons
        
        // Convert first ring to Web Mercator coordinates
        const auto& ring = feature.polygons[0];
        std::vector<double> polygon;
        
        for (const auto& point : ring) {
            auto [webX, webY] = CoordUtils::tileToWebMercator(point.first, point.second, z, x, y, 4096);
            polygon.push_back(webX);
            polygon.push_back(webY);
        }
        
        if (polygon.size() < 6) continue;
        
        std::cout << "\n=== Testing Polygon " << (tested + 1) << " ===\n";
        std::cout << "Vertices: " << polygon.size() / 2 << "\n";
        
        // Get polygon center for verification
        double centerX = 0, centerY = 0;
        for (size_t i = 0; i < polygon.size(); i += 2) {
            centerX += polygon[i];
            centerY += polygon[i + 1];
        }
        centerX /= (polygon.size() / 2);
        centerY /= (polygon.size() / 2);
        
        auto [lat, lon] = CoordUtils::webMercatorToLatLon(centerX, centerY);
        std::cout << "Center: " << std::fixed << std::setprecision(5) << lat << "°N, " << lon << "°E\n";
        std::cout << "Google Maps: https://www.google.com/maps/@" << lat << "," << lon << ",18z\n";
        
        // Print properties if available
        for (const auto& prop : feature.properties) {
            std::cout << "Property " << prop.first << ": " << prop.second << "\n";
        }
        
        // Test MeshCut vs Earcut
        meshcut::MeshCutOptions options;
        options.gridOriginX = tileMinX;
        options.gridOriginY = tileMinY;
        options.cellSize = (tileMaxX - tileMinX) / 32.0;
        options.gridWidth = 32;
        options.gridHeight = 32;
        
        auto meshcutResult = meshcut::meshcut_full(polygon, {}, options);
        
        // Earcut comparison
        using Point = std::array<double, 2>;
        std::vector<std::vector<Point>> earcutPolygon;
        earcutPolygon.push_back({});
        for (size_t i = 0; i < polygon.size(); i += 2) {
            earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
        }
        auto earcutIndices = mapbox::earcut<uint32_t>(earcutPolygon);
        
        std::cout << "Earcut: " << earcutIndices.size() / 3 << " triangles\n";
        std::cout << "MeshCut: " << meshcutResult.indices.size() / 3 << " triangles\n";
        std::cout << "Density improvement: " << std::fixed << std::setprecision(1) 
                  << ((double)meshcutResult.indices.size() / earcutIndices.size()) << "x\n";
        
        // Create visualization
        std::string filename = "../visualizations/real_tile_polygon_" + std::to_string(tested + 1) + ".svg";
        TileVisualization viz(filename, tileMinX, tileMinY, tileMaxX, tileMaxY);
        
        viz.drawPolygon(polygon, "#4CAF50", false);
        viz.drawTriangles(polygon, earcutIndices, "#FF5722", false);
        
        viz.drawPolygon(polygon, "#4CAF50", true);
        viz.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", true);
        
        std::cout << "Visualization saved to: " << filename << "\n";
        
        tested++;
    }
    
    std::cout << "\n=== REAL VECTOR TILE TEST COMPLETE ===\n";
    std::cout << "Tested " << tested << " real polygons from Innsbruck tile\n";
    std::cout << "All visualizations saved to ../visualizations/real_tile_polygon_*.svg\n";
    
    return 0;
}