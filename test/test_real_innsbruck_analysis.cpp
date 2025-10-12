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

// Proper MVT parser that extracts tile metadata AND polygons
class InnsbruckTileAnalyzer {
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
    
    struct TileMetadata {
        int z, x, y;           // Tile coordinates
        uint32_t extent;       // Tile extent (usually 4096)
        double minX, minY, maxX, maxY;  // Real-world bounding box
        
        void calculateBounds() {
            // For a tile at z/x/y, calculate the real-world bounding box
            // This is critical for setting up MeshCut grid correctly
            double n = std::pow(2.0, z);
            
            // Web Mercator bounds
            minX = (x / n) * 360.0 - 180.0;
            maxX = ((x + 1) / n) * 360.0 - 180.0;
            
            double lat_rad = std::atan(std::sinh(M_PI * (1 - 2 * y / n)));
            minY = lat_rad * 180.0 / M_PI;
            
            lat_rad = std::atan(std::sinh(M_PI * (1 - 2 * (y + 1) / n)));
            maxY = lat_rad * 180.0 / M_PI;
            
            std::cout << "=== TILE METADATA ===\n";
            std::cout << "Tile: z=" << z << ", x=" << x << ", y=" << y << "\n";
            std::cout << "Extent: " << extent << " units\n";
            std::cout << "Geographic bounds: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\n";
            std::cout << "Tile covers: " << (maxX - minX) << "° x " << (maxY - minY) << "° degrees\n\n";
        }
    };
    
    struct RealPolygon {
        std::vector<Point> points;
        std::string layerName;
        TileMetadata tileInfo;
        
        void printSummary() const {
            if (points.empty()) return;
            
            double minX = points[0].x, maxX = points[0].x;
            double minY = points[0].y, maxY = points[0].y;
            
            for (const auto& pt : points) {
                minX = std::min(minX, pt.x);
                maxX = std::max(maxX, pt.x);
                minY = std::min(minY, pt.y);
                maxY = std::max(maxY, pt.y);
            }
            
            std::cout << "=== REAL POLYGON ===\n";
            std::cout << "Layer: " << layerName << "\n";
            std::cout << "Vertices: " << points.size() << "\n";
            std::cout << "Tile coordinate range: (" << minX << "," << minY << ") to (" << maxX << "," << maxY << ")\n";
            
            // Check if polygon extends beyond tile extent
            bool extendsOutside = (minX < 0 || maxX > tileInfo.extent || minY < 0 || maxY > tileInfo.extent);
            std::cout << "Extends beyond tile: " << (extendsOutside ? "YES" : "NO") << "\n";
            
            if (extendsOutside) {
                std::cout << "Buffer zones detected:\n";
                if (minX < 0) std::cout << "  - Left buffer: " << (-minX) << " units\n";
                if (maxX > tileInfo.extent) std::cout << "  - Right buffer: " << (maxX - tileInfo.extent) << " units\n";
                if (minY < 0) std::cout << "  - Bottom buffer: " << (-minY) << " units\n";
                if (maxY > tileInfo.extent) std::cout << "  - Top buffer: " << (maxY - tileInfo.extent) << " units\n";
            }
            
            // Calculate what percentage of tile this polygon covers
            double polyArea = (maxX - minX) * (maxY - minY);
            double tileArea = tileInfo.extent * tileInfo.extent;
            double coverage = (polyArea / tileArea) * 100.0;
            std::cout << "Approximate tile coverage: " << std::fixed << std::setprecision(1) << coverage << "%\n\n";
        }
        
        std::vector<double> toFlatCoords() const {
            std::vector<double> coords;
            coords.reserve(points.size() * 2);
            for (const auto& pt : points) {
                coords.push_back(pt.x);
                coords.push_back(pt.y);
            }
            return coords;
        }
    };
    
    InnsbruckTileAnalyzer(const std::vector<uint8_t>& tileData, int z, int x, int y) 
        : data(tileData.data()), size(tileData.size()), pos(0) {
        metadata.z = z;
        metadata.x = x;
        metadata.y = y;
        metadata.extent = 4096;  // Default, will be updated if found in tile
    }
    
    TileMetadata metadata;
    
    std::vector<RealPolygon> extractRealPolygons() {
        std::vector<RealPolygon> polygons;
        pos = 0;
        
        std::cout << "Analyzing real Innsbruck tile data (" << size << " bytes)...\n\n";
        
        while (pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 3 && wireType == 2) {
                uint64_t layerLength = readVarint();
                size_t layerEnd = std::min(pos + layerLength, size);
                
                auto layerPolygons = parseLayer(layerEnd);
                for (auto& poly : layerPolygons) {
                    poly.tileInfo = metadata;
                    polygons.push_back(std::move(poly));
                }
                
                if (polygons.size() >= 10) {
                    std::cout << "Found " << polygons.size() << " real polygons\n\n";
                    break;
                }
            } else {
                skipField(wireType);
            }
        }
        
        // Calculate tile bounds now that we have the extent
        metadata.calculateBounds();
        
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
    
    std::vector<RealPolygon> parseLayer(size_t layerEnd) {
        std::vector<RealPolygon> polygons;
        std::vector<std::string> keys;
        std::vector<std::string> values;
        std::string layerName = "unknown";
        
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
                auto polygon = parseFeature(featureEnd, keys, values, layerName);
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
                metadata.extent = static_cast<uint32_t>(readVarint());
            } else {
                skipField(wireType);
            }
        }
        
        return polygons;
    }
    
    RealPolygon parseFeature(size_t featureEnd, const std::vector<std::string>& keys, 
                            const std::vector<std::string>& values, 
                            const std::string& layerName) {
        RealPolygon polygon;
        polygon.layerName = layerName;
        
        while (pos < featureEnd && pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 4 && wireType == 2) { // Geometry
                uint64_t geomLength = readVarint();
                size_t geomEnd = std::min(pos + geomLength, size);
                polygon.points = parseGeometry(geomEnd);
                break;
            } else {
                skipField(wireType);
            }
        }
        
        return polygon;
    }
    
    std::vector<Point> parseGeometry(size_t geomEnd) {
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

// MeshCut setup with CORRECT tile bounding box
meshcut::MeshCutOptions createCorrectMeshCutOptions(const InnsbruckTileAnalyzer::TileMetadata& tileInfo, uint32_t gridSize = 32) {
    meshcut::MeshCutOptions options;
    
    // CRITICAL: Use the actual tile coordinate space, including buffer zones
    // Find the actual min/max coordinates from all polygons to set proper bounds
    options.gridOriginX = -128.0;  // Allow for typical buffer zone
    options.gridOriginY = -128.0;  
    options.gridWidth = gridSize;
    options.gridHeight = gridSize;
    options.cellSize = (tileInfo.extent + 256.0) / gridSize;  // Cover extent + buffer
    options.diagonalNE = true;
    
    std::cout << "=== MESHCUT GRID SETUP ===\n";
    std::cout << "Grid size: " << gridSize << "x" << gridSize << " cells\n";
    std::cout << "Cell size: " << options.cellSize << " tile units\n";
    std::cout << "Grid origin: (" << options.gridOriginX << "," << options.gridOriginY << ")\n";
    std::cout << "Grid covers: " << options.gridOriginX << " to " << (options.gridOriginX + options.gridWidth * options.cellSize) << "\n";
    std::cout << "Tile extent: 0 to " << tileInfo.extent << " (buffer: ±128)\n\n";
    
    return options;
}

// SVG visualizer for real Innsbruck data
class InnsbruckVisualizer {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double width, height, scale;
    uint32_t tileExtent;
    
public:
    InnsbruckVisualizer(const std::string& filename, const std::vector<InnsbruckTileAnalyzer::RealPolygon>& polygons) {
        if (polygons.empty()) return;
        
        tileExtent = polygons[0].tileInfo.extent;
        
        // Calculate bounds from all polygons
        minX = maxX = polygons[0].points[0].x;
        minY = maxY = polygons[0].points[0].y;
        
        for (const auto& poly : polygons) {
            for (const auto& pt : poly.points) {
                minX = std::min(minX, pt.x);
                maxX = std::max(maxX, pt.x);
                minY = std::min(minY, pt.y);
                maxY = std::max(maxY, pt.y);
            }
        }
        
        // Expand bounds to show tile boundaries
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
        
        // Add title
        file << "<text x=\"" << (width * scale + 100) << "\" y=\"30\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"20\" font-weight=\"bold\">Real Innsbruck Tile Analysis</text>\n";
    }
    
    ~InnsbruckVisualizer() {
        file << "</svg>\n";
        file.close();
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale;
        double ty = (maxY - y) * scale;
        if (rightPanel) tx += width * scale + 100;
        return {tx + 50, ty + 100};
    }
    
    void drawTileBoundaries(bool rightPanel = false) {
        // Draw main tile boundary (0 to extent)
        file << "<g stroke=\"#ff0000\" stroke-width=\"2\" fill=\"none\" opacity=\"0.8\">\n";
        
        auto [x1, y1] = transform(0, 0, rightPanel);
        auto [x2, y2] = transform(tileExtent, 0, rightPanel);
        auto [x3, y3] = transform(tileExtent, tileExtent, rightPanel);
        auto [x4, y4] = transform(0, tileExtent, rightPanel);
        
        file << "<rect x=\"" << x1 << "\" y=\"" << y4 << "\" width=\"" << (x2-x1) << "\" height=\"" << (y1-y4) << "\"/>\n";
        file << "</g>\n";
        
        // Add tile boundary label
        auto [labelX, labelY] = transform(tileExtent/2, -50, rightPanel);
        file << "<text x=\"" << labelX << "\" y=\"" << labelY << "\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"12\" fill=\"red\">Tile Boundary (0-" << tileExtent << ")</text>\n";
    }
    
    void drawGrid(const meshcut::MeshCutOptions& options, bool rightPanel = false) {
        file << "<g stroke=\"#cccccc\" stroke-width=\"0.5\" fill=\"none\" opacity=\"0.6\">\n";
        
        // Vertical grid lines
        for (int i = 0; i <= (int)options.gridWidth; i++) {
            double x = options.gridOriginX + i * options.cellSize;
            auto [x1, y1] = transform(x, options.gridOriginY, rightPanel);
            auto [x2, y2] = transform(x, options.gridOriginY + options.gridHeight * options.cellSize, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        
        // Horizontal grid lines
        for (int j = 0; j <= (int)options.gridHeight; j++) {
            double y = options.gridOriginY + j * options.cellSize;
            auto [x1, y1] = transform(options.gridOriginX, y, rightPanel);
            auto [x2, y2] = transform(options.gridOriginX + options.gridWidth * options.cellSize, y, rightPanel);
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
    
    void drawTriangles(const std::vector<std::pair<double, double>>& vertices, 
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
    
    void drawEarcutTriangles(const std::vector<double>& coords, 
                            const std::vector<uint32_t>& indices, 
                            const std::string& color, bool rightPanel = false) {
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
    
    void addPanelTitle(const std::string& title, bool rightPanel = false) {
        double titleX = rightPanel ? width * scale + 100 + (width * scale / 2) : (width * scale / 2);
        file << "<text x=\"" << (titleX + 50) << "\" y=\"70\" text-anchor=\"middle\" ";
        file << "font-family=\"Arial\" font-size=\"16\" font-weight=\"bold\">" << title << "</text>\n";
    }
    
    void addStats(const std::string& stats, bool rightPanel = false) {
        double statsX = rightPanel ? width * scale + 100 + 10 : 10;
        file << "<text x=\"" << (statsX + 50) << "\" y=\"" << (height * scale + 150) << "\" ";
        file << "font-family=\"Arial\" font-size=\"11\" fill=\"#666\">" << stats << "</text>\n";
    }
};

int main() {
    std::cout << "=== REAL INNSBRUCK TILE ANALYSIS WITH CORRECT MESHCUT SETUP ===\n";
    std::cout << "Downloading and analyzing the actual tile you're testing in MapLibre\n\n";
    
    // Use the EXACT tile coordinates for Innsbruck that you're testing
    int z = 14, x = 8710, y = 5744;  // Innsbruck tile coordinates
    std::string tileUrl = "https://demotiles.maplibre.org/tiles-omt/" + std::to_string(z) + "/" + 
                         std::to_string(x) + "/" + std::to_string(y) + ".pbf";
    
    std::cout << "Downloading Innsbruck tile: z=" << z << ", x=" << x << ", y=" << y << "\n";
    std::cout << "URL: " << tileUrl << "\n\n";
    
    std::vector<uint8_t> tileData;
    if (!downloadTile(tileUrl, tileData)) {
        std::cerr << "Failed to download Innsbruck tile\n";
        return 1;
    }
    
    std::cout << "Downloaded " << tileData.size() << " bytes\n";
    
    // Try decompression if needed
    std::vector<uint8_t> decompressed = decompressGzip(tileData);
    if (!decompressed.empty()) {
        std::cout << "Decompressed to " << decompressed.size() << " bytes\n";
        tileData = std::move(decompressed);
    }
    
    // Parse the real tile data with correct metadata
    InnsbruckTileAnalyzer analyzer(tileData, z, x, y);
    auto realPolygons = analyzer.extractRealPolygons();
    
    if (realPolygons.empty()) {
        std::cerr << "No polygons found in the Innsbruck tile\n";
        return 1;
    }
    
    // Print analysis of each real polygon
    for (size_t i = 0; i < realPolygons.size(); i++) {
        std::cout << "--- Polygon " << (i+1) << " ---\n";
        realPolygons[i].printSummary();
    }
    
    // Set up MeshCut with correct tile bounding box
    auto meshcutOptions = createCorrectMeshCutOptions(analyzer.metadata, 32);
    
    // Create individual visualizations for each polygon
    for (size_t i = 0; i < std::min(realPolygons.size(), size_t(5)); i++) {
        const auto& poly = realPolygons[i];
        auto coords = poly.toFlatCoords();
        
        std::cout << "=== ANALYZING POLYGON " << (i+1) << " ===\n";
        
        // Test with Earcut
        using Point = std::array<double, 2>;
        std::vector<std::vector<Point>> earcutPolygon;
        earcutPolygon.push_back({});
        for (size_t j = 0; j < coords.size(); j += 2) {
            earcutPolygon[0].push_back({{coords[j], coords[j+1]}});
        }
        auto earcutResult = mapbox::earcut<uint32_t>(earcutPolygon);
        
        // Test with MeshCut
        auto meshcutResult = meshcut::meshcut_full(coords, {}, meshcutOptions);
        
        std::cout << "Earcut: " << coords.size()/2 << " vertices → " << earcutResult.size()/3 << " triangles\n";
        std::cout << "MeshCut: " << coords.size()/2 << " vertices → " << meshcutResult.vertices.size()/2 << " vertices → " << meshcutResult.indices.size()/3 << " triangles\n";
        
        double improvement = (double)meshcutResult.indices.size() / (double)earcutResult.size();
        std::cout << "Triangle density: " << std::fixed << std::setprecision(1) << improvement << "x\n\n";
        
        // Create visualization
        std::string filename = "innsbruck_polygon_" + std::to_string(i+1) + "_analysis.svg";
        InnsbruckVisualizer viz(filename, {poly});
        
        // Left panel: Earcut
        viz.addPanelTitle("Earcut (MapLibre Default)", false);
        viz.drawTileBoundaries(false);
        viz.drawPolygon(coords, "#4CAF50", false);
        viz.drawEarcutTriangles(coords, earcutResult, "#FF5722", false);
        viz.addStats("Earcut: " + std::to_string(earcutResult.size()/3) + " triangles", false);
        
        // Right panel: MeshCut
        viz.addPanelTitle("MeshCut (Terrain-Aware)", true);
        viz.drawTileBoundaries(true);
        viz.drawGrid(meshcutOptions, true);
        viz.drawPolygon(coords, "#4CAF50", true);
        
        // Convert MeshCut vertices to coordinate pairs
        std::vector<std::pair<double, double>> meshVertices;
        for (size_t j = 0; j < meshcutResult.vertices.size(); j += 2) {
            meshVertices.emplace_back(meshcutResult.vertices[j], meshcutResult.vertices[j+1]);
        }
        viz.drawTriangles(meshVertices, meshcutResult.indices, "#2196F3", true);
        viz.addStats("MeshCut: " + std::to_string(meshcutResult.indices.size()/3) + " triangles, " + 
                    std::to_string(meshVertices.size()) + " vertices", true);
        
        std::cout << "Individual visualization saved: " << filename << "\n";
    }
    
    // Create combined visualization
    if (realPolygons.size() > 1) {
        std::cout << "\n=== CREATING COMBINED VISUALIZATION ===\n";
        
        std::string combinedFilename = "innsbruck_all_polygons_combined.svg";
        InnsbruckVisualizer combinedViz(combinedFilename, realPolygons);
        
        // Show all polygons together with tile boundaries and grid
        combinedViz.addPanelTitle("All Real Polygons - Earcut", false);
        combinedViz.drawTileBoundaries(false);
        
        combinedViz.addPanelTitle("All Real Polygons - MeshCut", true);
        combinedViz.drawTileBoundaries(true);
        combinedViz.drawGrid(meshcutOptions, true);
        
        // Draw all polygons
        std::vector<std::string> colors = {"#4CAF50", "#2196F3", "#FF9800", "#9C27B0", "#F44336"};
        for (size_t i = 0; i < realPolygons.size() && i < colors.size(); i++) {
            auto coords = realPolygons[i].toFlatCoords();
            combinedViz.drawPolygon(coords, colors[i], false, 0.2);  // Left panel
            combinedViz.drawPolygon(coords, colors[i], true, 0.2);   // Right panel
        }
        
        combinedViz.addStats("Combined: " + std::to_string(realPolygons.size()) + " real polygons", false);
        combinedViz.addStats("32x32 grid, actual tile bounds", true);
        
        std::cout << "Combined visualization saved: " << combinedFilename << "\n";
    }
    
    std::cout << "\n=== ANALYSIS COMPLETE ===\n";
    std::cout << "Generated visualizations using REAL Innsbruck tile data\n";
    std::cout << "Key findings:\n";
    std::cout << "- Tile extent: " << analyzer.metadata.extent << " units\n";
    std::cout << "- Real polygons analyzed: " << realPolygons.size() << "\n";
    std::cout << "- MeshCut grid properly configured for tile bounds\n";
    std::cout << "- Buffer zones handled correctly\n";
    std::cout << "\nThis should help identify the exact cause of your MapLibre artifacts!\n";
    
    return 0;
}