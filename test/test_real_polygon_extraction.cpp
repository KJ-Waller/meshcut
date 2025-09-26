#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <curl/curl.h>
#include <zlib.h>
#include <cmath>
#include <iomanip>
#include <string>

struct HttpResponse {
    std::vector<uint8_t> data;
};

size_t WriteCallback(void* contents, size_t size, size_t nmemb, HttpResponse* response) {
    size_t totalSize = size * nmemb;
    const uint8_t* bytes = static_cast<const uint8_t*>(contents);
    response->data.insert(response->data.end(), bytes, bytes + totalSize);
    return totalSize;
}

bool downloadTile(int z, int x, int y, std::vector<uint8_t>& tileData) {
    CURL* curl = curl_easy_init();
    if (!curl) return false;
    
    std::string url = "https://demotiles.maplibre.org/tiles-omt/" + 
                      std::to_string(z) + "/" + 
                      std::to_string(x) + "/" + 
                      std::to_string(y) + ".pbf";

    // Print out the URL being downloaded
    std::cout << "Downloading tile from: " << url << std::endl;
    
    HttpResponse response;
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &response);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 10L);
    
    CURLcode res = curl_easy_perform(curl);
    curl_easy_cleanup(curl);
    
    if (res != CURLE_OK || response.data.empty()) {
        return false;
    }
    
    tileData = std::move(response.data);
    return true;
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

// Simple Protocol Buffer parser for MVT (Mapbox Vector Tiles)
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
            if (shift >= 64) break; // Prevent overflow
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
                
                if (polygons.size() >= 5) {
                    std::cout << "Found " << polygons.size() << " polygons, stopping for demo\n";
                    break; // Limit for demo
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
            case 0: readVarint(); break; // Varint
            case 1: skipBytes(8); break; // 64-bit
            case 2: skipBytes(readVarint()); break; // Length-delimited
            case 5: skipBytes(4); break; // 32-bit
            default: break;
        }
    }
    
    std::vector<Polygon> parseLayer(size_t layerEnd) {
        std::vector<Polygon> polygons;
        std::vector<std::string> keys;
        std::vector<std::string> values;
        std::string layerName = "unknown";
        uint32_t extent = 4096; // Default MVT extent
        
        size_t startPos = pos;
        
        while (pos < layerEnd && pos < size) {
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
                    layerName = readString(length);
                    break;
                }
                case 3: { // features
                    uint64_t featureLength = readVarint();
                    size_t featureEnd = std::min(pos + featureLength, size);
                    
                    auto polygon = parseFeature(featureEnd, keys, values, extent);
                    if (!polygon.points.empty()) {
                        polygon.layerName = layerName;
                        polygons.push_back(polygon);
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
                    size_t valueEnd = std::min(pos + length, size);
                    values.push_back(parseValue(valueEnd));
                    break;
                }
                case 15: { // extent
                    extent = readVarint();
                    break;
                }
                default:
                    skipField(wireType);
                    break;
            }
        }
        
        std::cout << "Layer '" << layerName << "': found " << polygons.size() << " polygons\n";
        return polygons;
    }
    
    std::string parseValue(size_t valueEnd) {
        // Simplified - just read string values for now
        while (pos < valueEnd && pos < size) {
            uint64_t tag = readVarint();
            uint8_t wireType = tag & 0x7;
            uint32_t fieldNumber = tag >> 3;
            
            if (fieldNumber == 5 && wireType == 2) { // string_value
                uint64_t length = readVarint();
                return readString(length);
            } else {
                skipField(wireType);
            }
        }
        return "";
    }
    
    Polygon parseFeature(size_t featureEnd, const std::vector<std::string>& keys, 
                        const std::vector<std::string>& values, uint32_t extent) {
        Polygon polygon;
        std::vector<uint32_t> tags;
        uint32_t geomType = 0;
        std::vector<uint32_t> geometry;
        
        while (pos < featureEnd && pos < size) {
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
                    size_t tagsEnd = std::min(pos + length, size);
                    while (pos < tagsEnd && pos < size) {
                        tags.push_back(readVarint());
                    }
                    break;
                }
                case 3: { // type
                    geomType = readVarint();
                    break;
                }
                case 4: { // geometry
                    uint64_t length = readVarint();
                    size_t geomEnd = std::min(pos + length, size);
                    while (pos < geomEnd && pos < size) {
                        geometry.push_back(readVarint());
                    }
                    break;
                }
                default:
                    skipField(wireType);
                    break;
            }
        }
        
        // Only process polygons (type 3)
        if (geomType == 3 && !geometry.empty()) {
            polygon.points = decodePolygonGeometry(geometry, extent);
        }
        
        return polygon;
    }
    
    std::vector<Point> decodePolygonGeometry(const std::vector<uint32_t>& geometry, uint32_t extent) {
        std::vector<Point> points;
        
        if (geometry.empty()) return points;
        
        int32_t x = 0, y = 0;
        size_t i = 0;
        
        while (i < geometry.size()) {
            if (i + 2 >= geometry.size()) break;
            
            uint32_t command = geometry[i++];
            uint32_t cmdId = command & 0x7;
            uint32_t cmdCount = command >> 3;
            
            if (cmdId == 1) { // MoveTo
                for (uint32_t j = 0; j < cmdCount && i + 1 < geometry.size(); j++) {
                    int32_t dx = ((geometry[i] >> 1) ^ (-(geometry[i] & 1)));
                    int32_t dy = ((geometry[i + 1] >> 1) ^ (-(geometry[i + 1] & 1)));
                    i += 2;
                    
                    x += dx;
                    y += dy;
                    
                    // Convert from tile coordinates to normalized coordinates
                    points.push_back({(double)x / extent, (double)y / extent});
                }
            } else if (cmdId == 2) { // LineTo
                for (uint32_t j = 0; j < cmdCount && i + 1 < geometry.size(); j++) {
                    int32_t dx = ((geometry[i] >> 1) ^ (-(geometry[i] & 1)));
                    int32_t dy = ((geometry[i + 1] >> 1) ^ (-(geometry[i + 1] & 1)));
                    i += 2;
                    
                    x += dx;
                    y += dy;
                    
                    points.push_back({(double)x / extent, (double)y / extent});
                }
            } else if (cmdId == 7) { // ClosePath
                // Close the polygon by connecting back to first point
                if (!points.empty()) {
                    points.push_back(points[0]);
                }
                break;
            }
        }
        
        return points;
    }
};

class RealPolygonVisualization {
private:
    std::ofstream file;
    double minX, minY, maxX, maxY;
    double scale;
    
public:
    RealPolygonVisualization(const std::string& filename, const std::vector<double>& polygon) {
        // Calculate bounds from actual polygon
        minX = maxX = polygon[0];
        minY = maxY = polygon[1];
        for (size_t i = 0; i < polygon.size(); i += 2) {
            minX = std::min(minX, polygon[i]);
            maxX = std::max(maxX, polygon[i]);
            minY = std::min(minY, polygon[i+1]);
            maxY = std::max(maxY, polygon[i+1]);
        }
        
        // Add padding
        double padding = std::max(maxX - minX, maxY - minY) * 0.1;
        minX -= padding; maxX += padding;
        minY -= padding; maxY += padding;
        
        scale = 400.0 / std::max(maxX - minX, maxY - minY);
        
        file.open(filename);
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
        file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"500\" viewBox=\"0 0 900 500\">\n";
        file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
        file << "<text x=\"225\" y=\"25\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">Earcut (Sparse Triangulation)</text>\n";
        file << "<text x=\"675\" y=\"25\" text-anchor=\"middle\" font-size=\"16\" font-weight=\"bold\">MeshCut (Dense Grid)</text>\n";
        file << "<line x1=\"450\" y1=\"40\" x2=\"450\" y2=\"480\" stroke=\"#ccc\" stroke-width=\"2\"/>\n";
    }
    
    ~RealPolygonVisualization() {
        file << "</svg>\n";
    }
    
    std::pair<double, double> transform(double x, double y, bool rightPanel = false) {
        double tx = (x - minX) * scale + 25;
        double ty = (y - minY) * scale + 50; // Don't flip Y - MVT coordinates are already correct
        if (rightPanel) tx += 450;
        return {tx, ty};
    }
    
    void drawPolygon(const std::vector<double>& polygon, const std::string& color, bool rightPanel = false) {
        file << "<polygon points=\"";
        for (size_t i = 0; i < polygon.size(); i += 2) {
            auto [x, y] = transform(polygon[i], polygon[i+1], rightPanel);
            file << x << "," << y << " ";
        }
        file << "\" fill=\"" << color << "\" fill-opacity=\"0.3\" stroke=\"" << color << "\" stroke-width=\"2\"/>\n";
    }
    
    void drawTriangles(const std::vector<double>& vertices, const std::vector<uint32_t>& indices, 
                       const std::string& color, bool rightPanel = false) {
        file << "<g fill=\"none\" stroke=\"" << color << "\" stroke-width=\"0.5\" opacity=\"0.8\">\n";
        for (size_t i = 0; i < indices.size(); i += 3) {
            auto [x0, y0] = transform(vertices[indices[i]*2], vertices[indices[i]*2+1], rightPanel);
            auto [x1, y1] = transform(vertices[indices[i+1]*2], vertices[indices[i+1]*2+1], rightPanel);
            auto [x2, y2] = transform(vertices[indices[i+2]*2], vertices[indices[i+2]*2+1], rightPanel);
            file << "<polygon points=\"" << x0 << "," << y0 << " " << x1 << "," << y1 << " " << x2 << "," << y2 << "\"/>\n";
        }
        file << "</g>\n";
    }
    
    void drawGrid(const meshcut::MeshCutOptions& options, bool rightPanel = false) {
        if (!rightPanel) return;
        
        file << "<g stroke=\"#e0e0e0\" stroke-width=\"0.3\" fill=\"none\" opacity=\"0.5\">\n";
        for (int i = 0; i <= options.gridWidth; i++) {
            double x = options.gridOriginX + i * options.cellSize;
            auto [x1, y1] = transform(x, options.gridOriginY, rightPanel);
            auto [x2, y2] = transform(x, options.gridOriginY + options.gridHeight * options.cellSize, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        for (int j = 0; j <= options.gridHeight; j++) {
            double y = options.gridOriginY + j * options.cellSize;
            auto [x1, y1] = transform(options.gridOriginX, y, rightPanel);
            auto [x2, y2] = transform(options.gridOriginX + options.gridWidth * options.cellSize, y, rightPanel);
            file << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\"" << x2 << "\" y2=\"" << y2 << "\"/>\n";
        }
        file << "</g>\n";
    }
};

int main() {
    std::cout << "MeshCut REAL Vector Tile Polygon Test\n";
    std::cout << "=====================================\n";
    
    // Download the actual Innsbruck tile
    int z = 14, x = 8710, y = 5744;
    std::vector<uint8_t> tileData;
    
    std::cout << "Downloading Innsbruck tile: z=" << z << ", x=" << x << ", y=" << y << "\n";
    if (!downloadTile(z, x, y, tileData)) {
        std::cerr << "Failed to download tile\n";
        return 1;
    }
    
    std::cout << "Downloaded " << tileData.size() << " bytes\n";
    
    // Try to decompress if gzipped
    if (tileData.size() > 2 && tileData[0] == 0x1f && tileData[1] == 0x8b) {
        std::cout << "Decompressing gzipped tile...\n";
        auto decompressed = decompressGzip(tileData);
        if (!decompressed.empty()) {
            tileData = std::move(decompressed);
            std::cout << "Decompressed to " << tileData.size() << " bytes\n";
        }
    }
    
    // Parse the MVT and extract real polygons
    MVTParser parser(tileData);
    auto polygons = parser.extractPolygons();
    
    std::cout << "\nExtracted " << polygons.size() << " real polygons from tile\n";
    
    if (polygons.empty()) {
        std::cout << "No polygons found in tile data\n";
        return 1;
    }
    
    // Use the first polygon with enough vertices
    MVTParser::Polygon* selectedPolygon = nullptr;
    for (auto& poly : polygons) {
        if (poly.points.size() >= 4) { // Need at least 4 points (triangle + close)
            selectedPolygon = &poly;
            break;
        }
    }
    
    if (!selectedPolygon) {
        std::cout << "No suitable polygons found (need at least 4 vertices)\n";
        return 1;
    }
    
    std::cout << "\nUsing polygon from layer '" << selectedPolygon->layerName << "'\n";
    std::cout << "Polygon has " << selectedPolygon->points.size() << " vertices\n";
    
    // Convert to flat array for meshcut
    std::vector<double> flatPolygon;
    for (const auto& point : selectedPolygon->points) {
        flatPolygon.push_back(point.x);
        flatPolygon.push_back(point.y);
    }
    
    // Print actual coordinates to verify they're real
    std::cout << "First few coordinates: ";
    for (size_t i = 0; i < std::min(6, (int)flatPolygon.size()); i += 2) {
        std::cout << "(" << std::fixed << std::setprecision(4) << flatPolygon[i] 
                  << "," << flatPolygon[i+1] << ") ";
    }
    std::cout << "\n";
    
    // Set up 32x32 grid as requested
    meshcut::MeshCutOptions options;
    
    // Calculate polygon bounds for grid
    double minX = flatPolygon[0], maxX = flatPolygon[0];
    double minY = flatPolygon[1], maxY = flatPolygon[1];
    for (size_t i = 0; i < flatPolygon.size(); i += 2) {
        minX = std::min(minX, flatPolygon[i]);
        maxX = std::max(maxX, flatPolygon[i]);
        minY = std::min(minY, flatPolygon[i+1]);
        maxY = std::max(maxY, flatPolygon[i+1]);
    }
    
    // Expand grid slightly beyond polygon bounds
    double padding = std::max(maxX - minX, maxY - minY) * 0.1;
    options.gridOriginX = minX - padding;
    options.gridOriginY = minY - padding;
    options.cellSize = (std::max(maxX - minX, maxY - minY) + 2*padding) / 32.0; // 32x32 grid
    options.gridWidth = 32;
    options.gridHeight = 32;
    options.diagonalNE = true;
    
    std::cout << "\nGrid setup: 32x32, cell size = " << std::setprecision(6) << options.cellSize << "\n";
    
    // Test with MeshCut
    std::cout << "Running MeshCut triangulation...\n";
    auto meshcutResult = meshcut::meshcut_full(flatPolygon, {}, options);
    
    std::cout << "MeshCut result: " << meshcutResult.vertices.size()/2 << " vertices, " 
              << meshcutResult.indices.size()/3 << " triangles\n";
    
    // Test with Earcut for comparison
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
    
    // Create visualization of REAL polygon
    std::string filename = "../visualizations/real_innsbruck_polygon.svg";
    RealPolygonVisualization viz(filename, flatPolygon);
    
    // Left panel - Earcut
    viz.drawPolygon(flatPolygon, "#4CAF50", false);
    viz.drawTriangles(flatPolygon, earcutIndices, "#FF5722", false);
    
    // Right panel - MeshCut with grid
    viz.drawGrid(options, true);
    viz.drawPolygon(flatPolygon, "#4CAF50", true);
    viz.drawTriangles(meshcutResult.vertices, meshcutResult.indices, "#2196F3", true);
    
    std::cout << "\nVisualization saved: " << filename << "\n";
    
    double improvement = (double)meshcutResult.indices.size() / (double)earcutIndices.size();
    std::cout << "Triangle density improvement: " << std::fixed << std::setprecision(1) << improvement << "x\n";
    
    std::cout << "\n=== REAL DATA SUCCESS ===\n";
    std::cout << "✓ Downloaded real vector tile: " << tileData.size() << " bytes\n";
    std::cout << "✓ Parsed " << polygons.size() << " real polygons from MVT data\n";
    std::cout << "✓ Used ACTUAL polygon from layer '" << selectedPolygon->layerName << "'\n";
    std::cout << "✓ Applied 32x32 MeshCut grid as requested\n";
    std::cout << "✓ Compared against Earcut baseline\n";
    std::cout << "✓ Generated visualization of real polygon triangulation\n";
    
    return 0;
}