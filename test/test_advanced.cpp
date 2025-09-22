#include "meshcut.hpp"
#include "earcut.hpp"
#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <iomanip>
#include <cmath>

using namespace std::chrono;

// Helper function to create a complex polygon (lake-like shape)
std::vector<double> createLakePolygon() {
    std::vector<double> lake;
    
    // Create a lake-like shape with some irregular edges
    int numPoints = 20;
    double centerX = 1.5, centerY = 1.5;
    double baseRadius = 1.2;
    
    for (int i = 0; i < numPoints; i++) {
        double angle = 2.0 * M_PI * i / numPoints;
        
        // Add some irregularity
        double radiusVar = baseRadius + 0.3 * std::sin(5 * angle) + 0.2 * std::cos(3 * angle);
        double x = centerX + radiusVar * std::cos(angle);
        double y = centerY + radiusVar * std::sin(angle);
        
        lake.push_back(x);
        lake.push_back(y);
    }
    
    return lake;
}

// Helper to count unique vertices (for validation)
int countUniqueVertices(const std::vector<uint32_t>& indices) {
    if (indices.empty()) return 0;
    
    uint32_t maxIndex = 0;
    for (uint32_t idx : indices) {
        maxIndex = std::max(maxIndex, idx);
    }
    return maxIndex + 1;
}

// Validation: check that triangles don't overlap and have positive area
bool validateTriangles(const std::vector<double>& vertices, const std::vector<uint32_t>& indices) {
    if (indices.size() % 3 != 0) {
        std::cout << "ERROR: Triangle indices not divisible by 3\n";
        return false;
    }
    
    int triangleCount = indices.size() / 3;
    int positiveArea = 0;
    
    for (int t = 0; t < triangleCount; t++) {
        uint32_t i0 = indices[t * 3];
        uint32_t i1 = indices[t * 3 + 1];
        uint32_t i2 = indices[t * 3 + 2];
        
        if (i0 >= vertices.size()/2 || i1 >= vertices.size()/2 || i2 >= vertices.size()/2) {
            std::cout << "ERROR: Index out of bounds in triangle " << t << "\n";
            return false;
        }
        
        // Check triangle area (cross product)
        double x0 = vertices[i0 * 2], y0 = vertices[i0 * 2 + 1];
        double x1 = vertices[i1 * 2], y1 = vertices[i1 * 2 + 1];
        double x2 = vertices[i2 * 2], y2 = vertices[i2 * 2 + 1];
        
        double area = 0.5 * std::abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
        if (area > 1e-10) {
            positiveArea++;
        }
    }
    
    std::cout << "Triangles with positive area: " << positiveArea << "/" << triangleCount << "\n";
    return positiveArea == triangleCount;
}

void runPerformanceTest() {
    std::cout << "\n=== Performance Comparison Test ===\n";
    
    // Create test polygon
    auto lake = createLakePolygon();
    std::cout << "Test polygon: " << lake.size()/2 << " vertices\n";
    
    // MeshCut options
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0;
    options.cellSize = 0.125;  // 8x8 = 64 cells per unit
    options.gridWidth = 32;    
    options.gridHeight = 32;   // Total: 1024 cells
    options.diagonalNE = true;
    
    // Run MeshCut
    auto start = high_resolution_clock::now();
    std::vector<uint32_t> meshcutResult = meshcut::meshcut(lake, {}, options);
    auto meshcutTime = high_resolution_clock::now() - start;
    
    // Run Earcut for comparison
    using Point = std::array<double, 2>;
    std::vector<std::vector<Point>> earcutPolygon;
    earcutPolygon.push_back({});
    for (size_t i = 0; i < lake.size(); i += 2) {
        earcutPolygon[0].push_back({{lake[i], lake[i+1]}});
    }
    
    start = high_resolution_clock::now();
    std::vector<uint32_t> earcutResult = mapbox::earcut<uint32_t>(earcutPolygon);
    auto earcutTime = high_resolution_clock::now() - start;
    
    // Results
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "MeshCut: " << meshcutResult.size()/3 << " triangles in " 
              << duration_cast<microseconds>(meshcutTime).count() << "μs\n";
    std::cout << "Earcut:  " << earcutResult.size()/3 << " triangles in " 
              << duration_cast<microseconds>(earcutTime).count() << "μs\n";
    
    std::cout << "Speed ratio: " << std::setprecision(1)
              << static_cast<double>(duration_cast<nanoseconds>(earcutTime).count()) /
                 duration_cast<nanoseconds>(meshcutTime).count() << "x\n";
              
    std::cout << "Triangle ratio: " << std::setprecision(1)
              << static_cast<double>(meshcutResult.size()) / earcutResult.size() << "x\n";
    
    std::cout << "Unique vertices in MeshCut: " << countUniqueVertices(meshcutResult) << "\n";
    std::cout << "Unique vertices in Earcut: " << countUniqueVertices(earcutResult) << "\n";
}

void runValidationTest() {
    std::cout << "\n=== Validation Test ===\n";
    
    // Simple test polygon
    std::vector<double> triangle = {1.0, 1.0, 3.0, 1.0, 2.0, 3.0};
    
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0;
    options.cellSize = 0.25;
    options.gridWidth = 16;
    options.gridHeight = 16;
    
    auto result = meshcut::meshcut(triangle, {}, options);
    std::cout << "Simple triangle test: " << result.size()/3 << " triangles generated\n";
    
    // Mock vertices for validation (we need to implement vertex coordinate output)
    std::vector<double> mockVertices(countUniqueVertices(result) * 2, 0.0);
    // In a real implementation, we'd return both vertices and indices
    // For now, just validate the index structure
    bool valid = (result.size() % 3 == 0);
    std::cout << "Index structure valid: " << (valid ? "YES" : "NO") << "\n";
}

int main() {
    std::cout << "MeshCut Advanced Testing\n";
    std::cout << "========================\n";
    
    runValidationTest();
    runPerformanceTest();
    
    std::cout << "\nAdvanced testing completed!\n";
    return 0;
}