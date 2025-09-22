#include "meshcut.hpp"
#include "earcut.hpp" 
#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <random>

using namespace std::chrono;

// Create various test polygons of different complexity
class PolygonGenerator {
public:
    // Simple regular polygon
    static std::vector<double> createRegularPolygon(int vertices, double radius = 1.5, double centerX = 2.0, double centerY = 2.0) {
        std::vector<double> polygon;
        polygon.reserve(vertices * 2);
        
        for (int i = 0; i < vertices; i++) {
            double angle = 2.0 * M_PI * i / vertices;
            polygon.push_back(centerX + radius * std::cos(angle));
            polygon.push_back(centerY + radius * std::sin(angle));
        }
        return polygon;
    }
    
    // Complex irregular polygon (like a lake)
    static std::vector<double> createComplexPolygon(int baseVertices = 50, double radius = 1.8, double centerX = 2.0, double centerY = 2.0) {
        std::vector<double> polygon;
        polygon.reserve(baseVertices * 2);
        
        std::mt19937 rng(42); // Fixed seed for reproducible results
        std::uniform_real_distribution<double> radiusVar(0.7, 1.3);
        std::uniform_real_distribution<double> angleVar(-0.1, 0.1);
        
        for (int i = 0; i < baseVertices; i++) {
            double baseAngle = 2.0 * M_PI * i / baseVertices;
            double angle = baseAngle + angleVar(rng);
            double r = radius * radiusVar(rng);
            
            // Add some higher-frequency variation
            r += 0.2 * std::sin(5 * baseAngle) + 0.1 * std::cos(7 * baseAngle);
            
            polygon.push_back(centerX + r * std::cos(angle));
            polygon.push_back(centerY + r * std::sin(angle));
        }
        return polygon;
    }
    
    // Very complex polygon with many vertices
    static std::vector<double> createVeryComplexPolygon(int vertices = 200) {
        return createComplexPolygon(vertices, 1.9, 2.0, 2.0);
    }
};

// Performance benchmark runner
class PerformanceBenchmark {
private:
    struct BenchmarkResult {
        std::string name;
        int polygonVertices;
        int gridCells;
        int trianglesGenerated;
        double timeMs;
        double trianglesPerMs;
    };
    
    std::vector<BenchmarkResult> results;
    
public:
    void runBenchmark(const std::string& name, const std::vector<double>& polygon, const meshcut::MeshCutOptions& options, int iterations = 100) {
        std::cout << "Running benchmark: " << name << " (" << iterations << " iterations)\n";
        
        // Warmup run
        meshcut::meshcut(polygon, {}, options);
        
        // Timed runs
        auto start = high_resolution_clock::now();
        std::vector<uint32_t> result;
        
        for (int i = 0; i < iterations; i++) {
            result = meshcut::meshcut(polygon, {}, options);
        }
        
        auto end = high_resolution_clock::now();
        double totalTimeMs = duration_cast<nanoseconds>(end - start).count() / 1e6;
        double avgTimeMs = totalTimeMs / iterations;
        
        BenchmarkResult bench;
        bench.name = name;
        bench.polygonVertices = polygon.size() / 2;
        bench.gridCells = options.gridWidth * options.gridHeight;
        bench.trianglesGenerated = result.size() / 3;
        bench.timeMs = avgTimeMs;
        bench.trianglesPerMs = bench.trianglesGenerated / avgTimeMs;
        
        results.push_back(bench);
        
        std::cout << "  " << std::fixed << std::setprecision(3) 
                  << avgTimeMs << "ms avg, " << bench.trianglesGenerated << " triangles\n";
    }
    
    void runEarcutComparison(const std::string& name, const std::vector<double>& polygon, int iterations = 100) {
        std::cout << "Running earcut comparison: " << name << "\n";
        
        // Convert to earcut format
        using Point = std::array<double, 2>;
        std::vector<std::vector<Point>> earcutPolygon;
        earcutPolygon.push_back({});
        for (size_t i = 0; i < polygon.size(); i += 2) {
            earcutPolygon[0].push_back({{polygon[i], polygon[i+1]}});
        }
        
        // Warmup
        mapbox::earcut<uint32_t>(earcutPolygon);
        
        // Timed runs
        auto start = high_resolution_clock::now();
        std::vector<uint32_t> result;
        
        for (int i = 0; i < iterations; i++) {
            result = mapbox::earcut<uint32_t>(earcutPolygon);
        }
        
        auto end = high_resolution_clock::now();
        double totalTimeMs = duration_cast<nanoseconds>(end - start).count() / 1e6;
        double avgTimeMs = totalTimeMs / iterations;
        
        std::cout << "  Earcut: " << std::fixed << std::setprecision(3) 
                  << avgTimeMs << "ms avg, " << result.size()/3 << " triangles\n";
    }
    
    void printSummary() {
        std::cout << "\n" << std::string(80, '=') << "\n";
        std::cout << "PERFORMANCE BENCHMARK SUMMARY\n";
        std::cout << std::string(80, '=') << "\n";
        
        std::cout << std::left << std::setw(20) << "Test" 
                  << std::setw(10) << "Vertices" 
                  << std::setw(10) << "Grid"
                  << std::setw(12) << "Triangles"
                  << std::setw(12) << "Time (ms)"
                  << std::setw(15) << "Tri/ms" << "\n";
        std::cout << std::string(80, '-') << "\n";
        
        for (const auto& result : results) {
            std::cout << std::left << std::setw(20) << result.name
                      << std::setw(10) << result.polygonVertices
                      << std::setw(10) << (std::to_string(result.gridCells) + " cells")
                      << std::setw(12) << result.trianglesGenerated
                      << std::setw(12) << std::fixed << std::setprecision(3) << result.timeMs
                      << std::setw(15) << std::fixed << std::setprecision(1) << result.trianglesPerMs << "\n";
        }
        
        if (results.size() > 1) {
            std::cout << "\nSCALING ANALYSIS:\n";
            for (size_t i = 1; i < results.size(); i++) {
                double speedupRatio = results[0].timeMs / results[i].timeMs;
                double complexityRatio = static_cast<double>(results[i].polygonVertices) / results[0].polygonVertices;
                std::cout << "  " << results[i].name << " vs " << results[0].name 
                          << ": " << std::setprecision(2) << complexityRatio << "x complexity"
                          << ", " << (speedupRatio > 1 ? speedupRatio : -1.0/speedupRatio) << "x " 
                          << (speedupRatio > 1 ? "slower" : "faster") << "\n";
            }
        }
    }
};

int main() {
    std::cout << "🚀 MeshCut Performance Benchmark Suite\n";
    std::cout << "=======================================\n\n";
    
    PerformanceBenchmark benchmark;
    
    // Standard grid options
    meshcut::MeshCutOptions options;
    options.gridOriginX = 0.0;
    options.gridOriginY = 0.0; 
    options.cellSize = 0.25;
    options.gridWidth = 16;
    options.gridHeight = 16;
    options.diagonalNE = true;
    
    // Test 1: Simple polygon (baseline)
    auto simplePolygon = PolygonGenerator::createRegularPolygon(6, 1.5, 2.0, 2.0);
    benchmark.runBenchmark("Simple (6-gon)", simplePolygon, options, 1000);
    benchmark.runEarcutComparison("Simple (6-gon)", simplePolygon, 1000);
    
    std::cout << "\n";
    
    // Test 2: Medium complexity
    auto mediumPolygon = PolygonGenerator::createComplexPolygon(20, 1.8, 2.0, 2.0);
    benchmark.runBenchmark("Medium (20 vertices)", mediumPolygon, options, 500);
    benchmark.runEarcutComparison("Medium (20 vertices)", mediumPolygon, 500);
    
    std::cout << "\n";
    
    // Test 3: High complexity
    auto complexPolygon = PolygonGenerator::createComplexPolygon(50, 1.8, 2.0, 2.0);
    benchmark.runBenchmark("Complex (50 vertices)", complexPolygon, options, 200);
    benchmark.runEarcutComparison("Complex (50 vertices)", complexPolygon, 200);
    
    std::cout << "\n";
    
    // Test 4: Very high complexity
    auto veryComplexPolygon = PolygonGenerator::createVeryComplexPolygon(100);
    benchmark.runBenchmark("Very Complex (100)", veryComplexPolygon, options, 100);
    benchmark.runEarcutComparison("Very Complex (100)", veryComplexPolygon, 100);
    
    std::cout << "\n";
    
    // Test 5: Different grid sizes
    meshcut::MeshCutOptions fineGrid = options;
    fineGrid.cellSize = 0.125;
    fineGrid.gridWidth = 32;  
    fineGrid.gridHeight = 32;
    
    benchmark.runBenchmark("Fine Grid (32x32)", mediumPolygon, fineGrid, 200);
    
    benchmark.printSummary();
    
    std::cout << "\n🎯 This benchmark tests the spatial indexing optimization.\n";
    std::cout << "   Look for consistent performance scaling with polygon complexity.\n";
    std::cout << "   Before spatial indexing: O(n*m) scaling was terrible for complex polygons.\n";
    std::cout << "   After spatial indexing: Should scale much better!\n\n";
    
    return 0;
}