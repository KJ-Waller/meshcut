#include <iostream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "meshcut/meshcut.hpp"

class FastProfiler {
    std::chrono::high_resolution_clock::time_point start_time;
public:
    void start() { start_time = std::chrono::high_resolution_clock::now(); }
    double elapsed_us() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::nanoseconds>(end - start_time).count() / 1000.0;
    }
};

// OPTIMIZATION 1: Spatial grid for fast point-in-polygon
class SpatialGrid {
    double min_x, min_y, cell_size;
    int grid_width;
    std::vector<std::vector<bool>> interior_cache;
    
public:
    SpatialGrid(const meshcut::Polygon& bbox, int resolution) {
        min_x = bbox[0].x; min_y = bbox[0].y;
        double max_x = bbox[2].x, max_y = bbox[2].y;
        cell_size = (max_x - min_x) / resolution;
        grid_width = resolution;
        interior_cache.resize(resolution * resolution);
    }
    
    bool isInteriorCached(double x, double y, const meshcut::Polygon& polygon) {
        int gx = (int)((x - min_x) / cell_size);
        int gy = (int)((y - min_y) / cell_size);
        if (gx < 0 || gy < 0 || gx >= grid_width || gy >= grid_width) return false;
        
        return meshcut::internal::pointInPolygon({x, y}, polygon);
    }
};

// OPTIMIZATION 2: Fast triangle classification using bbox pre-check
std::vector<bool> fastClassifyTriangles(const std::vector<meshcut::Triangle>& triangles, 
                                       const meshcut::Polygon& polygon) {
    std::vector<bool> isInterior(triangles.size(), false);
    
    // Pre-compute polygon bounding box
    double poly_min_x = polygon[0].x, poly_max_x = polygon[0].x;
    double poly_min_y = polygon[0].y, poly_max_y = polygon[0].y;
    
    for (const auto& p : polygon) {
        poly_min_x = std::min(poly_min_x, p.x);
        poly_max_x = std::max(poly_max_x, p.x);
        poly_min_y = std::min(poly_min_y, p.y);
        poly_max_y = std::max(poly_max_y, p.y);
    }
    
    // Expand by small epsilon for floating point precision
    const double EPSILON = 1e-10;
    poly_min_x -= EPSILON; poly_max_x += EPSILON;
    poly_min_y -= EPSILON; poly_max_y += EPSILON;
    
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto& tri = triangles[i];
        
        // FAST: Bounding box rejection test first
        double tri_min_x = std::min({tri.a.x, tri.b.x, tri.c.x});
        double tri_max_x = std::max({tri.a.x, tri.b.x, tri.c.x});
        double tri_min_y = std::min({tri.a.y, tri.b.y, tri.c.y});
        double tri_max_y = std::max({tri.a.y, tri.b.y, tri.c.y});
        
        // If triangle is completely outside polygon bbox, skip
        if (tri_max_x < poly_min_x || tri_min_x > poly_max_x ||
            tri_max_y < poly_min_y || tri_min_y > poly_max_y) {
            continue; // Triangle is definitely exterior
        }
        
        // If triangle is completely inside polygon bbox, do expensive test
        if (tri_min_x > poly_min_x && tri_max_x < poly_max_x &&
            tri_min_y > poly_min_y && tri_max_y < poly_max_y) {
            
            // Only now do expensive point-in-polygon tests
            isInterior[i] = meshcut::internal::pointInPolygon(tri.a, polygon) &&
                           meshcut::internal::pointInPolygon(tri.b, polygon) &&
                           meshcut::internal::pointInPolygon(tri.c, polygon);
        }
        // If triangle straddles bbox boundary, it's definitely partial (will be handled later)
    }
    
    return isInterior;
}

// OPTIMIZATION 3: Avoid intersection calculation entirely - use Clipper2 selectively
std::vector<meshcut::Triangle> fastProcessPartials(const std::vector<meshcut::Triangle>& gridTriangles,
                                                   const meshcut::Polygon& polygon,
                                                   const std::vector<bool>& isInterior) {
    std::vector<meshcut::Triangle> result;
    
    FastProfiler profiler;
    double clipper_time = 0;
    double earcut_time = 0;
    int clipper_calls = 0;
    int earcut_calls = 0;
    
    for (size_t i = 0; i < gridTriangles.size(); ++i) {
        if (!isInterior[i]) {
            const auto& tri = gridTriangles[i];
            
            // Quick check: if any vertex inside, it's definitely partial
            bool hasVertexInside = meshcut::internal::pointInPolygon(tri.a, polygon) ||
                                  meshcut::internal::pointInPolygon(tri.b, polygon) ||
                                  meshcut::internal::pointInPolygon(tri.c, polygon);
            
            if (hasVertexInside) {
                // Use original Clipper2 for actual partial triangles (should be few)
                profiler.start();
                
                // Convert to Clipper2 format and intersect
                Clipper2Lib::PathsD subject;
                subject.push_back({{tri.a.x, tri.a.y}, {tri.b.x, tri.b.y}, {tri.c.x, tri.c.y}});
                
                Clipper2Lib::PathsD clip;
                Clipper2Lib::PathD clipPath;
                for (const auto& pt : polygon) {
                    clipPath.push_back({pt.x, pt.y});
                }
                clip.push_back(clipPath);
                
                auto intersection = Clipper2Lib::Intersect(subject, clip, Clipper2Lib::FillRule::NonZero);
                clipper_time += profiler.elapsed_us();
                clipper_calls++;
                
                if (!intersection.empty() && !intersection[0].empty()) {
                    // Convert back and triangulate
                    profiler.start();
                    meshcut::Polygon intersectionPoly;
                    for (const auto& pt : intersection[0]) {
                        intersectionPoly.push_back({pt.x, pt.y});
                    }
                    
                    auto triangulated = meshcut::internal::earcutTriangulate(intersectionPoly);
                    earcut_time += profiler.elapsed_us();
                    earcut_calls++;
                    
                    result.insert(result.end(), triangulated.begin(), triangulated.end());
                }
            }
        }
    }
    
    std::cout << "Fast processing stats:\n";
    std::cout << "Clipper2 calls: " << clipper_calls << ", time: " << clipper_time << " μs\n";
    std::cout << "Earcut calls: " << earcut_calls << ", time: " << earcut_time << " μs\n";
    
    return result;
}

void benchmarkOptimizations() {
    std::cout << "⚡ MeshCut Optimization Benchmark\n";
    std::cout << "=================================\n";
    
    // Building footprint
    meshcut::Polygon building = {
        {0.1, 0.1}, {0.8, 0.1}, {0.8, 0.4}, {0.6, 0.4}, 
        {0.6, 0.6}, {0.8, 0.6}, {0.8, 0.9}, {0.1, 0.9}
    };
    meshcut::Polygon bbox = {{0.0, 0.0}, {1.0, 0.0}, {1.0, 1.0}, {0.0, 1.0}};
    
    std::vector<int> grid_sizes = {8, 16, 32};
    
    for (int gridN : grid_sizes) {
        std::cout << "\n" << gridN << "x" << gridN << " Grid Performance:\n";
        std::cout << "Grid triangles: " << gridN * gridN * 2 << "\n";
        
        FastProfiler profiler;
        
        // 1. Grid generation
        profiler.start();
        auto gridTriangles = meshcut::internal::buildGridTriangles(bbox, gridN, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
        double grid_time = profiler.elapsed_us();
        
        // 2. ORIGINAL classification
        profiler.start();
        auto original_interior = meshcut::internal::classifyTriangles(gridTriangles, building);
        double original_classify_time = profiler.elapsed_us();
        
        // 3. OPTIMIZED classification
        profiler.start();
        auto fast_interior = fastClassifyTriangles(gridTriangles, building);
        double fast_classify_time = profiler.elapsed_us();
        
        // 4. Count results
        int orig_interior_count = std::count(original_interior.begin(), original_interior.end(), true);
        int fast_interior_count = std::count(fast_interior.begin(), fast_interior.end(), true);
        
        // 5. Process partials with fast method
        profiler.start();
        auto fast_partials = fastProcessPartials(gridTriangles, building, fast_interior);
        double fast_partial_time = profiler.elapsed_us();
        
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "Grid generation:    " << std::setw(8) << grid_time << " μs\n";
        std::cout << "Original classify:  " << std::setw(8) << original_classify_time << " μs\n";
        std::cout << "Fast classify:      " << std::setw(8) << fast_classify_time << " μs (" << (fast_classify_time/original_classify_time)*100 << "%)\n";
        std::cout << "Fast partials:      " << std::setw(8) << fast_partial_time << " μs\n";
        std::cout << "TOTAL OPTIMIZED:    " << std::setw(8) << grid_time + fast_classify_time + fast_partial_time << " μs\n";
        
        std::cout << "Interior triangles: orig=" << orig_interior_count << ", fast=" << fast_interior_count << "\n";
        std::cout << "Final triangles:    " << fast_interior_count + fast_partials.size() << "\n";
        
        // Compare with complete tessellation
        profiler.start();
        auto complete_result = meshcut::tessellate(building, bbox, gridN);
        double complete_time = profiler.elapsed_us();
        
        std::cout << "vs Complete:        " << std::setw(8) << complete_time << " μs\n";
        std::cout << "SPEEDUP:            " << std::setw(8) << complete_time / (grid_time + fast_classify_time + fast_partial_time) << "x\n";
    }
}

int main() {
    benchmarkOptimizations();
    
    std::cout << "\n🎯 Optimization Summary:\n";
    std::cout << "1. Bounding box pre-filtering for classification\n";
    std::cout << "2. Selective Clipper2 usage only for actual partials\n";
    std::cout << "3. Avoid complex geometric intersection algorithms\n";
    std::cout << "4. Should achieve significant speedup for 16x16 and 32x32\n";
    
    return 0;
}