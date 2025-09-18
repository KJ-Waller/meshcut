#pragma once

#include "meshcut/meshcut.hpp"
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>

namespace meshcut {

/**
 * Earcut-compatible interface for MeshCut
 * Provides the same API as mapbox::earcut for drop-in replacement
 */

/**
 * Convert MeshCut Triangle output to Earcut-style indices
 * This is the key function that makes MeshCut compatible with Earcut
 */
template<typename N = uint32_t>
std::vector<N> triangles_to_indices(const std::vector<Triangle>& triangles, 
                                   std::vector<Point>& unique_vertices,
                                   double epsilon = 1e-10) {
    std::vector<N> indices;
    unique_vertices.clear();
    
    // Map from vertex coordinates to index
    std::unordered_map<uint64_t, N> vertex_map;
    
    auto hash_point = [epsilon](const Point& p) -> uint64_t {
        // Quantize coordinates to handle floating point precision
        int64_t x_quantized = static_cast<int64_t>(std::round(p.x / epsilon));
        int64_t y_quantized = static_cast<int64_t>(std::round(p.y / epsilon));
        return (static_cast<uint64_t>(x_quantized) << 32) | 
               (static_cast<uint64_t>(y_quantized) & 0xFFFFFFFF);
    };
    
    auto find_or_add_vertex = [&](const Point& p) -> N {
        uint64_t hash = hash_point(p);
        auto it = vertex_map.find(hash);
        if (it != vertex_map.end()) {
            return it->second;
        }
        
        N index = static_cast<N>(unique_vertices.size());
        unique_vertices.push_back(p);
        vertex_map[hash] = index;
        return index;
    };
    
    // Convert each triangle
    indices.reserve(triangles.size() * 3);
    for (const auto& triangle : triangles) {
        indices.push_back(find_or_add_vertex(triangle.a));
        indices.push_back(find_or_add_vertex(triangle.b));
        indices.push_back(find_or_add_vertex(triangle.c));
    }
    
    return indices;
}

/**
 * Auto-compute bounding box from polygon
 */
inline Polygon compute_bbox(const Polygon& polygon) {
    if (polygon.empty()) {
        return {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    }
    
    double min_x = polygon[0].x, max_x = polygon[0].x;
    double min_y = polygon[0].y, max_y = polygon[0].y;
    
    for (const auto& point : polygon) {
        min_x = std::min(min_x, point.x);
        max_x = std::max(max_x, point.x);
        min_y = std::min(min_y, point.y);
        max_y = std::max(max_y, point.y);
    }
    
    // Add small padding to avoid edge cases
    double padding = std::max(max_x - min_x, max_y - min_y) * 0.01;
    min_x -= padding;
    max_x += padding;
    min_y -= padding;
    max_y += padding;
    
    return {
        {min_x, min_y},  // bottom-left
        {max_x, min_y},  // bottom-right
        {max_x, max_y},  // top-right
        {min_x, max_y}   // top-left
    };
}

/**
 * Drop-in replacement for mapbox::earcut
 * 
 * @param polygon Input polygon (same format as Earcut)
 * @param grid_size Grid resolution (default 32 for 32x32 grid)
 * @param orientation Triangle orientation within grid cells
 * @return Vector of vertex indices (same format as Earcut)
 */
template<typename N = uint32_t, typename InputPolygon>
std::vector<N> earcut_compatible(const InputPolygon& polygon, 
                                int grid_size = 32,
                                TriangleOrientation orientation = TriangleOrientation::TOPLEFT_BOTTOMRIGHT) {
    // Convert input polygon to MeshCut format
    Polygon meshcut_polygon;
    meshcut_polygon.reserve(polygon.size());
    
    // Handle different input formats
    for (const auto& point : polygon) {
        // Support both Point objects and coordinate pairs
        if constexpr (std::is_same_v<typename InputPolygon::value_type, Point>) {
            meshcut_polygon.push_back(point);
        } else {
            // Assume it's a coordinate pair (array, tuple, etc.)
            meshcut_polygon.push_back({static_cast<double>(point[0]), static_cast<double>(point[1])});
        }
    }
    
    // Auto-compute bounding box
    auto bbox = compute_bbox(meshcut_polygon);
    
    // Call MeshCut
    auto triangles = tessellate(meshcut_polygon, bbox, grid_size, orientation);
    
    // Convert to Earcut-compatible indices
    std::vector<Point> vertices;
    return triangles_to_indices<N>(triangles, vertices);
}

/**
 * Alternative interface that also returns the vertex list
 * Useful when you need both the indices and the actual coordinates
 */
template<typename N = uint32_t, typename InputPolygon>
std::pair<std::vector<N>, std::vector<Point>> earcut_compatible_with_vertices(
    const InputPolygon& polygon,
    int grid_size = 32,
    TriangleOrientation orientation = TriangleOrientation::TOPLEFT_BOTTOMRIGHT) {
    
    // Convert input polygon to MeshCut format
    Polygon meshcut_polygon;
    meshcut_polygon.reserve(polygon.size());
    
    for (const auto& point : polygon) {
        if constexpr (std::is_same_v<typename InputPolygon::value_type, Point>) {
            meshcut_polygon.push_back(point);
        } else {
            meshcut_polygon.push_back({static_cast<double>(point[0]), static_cast<double>(point[1])});
        }
    }
    
    auto bbox = compute_bbox(meshcut_polygon);
    auto triangles = tessellate(meshcut_polygon, bbox, grid_size, orientation);
    
    std::vector<Point> vertices;
    auto indices = triangles_to_indices<N>(triangles, vertices);
    
    return {std::move(indices), std::move(vertices)};
}

/**
 * Exact mapbox::earcut signature for drop-in replacement
 * This allows you to replace:
 *   auto indices = mapbox::earcut<uint32_t>(polygon);
 * With:
 *   auto indices = meshcut::earcut<uint32_t>(polygon);
 */
template<typename N = uint32_t, typename Polygon>
std::vector<N> earcut(const Polygon& polygon, int grid_size = 32) {
    return earcut_compatible<N>(polygon, grid_size);
}

} // namespace meshcut