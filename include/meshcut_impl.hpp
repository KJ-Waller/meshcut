#pragma once

#include "earcut.hpp"

#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <tuple>
#include <type_traits>

namespace meshcut {

// Forward declarations for internal functions
namespace detail {

/**
 * World coordinate (for precise vertex deduplication)
 */
struct WorldPoint {
    double x, y;
    
    bool operator==(const WorldPoint& other) const {
        const double eps = 1e-9;
        return std::abs(x - other.x) < eps && std::abs(y - other.y) < eps;
    }
};

/**
 * Hash function for WorldPoint (for vertex deduplication)
 */
struct WorldPointHash {
    std::size_t operator()(const WorldPoint& p) const {
        // Quantize to avoid floating-point precision issues
        const double scale = 1e6;  // 1 micrometer precision
        long long ix = static_cast<long long>(std::round(p.x * scale));
        long long iy = static_cast<long long>(std::round(p.y * scale));
        return std::hash<long long>()(ix) ^ (std::hash<long long>()(iy) << 1);
    }
};

/**
 * Grid coordinate (integer grid space)
 */
struct GridPoint {
    int x, y;
    
    bool operator==(const GridPoint& other) const {
        return x == other.x && y == other.y;
    }
};

/**
 * Hash function for GridPoint (for vertex deduplication)
 */
struct GridPointHash {
    std::size_t operator()(const GridPoint& p) const {
        // Combine hash of x and y coordinates
        return std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1);
    }
};

/**
 * Grid cell classification
 */
enum CellState {
    FULL_IN,    // All corners inside polygon - emit grid triangles directly
    FULL_OUT,   // No intersection with polygon - skip
    BOUNDARY    // Partial intersection - needs clipping and triangulation
};

/**
 * Coordinate system transformation functions
 */
class CoordinateTransform {
public:
    CoordinateTransform(const MeshCutOptions& options)
        : originX(options.gridOriginX)
        , originY(options.gridOriginY) 
        , cellSize(options.cellSize)
        , invCellSize(1.0 / options.cellSize)
    {}
    
    // World to grid coordinate conversion
    GridPoint worldToGrid(double x, double y) const {
        return {
            static_cast<int>(std::floor((x - originX) * invCellSize)),
            static_cast<int>(std::floor((y - originY) * invCellSize))
        };
    }
    
    // Grid to world coordinate conversion  
    std::pair<double, double> gridToWorld(int gx, int gy) const {
        return {
            originX + gx * cellSize,
            originY + gy * cellSize
        };
    }
    
    // Grid to world coordinate conversion (for grid corners)
    std::pair<double, double> gridToWorld(const GridPoint& p) const {
        return gridToWorld(p.x, p.y);
    }
    
private:
    double originX, originY, cellSize, invCellSize;
};

/**
 * Edge representation for scanline algorithm
 */
struct PolygonEdge {
    double x1, y1, x2, y2;  // Edge endpoints
    double minY, maxY;      // Y bounds for quick culling
    
    PolygonEdge(double x1_, double y1_, double x2_, double y2_) 
        : x1(x1_), y1(y1_), x2(x2_), y2(y2_)
        , minY(std::min(y1_, y2_)), maxY(std::max(y1_, y2_)) {}
    
    // Compute X intersection at given Y
    double getXAtY(double y) const {
        if (std::abs(y2 - y1) < 1e-10) return x1; // Horizontal edge
        return x1 + (x2 - x1) * (y - y1) / (y2 - y1);
    }
    
    // Check if edge intersects horizontal ray at (x, y) going right
    bool intersectsRay(double x, double y) const {
        // Quick Y bounds check
        if (y < minY || y > maxY) return false;
        
        // Handle edge cases - endpoint exactly on ray
        if (std::abs(y - y1) < 1e-10) {
            // Only count if other endpoint is above ray
            return y2 > y;
        }
        if (std::abs(y - y2) < 1e-10) {
            // Only count if other endpoint is above ray  
            return y1 > y;
        }
        
        // Compute intersection X
        double intersectX = getXAtY(y);
        return intersectX > x;
    }
};

/**
 * Optimized point-in-polygon using precomputed edge list
 */
class PolygonTester {
private:
    std::vector<PolygonEdge> edges;
    
public:
    PolygonTester(const std::vector<double>& polygon) {
        if (polygon.size() < 6) return;
        
        int n = polygon.size() / 2;
        edges.reserve(n);
        
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            double x1 = polygon[i * 2];
            double y1 = polygon[i * 2 + 1];
            double x2 = polygon[j * 2];  
            double y2 = polygon[j * 2 + 1];
            
            // Skip degenerate edges
            if (std::abs(x1 - x2) < 1e-10 && std::abs(y1 - y2) < 1e-10) continue;
            
            edges.emplace_back(x1, y1, x2, y2);
        }
    }
    
    bool isInside(double x, double y) const {
        int intersections = 0;
        
        for (const auto& edge : edges) {
            if (edge.intersectsRay(x, y)) {
                intersections++;
            }
        }
        
        return (intersections % 2) == 1;
    }
    
    // Test if any polygon edges intersect the axis-aligned rectangle
    bool intersectsRect(double left, double right, double bottom, double top) const {
        for (const auto& edge : edges) {
            // Quick bounds check
            double edgeMinX = std::min(edge.x1, edge.x2);
            double edgeMaxX = std::max(edge.x1, edge.x2);
            if (edgeMaxX < left || edgeMinX > right || edge.maxY < bottom || edge.minY > top) {
                continue; // No intersection possible
            }
            
            // Check if edge intersects any of the 4 rectangle edges
            if (lineSegmentIntersectsRect(edge.x1, edge.y1, edge.x2, edge.y2, left, right, bottom, top)) {
                return true;
            }
        }
        return false;
    }

private:
    // Helper: Check if line segment intersects axis-aligned rectangle
    bool lineSegmentIntersectsRect(double x1, double y1, double x2, double y2,
                                  double left, double right, double bottom, double top) const {
        // Cohen-Sutherland-like approach
        auto outcode = [&](double x, double y) {
            int code = 0;
            if (x < left) code |= 1;      // LEFT
            if (x > right) code |= 2;     // RIGHT  
            if (y < bottom) code |= 4;    // BOTTOM
            if (y > top) code |= 8;       // TOP
            return code;
        };
        
        int code1 = outcode(x1, y1);
        int code2 = outcode(x2, y2);
        
        // Both endpoints inside rectangle
        if (code1 == 0 || code2 == 0) return true;
        
        // Both endpoints on same side of rectangle
        if (code1 & code2) return false;
        
        // Need more detailed intersection test
        // Check intersection with each rectangle edge
        
        // Left edge
        if ((code1 & 1) != (code2 & 1)) {
            double y = y1 + (y2 - y1) * (left - x1) / (x2 - x1);
            if (y >= bottom && y <= top) return true;
        }
        
        // Right edge  
        if ((code1 & 2) != (code2 & 2)) {
            double y = y1 + (y2 - y1) * (right - x1) / (x2 - x1);
            if (y >= bottom && y <= top) return true;
        }
        
        // Bottom edge
        if ((code1 & 4) != (code2 & 4)) {
            double x = x1 + (x2 - x1) * (bottom - y1) / (y2 - y1);
            if (x >= left && x <= right) return true;
        }
        
        // Top edge
        if ((code1 & 8) != (code2 & 8)) {
            double x = x1 + (x2 - x1) * (top - y1) / (y2 - y1);
            if (x >= left && x <= right) return true;
        }
        
        return false;
    }
};

/**
 * Simple point-in-polygon test (placeholder - will optimize with scanline later)
 */

/**
 * Classify a grid cell based on polygon intersection
 */
CellState classifyCell(int i, int j, const PolygonTester& tester, const CoordinateTransform& transform) {
    // Get the 4 corners of the cell in world coordinates
    auto [x0, y0] = transform.gridToWorld(i, j);
    auto [x1, y1] = transform.gridToWorld(i + 1, j);
    auto [x2, y2] = transform.gridToWorld(i + 1, j + 1);
    auto [x3, y3] = transform.gridToWorld(i, j + 1);
    
    // Test all 4 corners
    bool in0 = tester.isInside(x0, y0);
    bool in1 = tester.isInside(x1, y1);
    bool in2 = tester.isInside(x2, y2);
    bool in3 = tester.isInside(x3, y3);
    
    // All corners inside = fully inside
    if (in0 && in1 && in2 && in3) {
        return FULL_IN;
    }
    
    // No corners inside - check if polygon edges intersect cell
    if (!in0 && !in1 && !in2 && !in3) {
        // Use improved edge-cell intersection test
        if (tester.intersectsRect(x0, x1, y0, y2)) {
            return BOUNDARY;
        } else {
            return FULL_OUT;
        }
    }
    
    // Mixed corners = boundary
    return BOUNDARY;
}

/**
 * Generate two triangles for a grid cell (FULL_IN case)
 */
template<typename N>
void emitGridTriangles(int i, int j, bool diagonalNE, 
                      std::vector<N>& indices,
                      std::unordered_map<WorldPoint, N, WorldPointHash>& vertexMap,
                      N& nextVertexIndex,
                      const CoordinateTransform& transform) {
    
    // Get world coordinates for the 4 corners
    auto [x0, y0] = transform.gridToWorld(i, j);
    auto [x1, y1] = transform.gridToWorld(i + 1, j);
    auto [x2, y2] = transform.gridToWorld(i + 1, j + 1);
    auto [x3, y3] = transform.gridToWorld(i, j + 1);
    
    WorldPoint corners[4] = {{x0, y0}, {x1, y1}, {x2, y2}, {x3, y3}};
    N vertexIndices[4];
    
    // Get or create vertex indices for the 4 corners
    for (int c = 0; c < 4; c++) {
        auto it = vertexMap.find(corners[c]);
        if (it != vertexMap.end()) {
            vertexIndices[c] = it->second;
        } else {
            vertexIndices[c] = nextVertexIndex++;
            vertexMap[corners[c]] = vertexIndices[c];
        }
    }
    
    // Emit triangles based on diagonal direction
    if (diagonalNE) {
        // NE diagonal: (0,1,3) and (1,2,3)
        indices.push_back(vertexIndices[0]);
        indices.push_back(vertexIndices[1]);
        indices.push_back(vertexIndices[3]);
        
        indices.push_back(vertexIndices[1]);
        indices.push_back(vertexIndices[2]);
        indices.push_back(vertexIndices[3]);
    } else {
        // NW diagonal: (0,1,2) and (0,2,3)
        indices.push_back(vertexIndices[0]);
        indices.push_back(vertexIndices[1]);
        indices.push_back(vertexIndices[2]);
        
        indices.push_back(vertexIndices[0]);
        indices.push_back(vertexIndices[2]);
        indices.push_back(vertexIndices[3]);
    }
}

/**
 * Sutherland-Hodgman polygon clipping against axis-aligned rectangle
 * Clips input polygon to the rectangle [left, right] x [bottom, top]
 */
std::vector<std::pair<double, double>> clipPolygonToRect(
    const std::vector<double>& polygon, 
    double left, double right, double bottom, double top) {
    
    if (polygon.size() < 6) return {}; // Need at least 3 vertices
    
    // Convert to point pairs for easier processing
    std::vector<std::pair<double, double>> points;
    for (size_t i = 0; i < polygon.size(); i += 2) {
        points.emplace_back(polygon[i], polygon[i + 1]);
    }
    
    // Clip against each edge of the rectangle
    // Order: left, right, bottom, top
    auto clipLeft = [&](const std::vector<std::pair<double, double>>& input) {
        std::vector<std::pair<double, double>> output;
        if (input.empty()) return output;
        
        auto prev = input.back();
        for (const auto& curr : input) {
            if (curr.first >= left) { // Current point is inside
                if (prev.first < left) { // Previous point was outside
                    // Compute intersection with left edge
                    double t = (left - prev.first) / (curr.first - prev.first);
                    double y = prev.second + t * (curr.second - prev.second);
                    output.emplace_back(left, y);
                }
                output.push_back(curr);
            } else if (prev.first >= left) { // Current outside, previous inside
                // Compute intersection with left edge
                double t = (left - prev.first) / (curr.first - prev.first);
                double y = prev.second + t * (curr.second - prev.second);
                output.emplace_back(left, y);
            }
            prev = curr;
        }
        return output;
    };
    
    auto clipRight = [&](const std::vector<std::pair<double, double>>& input) {
        std::vector<std::pair<double, double>> output;
        if (input.empty()) return output;
        
        auto prev = input.back();
        for (const auto& curr : input) {
            if (curr.first <= right) { // Current point is inside
                if (prev.first > right) { // Previous point was outside
                    double t = (right - prev.first) / (curr.first - prev.first);
                    double y = prev.second + t * (curr.second - prev.second);
                    output.emplace_back(right, y);
                }
                output.push_back(curr);
            } else if (prev.first <= right) { // Current outside, previous inside
                double t = (right - prev.first) / (curr.first - prev.first);
                double y = prev.second + t * (curr.second - prev.second);
                output.emplace_back(right, y);
            }
            prev = curr;
        }
        return output;
    };
    
    auto clipBottom = [&](const std::vector<std::pair<double, double>>& input) {
        std::vector<std::pair<double, double>> output;
        if (input.empty()) return output;
        
        auto prev = input.back();
        for (const auto& curr : input) {
            if (curr.second >= bottom) { // Current point is inside
                if (prev.second < bottom) { // Previous point was outside
                    double t = (bottom - prev.second) / (curr.second - prev.second);
                    double x = prev.first + t * (curr.first - prev.first);
                    output.emplace_back(x, bottom);
                }
                output.push_back(curr);
            } else if (prev.second >= bottom) { // Current outside, previous inside
                double t = (bottom - prev.second) / (curr.second - prev.second);
                double x = prev.first + t * (curr.first - prev.first);
                output.emplace_back(x, bottom);
            }
            prev = curr;
        }
        return output;
    };
    
    auto clipTop = [&](const std::vector<std::pair<double, double>>& input) {
        std::vector<std::pair<double, double>> output;
        if (input.empty()) return output;
        
        auto prev = input.back();
        for (const auto& curr : input) {
            if (curr.second <= top) { // Current point is inside
                if (prev.second > top) { // Previous point was outside
                    double t = (top - prev.second) / (curr.second - prev.second);
                    double x = prev.first + t * (curr.first - prev.first);
                    output.emplace_back(x, top);
                }
                output.push_back(curr);
            } else if (prev.second <= top) { // Current outside, previous inside
                double t = (top - prev.second) / (curr.second - prev.second);
                double x = prev.first + t * (curr.first - prev.first);
                output.emplace_back(x, top);
            }
            prev = curr;
        }
        return output;
    };
    
    // Apply clipping stages sequentially
    auto clipped = clipLeft(points);
    clipped = clipRight(clipped);
    clipped = clipBottom(clipped);
    clipped = clipTop(clipped);
    
    return clipped;
}

/**
 * Process a boundary cell - clip polygon and triangulate
 */
template<typename N>
void processBoundaryCell(int i, int j, const std::vector<double>& polygon, 
                        const CoordinateTransform& transform,
                        std::vector<N>& indices,
                        std::unordered_map<WorldPoint, N, WorldPointHash>& vertexMap,
                        N& nextVertexIndex) {
    
    // Get cell bounds in world coordinates
    auto [left, bottom] = transform.gridToWorld(i, j);
    auto [right, top] = transform.gridToWorld(i + 1, j + 1);
    
    // Clip polygon to cell rectangle
    auto clippedPoints = clipPolygonToRect(polygon, left, right, bottom, top);
    
    if (clippedPoints.size() < 3) {
        return; // Clipped polygon is degenerate
    }
    
    // Convert clipped polygon to earcut format
    std::vector<std::vector<std::pair<double, double>>> earcutPolygon;
    earcutPolygon.push_back(clippedPoints);
    
    // Triangulate with earcut
    std::vector<N> triangleIndices = mapbox::earcut<N>(earcutPolygon);
    
    // Map earcut indices to global vertex indices with deduplication
    std::vector<N> globalIndices;
    globalIndices.reserve(clippedPoints.size());
    
    for (const auto& point : clippedPoints) {
        WorldPoint worldPoint{point.first, point.second};
        
        // Check if we already have this vertex
        auto it = vertexMap.find(worldPoint);
        if (it != vertexMap.end()) {
            globalIndices.push_back(it->second);
        } else {
            // Create new vertex index
            N newIndex = nextVertexIndex++;
            vertexMap[worldPoint] = newIndex;
            globalIndices.push_back(newIndex);
        }
    }
    
    // Add triangles to output using global indices
    for (N triangleIndex : triangleIndices) {
        indices.push_back(globalIndices[triangleIndex]);
    }
}

} // namespace detail

/**
 * Main meshcut implementation (full result with vertices)
 */
template <typename N>
MeshCutResult meshcut_full(const std::vector<double>& polygon, const std::vector<N>& holes, const MeshCutOptions& options) {
    
    // Handle holes later - for now just process outer ring
    if (!holes.empty()) {
        // TODO: Implement hole support
        // For now, just process the outer ring
    }
    
    if (polygon.size() < 6) {
        return {}; // Need at least 3 vertices (6 coordinates)
    }
    
    std::vector<N> indices;
    detail::CoordinateTransform transform(options);
    detail::PolygonTester tester(polygon);
    
    // Compute polygon bounding box in grid coordinates
    double minX = polygon[0], maxX = polygon[0];
    double minY = polygon[1], maxY = polygon[1];
    for (size_t i = 2; i < polygon.size(); i += 2) {
        minX = std::min(minX, polygon[i]);
        maxX = std::max(maxX, polygon[i]);
        minY = std::min(minY, polygon[i + 1]);
        maxY = std::max(maxY, polygon[i + 1]);
    }
    
    auto minGrid = transform.worldToGrid(minX, minY);
    auto maxGrid = transform.worldToGrid(maxX, maxY);
    
    // Clamp to actual grid bounds
    int minI = std::max(0, minGrid.x);
    int maxI = std::min(options.gridWidth - 1, maxGrid.x);
    int minJ = std::max(0, minGrid.y);
    int maxJ = std::min(options.gridHeight - 1, maxGrid.y);
    
    // Vertex deduplication map and counter
    std::unordered_map<detail::WorldPoint, N, detail::WorldPointHash> vertexMap;
    N nextVertexIndex = 0;
    
    // Process each grid cell
    for (int j = minJ; j <= maxJ; j++) {
        for (int i = minI; i <= maxI; i++) {
            detail::CellState state = detail::classifyCell(i, j, tester, transform);
            
            switch (state) {
                case detail::FULL_IN:
                    detail::emitGridTriangles(i, j, options.diagonalNE, indices, vertexMap, nextVertexIndex, transform);
                    break;
                    
                case detail::FULL_OUT:
                    // Skip - no triangles to emit
                    break;
                    
                case detail::BOUNDARY:
                    detail::processBoundaryCell(i, j, polygon, transform, indices, vertexMap, nextVertexIndex);
                    break;
            }
        }
    }
    
    // Build vertex array from vertex map
    std::vector<double> vertices(nextVertexIndex * 2);
    for (const auto& [worldPoint, index] : vertexMap) {
        vertices[index * 2] = worldPoint.x;
        vertices[index * 2 + 1] = worldPoint.y;
    }
    
    // Convert indices to uint32_t if needed
    std::vector<uint32_t> finalIndices;
    finalIndices.reserve(indices.size());
    for (N idx : indices) {
        finalIndices.push_back(static_cast<uint32_t>(idx));
    }
    
    return {std::move(vertices), std::move(finalIndices)};
}

/**
 * Main meshcut implementation (indices only, earcut-compatible)
 */
template <typename N>
std::vector<N> meshcut(const std::vector<double>& polygon, const std::vector<N>& holes, const MeshCutOptions& options) {
    auto result = meshcut_full(polygon, holes, options);
    
    // Convert back to requested index type
    std::vector<N> indices;
    indices.reserve(result.indices.size());
    for (uint32_t idx : result.indices) {
        indices.push_back(static_cast<N>(idx));
    }
    return indices;
}

/**
 * Template specialization for custom polygon types (full result)
 */
template <typename N, typename Polygon>
MeshCutResult meshcut_full(const Polygon& polygon, const MeshCutOptions& options) {
    // Convert custom polygon format to vector<double> format
    std::vector<double> coords;
    
    // Assume polygon is vector of vector of points (like earcut)
    // Extract outer ring (first element)
    if (!polygon.empty()) {
        const auto& ring = polygon[0];
        coords.reserve(ring.size() * 2);
        
        for (const auto& point : ring) {
            // Handle different point types (array, pair, etc.)
            if constexpr (std::is_same_v<std::decay_t<decltype(point)>, std::array<double, 2>>) {
                coords.push_back(point[0]);
                coords.push_back(point[1]);
            } else {
                // Use mapbox::util::nth for other types (earcut compatibility)
                coords.push_back(mapbox::util::nth<0, decltype(point)>::get(point));
                coords.push_back(mapbox::util::nth<1, decltype(point)>::get(point));
            }
        }
    }
    
    // TODO: Handle holes (additional rings in polygon vector)
    std::vector<N> holes; // Empty for now
    
    return meshcut_full<N>(coords, holes, options);
}

/**
 * Template specialization for custom polygon types (indices only, earcut compatibility)
 */
template <typename N, typename Polygon>
std::vector<N> meshcut(const Polygon& polygon, const MeshCutOptions& options) {
    // Convert custom polygon format to vector<double> format
    std::vector<double> coords;
    
    // Assume polygon is vector of vector of points (like earcut)
    // Extract outer ring (first element)
    if (!polygon.empty()) {
        const auto& ring = polygon[0];
        coords.reserve(ring.size() * 2);
        
        for (const auto& point : ring) {
            // Handle different point types (array, pair, etc.)
            if constexpr (std::is_same_v<std::decay_t<decltype(point)>, std::array<double, 2>>) {
                coords.push_back(point[0]);
                coords.push_back(point[1]);
            } else {
                // Use mapbox::util::nth for other types (earcut compatibility)
                coords.push_back(mapbox::util::nth<0, decltype(point)>::get(point));
                coords.push_back(mapbox::util::nth<1, decltype(point)>::get(point));
            }
        }
    }
    
    // TODO: Handle holes (additional rings in polygon vector)
    std::vector<N> holes; // Empty for now
    
    return meshcut<N>(coords, holes, options);
}

// Explicit template instantiations for common types
template std::vector<uint32_t> meshcut(const std::vector<double>&, const std::vector<uint32_t>&, const MeshCutOptions&);
template std::vector<uint16_t> meshcut(const std::vector<double>&, const std::vector<uint16_t>&, const MeshCutOptions&);

template MeshCutResult meshcut_full<uint32_t>(const std::vector<double>&, const std::vector<uint32_t>&, const MeshCutOptions&);
template MeshCutResult meshcut_full<uint16_t>(const std::vector<double>&, const std::vector<uint16_t>&, const MeshCutOptions&);

} // namespace meshcut