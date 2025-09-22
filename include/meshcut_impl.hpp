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
 * Simple point-in-polygon test (placeholder - will optimize with scanline later)
 */
bool pointInPolygon(double x, double y, const std::vector<double>& polygon) {
    // Ray casting algorithm - count intersections with polygon edges
    int intersections = 0;
    int n = polygon.size() / 2;
    
    for (int i = 0; i < n; i++) {
        int j = (i + 1) % n;
        double x1 = polygon[i * 2];
        double y1 = polygon[i * 2 + 1];  
        double x2 = polygon[j * 2];
        double y2 = polygon[j * 2 + 1];
        
        // Check if ray crosses edge
        if (((y1 > y) != (y2 > y)) && 
            (x < (x2 - x1) * (y - y1) / (y2 - y1) + x1)) {
            intersections++;
        }
    }
    
    return (intersections % 2) == 1;
}

/**
 * Classify a grid cell based on polygon intersection
 */
CellState classifyCell(int i, int j, const std::vector<double>& polygon, const CoordinateTransform& transform) {
    // Get the 4 corners of the cell in world coordinates
    auto [x0, y0] = transform.gridToWorld(i, j);
    auto [x1, y1] = transform.gridToWorld(i + 1, j);
    auto [x2, y2] = transform.gridToWorld(i + 1, j + 1);
    auto [x3, y3] = transform.gridToWorld(i, j + 1);
    
    // Test all 4 corners
    bool in0 = pointInPolygon(x0, y0, polygon);
    bool in1 = pointInPolygon(x1, y1, polygon);
    bool in2 = pointInPolygon(x2, y2, polygon);
    bool in3 = pointInPolygon(x3, y3, polygon);
    
    // All corners inside = fully inside
    if (in0 && in1 && in2 && in3) {
        return FULL_IN;
    }
    
    // No corners inside - could still intersect, but for now assume fully outside
    // TODO: Add proper edge-cell intersection test
    if (!in0 && !in1 && !in2 && !in3) {
        return FULL_OUT;
    }
    
    // Mixed or edge case = boundary
    return BOUNDARY;
}

/**
 * Generate two triangles for a grid cell (FULL_IN case)
 */
template<typename N>
void emitGridTriangles(int i, int j, bool diagonalNE, 
                      std::vector<N>& indices,
                      std::unordered_map<GridPoint, N, GridPointHash>& vertexMap,
                      N& nextVertexIndex) {
    // Get or create vertex indices for the 4 corners
    GridPoint corners[4] = {{i, j}, {i+1, j}, {i+1, j+1}, {i, j+1}};
    N vertexIndices[4];
    
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
 * Process a boundary cell - clip polygon and triangulate
 * (Placeholder - will implement Sutherland-Hodgman clipping later)
 */
template<typename N>
void processBoundaryCell(int i, int j, const std::vector<double>& polygon, 
                        const CoordinateTransform& transform,
                        std::vector<N>& indices,
                        std::unordered_map<GridPoint, N, GridPointHash>& vertexMap,
                        N& nextVertexIndex) {
    // TODO: Implement Sutherland-Hodgman clipping to cell rectangle
    // TODO: Run earcut on clipped polygon
    // TODO: Add resulting triangles to indices
    
    // For now, just skip boundary cells (placeholder)
    // This will be implemented in the next phase
}

} // namespace detail

/**
 * Main meshcut implementation
 */
template <typename N>
std::vector<N> meshcut(const std::vector<double>& polygon, const std::vector<N>& holes, const MeshCutOptions& options) {
    
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
    std::unordered_map<detail::GridPoint, N, detail::GridPointHash> vertexMap;
    N nextVertexIndex = 0;
    
    // Process each grid cell
    for (int j = minJ; j <= maxJ; j++) {
        for (int i = minI; i <= maxI; i++) {
            detail::CellState state = detail::classifyCell(i, j, polygon, transform);
            
            switch (state) {
                case detail::FULL_IN:
                    detail::emitGridTriangles(i, j, options.diagonalNE, indices, vertexMap, nextVertexIndex);
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
    
    return indices;
}

/**
 * Template specialization for custom polygon types (earcut compatibility)
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

} // namespace meshcut