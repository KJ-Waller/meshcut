#pragma once

#include <vector>
#include <cstdint>
#include <array>

namespace meshcut {

/**
 * MeshCut Options - specifies the grid for mesh-based triangulation
 */
struct MeshCutOptions {
    // Grid specification in world coordinates
    double gridOriginX = 0.0;      // Bottom-left corner X of grid in world coordinates
    double gridOriginY = 0.0;      // Bottom-left corner Y of grid in world coordinates
    double cellSize = 1.0;         // Size of each grid cell in world units
    int gridWidth = 32;            // Number of cells in X direction
    int gridHeight = 32;           // Number of cells in Y direction
    bool diagonalNE = true;        // true for / diagonal, false for \ diagonal
};

/**
 * MeshCut - Main triangulation function
 * 
 * Compatible with earcut API but with added mesh constraints.
 * Takes the same polygon format as earcut: vector of coordinate pairs.
 * 
 * @param polygon Vector of coordinates [x1,y1, x2,y2, ...] defining the polygon
 * @param holes Vector of indices where holes start (earcut format)
 * @param options Grid specification for mesh-constrained triangulation  
 * @return Vector of triangle indices referring to input vertices (+ any new intersection vertices)
 */
template <typename N = uint32_t>
std::vector<N> meshcut(
    const std::vector<double>& polygon,
    const std::vector<N>& holes = {},
    const MeshCutOptions& options = {}
);

/**
 * MeshCut with custom polygon format (compatible with earcut template usage)
 * 
 * @param polygon Custom polygon format (same as earcut - vector of rings)
 * @param options Grid specification for mesh-constrained triangulation
 * @return Vector of triangle indices
 */
template <typename N = uint32_t, typename Polygon>
std::vector<N> meshcut(
    const Polygon& polygon,
    const MeshCutOptions& options = {}
);

} // namespace meshcut

// Include implementation (header-only like earcut)
#include "meshcut_impl.hpp"