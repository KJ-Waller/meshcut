#pragma once

#include <vector>
#include <array>
#include <cstdint>

namespace meshcut {

enum class TriangleOrientation {
    TOPLEFT_BOTTOMRIGHT,  // Default: diagonal from top-left to bottom-right (\)
    TOPRIGHT_BOTTOMLEFT   // Alternative: diagonal from top-right to bottom-left (/)
};

struct Point {
    double x, y;
    
    Point() : x(0.0), y(0.0) {}
    Point(double x_, double y_) : x(x_), y(y_) {}
};

struct Triangle {
    Point a, b, c;
    
    Triangle() = default;
    Triangle(const Point& a_, const Point& b_, const Point& c_) : a(a_), b(b_), c(c_) {}
};

using Polygon = std::vector<Point>;

/**
 * Main API function: Generate mesh of triangles covering a polygon using grid-based approach
 * 
 * @param polygon The polygon to fill (assumes simple polygon for now)
 * @param bbox Bounding box as 4 points [bottom-left, bottom-right, top-right, top-left]
 * @param gridN Grid resolution (creates gridN x gridN cells)
 * @param orientation Triangle orientation within grid cells
 * @return Vector of triangles covering the polygon
 */
std::vector<Triangle> tessellate(const Polygon& polygon,
                                const Polygon& bbox,
                                int gridN,
                                TriangleOrientation orientation = TriangleOrientation::TOPLEFT_BOTTOMRIGHT);

/**
 * Utility function for MapLibre tile coordinates
 * Tessellate polygon on a tile grid with specified subdivisions
 * 
 * @param polygon The polygon to tessellate (in tile coordinates 0-1)
 * @param tileSubdivisions Number of subdivisions per tile edge (typically 32 for 32x32)
 * @param orientation Triangle orientation within grid cells  
 * @return Vector of triangles covering the polygon
 */
std::vector<Triangle> tessellate_tile(const Polygon& polygon,
                                     uint32_t tileSubdivisions = 32,
                                     TriangleOrientation orientation = TriangleOrientation::TOPLEFT_BOTTOMRIGHT);

// Internal functions (exposed for advanced usage)
namespace internal {
    /**
     * Build regular grid triangles within bounding box
     */
    std::vector<Triangle> buildGridTriangles(const Polygon& bbox, int gridN, 
                                           TriangleOrientation orientation = TriangleOrientation::TOPLEFT_BOTTOMRIGHT);
    
    /**
     * Classify triangles as fully inside polygon or not
     */
    std::vector<bool> classifyTriangles(const std::vector<Triangle>& triangles, 
                                       const Polygon& polygon);
    
    /**
     * Check if point is inside polygon (ray casting algorithm)
     */
    bool pointInPolygon(const Point& point, const Polygon& polygon);
    
    /**
     * Check if triangle intersects with any polygon edge
     */
    bool triangleIntersectsPolygonEdge(const Triangle& triangle, const Polygon& polygon);
    
    /**
     * Find triangles that partially intersect with polygon (not fully inside, but do intersect)
     */
    std::vector<bool> findPartiallyIntersectedTriangles(const std::vector<Triangle>& triangles,
                                                        const Polygon& polygon,
                                                        const std::vector<bool>& isInterior);
    
    /**
     * Compute intersection of a triangle with a polygon
     * Returns the polygon representing the intersection (empty if no intersection)
     */
    Polygon computeTrianglePolygonIntersection(const Triangle& triangle, const Polygon& polygon);
    
    /**
     * Check if a point is inside a triangle using barycentric coordinates
     */
    bool pointInTriangle(const Point& point, const Triangle& triangle);
    
    /**
     * Check if two line segments intersect
     */
    bool lineSegmentsIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2);
    
    /**
     * Triangulate polygon using Earcut algorithm
     */
    std::vector<Triangle> earcutTriangulate(const Polygon& polygon);
}

} // namespace meshcut

// C wrapper for easier Python integration
#ifdef __cplusplus
extern "C" {
#endif

/**
 * Simple C interface for tessellation
 * Returns 0 on success, -1 on error
 */
int tessellate_simple(
    const double* polygon_points,  // x,y pairs
    int num_points,
    double bbox_min_x, double bbox_min_y, double bbox_max_x, double bbox_max_y,
    int grid_size,
    int orientation,  // 0 = TOPLEFT_BOTTOMRIGHT (\), 1 = TOPRIGHT_BOTTOMLEFT (/)
    double* output_triangles,  // x1,y1,x2,y2,x3,y3 per triangle
    int* triangle_count
);

#ifdef __cplusplus
}
#endif