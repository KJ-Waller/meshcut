#include "meshcut/meshcut.hpp"
#include <cmath>
#include <algorithm>
#include <limits>
#include <set>
#include <map>
#include <unordered_set>
#include <cstdio>
#include "clipper2/clipper.h"
#include "mapbox/earcut.hpp"

namespace meshcut {

std::vector<Triangle> tessellate(const Polygon& polygon,
                                const Polygon& bbox,
                                int gridN,
                                TriangleOrientation orientation) {
    // Step 1: Build regular grid triangles
    auto gridTriangles = internal::buildGridTriangles(bbox, gridN, orientation);
    
    // Step 2: Classify triangles as fully inside or not
    auto isInterior = internal::classifyTriangles(gridTriangles, polygon);
    
    // Step 3: Collect interior triangles
    std::vector<Triangle> interiorTriangles;
    for (size_t i = 0; i < gridTriangles.size(); ++i) {
        if (isInterior[i]) {
            interiorTriangles.push_back(gridTriangles[i]);
        }
    }
    
    // Step 4: Find triangles that partially intersect with polygon
    auto isPartiallyIntersected = internal::findPartiallyIntersectedTriangles(gridTriangles, polygon, isInterior);
    
    // Step 5: For each partially intersected triangle, compute its intersection with polygon
    // and triangulate that intersection individually
    std::vector<Triangle> partialTriangles;
    for (size_t i = 0; i < gridTriangles.size(); ++i) {
        if (isPartiallyIntersected[i]) {
            // Compute intersection of this triangle with the polygon
            auto intersection = internal::computeTrianglePolygonIntersection(gridTriangles[i], polygon);
            
            // If intersection is non-empty, triangulate it
            if (intersection.size() >= 3) {
                auto triangulated = internal::earcutTriangulate(intersection);
                partialTriangles.insert(partialTriangles.end(), triangulated.begin(), triangulated.end());
            }
        }
    }
    
    // Step 6: Combine interior + partial triangles
    std::vector<Triangle> finalMesh;
    finalMesh.reserve(interiorTriangles.size() + partialTriangles.size());
    finalMesh.insert(finalMesh.end(), interiorTriangles.begin(), interiorTriangles.end());
    finalMesh.insert(finalMesh.end(), partialTriangles.begin(), partialTriangles.end());
    
    return finalMesh;
}

namespace internal {

// Helper function to check if two points are approximately equal
bool pointsEqual(const Point& a, const Point& b, double epsilon = 1e-9) {
    return std::abs(a.x - b.x) < epsilon && std::abs(a.y - b.y) < epsilon;
}

// Helper function to find shared vertices between two triangles
std::vector<Point> findSharedVertices(const Triangle& t1, const Triangle& t2) {
    std::vector<Point> shared;
    std::vector<Point> v1 = {t1.a, t1.b, t1.c};
    std::vector<Point> v2 = {t2.a, t2.b, t2.c};
    
    for (const auto& p1 : v1) {
        for (const auto& p2 : v2) {
            if (pointsEqual(p1, p2)) {
                shared.push_back(p1);
                break;
            }
        }
    }
    return shared;
}

// Convert triangles to grid cells (quads) by grouping adjacent triangles
std::vector<Polygon> trianglesToGridCells(const std::vector<Triangle>& triangles) {
    std::vector<Polygon> gridCells;
    std::unordered_set<size_t> usedTriangles;
    
    for (size_t i = 0; i < triangles.size(); ++i) {
        if (usedTriangles.count(i)) continue;
        
        const Triangle& tri1 = triangles[i];
        
        // Try to find a matching triangle to form a quad
        for (size_t j = i + 1; j < triangles.size(); ++j) {
            if (usedTriangles.count(j)) continue;
            
            const Triangle& tri2 = triangles[j];
            auto shared = findSharedVertices(tri1, tri2);
            
            if (shared.size() == 2) {
                // Found a pair! Create a quad from the 4 unique vertices
                std::set<std::pair<double, double>> uniqueVertices;
                uniqueVertices.insert({tri1.a.x, tri1.a.y});
                uniqueVertices.insert({tri1.b.x, tri1.b.y});
                uniqueVertices.insert({tri1.c.x, tri1.c.y});
                uniqueVertices.insert({tri2.a.x, tri2.a.y});
                uniqueVertices.insert({tri2.b.x, tri2.b.y});
                uniqueVertices.insert({tri2.c.x, tri2.c.y});
                
                if (uniqueVertices.size() == 4) {
                    // Convert to Polygon and sort vertices
                    Polygon quad;
                    for (const auto& v : uniqueVertices) {
                        quad.emplace_back(v.first, v.second);
                    }
                    
                    // Sort vertices clockwise around centroid
                    Point centroid(0, 0);
                    for (const auto& p : quad) {
                        centroid.x += p.x;
                        centroid.y += p.y;
                    }
                    centroid.x /= quad.size();
                    centroid.y /= quad.size();
                    
                    std::sort(quad.begin(), quad.end(), [&](const Point& a, const Point& b) {
                        double angleA = std::atan2(a.y - centroid.y, a.x - centroid.x);
                        double angleB = std::atan2(b.y - centroid.y, b.x - centroid.x);
                        return angleA < angleB;
                    });
                    
                    gridCells.push_back(quad);
                    usedTriangles.insert(i);
                    usedTriangles.insert(j);
                    break;
                }
            }
        }
    }
    
    // Add any remaining individual triangles that couldn't be paired
    for (size_t i = 0; i < triangles.size(); ++i) {
        if (!usedTriangles.count(i)) {
            const Triangle& tri = triangles[i];
            Polygon triPoly = {tri.a, tri.b, tri.c};
            gridCells.push_back(triPoly);
        }
    }
    
    return gridCells;
}

// Simple polygon union using axis-aligned bounding box merging
// This is a simplified approach - in a full implementation you'd use Clipper
Polygon computeSimplifiedUnion(const std::vector<Polygon>& polygons) {
    if (polygons.empty()) return {};
    
    // For now, just return the bounding box of all polygons
    // This is a placeholder until proper boolean operations are implemented
    double minX = std::numeric_limits<double>::max();
    double minY = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double maxY = std::numeric_limits<double>::lowest();
    
    for (const auto& poly : polygons) {
        for (const auto& point : poly) {
            minX = std::min(minX, point.x);
            minY = std::min(minY, point.y);
            maxX = std::max(maxX, point.x);
            maxY = std::max(maxY, point.y);
        }
    }
    
    // Return bounding rectangle
    return {
        {minX, minY},
        {maxX, minY}, 
        {maxX, maxY},
        {minX, maxY}
    };
}

std::vector<Triangle> buildGridTriangles(const Polygon& bbox, int gridN, TriangleOrientation orientation) {
    if (bbox.size() < 4) {
        return {}; // Invalid bounding box
    }
    
    // Extract bounding box coordinates
    // Assume bbox = [bottom-left, bottom-right, top-right, top-left]
    double minX = bbox[0].x;
    double minY = bbox[0].y;
    double maxX = bbox[2].x;
    double maxY = bbox[2].y;
    
    double cellWidth = (maxX - minX) / gridN;
    double cellHeight = (maxY - minY) / gridN;
    
    std::vector<Triangle> triangles;
    triangles.reserve(gridN * gridN * 2); // 2 triangles per cell
    
    for (int i = 0; i < gridN; ++i) {
        for (int j = 0; j < gridN; ++j) {
            double x0 = minX + i * cellWidth;
            double y0 = minY + j * cellHeight;
            double x1 = x0 + cellWidth;
            double y1 = y0 + cellHeight;
            
            // Create 4 corners of the cell
            Point p00(x0, y0); // bottom-left
            Point p10(x1, y0); // bottom-right
            Point p11(x1, y1); // top-right
            Point p01(x0, y1); // top-left
            
            // Create 2 triangles per cell based on orientation
            if (orientation == TriangleOrientation::TOPLEFT_BOTTOMRIGHT) {
                // Diagonal from top-left to bottom-right (/)
                // Triangle 1: bottom-left, bottom-right, top-left
                triangles.emplace_back(p00, p10, p01);
                // Triangle 2: bottom-right, top-right, top-left
                triangles.emplace_back(p10, p11, p01);
            } else { // TOPRIGHT_BOTTOMLEFT
                // Diagonal from top-right to bottom-left (\)
                // Triangle 1: bottom-left, top-right, bottom-right
                triangles.emplace_back(p00, p11, p10);
                // Triangle 2: bottom-left, top-left, top-right
                triangles.emplace_back(p00, p01, p11);
            }
        }
    }
    
    return triangles;
}

std::vector<bool> classifyTriangles(const std::vector<Triangle>& triangles, 
                                   const Polygon& polygon) {
    std::vector<bool> isInterior(triangles.size(), false);
    
    for (size_t i = 0; i < triangles.size(); ++i) {
        const auto& tri = triangles[i];
        
        // Check if all vertices are inside polygon
        // Use the original approach: ALL vertices must be inside
        bool allVerticesInside = pointInPolygon(tri.a, polygon) &&
                               pointInPolygon(tri.b, polygon) &&
                               pointInPolygon(tri.c, polygon);
        
        // Check if triangle doesn't intersect polygon edges
        bool noEdgeIntersection = !triangleIntersectsPolygonEdge(tri, polygon);
        
        // Triangle is fully interior if all vertices are inside AND no edge intersections
        isInterior[i] = allVerticesInside && noEdgeIntersection;
    }
    
    return isInterior;
}

bool pointInPolygon(const Point& point, const Polygon& polygon) {
    if (polygon.size() < 3) return false;
    
    const double EPSILON = 1e-6; // Looser tolerance for boundary detection
    size_t n = polygon.size();
    
    // First check if point is exactly on a vertex
    for (const auto& vertex : polygon) {
        if (std::abs(point.x - vertex.x) < EPSILON && std::abs(point.y - vertex.y) < EPSILON) {
            return true; // Point is on vertex, consider it inside
        }
    }
    
    // Check if point is on any edge
    for (size_t i = 0; i < n; ++i) {
        size_t j = (i + 1) % n;
        const Point& p1 = polygon[i];
        const Point& p2 = polygon[j];
        
        // Vector from p1 to p2
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        
        // Vector from p1 to point
        double px = point.x - p1.x;
        double py = point.y - p1.y;
        
        // Check if point lies on the line segment
        // First check if vectors are parallel (cross product = 0)
        double cross = dx * py - dy * px;
        if (std::abs(cross) < EPSILON) {
            // Point is on the line, now check if it's within the segment
            double dot = px * dx + py * dy;
            double len_sq = dx * dx + dy * dy;
            
            if (dot >= -EPSILON && dot <= len_sq + EPSILON) {
                return true; // Point is on edge, consider it inside
            }
        }
    }
    
    // Standard ray casting algorithm for interior points
    bool inside = false;
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        const Point& pi = polygon[i];
        const Point& pj = polygon[j];
        
        if (((pi.y > point.y) != (pj.y > point.y)) &&
            (point.x < (pj.x - pi.x) * (point.y - pi.y) / (pj.y - pi.y) + pi.x)) {
            inside = !inside;
        }
    }
    
    return inside;
}

bool triangleIntersectsPolygonEdge(const Triangle& triangle, const Polygon& polygon) {
    if (polygon.size() < 3) return false;
    
    // Get triangle edges
    std::vector<std::pair<Point, Point>> triangleEdges = {
        {triangle.a, triangle.b},
        {triangle.b, triangle.c},
        {triangle.c, triangle.a}
    };
    
    // Check intersection with each polygon edge
    size_t n = polygon.size();
    for (size_t i = 0; i < n; ++i) {
        Point polygonStart = polygon[i];
        Point polygonEnd = polygon[(i + 1) % n];
        
        // Check if this polygon edge intersects with any triangle edge
        for (const auto& triangleEdge : triangleEdges) {
            if (lineSegmentsIntersect(triangleEdge.first, triangleEdge.second, 
                                    polygonStart, polygonEnd)) {
                return true;
            }
        }
    }
    
    return false;
}

// Helper function to check if two line segments intersect
bool lineSegmentsIntersect(const Point& p1, const Point& q1, const Point& p2, const Point& q2) {
    auto orientation = [](const Point& p, const Point& q, const Point& r) -> int {
        double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
        if (std::abs(val) < 1e-10) return 0; // collinear
        return (val > 0) ? 1 : 2; // clock or counterclock wise
    };
    
    auto onSegment = [](const Point& p, const Point& q, const Point& r) -> bool {
        return q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x) &&
               q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y);
    };
    
    int o1 = orientation(p1, q1, p2);
    int o2 = orientation(p1, q1, q2);
    int o3 = orientation(p2, q2, p1);
    int o4 = orientation(p2, q2, q1);
    
    // General case
    if (o1 != o2 && o3 != o4) return true;
    
    // Special cases
    // p1, q1 and p2 are collinear and p2 lies on segment p1q1
    if (o1 == 0 && onSegment(p1, p2, q1)) return true;
    // p1, q1 and q2 are collinear and q2 lies on segment p1q1
    if (o2 == 0 && onSegment(p1, q2, q1)) return true;
    // p2, q2 and p1 are collinear and p1 lies on segment p2q2
    if (o3 == 0 && onSegment(p2, p1, q2)) return true;
    // p2, q2 and q1 are collinear and q1 lies on segment p2q2
    if (o4 == 0 && onSegment(p2, q1, q2)) return true;
    
    return false;
}

Polygon computeRemainder(const Polygon& originalPolygon,
                        const std::vector<Triangle>& interiorTriangles) {
    if (interiorTriangles.empty()) {
        return originalPolygon;
    }
    
    // Convert original polygon to Clipper2 format (subject)
    Clipper2Lib::PathsD subject;
    Clipper2Lib::PathD subjectPath;
    for (const auto& point : originalPolygon) {
        subjectPath.push_back(Clipper2Lib::PointD(point.x, point.y));
    }
    subject.push_back(subjectPath);
    
    // Convert interior triangles to Clipper2 format (clips)
    Clipper2Lib::PathsD clips;
    for (const auto& triangle : interiorTriangles) {
        Clipper2Lib::PathD trianglePath;
        trianglePath.push_back(Clipper2Lib::PointD(triangle.a.x, triangle.a.y));
        trianglePath.push_back(Clipper2Lib::PointD(triangle.b.x, triangle.b.y));
        trianglePath.push_back(Clipper2Lib::PointD(triangle.c.x, triangle.c.y));
        clips.push_back(trianglePath);
    }
    
    // Use the high-level Difference function - this is the correct way!
    Clipper2Lib::PathsD solution = Clipper2Lib::Difference(subject, clips, Clipper2Lib::FillRule::NonZero);
    
    if (solution.empty()) {
        // Empty solution means all interior triangles covered the polygon completely
        // Return empty polygon (no remainder)
        return {};
    }
    
    // Convert the largest result polygon back to our format
    double maxArea = 0.0;
    size_t maxIndex = 0;
    for (size_t i = 0; i < solution.size(); ++i) {
        double area = std::abs(Clipper2Lib::Area(solution[i]));
        if (area > maxArea) {
            maxArea = area;
            maxIndex = i;
        }
    }
    
    // If the largest remaining area is very small, consider it numerical error
    // and return empty polygon (complete coverage)
    const double AREA_TOLERANCE = 1e-14;
    if (maxArea < AREA_TOLERANCE) {
        return {};
    }
    
    // Convert back to our Polygon format
    Polygon result;
    for (const auto& point : solution[maxIndex]) {
        result.push_back(Point(point.x, point.y));
    }
    
    return result;
}

std::vector<Triangle> earcutTriangulate(const Polygon& polygon) {
    std::vector<Triangle> triangles;
    
    if (polygon.size() < 3) return triangles;
    
    // Convert polygon to Earcut format (currently only supports outer ring)
    // TODO: Add support for holes by modifying this function to accept holes parameter
    using Coord = double;
    using N = uint32_t;
    std::vector<std::vector<std::array<Coord, 2>>> earcut_polygon;
    
    // Outer ring
    std::vector<std::array<Coord, 2>> outer_ring;
    for (const auto& point : polygon) {
        outer_ring.push_back({point.x, point.y});
    }
    earcut_polygon.push_back(outer_ring);
    
    // Run Earcut triangulation
    std::vector<N> indices = mapbox::earcut<N>(earcut_polygon);
    
    // Convert indices back to triangles
    for (size_t i = 0; i < indices.size(); i += 3) {
        if (i + 2 < indices.size()) {
            Point a(polygon[indices[i]].x, polygon[indices[i]].y);
            Point b(polygon[indices[i + 1]].x, polygon[indices[i + 1]].y);
            Point c(polygon[indices[i + 2]].x, polygon[indices[i + 2]].y);
            triangles.emplace_back(a, b, c);
        }
    }
    
    return triangles;
}

// Add new function that supports holes
std::vector<Triangle> earcutTriangulateWithHoles(const Polygon& exterior, 
                                                 const std::vector<Polygon>& holes) {
    std::vector<Triangle> triangles;
    
    if (exterior.size() < 3) return triangles;
    
    using Coord = double;
    using N = uint32_t;
    std::vector<std::vector<std::array<Coord, 2>>> earcut_polygon;
    
    // Add exterior ring
    std::vector<std::array<Coord, 2>> outer_ring;
    for (const auto& point : exterior) {
        outer_ring.push_back({point.x, point.y});
    }
    earcut_polygon.push_back(outer_ring);
    
    // Add hole rings
    for (const auto& hole : holes) {
        std::vector<std::array<Coord, 2>> hole_ring;
        for (const auto& point : hole) {
            hole_ring.push_back({point.x, point.y});
        }
        earcut_polygon.push_back(hole_ring);
    }
    
    // Run Earcut triangulation with holes
    std::vector<N> indices = mapbox::earcut<N>(earcut_polygon);
    
    // Convert indices back to triangles
    // Note: indices reference all vertices in order: exterior first, then holes
    std::vector<Point> all_vertices;
    for (const auto& point : exterior) {
        all_vertices.emplace_back(point.x, point.y);
    }
    for (const auto& hole : holes) {
        for (const auto& point : hole) {
            all_vertices.emplace_back(point.x, point.y);
        }
    }
    
    for (size_t i = 0; i < indices.size(); i += 3) {
        if (i + 2 < indices.size()) {
            if (indices[i] < all_vertices.size() && 
                indices[i + 1] < all_vertices.size() && 
                indices[i + 2] < all_vertices.size()) {
                triangles.emplace_back(all_vertices[indices[i]], 
                                     all_vertices[indices[i + 1]], 
                                     all_vertices[indices[i + 2]]);
            }
        }
    }
    
    return triangles;
}

std::vector<bool> findPartiallyIntersectedTriangles(const std::vector<Triangle>& triangles,
                                                   const Polygon& polygon,
                                                   const std::vector<bool>& isInterior) {
    std::vector<bool> isPartiallyIntersected(triangles.size(), false);
    
    for (size_t i = 0; i < triangles.size(); ++i) {
        // Skip if triangle is fully interior
        if (isInterior[i]) {
            continue;
        }
        
        const auto& tri = triangles[i];
        
        // Check if triangle has any intersection with polygon
        // A triangle partially intersects if:
        // 1. At least one vertex is inside the polygon, OR
        // 2. The triangle intersects with polygon edges
        
        bool hasVertexInside = pointInPolygon(tri.a, polygon) ||
                              pointInPolygon(tri.b, polygon) ||
                              pointInPolygon(tri.c, polygon);
        
        bool intersectsEdges = triangleIntersectsPolygonEdge(tri, polygon);
        
        // Also check if any polygon vertex is inside the triangle
        bool hasPolygonVertexInside = false;
        for (const auto& polyPoint : polygon) {
            if (pointInTriangle(polyPoint, tri)) {
                hasPolygonVertexInside = true;
                break;
            }
        }
        
        isPartiallyIntersected[i] = hasVertexInside || intersectsEdges || hasPolygonVertexInside;
    }
    
    return isPartiallyIntersected;
}

Polygon computeTrianglePolygonIntersection(const Triangle& triangle, const Polygon& polygon) {
    // Convert triangle to Clipper2 format
    Clipper2Lib::PathD trianglePath;
    trianglePath.push_back(Clipper2Lib::PointD(triangle.a.x, triangle.a.y));
    trianglePath.push_back(Clipper2Lib::PointD(triangle.b.x, triangle.b.y));
    trianglePath.push_back(Clipper2Lib::PointD(triangle.c.x, triangle.c.y));
    
    Clipper2Lib::PathsD trianglePaths;
    trianglePaths.push_back(trianglePath);
    
    // Convert polygon to Clipper2 format
    Clipper2Lib::PathD polygonPath;
    for (const auto& point : polygon) {
        polygonPath.push_back(Clipper2Lib::PointD(point.x, point.y));
    }
    
    Clipper2Lib::PathsD polygonPaths;
    polygonPaths.push_back(polygonPath);
    
    // Compute intersection using Clipper2
    Clipper2Lib::PathsD solution = Clipper2Lib::Intersect(trianglePaths, polygonPaths, Clipper2Lib::FillRule::NonZero);
    
    if (solution.empty()) {
        return {}; // No intersection
    }
    
    // Return the largest intersection polygon (in case there are multiple)
    double maxArea = 0.0;
    size_t maxIndex = 0;
    for (size_t i = 0; i < solution.size(); ++i) {
        double area = std::abs(Clipper2Lib::Area(solution[i]));
        if (area > maxArea) {
            maxArea = area;
            maxIndex = i;
        }
    }
    
    // Convert back to our Polygon format
    Polygon result;
    for (const auto& point : solution[maxIndex]) {
        result.push_back(Point(point.x, point.y));
    }
    
    return result;
}

// Helper function to check if a point is inside a triangle
bool pointInTriangle(const Point& point, const Triangle& triangle) {
    // Use barycentric coordinates to check if point is inside triangle
    const Point& a = triangle.a;
    const Point& b = triangle.b;
    const Point& c = triangle.c;
    const Point& p = point;
    
    double denom = (b.y - c.y) * (a.x - c.x) + (c.x - b.x) * (a.y - c.y);
    if (std::abs(denom) < 1e-12) return false; // Degenerate triangle
    
    double u = ((b.y - c.y) * (p.x - c.x) + (c.x - b.x) * (p.y - c.y)) / denom;
    double v = ((c.y - a.y) * (p.x - c.x) + (a.x - c.x) * (p.y - c.y)) / denom;
    double w = 1.0 - u - v;
    
    const double EPSILON = 1e-9;
    return (u >= -EPSILON) && (v >= -EPSILON) && (w >= -EPSILON);
}

} // namespace internal

std::vector<Triangle> tessellate_tile(const Polygon& polygon,
                                     uint32_t tileSubdivisions,
                                     TriangleOrientation orientation) {
    // Create tile bounding box (0,0) to (1,1)
    Polygon tileBbox = {
        {0.0, 0.0},  // bottom-left
        {1.0, 0.0},  // bottom-right
        {1.0, 1.0},  // top-right
        {0.0, 1.0}   // top-left
    };
    
    return tessellate(polygon, tileBbox, static_cast<int>(tileSubdivisions), orientation);
}

} // namespace meshcut

// C wrapper for easier Python integration
extern "C" {

int tessellate_simple(
    const double* polygon_points,  // x,y pairs
    int num_points,
    double bbox_min_x, double bbox_min_y, double bbox_max_x, double bbox_max_y,
    int grid_size,
    int orientation,  // 0 = TL-BR, 1 = TR-BL
    double* output_triangles,  // x1,y1,x2,y2,x3,y3 per triangle
    int* triangle_count
) {
    try {
        // Convert input polygon
        meshcut::Polygon polygon;
        for (int i = 0; i < num_points; ++i) {
            polygon.push_back({polygon_points[i*2], polygon_points[i*2+1]});
        }
        
        // Create bounding box
        meshcut::Polygon bbox = {
            {bbox_min_x, bbox_min_y},
            {bbox_max_x, bbox_min_y}, 
            {bbox_max_x, bbox_max_y},
            {bbox_min_x, bbox_max_y}
        };
        
        // Convert orientation
        meshcut::TriangleOrientation orient = (orientation == 0) ? 
            meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT :
            meshcut::TriangleOrientation::TOPRIGHT_BOTTOMLEFT;
        
        // Tessellate
        auto triangles = meshcut::tessellate(polygon, bbox, grid_size, orient);
        
        // Convert output
        *triangle_count = static_cast<int>(triangles.size());
        for (size_t i = 0; i < triangles.size(); ++i) {
            int base = i * 6;
            output_triangles[base + 0] = triangles[i].a.x;
            output_triangles[base + 1] = triangles[i].a.y;
            output_triangles[base + 2] = triangles[i].b.x;
            output_triangles[base + 3] = triangles[i].b.y;
            output_triangles[base + 4] = triangles[i].c.x;
            output_triangles[base + 5] = triangles[i].c.y;
        }
        
        return 0;  // Success
    } catch (...) {
        return -1;  // Error
    }
}

}
