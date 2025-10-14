#pragma once

#include "earcut.hpp"

#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <tuple>
#include <type_traits>
#include <map>

namespace meshcut {

// Forward declarations for internal functions
namespace detail {

/**
 * Thread-local buffer pool for reusing memory allocations
 * Eliminates per-cell allocation overhead in boundary processing
 */
class BufferPool {
private:
    // Clipping operation buffers
    mutable std::vector<std::pair<double, double>> pointBuffer1;
    mutable std::vector<std::pair<double, double>> pointBuffer2;
    mutable std::vector<std::pair<double, double>> pointBuffer3;
    mutable std::vector<std::pair<double, double>> pointBuffer4;
    mutable std::vector<std::pair<double, double>> pointBuffer5;
    
    // Earcut triangulation buffers
    mutable std::vector<std::vector<std::pair<double, double>>> earcutPolygonBuffer;
    mutable std::vector<uint32_t> triangleBuffer;
    mutable std::vector<uint32_t> globalIndicesBuffer;
    
    // Grid processing buffers
    mutable std::vector<int> gridRowBuffer1;
    mutable std::vector<int> gridRowBuffer2;
    mutable std::vector<double> intersectionBuffer;
    
public:
    BufferPool() {
        // Pre-allocate reasonable sizes to avoid small reallocations
        pointBuffer1.reserve(32);
        pointBuffer2.reserve(32);
        pointBuffer3.reserve(32);
        pointBuffer4.reserve(32);
        pointBuffer5.reserve(32);
        
        earcutPolygonBuffer.reserve(2);
        triangleBuffer.reserve(96);
        globalIndicesBuffer.reserve(32);
        
        gridRowBuffer1.reserve(256);
        gridRowBuffer2.reserve(256);
        intersectionBuffer.reserve(64);
    }
    
    // Get buffers for clipping operations (5 stages: input + 4 clipping stages)
    std::vector<std::pair<double, double>>& getPointBuffer(int stage) const {
        switch (stage) {
            case 0: return pointBuffer1;
            case 1: return pointBuffer2;
            case 2: return pointBuffer3;
            case 3: return pointBuffer4;
            case 4: return pointBuffer5;
            default: return pointBuffer1;
        }
    }
    
    // Get buffers for earcut triangulation
    std::vector<std::vector<std::pair<double, double>>>& getEarcutPolygonBuffer() const {
        earcutPolygonBuffer.clear();
        return earcutPolygonBuffer;
    }
    
    std::vector<uint32_t>& getTriangleBuffer() const {
        triangleBuffer.clear();
        return triangleBuffer;
    }
    
    std::vector<uint32_t>& getGlobalIndicesBuffer() const {
        globalIndicesBuffer.clear();
        return globalIndicesBuffer;
    }
    
    // Get buffers for grid row processing
    std::vector<int>& getGridRowBuffer1() const {
        gridRowBuffer1.clear();
        return gridRowBuffer1;
    }
    
    std::vector<int>& getGridRowBuffer2() const {
        gridRowBuffer2.clear();
        return gridRowBuffer2;
    }
    
    std::vector<double>& getIntersectionBuffer() const {
        intersectionBuffer.clear();
        return intersectionBuffer;
    }
};

/**
 * Thread-local buffer pool instance
 */
thread_local BufferPool g_bufferPool;

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
    std::size_t operator()(const WorldPoint& p) const noexcept {
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
 * Integer coordinate system for fast arithmetic operations
 * All coordinates are scaled to integer space for better performance
 */
class IntegerCoordinateTransform {
public:
    // Constructor from MeshCutOptions
    IntegerCoordinateTransform(const MeshCutOptions& options)
        : originX(options.gridOriginX)
        , originY(options.gridOriginY) 
        , cellSize(options.cellSize)
        , gridWidth(options.gridWidth)
        , gridHeight(options.gridHeight)
    {
        // Dynamic scale factor: ~100 integer units per grid cell regardless of grid size
        // This provides good precision for any grid resolution
        int maxDimension = std::max(gridWidth, gridHeight);
        scaleFactor = maxDimension * 100; // 100 units per cell
        invScaleFactor = 1.0 / scaleFactor;
    }
    
    // Constructor with explicit grid dimensions and bounds
    IntegerCoordinateTransform(int gw, int gh, double minX, double minY, double maxX, double maxY)
        : originX(minX)
        , originY(minY)
        , gridWidth(gw)
        , gridHeight(gh)
    {
        cellSize = (maxX - minX) / gw;
        
        // Dynamic scale factor: ~100 integer units per grid cell regardless of grid size
        int maxDimension = std::max(gridWidth, gridHeight);
        scaleFactor = maxDimension * 100; // 100 units per cell
        invScaleFactor = 1.0 / scaleFactor;
    }
    
    // Convert world coordinates to integer grid space
    std::pair<int64_t, int64_t> worldToInteger(double x, double y) const {
        // Scale relative to grid origin, then multiply by scale factor
        double relativeX = (x - originX) / cellSize;
        double relativeY = (y - originY) / cellSize;
        
        return {
            static_cast<int64_t>(std::round(relativeX * scaleFactor)),
            static_cast<int64_t>(std::round(relativeY * scaleFactor))
        };
    }
    
    // Convert integer coordinates back to world space
    std::pair<double, double> integerToWorld(int64_t ix, int64_t iy) const {
        double relativeX = ix * invScaleFactor;
        double relativeY = iy * invScaleFactor;
        
        return {
            originX + relativeX * cellSize,
            originY + relativeY * cellSize
        };
    }
    
    // Get integer coordinates for grid cell corner
    std::pair<int64_t, int64_t> gridCellToInteger(int gx, int gy) const {
        return {
            static_cast<int64_t>(gx) * scaleFactor,
            static_cast<int64_t>(gy) * scaleFactor
        };
    }
    
    // Get scale factor for calculations
    int getScaleFactor() const {
        return scaleFactor;
    }
    
    // Convert back to floating-point coordinate transform for compatibility
    std::pair<double, double> gridToWorld(int gx, int gy) const {
        return {
            originX + gx * cellSize,
            originY + gy * cellSize
        };
    }
    
private:
    double originX, originY, cellSize;
    int gridWidth, gridHeight;
    int scaleFactor;
    double invScaleFactor;
};

/**
 * Integer-space edge representation for fast scanline processing
 */
struct IntegerPolygonEdge {
    int64_t x1, y1, x2, y2;  // Edge endpoints in integer space
    int64_t minY, maxY;      // Y bounds for quick culling
    
    IntegerPolygonEdge(int64_t x1_, int64_t y1_, int64_t x2_, int64_t y2_) 
        : x1(x1_), y1(y1_), x2(x2_), y2(y2_)
        , minY(std::min(y1_, y2_)), maxY(std::max(y1_, y2_)) {}
    
    // Compute X intersection at given Y using integer arithmetic
    int64_t getXAtY(int64_t y) const {
        if (y2 == y1) return x1; // Horizontal edge
        
        // Integer-based intersection calculation
        // x = x1 + (x2 - x1) * (y - y1) / (y2 - y1)
        int64_t numerator = (x2 - x1) * (y - y1);
        int64_t denominator = (y2 - y1);
        
        return x1 + numerator / denominator;
    }
    
    // Check if edge intersects horizontal ray at (x, y) going right
    bool intersectsRay(int64_t x, int64_t y) const {
        // Quick Y bounds check
        if (y < minY || y > maxY) return false;
        
        // Handle edge cases - endpoint exactly on ray
        if (y == y1) {
            // Only count if other endpoint is above ray
            return y2 > y;
        }
        if (y == y2) {
            // Only count if other endpoint is above ray  
            return y1 > y;
        }
        
        // Compute intersection X using integer arithmetic
        int64_t intersectX = getXAtY(y);
        return intersectX > x;
    }
    
    // Check if point is exactly on this edge (including endpoints)
    bool isPointOnEdge(int64_t x, int64_t y, int64_t tolerance = 100) const {
        // Check if point is exactly on an endpoint
        if ((std::abs(x - x1) <= tolerance && std::abs(y - y1) <= tolerance) ||
            (std::abs(x - x2) <= tolerance && std::abs(y - y2) <= tolerance)) {
            return true;
        }
        
        // Check if point is on the edge line segment
        // First check if point is within bounding box
        int64_t minX = std::min(x1, x2), maxX = std::max(x1, x2);
        if (x < minX - tolerance || x > maxX + tolerance || y < minY - tolerance || y > maxY + tolerance) {
            return false;
        }
        
        // Check if point lies on the line segment using cross product
        int64_t dx1 = x - x1, dy1 = y - y1;
        int64_t dx2 = x2 - x1, dy2 = y2 - y1;
        
        // Cross product should be near zero if point is on line
        int64_t cross = dx1 * dy2 - dy1 * dx2;
        if (std::abs(cross) > tolerance * std::max(std::abs(dx2), std::abs(dy2))) return false;
        
        // Check if point is within the segment bounds using dot product
        int64_t dot = dx1 * dx2 + dy1 * dy2;
        int64_t lenSq = dx2 * dx2 + dy2 * dy2;
        
        return dot >= -tolerance && dot <= lenSq + tolerance;
    }
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
    
    // Get cell size for calculations
    double getCellSize() const {
        return cellSize;
    }
    
    // Get coordinate system bounds (needed for integer coordinate transforms)
    struct Bounds {
        double minX, minY, maxX, maxY;
    };
    
    Bounds getBounds() const {
        // For now, return a reasonable default - could be enhanced to track actual bounds
        return {originX, originY, originX + 1000 * cellSize, originY + 1000 * cellSize};
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
    
    // Check if point is exactly on this edge (including endpoints)
    bool isPointOnEdge(double x, double y, double tolerance = 1e-10) const {
        // Check if point is exactly on an endpoint
        if ((std::abs(x - x1) < tolerance && std::abs(y - y1) < tolerance) ||
            (std::abs(x - x2) < tolerance && std::abs(y - y2) < tolerance)) {
            return true;
        }
        
        // Check if point is on the edge line segment
        // First check if point is within bounding box
        double minX = std::min(x1, x2), maxX = std::max(x1, x2);
        if (x < minX - tolerance || x > maxX + tolerance || y < minY - tolerance || y > maxY + tolerance) {
            return false;
        }
        
        // Check if point lies on the line segment using cross product
        double dx1 = x - x1, dy1 = y - y1;
        double dx2 = x2 - x1, dy2 = y2 - y1;
        
        // Cross product should be near zero if point is on line
        double cross = dx1 * dy2 - dy1 * dx2;
        if (std::abs(cross) > tolerance) return false;
        
        // Check if point is within the segment bounds using dot product
        double dot = dx1 * dx2 + dy1 * dy2;
        double lenSq = dx2 * dx2 + dy2 * dy2;
        
        return dot >= -tolerance && dot <= lenSq + tolerance;
    }
};

/**
 * Integer-space scanline-based point-in-polygon tester
 * All operations in integer arithmetic for maximum performance
 */
class IntegerScanlinePolygonTester {
private:
    std::vector<IntegerPolygonEdge> edges;
    
    // Y-bucketed edge tables for scanline algorithm
    struct IntegerActiveEdge {
        int64_t x;           // Current X intersection
        int64_t deltaX;      // X increment per Y step (scaled)
        int64_t deltaY;      // Y step size
        int64_t maxY;        // Y coordinate where edge ends
        
        IntegerActiveEdge(const IntegerPolygonEdge& edge, int64_t y) {
            if (edge.y2 == edge.y1) {
                // Horizontal edge - skip
                deltaX = 0;
                deltaY = 1;
                x = edge.x1;
                maxY = edge.maxY;
            } else {
                deltaY = edge.y2 - edge.y1;
                deltaX = edge.x2 - edge.x1;  // Will divide by deltaY when stepping
                x = edge.getXAtY(y);
                maxY = edge.maxY;
            }
        }
        
        void stepY(int64_t dy) {
            if (deltaY != 0) {
                x += (deltaX * dy) / deltaY;
            }
        }
    };
    
    // Y-coordinate to edges that start at that Y
    std::map<int64_t, std::vector<int>> edgeStartTable;
    int64_t minY, maxY;
    
    // Cache for grid row processing
    mutable std::vector<IntegerActiveEdge> activeEdges;
    
    void buildYBuckets() {
        if (edges.empty()) return;
        
        minY = edges[0].minY;
        maxY = edges[0].maxY;
        
        for (size_t i = 0; i < edges.size(); i++) {
            const auto& edge = edges[i];
            minY = std::min(minY, edge.minY);
            maxY = std::max(maxY, edge.maxY);
            
            // Add edge to bucket at its minimum Y coordinate
            edgeStartTable[edge.minY].push_back(i);
        }
    }
    
public:
    IntegerScanlinePolygonTester(const std::vector<double>& polygon, const IntegerCoordinateTransform& transform) {
        if (polygon.size() < 6) return;
        
        int n = polygon.size() / 2;
        edges.reserve(n);
        
        for (int i = 0; i < n; i++) {
            int j = (i + 1) % n;
            
            // Convert to integer coordinates
            auto [ix1, iy1] = transform.worldToInteger(polygon[i * 2], polygon[i * 2 + 1]);
            auto [ix2, iy2] = transform.worldToInteger(polygon[j * 2], polygon[j * 2 + 1]);
            
            // Skip degenerate edges
            if (ix1 == ix2 && iy1 == iy2) continue;
            
            edges.emplace_back(ix1, iy1, ix2, iy2);
        }
        
        buildYBuckets();
    }
    
    // Process an entire row of grid points at Y coordinate (in integer space)
    void classifyGridRow(int64_t intY, int64_t leftIntX, int64_t rightIntX, int64_t stepIntX, std::vector<int>& results) const {
        results.clear();
        
        if (intY < minY || intY > maxY) {
            // Y is outside polygon bounds - all points are outside
            int numPoints = static_cast<int>((rightIntX - leftIntX) / stepIntX) + 1;
            results.resize(numPoints, 0);
            return;
        }
        
        // Update active edges for this Y coordinate
        activeEdges.clear();
        
        // Add edges that start at or before this Y
        for (const auto& bucket : edgeStartTable) {
            if (bucket.first > intY) break; // No more edges can be active
            
            for (int edgeIdx : bucket.second) {
                const auto& edge = edges[edgeIdx];
                
                // Check if edge is active at this Y
                if (edge.minY <= intY && edge.maxY > intY) {
                    // Skip horizontal edges
                    if (edge.y1 != edge.y2) {
                        activeEdges.emplace_back(edge, intY);
                    }
                }
            }
        }
        
        // Collect intersection X coordinates using integer arithmetic
        auto& intersectionBuffer = g_bufferPool.getIntersectionBuffer();
        intersectionBuffer.clear();
        for (const auto& activeEdge : activeEdges) {
            intersectionBuffer.push_back(static_cast<double>(activeEdge.x));
        }
        std::sort(intersectionBuffer.begin(), intersectionBuffer.end());
        
        // Classify each point in the row using integer comparisons
        int numPoints = static_cast<int>((rightIntX - leftIntX) / stepIntX) + 1;
        results.resize(numPoints);
        
        for (int i = 0; i < numPoints; i++) {
            int64_t intX = leftIntX + i * stepIntX;
            
            // Count intersections to the right of point using integer comparison
            int intersectionCount = 0;
            for (double intersectX : intersectionBuffer) {
                if (static_cast<int64_t>(intersectX) > intX) {
                    intersectionCount++;
                }
            }
            
            results[i] = (intersectionCount % 2 == 1) ? 1 : 0;
        }
    }
    
    // Legacy single-point interface (converts to integer internally)
    bool isInside(double x, double y, const IntegerCoordinateTransform& transform) const {
        auto [intX, intY] = transform.worldToInteger(x, y);
        
        if (intY < minY || intY > maxY) return false;
        
        int intersections = 0;
        for (const auto& edge : edges) {
            if (edge.intersectsRay(intX, intY)) {
                intersections++;
            }
        }
        
        return (intersections % 2) == 1;
    }
    
    // Check if point is exactly on the polygon boundary (integer space)
    bool isOnBoundary(double x, double y, const IntegerCoordinateTransform& transform) const {
        auto [intX, intY] = transform.worldToInteger(x, y);
        
        for (const auto& edge : edges) {
            if (edge.isPointOnEdge(intX, intY)) {
                return true;
            }
        }
        return false;
    }
    
    // Get point classification: 0 = outside, 1 = inside, 2 = on boundary
    int classifyPoint(double x, double y, const IntegerCoordinateTransform& transform) const {
        // First check if point is on boundary
        if (isOnBoundary(x, y, transform)) {
            return 2; // On boundary
        }
        
        // If not on boundary, use scanline for inside/outside
        return isInside(x, y, transform) ? 1 : 0;
    }
};

/**
 * Scanline-based point-in-polygon tester with Y-bucketed edge processing
 */
class ScanlinePolygonTester {
private:
    std::vector<PolygonEdge> edges;
    
    // Y-bucketed edge tables for scanline algorithm
    struct ActiveEdge {
        double x;           // Current X intersection
        double deltaX;      // X increment per Y step
        int edgeIndex;      // Original edge index
        double maxY;        // Y coordinate where edge ends
        
        ActiveEdge(const PolygonEdge& edge, double y) : edgeIndex(-1) {
            if (std::abs(edge.y2 - edge.y1) < 1e-10) {
                // Horizontal edge - skip
                deltaX = 0;
                x = edge.x1;
                maxY = edge.maxY;
            } else {
                deltaX = (edge.x2 - edge.x1) / (edge.y2 - edge.y1);
                x = edge.getXAtY(y);
                maxY = edge.maxY;
            }
        }
        
        void stepY(double dy) {
            x += deltaX * dy;
        }
    };
    
    // Y-coordinate to edges that start at that Y
    std::map<double, std::vector<int>> edgeStartTable;
    double minY, maxY;
    
    // Cache for grid row processing - using buffer pool now
    mutable std::vector<ActiveEdge> activeEdges;
    
    void buildYBuckets() {
        if (edges.empty()) return;
        
        minY = edges[0].minY;
        maxY = edges[0].maxY;
        
        for (size_t i = 0; i < edges.size(); i++) {
            const auto& edge = edges[i];
            minY = std::min(minY, edge.minY);
            maxY = std::max(maxY, edge.maxY);
            
            // Add edge to bucket at its minimum Y coordinate
            edgeStartTable[edge.minY].push_back(i);
        }
    }
    
public:
    ScanlinePolygonTester(const std::vector<double>& polygon) {
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
        
        buildYBuckets();
    }
    
    // Process an entire row of grid points at Y coordinate
    void classifyGridRow(double y, double leftX, double rightX, double stepX, std::vector<int>& results) const {
        results.clear();
        
        if (y < minY || y > maxY) {
            // Y is outside polygon bounds - all points are outside
            int numPoints = static_cast<int>((rightX - leftX) / stepX) + 1;
            results.resize(numPoints, 0);
            return;
        }
        
        // Update active edges for this Y coordinate
        activeEdges.clear();
        
        // Add edges that start at or before this Y
        for (const auto& bucket : edgeStartTable) {
            if (bucket.first > y) break; // No more edges can be active
            
            for (int edgeIdx : bucket.second) {
                const auto& edge = edges[edgeIdx];
                
                // Check if edge is active at this Y
                if (edge.minY <= y && edge.maxY > y) {
                    // Skip horizontal edges
                    if (std::abs(edge.y2 - edge.y1) > 1e-10) {
                        activeEdges.emplace_back(edge, y);
                    }
                }
            }
        }
        
        // Use buffer pool for intersections to avoid allocation
        auto& intersectionBuffer = g_bufferPool.getIntersectionBuffer();
        for (const auto& activeEdge : activeEdges) {
            intersectionBuffer.push_back(activeEdge.x);
        }
        std::sort(intersectionBuffer.begin(), intersectionBuffer.end());
        
        // Classify each point in the row
        int numPoints = static_cast<int>((rightX - leftX) / stepX) + 1;
        results.resize(numPoints);
        
        for (int i = 0; i < numPoints; i++) {
            double x = leftX + i * stepX;
            
            // Count intersections to the right of point
            int intersectionCount = 0;
            for (double intersectX : intersectionBuffer) {
                if (intersectX > x) {
                    intersectionCount++;
                }
            }
            
            results[i] = (intersectionCount % 2 == 1) ? 1 : 0;
        }
    }
    
    // Legacy single-point interface for boundary detection
    bool isInside(double x, double y) const {
        if (y < minY || y > maxY) return false;
        
        int intersections = 0;
        for (const auto& edge : edges) {
            if (edge.intersectsRay(x, y)) {
                intersections++;
            }
        }
        
        return (intersections % 2) == 1;
    }
    
    // Check if point is exactly on the polygon boundary
    bool isOnBoundary(double x, double y, double tolerance = 1e-10) const {
        for (const auto& edge : edges) {
            if (edge.isPointOnEdge(x, y, tolerance)) {
                return true;
            }
        }
        return false;
    }
    
    // Get point classification: 0 = outside, 1 = inside, 2 = on boundary
    int classifyPoint(double x, double y, double tolerance = 1e-10) const {
        // First check if point is on boundary
        if (isOnBoundary(x, y, tolerance)) {
            return 2; // On boundary
        }
        
        // If not on boundary, use scanline for inside/outside
        return isInside(x, y) ? 1 : 0;
    }
};

/**
 * Original spatial-indexed polygon tester (kept for spatial intersection queries)
 */
class PolygonTester {
private:
    std::vector<PolygonEdge> edges;
    
    // Spatial index: grid cells -> list of edge indices that might intersect them
    std::vector<std::vector<std::vector<int>>> spatialIndex;
    double gridOriginX, gridOriginY, cellSize;
    int gridWidth, gridHeight;
    
    // Build spatial index for fast edge-cell queries
    void buildSpatialIndex(const MeshCutOptions& options) {
        gridOriginX = options.gridOriginX;
        gridOriginY = options.gridOriginY; 
        cellSize = options.cellSize;
        gridWidth = options.gridWidth;
        gridHeight = options.gridHeight;
        
        // Initialize spatial index grid
        spatialIndex.resize(gridHeight);
        for (int j = 0; j < gridHeight; j++) {
            spatialIndex[j].resize(gridWidth);
        }
        
        // For each edge, add it to all cells it might intersect
        for (size_t edgeIdx = 0; edgeIdx < edges.size(); edgeIdx++) {
            const auto& edge = edges[edgeIdx];
            
            // Get edge bounding box in grid coordinates
            double edgeMinX = std::min(edge.x1, edge.x2);
            double edgeMaxX = std::max(edge.x1, edge.x2);
            double edgeMinY = std::min(edge.y1, edge.y2);
            double edgeMaxY = std::max(edge.y1, edge.y2);
            
            // Convert to grid cell indices (with some padding for safety)
            int minI = std::max(0, static_cast<int>(std::floor((edgeMinX - gridOriginX) / cellSize)) - 1);
            int maxI = std::min(gridWidth - 1, static_cast<int>(std::ceil((edgeMaxX - gridOriginX) / cellSize)) + 1);
            int minJ = std::max(0, static_cast<int>(std::floor((edgeMinY - gridOriginY) / cellSize)) - 1);
            int maxJ = std::min(gridHeight - 1, static_cast<int>(std::ceil((edgeMaxY - gridOriginY) / cellSize)) + 1);
            
            // Add edge to all potentially intersecting cells
            for (int j = minJ; j <= maxJ; j++) {
                for (int i = minI; i <= maxI; i++) {
                    spatialIndex[j][i].push_back(edgeIdx);
                }
            }
        }
    }
    
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
    
    // Initialize spatial index (call this once after construction)
    void initSpatialIndex(const MeshCutOptions& options) {
        buildSpatialIndex(options);
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
    
    // Check if point is exactly on the polygon boundary
    bool isOnBoundary(double x, double y, double tolerance = 1e-10) const {
        for (const auto& edge : edges) {
            if (edge.isPointOnEdge(x, y, tolerance)) {
                return true;
            }
        }
        return false;
    }
    
    // Get point classification: 0 = outside, 1 = inside, 2 = on boundary
    int classifyPoint(double x, double y, double tolerance = 1e-10) const {
        // First check if point is on boundary
        if (isOnBoundary(x, y, tolerance)) {
            return 2; // On boundary
        }
        
        // If not on boundary, use ray casting for inside/outside
        return isInside(x, y) ? 1 : 0;
    }
    
    // Fast spatial-index-based rectangle intersection test
    bool intersectsRectFast(int cellI, int cellJ) const {
        if (spatialIndex.empty() || cellJ < 0 || cellJ >= gridHeight || cellI < 0 || cellI >= gridWidth) {
            return false; // Spatial index not built or cell out of bounds
        }
        
        // Get cell bounds in world coordinates
        double left = gridOriginX + cellI * cellSize;
        double right = gridOriginX + (cellI + 1) * cellSize;
        double bottom = gridOriginY + cellJ * cellSize;
        double top = gridOriginY + (cellJ + 1) * cellSize;
        
        // Only test edges in spatial index for this cell
        const auto& cellEdges = spatialIndex[cellJ][cellI];
        for (int edgeIdx : cellEdges) {
            const auto& edge = edges[edgeIdx];
            
            // Quick bounds check (redundant but fast)
            double edgeMinX = std::min(edge.x1, edge.x2);
            double edgeMaxX = std::max(edge.x1, edge.x2);
            if (edgeMaxX < left || edgeMinX > right || edge.maxY < bottom || edge.minY > top) {
                continue;
            }
            
            // Precise intersection test
            if (lineSegmentIntersectsRect(edge.x1, edge.y1, edge.x2, edge.y2, left, right, bottom, top)) {
                return true;
            }
        }
        return false;
    }
    
    // Legacy method for backward compatibility (slower)
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
 * Batch grid rasterization system using integer arithmetic for maximum performance
 * Eliminates floating-point overhead by working entirely in integer coordinate space
 */
class BatchGridRasterizer {
private:
    // Grid of point classifications: 0=out, 1=in, 2=boundary
    std::vector<std::vector<int>> gridPoints;
    int gridWidth, gridHeight;
    CoordinateTransform transform;
    std::unique_ptr<IntegerCoordinateTransform> intTransform;
    
public:
    BatchGridRasterizer(int width, int height, const CoordinateTransform& trans) 
        : gridWidth(width), gridHeight(height), transform(trans) {
        
        // Create integer coordinate transform for this grid
        auto bounds = transform.getBounds();
        intTransform = std::make_unique<IntegerCoordinateTransform>(
            width, height, bounds.minX, bounds.minY, bounds.maxX, bounds.maxY);
        
        // Allocate grid for (width+1) x (height+1) grid points
        gridPoints.resize(gridHeight + 1);
        for (int j = 0; j <= gridHeight; ++j) {
            gridPoints[j].resize(gridWidth + 1, 0); // Initialize all as "out"
        }
    }
    
    /**
     * Rasterize polygon to grid using integer scanline fill algorithm
     * This replaces individual point-in-polygon tests with batch processing
     */
    void rasterizePolygon(const IntegerScanlinePolygonTester& intScanlineTester, int minI, int maxI, int minJ, int maxJ) {
        // Process each row of grid points using integer scanline
        for (int j = minJ; j <= maxJ + 1; ++j) { // +1 to include top edge of top cells
            if (j < 0 || j > gridHeight) continue;
            
            // Get integer coordinates for this grid row
            auto [leftX, rowY] = transform.gridToWorld(minI, j);
            auto [rightX, _] = transform.gridToWorld(maxI + 1, j);
            
            // Convert to integer coordinate space
            auto [intLeftX, intRowY] = intTransform->worldToInteger(leftX, rowY);
            auto [intRightX, intRowY2] = intTransform->worldToInteger(rightX, rowY);
            
            double cellSize = transform.getCellSize();
            auto [intStepX, intStepY] = intTransform->worldToInteger(cellSize, 0);
            intStepX = intStepX - intTransform->worldToInteger(0, 0).first; // Get relative step
            
            // Use buffer pool for row results
            auto& rowResults = g_bufferPool.getGridRowBuffer1();
            
            // Classify entire row of grid points using integer arithmetic
            intScanlineTester.classifyGridRow(intRowY, intLeftX, intRightX, intStepX, rowResults);
            
            // Store results in our grid
            for (int i = minI; i <= maxI + 1; ++i) { // +1 to include right edge of rightmost cells
                if (i < 0 || i > gridWidth) continue;
                int idx = i - minI;
                if (idx >= 0 && idx < (int)rowResults.size()) {
                    gridPoints[j][i] = rowResults[idx];
                }
            }
        }
    }
    
    /**
     * Fast cell classification using pre-computed grid points with integer arithmetic
     * Replaces the expensive 4-corner point-in-polygon tests
     */
    CellState classifyCellFast(int i, int j, const PolygonTester& spatialTester) {
        // Check bounds
        if (i < 0 || i >= gridWidth || j < 0 || j >= gridHeight) {
            return FULL_OUT;
        }
        
        // Get pre-computed classifications for 4 corners
        int class0 = gridPoints[j][i];         // bottom-left
        int class1 = gridPoints[j][i + 1];     // bottom-right  
        int class2 = gridPoints[j + 1][i + 1]; // top-right
        int class3 = gridPoints[j + 1][i];     // top-left
        
        // Count corners that are strictly inside (not on boundary)
        int insideCount = (class0 == 1 ? 1 : 0) + (class1 == 1 ? 1 : 0) + 
                          (class2 == 1 ? 1 : 0) + (class3 == 1 ? 1 : 0);
        
        // Count corners that are on boundary 
        int boundaryCount = (class0 == 2 ? 1 : 0) + (class1 == 2 ? 1 : 0) + 
                            (class2 == 2 ? 1 : 0) + (class3 == 2 ? 1 : 0);
        
        // All corners strictly inside = fully inside
        if (insideCount == 4) {
            return FULL_IN;
        }
        
        // No corners inside or on boundary - check if polygon edges intersect cell
        if (insideCount == 0 && boundaryCount == 0) {
            // Use fast spatial-indexed intersection test
            if (spatialTester.intersectsRectFast(i, j)) {
                return BOUNDARY;
            } else {
                return FULL_OUT;
            }
        }
        
        // Any corners inside or on boundary, or mixed states = boundary
        return BOUNDARY;
    }
    
    /**
     * Get pre-computed classification for a specific grid point
     */
    int getGridPointClassification(int i, int j) const {
        if (i < 0 || i > gridWidth || j < 0 || j > gridHeight) {
            return 0; // out of bounds = outside
        }
        return gridPoints[j][i];
    }
};

/**
 * Optimized grid cell classification using scanline algorithm
 * Processes entire rows of cells efficiently using Y-bucketed edges
 */
CellState classifyCell(int i, int j, const ScanlinePolygonTester& scanlineTester, const PolygonTester& spatialTester, const CoordinateTransform& transform) {
    // Get the 4 corners of the cell in world coordinates
    auto [x0, y0] = transform.gridToWorld(i, j);
    auto [x1, y1] = transform.gridToWorld(i + 1, j);
    auto [x2, y2] = transform.gridToWorld(i + 1, j + 1);
    auto [x3, y3] = transform.gridToWorld(i, j + 1);
    
    // Test all 4 corners using boundary-aware classification
    int class0 = scanlineTester.classifyPoint(x0, y0);  // 0=out, 1=in, 2=boundary
    int class1 = scanlineTester.classifyPoint(x1, y1);
    int class2 = scanlineTester.classifyPoint(x2, y2);
    int class3 = scanlineTester.classifyPoint(x3, y3);
    
    // Count corners that are strictly inside (not on boundary)
    int insideCount = (class0 == 1 ? 1 : 0) + (class1 == 1 ? 1 : 0) + 
                      (class2 == 1 ? 1 : 0) + (class3 == 1 ? 1 : 0);
    
    // Count corners that are on boundary 
    int boundaryCount = (class0 == 2 ? 1 : 0) + (class1 == 2 ? 1 : 0) + 
                        (class2 == 2 ? 1 : 0) + (class3 == 2 ? 1 : 0);
    
    // All corners strictly inside = fully inside
    if (insideCount == 4) {
        return FULL_IN;
    }
    
    // No corners inside or on boundary - check if polygon edges intersect cell
    if (insideCount == 0 && boundaryCount == 0) {
        // Use fast spatial-indexed intersection test from the original tester
        if (spatialTester.intersectsRectFast(i, j)) {
            return BOUNDARY;
        } else {
            return FULL_OUT;
        }
    }
    
    // Any corners inside or on boundary, or mixed states = boundary
    return BOUNDARY;
}

/**
 * Batch process an entire row of grid cells using scanline algorithm
 * This is the core optimization - O(n+m) instead of O(n*m) for an entire row
 */
void classifyGridRow(int j, int minI, int maxI, const ScanlinePolygonTester& scanlineTester, 
                    const PolygonTester& spatialTester, const CoordinateTransform& transform,
                    std::vector<CellState>& rowResults) {
    rowResults.clear();
    rowResults.resize(maxI - minI + 1);
    
    // Process the bottom edge of all cells in this row (Y = j)
    auto [leftX, bottomY] = transform.gridToWorld(minI, j);
    auto [rightX, _] = transform.gridToWorld(maxI + 1, j);
    double cellSize = (rightX - leftX) / (maxI - minI + 1);
    
    // Use buffer pool to avoid allocations
    auto& bottomResults = g_bufferPool.getGridRowBuffer1();
    scanlineTester.classifyGridRow(bottomY, leftX, rightX, cellSize, bottomResults);
    
    // Process the top edge of all cells in this row (Y = j + 1)
    auto [leftX2, topY] = transform.gridToWorld(minI, j + 1);
    auto& topResults = g_bufferPool.getGridRowBuffer2();
    scanlineTester.classifyGridRow(topY, leftX2, rightX, cellSize, topResults);
    
    // Classify each cell based on its corners
    for (int i = minI; i <= maxI; i++) {
        int cellIdx = i - minI;
        
        // Get corner classifications (we have bottom and top edges, need to extrapolate for all 4 corners)
        // Bottom-left and bottom-right from bottom edge
        int class0 = cellIdx < static_cast<int>(bottomResults.size()) ? bottomResults[cellIdx] : 0;      // bottom-left
        int class1 = (cellIdx + 1) < static_cast<int>(bottomResults.size()) ? bottomResults[cellIdx + 1] : 0;  // bottom-right
        
        // Top-left and top-right from top edge  
        int class2 = (cellIdx + 1) < static_cast<int>(topResults.size()) ? topResults[cellIdx + 1] : 0;    // top-right
        int class3 = cellIdx < static_cast<int>(topResults.size()) ? topResults[cellIdx] : 0;      // top-left
        
        // Count corners that are strictly inside
        int insideCount = class0 + class1 + class2 + class3;
        
        // All corners inside = fully inside
        if (insideCount == 4) {
            rowResults[cellIdx] = FULL_IN;
        }
        // No corners inside - check spatial intersection
        else if (insideCount == 0) {
            if (spatialTester.intersectsRectFast(i, j)) {
                rowResults[cellIdx] = BOUNDARY;
            } else {
                rowResults[cellIdx] = FULL_OUT;
            }
        }
        // Mixed corners = boundary
        else {
            rowResults[cellIdx] = BOUNDARY;
        }
    }
}

/**
 * Legacy single-cell classification (for compatibility)
 */
CellState classifyCell(int i, int j, const PolygonTester& tester, const CoordinateTransform& transform) {
    // Get the 4 corners of the cell in world coordinates
    auto [x0, y0] = transform.gridToWorld(i, j);
    auto [x1, y1] = transform.gridToWorld(i + 1, j);
    auto [x2, y2] = transform.gridToWorld(i + 1, j + 1);
    auto [x3, y3] = transform.gridToWorld(i, j + 1);
    
    // Test all 4 corners using boundary-aware classification
    int class0 = tester.classifyPoint(x0, y0);  // 0=out, 1=in, 2=boundary
    int class1 = tester.classifyPoint(x1, y1);
    int class2 = tester.classifyPoint(x2, y2);
    int class3 = tester.classifyPoint(x3, y3);
    
    // Count corners that are strictly inside (not on boundary)
    int insideCount = (class0 == 1 ? 1 : 0) + (class1 == 1 ? 1 : 0) + 
                      (class2 == 1 ? 1 : 0) + (class3 == 1 ? 1 : 0);
    
    // Count corners that are on boundary 
    int boundaryCount = (class0 == 2 ? 1 : 0) + (class1 == 2 ? 1 : 0) + 
                        (class2 == 2 ? 1 : 0) + (class3 == 2 ? 1 : 0);
    
    // All corners strictly inside = fully inside
    if (insideCount == 4) {
        return FULL_IN;
    }
    
    // No corners inside or on boundary - check if polygon edges intersect cell
    if (insideCount == 0 && boundaryCount == 0) {
        // Use fast spatial-indexed intersection test
        if (tester.intersectsRectFast(i, j)) {
            return BOUNDARY;
        } else {
            return FULL_OUT;
        }
    }
    
    // Any corners inside or on boundary, or mixed states = boundary
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
 * Advanced Sutherland-Hodgman polygon clipping against axis-aligned rectangle
 * Optimized with early bounds checking, unified clipping logic, and ping-pong buffers
 */
std::vector<std::pair<double, double>> clipPolygonToRect(
    const std::vector<double>& polygon, 
    double left, double right, double bottom, double top) {
    
    if (polygon.size() < 6) return {}; // Need at least 3 vertices
    
    // Convert to point vector and compute bounding box for early exit
    auto& inputBuffer = g_bufferPool.getPointBuffer(0);
    inputBuffer.clear();
    double minX = std::numeric_limits<double>::max();
    double maxX = std::numeric_limits<double>::lowest();
    double minY = std::numeric_limits<double>::max(); 
    double maxY = std::numeric_limits<double>::lowest();
    
    for (size_t i = 0; i < polygon.size(); i += 2) {
        double x = polygon[i];
        double y = polygon[i + 1];
        inputBuffer.emplace_back(x, y);
        minX = std::min(minX, x);
        maxX = std::max(maxX, x);
        minY = std::min(minY, y);
        maxY = std::max(maxY, y);
    }
    
    // Early exit: if polygon bounding box doesn't overlap rectangle
    if (maxX < left || minX > right || maxY < bottom || minY > top) {
        return {}; // No intersection
    }
    
    // Early exit: if polygon is entirely inside rectangle
    if (minX >= left && maxX <= right && minY >= bottom && maxY <= top) {
        return inputBuffer; // Return copy of input
    }
    
    // Edge types for unified clipping function
    enum class EdgeType { LEFT, RIGHT, BOTTOM, TOP };
    
    // Unified clipping function with optimized intersection calculations
    auto clipAgainstEdge = [](const std::vector<std::pair<double, double>>& input,
                             std::vector<std::pair<double, double>>& output,
                             EdgeType edgeType, double edgeValue) {
        output.clear();
        if (input.empty()) return;
        
        auto prev = input.back();
        for (const auto& curr : input) {
            bool currInside, prevInside;
            double intersectionX, intersectionY;
            
            // Optimized inside/outside tests and intersection calculation per edge type
            switch (edgeType) {
                case EdgeType::LEFT:
                    currInside = curr.first >= edgeValue;
                    prevInside = prev.first >= edgeValue;
                    if (currInside != prevInside) {
                        // Optimized intersection: avoid division when possible
                        double dx = curr.first - prev.first;
                        if (std::abs(dx) > 1e-12) {
                            double t = (edgeValue - prev.first) / dx;
                            intersectionX = edgeValue;
                            intersectionY = prev.second + t * (curr.second - prev.second);
                        } else {
                            intersectionX = edgeValue;
                            intersectionY = prev.second;
                        }
                    }
                    break;
                    
                case EdgeType::RIGHT:
                    currInside = curr.first <= edgeValue;
                    prevInside = prev.first <= edgeValue;
                    if (currInside != prevInside) {
                        double dx = curr.first - prev.first;
                        if (std::abs(dx) > 1e-12) {
                            double t = (edgeValue - prev.first) / dx;
                            intersectionX = edgeValue;
                            intersectionY = prev.second + t * (curr.second - prev.second);
                        } else {
                            intersectionX = edgeValue;
                            intersectionY = prev.second;
                        }
                    }
                    break;
                    
                case EdgeType::BOTTOM:
                    currInside = curr.second >= edgeValue;
                    prevInside = prev.second >= edgeValue;
                    if (currInside != prevInside) {
                        double dy = curr.second - prev.second;
                        if (std::abs(dy) > 1e-12) {
                            double t = (edgeValue - prev.second) / dy;
                            intersectionX = prev.first + t * (curr.first - prev.first);
                            intersectionY = edgeValue;
                        } else {
                            intersectionX = prev.first;
                            intersectionY = edgeValue;
                        }
                    }
                    break;
                    
                case EdgeType::TOP:
                    currInside = curr.second <= edgeValue;
                    prevInside = prev.second <= edgeValue;
                    if (currInside != prevInside) {
                        double dy = curr.second - prev.second;
                        if (std::abs(dy) > 1e-12) {
                            double t = (edgeValue - prev.second) / dy;
                            intersectionX = prev.first + t * (curr.first - prev.first);
                            intersectionY = edgeValue;
                        } else {
                            intersectionX = prev.first;
                            intersectionY = edgeValue;
                        }
                    }
                    break;
            }
            
            // Add intersection point if crossing from outside to inside
            if (currInside && !prevInside) {
                output.emplace_back(intersectionX, intersectionY);
            }
            
            // Add current vertex if inside
            if (currInside) {
                output.push_back(curr);
            }
            
            // Add intersection point if crossing from inside to outside
            if (!currInside && prevInside) {
                output.emplace_back(intersectionX, intersectionY);
            }
            
            prev = curr;
        }
    };
    
    // Ping-pong between two buffers instead of using 4 separate ones
    auto& buffer1 = g_bufferPool.getPointBuffer(1);
    auto& buffer2 = g_bufferPool.getPointBuffer(2);
    
    // Apply 4 clipping stages with ping-pong buffering
    clipAgainstEdge(inputBuffer, buffer1, EdgeType::LEFT, left);
    if (buffer1.empty()) return {};
    
    clipAgainstEdge(buffer1, buffer2, EdgeType::RIGHT, right);
    if (buffer2.empty()) return {};
    
    clipAgainstEdge(buffer2, buffer1, EdgeType::BOTTOM, bottom);
    if (buffer1.empty()) return {};
    
    clipAgainstEdge(buffer1, buffer2, EdgeType::TOP, top);
    
    return buffer2;  // Return copy of final result
}

/**
 * Process a boundary cell - clip polygon and triangulate
 * Uses buffer pool to avoid allocations
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
    
    // Clip polygon to cell rectangle using buffer pool
    auto clippedPoints = clipPolygonToRect(polygon, left, right, bottom, top);
    
    if (clippedPoints.size() < 3) {
        return; // Clipped polygon is degenerate
    }
    
    // Use buffer pool for earcut triangulation
    auto& earcutPolygon = g_bufferPool.getEarcutPolygonBuffer();
    earcutPolygon.push_back(clippedPoints);
    
    // Triangulate with earcut using buffer pool
    auto& triangleIndices = g_bufferPool.getTriangleBuffer();
    triangleIndices = mapbox::earcut<uint32_t>(earcutPolygon);
    
    // Use buffer pool for global vertex mapping
    auto& globalIndices = g_bufferPool.getGlobalIndicesBuffer();
    globalIndices.reserve(clippedPoints.size());
    
    for (const auto& point : clippedPoints) {
        WorldPoint worldPoint{point.first, point.second};
        
        // Check if we already have this vertex
        auto it = vertexMap.find(worldPoint);
        if (it != vertexMap.end()) {
            globalIndices.push_back(static_cast<uint32_t>(it->second));
        } else {
            // Create new vertex index
            N newIndex = nextVertexIndex++;
            vertexMap[worldPoint] = newIndex;
            globalIndices.push_back(static_cast<uint32_t>(newIndex));
        }
    }
    
    // Add triangles to output using global indices
    for (uint32_t triangleIndex : triangleIndices) {
        indices.push_back(static_cast<N>(globalIndices[triangleIndex]));
    }
}

} // namespace detail

/**
 * Main meshcut implementation (full result with vertices)
 */
template <typename N>
MeshCutResult meshcut_full(const std::vector<double>& polygon, const std::vector<N>& holes, const MeshCutOptions& options) {
    
    if (polygon.size() < 6) {
        return {}; // Need at least 3 vertices (6 coordinates)
    }
    
    std::vector<N> indices;
    detail::CoordinateTransform transform(options);
    
    // Extract outer boundary and holes from the input format
    std::vector<double> outerRing;
    std::vector<std::vector<double>> holeRings;
    
    if (holes.empty()) {
        // No holes - use entire polygon as outer ring
        outerRing = polygon;
    } else {
        // Extract outer ring (from start to first hole)
        size_t outerEnd = holes[0] * 2; // holes contains vertex indices, not coordinate indices
        outerRing.assign(polygon.begin(), polygon.begin() + outerEnd);
        
        // Extract each hole
        for (size_t h = 0; h < holes.size(); h++) {
            size_t holeStart = holes[h] * 2;
            size_t holeEnd = (h + 1 < holes.size()) ? holes[h + 1] * 2 : polygon.size();
            
            std::vector<double> hole(polygon.begin() + holeStart, polygon.begin() + holeEnd);
            if (!hole.empty()) {
                holeRings.push_back(std::move(hole));
            }
        }
    }
    
    // Create both testers: integer scanline for point classification, spatial for edge-cell intersection
    // Use integer coordinate transform for maximum performance
    auto bounds = transform.getBounds();
    detail::IntegerCoordinateTransform intTransform(options.gridWidth, options.gridHeight, 
                                                   bounds.minX, bounds.minY, bounds.maxX, bounds.maxY);
    detail::IntegerScanlinePolygonTester intScanlineTester(outerRing, intTransform);
    detail::PolygonTester spatialTester(outerRing);
    
    // Create hole testers for point-in-hole detection
    std::vector<std::unique_ptr<detail::IntegerScanlinePolygonTester>> holeTesters;
    std::vector<std::unique_ptr<detail::PolygonTester>> holeSpatialTesters;
    
    for (const auto& hole : holeRings) {
        holeTesters.emplace_back(std::make_unique<detail::IntegerScanlinePolygonTester>(hole, intTransform));
        holeSpatialTesters.emplace_back(std::make_unique<detail::PolygonTester>(hole));
        holeSpatialTesters.back()->initSpatialIndex(options);
    }
    
    // Initialize spatial index for fast edge-cell intersection queries
    spatialTester.initSpatialIndex(options);
    
    // Compute polygon bounding box in grid coordinates using outer ring
    double minX = outerRing[0], maxX = outerRing[0];
    double minY = outerRing[1], maxY = outerRing[1];
    for (size_t i = 2; i < outerRing.size(); i += 2) {
        minX = std::min(minX, outerRing[i]);
        maxX = std::max(maxX, outerRing[i]);
        minY = std::min(minY, outerRing[i + 1]);
        maxY = std::max(maxY, outerRing[i + 1]);
    }
    
    auto minGrid = transform.worldToGrid(minX, minY);
    auto maxGrid = transform.worldToGrid(maxX, maxY);
    
    // Clamp to actual grid bounds
    int minI = std::max(0, minGrid.x);
    int maxI = std::min(options.gridWidth - 1, maxGrid.x);
    int minJ = std::max(0, minGrid.y);
    int maxJ = std::min(options.gridHeight - 1, maxGrid.y);
    
    // Create batch grid rasterizer and pre-compute all grid point classifications using integer arithmetic
    detail::BatchGridRasterizer gridRasterizer(options.gridWidth, options.gridHeight, transform);
    gridRasterizer.rasterizePolygon(intScanlineTester, minI, maxI, minJ, maxJ);
    
    // Vertex deduplication map and counter
    std::unordered_map<detail::WorldPoint, N, detail::WorldPointHash> vertexMap;
    N nextVertexIndex = 0;
    
    // Process grid cells using hole-aware classification
    for (int j = minJ; j <= maxJ; j++) {
        for (int i = minI; i <= maxI; i++) {
            // Use fast cell classification based on pre-computed grid points for outer polygon
            detail::CellState outerState = gridRasterizer.classifyCellFast(i, j, spatialTester);
            
            // If cell is outside the outer polygon, skip it entirely
            if (outerState == detail::FULL_OUT) {
                continue;
            }
            
            // Check if cell is inside any hole
            bool insideHole = false;
            for (size_t h = 0; h < holeSpatialTesters.size(); h++) {
                // Check if this cell intersects or is inside the hole
                // For simplicity, check the cell center
                auto [centerX, centerY] = transform.gridToWorld(i + 0.5, j + 0.5);
                if (holeTesters[h]->isInside(centerX, centerY, intTransform)) {
                    insideHole = true;
                    break;
                }
            }
            
            // If inside a hole, skip this cell
            if (insideHole) {
                continue;
            }
            
            // For cells that might intersect hole boundaries, we need more careful processing
            detail::CellState finalState = outerState;
            if (outerState == detail::FULL_IN) {
                // Check if any hole intersects this cell
                for (const auto& holeSpatialTester : holeSpatialTesters) {
                    if (holeSpatialTester->intersectsRectFast(i, j)) {
                        finalState = detail::BOUNDARY; // Need clipping
                        break;
                    }
                }
            }
            
            switch (finalState) {
                case detail::FULL_IN:
                    detail::emitGridTriangles(i, j, options.diagonalNE, indices, vertexMap, nextVertexIndex, transform);
                    break;
                    
                case detail::FULL_OUT:
                    // Skip - no triangles to emit
                    break;
                    
                case detail::BOUNDARY:
                    // For boundary cells, we need to clip against both outer polygon and holes
                    detail::processBoundaryCell(i, j, outerRing, transform, indices, vertexMap, nextVertexIndex);
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