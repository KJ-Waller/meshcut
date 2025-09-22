# MeshCut Implementation Blueprint

**Goal**: Create a blazingly fast C++ alternative to earcut that triangulates polygons while adhering to a regular grid mesh, for use in MapLibre Native.

## API Specification

### Input Format (Earcut-compatible + mesh parameters)
```cpp
struct MeshCutOptions {
    // Grid specification
    double gridOriginX, gridOriginY;  // Bottom-left corner of grid in world coordinates
    double cellSize;                  // Size of each grid cell in world units
    int gridWidth, gridHeight;        // Number of cells in each dimension
    bool diagonalNE = true;          // true for / diagonal, false for \ diagonal
};

// Main API - drop-in replacement for earcut
std::vector<uint32_t> meshcut(
    const std::vector<double>& coords,     // [x,y, x,y, ...] polygon vertices
    const std::vector<uint32_t>& holes,   // hole start indices (earcut format)
    const MeshCutOptions& options
);
```

### Output Format (Earcut-compatible)
- `std::vector<uint32_t>` of triangle indices referencing input vertices
- Vertices are NOT duplicated - indices reference original `coords` array plus any new intersection vertices
- New intersection vertices appended to a separate output array if needed

## Core Algorithm Overview

### Phase 1: Coordinate System Transformation
**Purpose**: Convert world coordinates to grid space for integer-based fast operations
- Transform polygon vertices: `gridX = (worldX - gridOriginX) / cellSize`
- All subsequent operations in grid space where cells are unit squares [i,i+1] × [j,j+1]
- Store transformation parameters for final coordinate conversion

### Phase 2: Grid Cell Classification
**Purpose**: Classify each grid cell to minimize expensive operations
- Compute polygon AABB in grid coordinates
- For each cell (i,j) in AABB:
  - Test 4 corner points for inside/outside polygon
  - If all 4 corners inside → **FULL_IN**
  - If no corners inside AND no polygon edges intersect cell → **FULL_OUT** 
  - Otherwise → **BOUNDARY**

### Phase 3: Triangle Generation

#### FULL_IN Cells
- Emit 2 triangles per cell directly from grid
- Respect diagonal direction (NE or NW)
- No clipping or additional processing needed

#### FULL_OUT Cells
- Skip entirely - no triangles generated

#### BOUNDARY Cells
1. Clip polygon to cell rectangle using Sutherland-Hodgman
2. If clipped polygon non-empty, triangulate with earcut
3. Transform triangle vertices back to world coordinates

## Detailed Component Requirements

### 1. Point-in-Polygon Testing
**Performance Critical**: This runs for every grid corner point

**Implementation**:
- Use scanline algorithm with edge buckets by integer Y coordinate
- Pre-process polygon edges into Y-bucketed active edge table
- Single sweep across all grid points gives inside/outside for all corners
- **Optimization**: Process entire rows of grid points at once

**Edge Cases**:
- Handle points exactly on polygon edges (consistent rule needed)
- Degenerate polygon cases (zero area, self-intersecting)

### 2. Polygon Edge vs Cell Intersection
**Purpose**: Detect BOUNDARY cells when no corners are inside

**Implementation**:
- Segment-rectangle intersection test for each polygon edge vs each cell
- **Optimization**: Use spatial indexing - build coarse grid of edges
- **Optimization**: Early rejection using edge bounding boxes
- Only test edges whose bbox overlaps cell bbox

### 3. Sutherland-Hodgman Polygon Clipping
**Purpose**: Clip polygon to axis-aligned cell rectangle

**Requirements**:
- Clip against 4 half-planes: x≥i, x≤i+1, y≥j, y≤j+1
- Handle multiple polygon rings (outer + holes) separately
- Preserve winding order for proper earcut input
- **Optimization**: Reuse clipping buffers, avoid allocations

**Edge Cases**:
- Clipped polygon becomes degenerate (area < epsilon)
- Multiple disconnected pieces after clipping
- Holes that disappear or merge with outer ring

### 4. Grid Triangle Generation
**Purpose**: Generate triangles for FULL_IN cells

**Requirements**:
- Two triangulation patterns based on diagonal direction:
  - NE diagonal (`diagonalNE=true`): (i,j)→(i+1,j)→(i,j+1) and (i+1,j)→(i+1,j+1)→(i,j+1)
  - NW diagonal (`diagonalNE=false`): (i,j)→(i+1,j)→(i+1,j+1) and (i,j)→(i+1,j+1)→(i,j+1)
- Generate vertex indices that reference shared vertex pool
- Transform grid coordinates back to world coordinates

### 5. Vertex Management
**Purpose**: Maintain coherent mesh with shared vertices

**Requirements**:
- Build vertex deduplication map: `(gridX, gridY) → vertexIndex`
- Hash grid coordinates to detect identical vertices across cell boundaries
- **Precision**: Use consistent coordinate quantization (e.g., snap to 1/1024 grid units)
- Append new intersection vertices to vertex array as needed

## Performance Optimizations

### Memory Management
- Pre-allocate vectors with estimated capacity
- Reuse temporary buffers for clipping operations
- Use object pools for frequently allocated temporary objects

### Spatial Indexing
- Build coarse spatial index of polygon edges for intersection testing
- Use simple grid-based indexing: `edgeGrid[cellY][cellX] = list<edgeIndices>`

### Parallelization Opportunities
- BOUNDARY cell processing is embarrassingly parallel
- Point-in-polygon testing can be parallelized by grid rows
- Consider thread pool for boundary cell clipping/triangulation

### Numerical Robustness
- Work in grid coordinate system to minimize floating-point errors
- Use consistent epsilon for geometric tests
- Handle edge cases with consistent tie-breaking rules

## Data Structures

### Core Types
```cpp
struct GridPoint { int x, y; };
struct Cell { int i, j; enum State {FULL_IN, FULL_OUT, BOUNDARY} state; };
struct Edge { GridPoint start, end; double minY, maxY; };
```

### Temporary Buffers (Reused)
```cpp
std::vector<GridPoint> clippingBuffer;
std::vector<uint32_t> triangulationBuffer;
std::unordered_map<GridPoint, uint32_t> vertexMap;
```

## Testing & Validation Requirements

### Correctness Tests
- Compare output triangle count with earcut for simple polygons
- Verify no overlapping triangles in output
- Ensure all triangles have positive area
- Test edge cases: polygon at grid boundaries, very small polygons

### Performance Benchmarks
- Compare speed vs earcut on realistic MapLibre polygons
- Test scalability with different grid resolutions (16x16 to 128x128)
- Memory usage profiling
- Target: ≤2x earcut runtime for equivalent triangle quality

### Integration Tests
- Drop-in replacement testing in MapLibre Native
- Visual validation of rendered triangulations
- Regression testing against known good outputs

## Implementation Phases

1. **Core Pipeline**: Basic algorithm with simple optimizations
2. **Performance Pass**: Add spatial indexing, memory optimizations
3. **Robustness Pass**: Handle edge cases, numerical stability
4. **Integration Pass**: MapLibre Native compatibility, final API polish

## Success Criteria
- **Speed**: Comparable to or faster than earcut for typical MapLibre use cases
- **Quality**: Produces triangulations suitable for GPU rendering
- **Compatibility**: Drop-in replacement for earcut in MapLibre Native
- **Robustness**: Handles real-world polygon data without crashes or artifacts