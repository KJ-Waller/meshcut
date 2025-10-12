# MapLibre MeshCut Integration Analysis Report

## Summary

I've successfully analyzed your MapLibre Native integration with MeshCut and identified several potential causes for the edge artifacts you're experiencing. The analysis replicates your exact integration approach and tests edge-crossing scenarios.

## Key Findings

### 1. **Grid Configuration Analysis**

Your MapLibre integration uses:
```cpp
const double extentDouble = static_cast<double>(EXTENT);  // 8192.0
options.gridWidth = 32;
options.gridHeight = 32; 
options.cellSize = (2.0 * extentDouble) / gridSize;      // 512.0 units per cell
options.gridOriginX = -extentDouble;                      // -8192
options.gridOriginY = -extentDouble;                      // -8192
```

**Analysis**: This creates a 32x32 grid with 512-unit cells covering -8192 to +8192, which is much larger than the MVT tile coordinate space (0 to 4096). This means:
- Most polygons only use a small portion of the grid
- Grid intersections may be creating vertices outside the expected tile boundaries

### 2. **Edge-Crossing Behavior**

Testing revealed significant differences between MeshCut and Earcut for edge-crossing polygons:

| Test Case | Earcut Triangles | MeshCut Triangles | Density Increase |
|-----------|------------------|-------------------|------------------|
| Left Edge Crossing | 3 | 17 | 5.67x |
| Right Edge Crossing | 3 | 17 | 5.67x |
| Large Spanning | 2 | 200 | 100.0x |
| Edge Touching | 4 | 25 | 6.25x |

### 3. **Coordinate System Issues**

**Good News**: No int16_t overflow detected in test cases (all coordinates within -32768 to 32767 range).

**Potential Issue**: Grid extends far beyond tile boundaries, creating unnecessary vertices outside the visible area.

### 4. **Integration Pattern Analysis**

Your integration correctly:
✅ Uses `meshcut_full()` to get both vertices and indices  
✅ Handles vertex array expansion properly  
✅ Falls back to earcut on errors  
✅ Manages coordinate conversion to int16_t  

**Potential Issues Identified**:

#### Issue #1: Grid Scale Mismatch
- Your grid covers -8192 to +8192 (16,384 unit range)
- MVT tiles use 0 to 4096 (4,096 unit range)  
- This means 75% of your grid is outside the tile boundaries

#### Issue #2: Edge Vertex Creation
- MeshCut creates grid intersection vertices that may be outside tile boundaries
- When tiles are rendered adjacent to each other, these "phantom" vertices could cause visual seams
- Polygon edges that cross tile boundaries get additional vertices that neighboring tiles don't have

#### Issue #3: Triangle Density Explosion
- Large polygons spanning multiple grid cells create 50-100x more triangles
- This could impact rendering performance, especially when zoomed out (where you noticed issues)

## Recommended Fixes

### Fix #1: Adjust Grid to Tile Boundaries
```cpp
// Instead of covering -8192 to +8192, cover actual tile space
options.gridOriginX = 0.0;              // Tile starts at 0
options.gridOriginY = 0.0;              // Tile starts at 0
options.cellSize = 4096.0 / gridSize;   // Cover 0-4096 tile space
options.gridWidth = gridSize;
options.gridHeight = gridSize;
```

### Fix #2: Polygon Clipping
Consider clipping polygons to tile boundaries before triangulation to prevent edge artifacts:
```cpp
// Clip polygon to tile boundaries (0,0) to (4096,4096) before MeshCut
```

### Fix #3: Adaptive Grid Size
Use smaller grid sizes for large polygons to reduce triangle explosion:
```cpp
// Reduce grid size for large polygons that span many tiles
uint32_t adaptiveGridSize = std::min(gridSize, 16u); // Cap at 16x16 for large polygons
```

### Fix #4: Edge Detection
Only use MeshCut for polygons that don't cross tile boundaries:
```cpp
bool crossesEdge = (minX < 10 || maxX > 4086 || minY < 10 || maxY > 4086);
if (crossesEdge) {
    // Use earcut for edge-crossing polygons to avoid artifacts
    return TriangulationResult(mapbox::earcut<uint32_t>(polygon));
}
```

## Generated Test Files

The analysis generated these visualization files:
- `edge_artifact_case_1.svg` - Left edge crossing polygon
- `edge_artifact_case_2.svg` - Right edge crossing polygon  
- `edge_artifact_case_3.svg` - Large spanning polygon (shows extreme triangle density)
- `edge_artifact_case_4.svg` - Edge touching polygon

Each shows side-by-side comparison of Earcut vs MeshCut triangulation with tile boundaries marked in red.

## Next Steps

1. **Try Fix #1 first**: Adjust grid to tile coordinate space (0-4096)
2. **Implement Fix #4**: Disable MeshCut for edge-crossing polygons
3. **Test on Android**: See if artifacts disappear, especially when zoomed out
4. **Performance test**: Monitor triangle count impact on rendering performance

The visualizations should help you see exactly how the triangulation differs near tile edges, which is likely the root cause of your zoom/rotation artifacts.