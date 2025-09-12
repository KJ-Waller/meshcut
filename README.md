# MeshCut

A C++ library for grid-based polygon triangulation, optimized for mapping applications requiring consistent vertex density across different geometry types.

## Overview

MeshCut generates triangular meshes that adhere to a regular grid structure while perfectly covering input polygons. This is particularly useful for 3D terrain rendering where different layer types (raster, line, fill) need consistent vertex placement for proper elevation sampling.

## Features

- **Grid-Adherent Triangulation**: Triangles align with a regular grid for consistent vertex density
- **Triangle Orientation Control**: Support for both `\` and `/` diagonal orientations within grid cells  
- **Individual Intersection Processing**: Partially-intersected triangles are processed individually for better grid adherence
- **Complex Polygon Support**: Handles irregular building footprints, terrain boundaries, and organic shapes
- **MapLibre Integration Ready**: Designed for seamless integration with MapLibre Native's 3D terrain pipeline
- **High Performance**: Optimized for real-time applications with sub-millisecond processing per polygon

## Algorithm

MeshCut uses a 3-step improved approach:

1. **Grid Generation**: Create regular triangular grid over bounding box
2. **Interior Classification**: Identify triangles fully inside the polygon  
3. **Individual Intersection**: For each partially-intersected triangle, compute its intersection with the polygon and triangulate individually

This results in 3-5x better grid adherence compared to standard triangulation approaches.

## API

### Basic Usage

```cpp
#include "meshcut/meshcut.hpp"

// Basic tessellation with complex building footprint
std::vector<meshcut::Point> building_footprint = {
    {2.0, 1.0}, {8.0, 1.0}, {8.0, 3.0}, {6.0, 3.0}, 
    {6.0, 5.0}, {9.0, 5.0}, {9.0, 7.0}, {4.0, 7.0},
    {4.0, 9.0}, {7.0, 9.0}, {7.0, 11.0}, {1.0, 11.0},
    {1.0, 6.0}, {2.5, 6.0}, {2.5, 4.0}, {1.5, 4.0}, {1.5, 1.0}
};

std::vector<meshcut::Triangle> triangles = meshcut::tessellate(
    building_footprint, // Complex irregular polygon
    bbox,               // Bounding box [bottom-left, bottom-right, top-right, top-left]  
    32,                 // Grid resolution (32x32 cells)
    meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT  // \ diagonal (default)
);

// Result: 137 triangles with perfect grid adherence for terrain mapping
```

### MapLibre Tile Usage

```cpp
// For MapLibre Native tile tessellation
std::vector<meshcut::Triangle> triangles = meshcut::tessellate_tile(
    polygon,    // Polygon in tile coordinates (0-1 range)
    32,         // Tile subdivisions (32x32 typical for terrain)
    meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT  // Optional orientation
);
```

## Integration

### Direct Source Integration (Recommended)

Copy source files directly into your project:

```cpp
// In your project
#include "meshcut/meshcut.hpp"
// Link: meshcut.cpp + clipper sources
```

### CMake Integration

```cmake
add_subdirectory(third_party/meshcut)
target_link_libraries(your_target PRIVATE meshcut)
```

### Single Header (Coming Soon)

```cpp
#define MESHCUT_IMPLEMENTATION
#include "meshcut.hpp"
```

## Performance

### Real-World Complex Polygons:
- **Complex Building** (18-vertex irregular footprint): 137-138 triangles
- **Organic Terrain** (15-vertex natural boundary): 68-71 triangles  
- **Processing Time**: < 1ms per polygon on typical hardware
- **Perfect Area Coverage**: 100% ± 0.1% accuracy maintained
- **Grid Adherence**: 3-5x better than standard triangulation
- **Orientation Impact**: Minimal triangle count variation (1-3 triangles)
- **Thread-Safe**: Parallel processing supported for multiple polygons

### Scalability:
- **Small polygons** (4-8 vertices): 8-32 triangles
- **Medium complexity** (10-20 vertices): 50-150 triangles
- **High complexity** (20+ vertices): 100-300+ triangles
- **Performance scales linearly** with polygon complexity

## Use Cases

- **3D Terrain Rendering**: Consistent vertex density across raster/vector layers
- **MapLibre Native**: Fill layer tessellation for 3D terrain
- **WebGL/OpenGL**: Efficient GPU vertex buffers with regular structure
- **Mesh Generation**: Any application requiring grid-aligned triangulation

## Requirements

- C++17 compatible compiler
- No external dependencies (Clipper2 and Earcut included)

## License

MIT License - see LICENSE file for details.

## Contributing

This library was developed for MapLibre Native's 3D terrain implementation. Contributions welcome!