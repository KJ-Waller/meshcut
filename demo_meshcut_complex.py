"""
MeshCut Complex Polygon Visualization
=====================================
Demonstrates triangle orientation control with complex realistic shapes.
Uses ONLY the MeshCut library from ./meshcut directory.

Triangle Orientations:
- \\ diagonals = TOPLEFT_BOTTOMRIGHT 
- / diagonals = TOPRIGHT_BOTTOMLEFT
"""

import sys
import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import LineCollection

# Add the MeshCut build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build'))

# Load the MeshCut library directly
meshcut_lib = ctypes.CDLL('./build/libmeshcut.so')

class Point:
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

class BoundingBox:
    def __init__(self, min_x, min_y, max_x, max_y):
        self.min_x = float(min_x)
        self.min_y = float(min_y) 
        self.max_x = float(max_x)
        self.max_y = float(max_y)

class Triangle:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

# Define C++ function signatures
meshcut_lib.tessellate_simple.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # polygon points (x,y pairs)
    ctypes.c_int,                     # number of points
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,  # bbox
    ctypes.c_int,                     # grid_size
    ctypes.c_int,                     # orientation (0=TL-BR, 1=TR-BL)
    ctypes.POINTER(ctypes.c_double),  # output triangles
    ctypes.POINTER(ctypes.c_int)      # output count
]
meshcut_lib.tessellate_simple.restype = ctypes.c_int

def meshcut_tessellate(polygon, bbox, grid_size, orientation):
    """
    Tessellate polygon using MeshCut library
    orientation: 0 = TOPLEFT_BOTTOMRIGHT (\\), 1 = TOPRIGHT_BOTTOMLEFT (/)
    """
    # Convert polygon to C array
    n_points = len(polygon)
    points_array = (ctypes.c_double * (n_points * 2))()
    for i, (x, y) in enumerate(polygon):
        points_array[i * 2] = x
        points_array[i * 2 + 1] = y
    
    # Prepare output arrays
    max_triangles = 1000
    triangles_array = (ctypes.c_double * (max_triangles * 6))()
    count = ctypes.c_int()
    
    # Call MeshCut
    result = meshcut_lib.tessellate_simple(
        points_array, n_points,
        bbox.min_x, bbox.min_y, bbox.max_x, bbox.max_y,
        grid_size, orientation,
        triangles_array, ctypes.byref(count)
    )
    
    if result != 0:
        raise RuntimeError(f"MeshCut tessellation failed with code {result}")
    
    # Convert result to triangles
    triangles = []
    for i in range(count.value):
        base = i * 6
        a = Point(triangles_array[base], triangles_array[base + 1])
        b = Point(triangles_array[base + 2], triangles_array[base + 3])
        c = Point(triangles_array[base + 4], triangles_array[base + 5])
        triangles.append(Triangle(a, b, c))
    
    return triangles

def visualize_complex_polygons():
    """Create comprehensive visualizations of complex polygon tessellation."""
    
    # Complex building footprint (L-shaped office complex)
    building = [
        (0, 0), (8, 0), (8, 3), (5, 3), (5, 5), (8, 5), 
        (8, 8), (3, 8), (3, 6), (0, 6), (0, 0)
    ]
    
    # Organic terrain boundary (irregular natural boundary)
    terrain = [
        (10, 1), (12, 0.5), (14, 1.5), (15, 3), (14.5, 4.5),
        (13, 5.5), (11.5, 6), (10.5, 5), (9.5, 3.5), (9, 2), (10, 1)
    ]
    
    # Complex star polygon (architectural feature)
    star_center = (20, 4)
    star_points = []
    for i in range(10):
        angle = i * np.pi / 5
        radius = 3 if i % 2 == 0 else 1.5
        x = star_center[0] + radius * np.cos(angle)
        y = star_center[1] + radius * np.sin(angle)
        star_points.append((x, y))
    
    polygons = [
        ("Complex Building Footprint", building, (0, 0, 8, 8)),
        ("Organic Terrain Boundary", terrain, (9, 0, 16, 7)), 
        ("Architectural Star Feature", star_points, (17, 1, 23, 7))
    ]
    
    fig, axes = plt.subplots(3, 2, figsize=(16, 20))
    fig.suptitle('MeshCut Complex Polygon Tessellation\nTriangle Orientation Comparison', 
                 fontsize=16, fontweight='bold')
    
    orientations = [
        ("Backslash Diagonals (\\)", 0),
        ("Forward Slash Diagonals (/)", 1)
    ]
    
    for poly_idx, (poly_name, polygon, bbox_coords) in enumerate(polygons):
        bbox = BoundingBox(*bbox_coords)
        
        for orient_idx, (orient_name, orient_code) in enumerate(orientations):
            ax = axes[poly_idx, orient_idx]
            
            # Tessellate with MeshCut
            try:
                triangles = meshcut_tessellate(polygon, bbox, 8, orient_code)
                triangle_count = len(triangles)
                
                # Draw triangles
                for triangle in triangles:
                    tri_points = [
                        [triangle.a.x, triangle.a.y],
                        [triangle.b.x, triangle.b.y], 
                        [triangle.c.x, triangle.c.y],
                        [triangle.a.x, triangle.a.y]
                    ]
                    tri_patch = patches.Polygon(
                        tri_points[:-1], 
                        fill=False, 
                        edgecolor='lightblue', 
                        linewidth=0.5
                    )
                    ax.add_patch(tri_patch)
                
                # Draw polygon boundary
                poly_points = list(polygon) + [polygon[0]]
                poly_x, poly_y = zip(*poly_points)
                ax.plot(poly_x, poly_y, 'r-', linewidth=2, label='Polygon boundary')
                
                # Draw grid
                grid_x = np.linspace(bbox.min_x, bbox.max_x, 9)
                grid_y = np.linspace(bbox.min_y, bbox.max_y, 9)
                
                for x in grid_x:
                    ax.axvline(x, color='gray', alpha=0.3, linewidth=0.5)
                for y in grid_y:
                    ax.axhline(y, color='gray', alpha=0.3, linewidth=0.5)
                
                # Configure plot
                ax.set_xlim(bbox.min_x - 0.5, bbox.max_x + 0.5)
                ax.set_ylim(bbox.min_y - 0.5, bbox.max_y + 0.5)
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.2)
                ax.set_title(f'{poly_name}\\n{orient_name}\\n{triangle_count} triangles')
                
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {str(e)}', 
                       transform=ax.transAxes, ha='center', va='center')
                ax.set_title(f'{poly_name}\\n{orient_name}\\nFailed')
    
    plt.tight_layout()
    plt.savefig('/home/kevin/gridcut/meshcut/meshcut_complex_demo.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    print("🎉 MeshCut complex polygon visualization completed!")
    print("📁 Saved as: meshcut_complex_demo.png")

if __name__ == "__main__":
    print("🔄 MeshCut Complex Polygon Tessellation Demo")
    print("============================================")
    print("Using ONLY MeshCut library from ./meshcut directory")
    print("")
    
    visualize_complex_polygons()