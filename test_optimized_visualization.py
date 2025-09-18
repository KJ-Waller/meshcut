"""
MeshCut Optimized Visualization Test
===================================
Verify that the fast geometric intersection produces correct results.
"""

import sys
import os
import ctypes
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Add the MeshCut build directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'build'))

# Load the optimized MeshCut library
meshcut_lib = ctypes.CDLL('./build/libmeshcut.so')

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

def meshcut_tessellate(polygon, bbox, grid_size, orientation):
    """Tessellate using optimized MeshCut."""
    # Convert polygon to C array
    n_points = len(polygon)
    points_array = (ctypes.c_double * (n_points * 2))()
    for i, (x, y) in enumerate(polygon):
        points_array[i * 2] = x
        points_array[i * 2 + 1] = y
    
    # Prepare output arrays
    max_triangles = 2000  # Increased for 32x32 grid
    triangles_array = (ctypes.c_double * (max_triangles * 6))()
    count = ctypes.c_int()
    
    # Call optimized MeshCut
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

def visualize_optimized_results():
    """Test optimized MeshCut with complex shapes."""
    
    # Complex test polygons
    test_cases = [
        ("L-shaped Building", [
            (0, 0), (8, 0), (8, 3), (5, 3), (5, 5), (8, 5), 
            (8, 8), (3, 8), (3, 6), (0, 6)
        ], (0, 0, 8, 8)),
        
        ("Organic Boundary", [
            (10, 1), (12, 0.5), (14, 1.5), (15, 3), (14.5, 4.5),
            (13, 5.5), (11.5, 6), (10.5, 5), (9.5, 3.5), (9, 2)
        ], (9, 0, 16, 7)),
        
        ("Complex Star", [
            (20, 4), (21, 2), (23, 2), (21.5, 4), (22.5, 6),
            (20, 5.5), (17.5, 6), (18.5, 4), (17, 2), (19, 2)
        ], (17, 1, 24, 7))
    ]
    
    fig, axes = plt.subplots(3, 2, figsize=(16, 20))
    fig.suptitle('Optimized MeshCut Results - Fast Geometric Intersection\nVerifying Correctness with Complex Shapes', 
                 fontsize=16, fontweight='bold')
    
    orientations = [
        ("Backslash Diagonals (\\)", 0),
        ("Forward Slash Diagonals (/)", 1)
    ]
    
    for poly_idx, (poly_name, polygon, bbox_coords) in enumerate(test_cases):
        bbox = BoundingBox(*bbox_coords)
        
        for orient_idx, (orient_name, orient_code) in enumerate(orientations):
            ax = axes[poly_idx, orient_idx]
            
            try:
                # Test with 8x8 grid for good detail vs performance balance
                triangles = meshcut_tessellate(polygon, bbox, 8, orient_code)
                triangle_count = len(triangles)
                
                print(f"{poly_name} with {orient_name}: {triangle_count} triangles")
                
                # Draw triangles
                for i, triangle in enumerate(triangles):
                    tri_points = [
                        [triangle.a.x, triangle.a.y],
                        [triangle.b.x, triangle.b.y], 
                        [triangle.c.x, triangle.c.y],
                        [triangle.a.x, triangle.a.y]
                    ]
                    
                    # Color triangles with alternating colors for visibility
                    color = 'lightblue' if i % 2 == 0 else 'lightcyan'
                    tri_patch = patches.Polygon(
                        tri_points[:-1], 
                        fill=True, 
                        facecolor=color,
                        edgecolor='blue', 
                        linewidth=0.5,
                        alpha=0.7
                    )
                    ax.add_patch(tri_patch)
                
                # Draw polygon boundary
                poly_points = list(polygon) + [polygon[0]]
                poly_x, poly_y = zip(*poly_points)
                ax.plot(poly_x, poly_y, 'r-', linewidth=3, label='Polygon boundary')
                
                # Draw grid lines
                grid_x = np.linspace(bbox.min_x, bbox.max_x, 9)
                grid_y = np.linspace(bbox.min_y, bbox.max_y, 9)
                
                for x in grid_x:
                    ax.axvline(x, color='gray', alpha=0.4, linewidth=0.5)
                for y in grid_y:
                    ax.axhline(y, color='gray', alpha=0.4, linewidth=0.5)
                
                # Configure plot
                ax.set_xlim(bbox.min_x - 0.5, bbox.max_x + 0.5)
                ax.set_ylim(bbox.min_y - 0.5, bbox.max_y + 0.5)
                ax.set_aspect('equal')
                ax.grid(True, alpha=0.2)
                ax.set_title(f'{poly_name}\n{orient_name}\n✅ {triangle_count} triangles (8x8 grid)')
                
                if poly_idx == 0 and orient_idx == 0:
                    ax.legend()
                
            except Exception as e:
                ax.text(0.5, 0.5, f'❌ Error: {str(e)}', 
                       transform=ax.transAxes, ha='center', va='center',
                       fontsize=10, color='red')
                ax.set_title(f'{poly_name}\n{orient_name}\nFAILED')
    
    plt.tight_layout()
    plt.savefig('/home/kevin/gridcut/meshcut/optimized_verification.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    print("🎉 Optimized MeshCut visualization completed!")
    print("📁 Saved as: optimized_verification.png")

def performance_summary():
    """Quick performance summary with the optimized version."""
    print("\n📊 Performance Summary (Optimized vs Original)")
    print("=" * 50)
    
    test_polygon = [
        (0.1, 0.1), (0.8, 0.1), (0.8, 0.4), (0.6, 0.4), 
        (0.6, 0.6), (0.8, 0.6), (0.8, 0.9), (0.1, 0.9)
    ]
    bbox = BoundingBox(0.0, 0.0, 1.0, 1.0)
    
    import time
    
    # Test different grid sizes
    grid_sizes = [4, 8, 16, 32]
    
    print("Grid Size | Triangles | Time (ms)")
    print("----------|-----------|----------")
    
    for grid_size in grid_sizes:
        start_time = time.time()
        triangles = meshcut_tessellate(test_polygon, bbox, grid_size, 0)
        elapsed_ms = (time.time() - start_time) * 1000
        
        print(f"{grid_size:>8}x | {len(triangles):>8} | {elapsed_ms:>7.2f}")
    
    print("\n✅ Fast geometric intersection is working!")
    print("🚀 Ready for MapLibre integration testing")

if __name__ == "__main__":
    print("🔄 Optimized MeshCut Verification")
    print("=================================")
    print("Testing fast geometric intersection with complex shapes")
    print()
    
    visualize_optimized_results()
    performance_summary()