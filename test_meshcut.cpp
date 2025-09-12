#include "meshcut/meshcut.hpp"
#include <iostream>
#include <cassert>

int main() {
    std::cout << "🔄 Testing MeshCut Library\n";
    std::cout << "=" << std::string(30, '=') << "\n";
    
    // Create a simple test polygon (square)
    meshcut::Polygon polygon = {
        {2.0, 2.0},
        {6.0, 2.0},
        {6.0, 6.0}, 
        {2.0, 6.0}
    };
    
    // Create bounding box
    meshcut::Polygon bbox = {
        {0.0, 0.0},   // bottom-left
        {8.0, 0.0},   // bottom-right  
        {8.0, 8.0},   // top-right
        {0.0, 8.0}    // top-left
    };
    
    // Test basic tessellation
    std::cout << "🔍 Testing basic tessellation...\n";
    auto triangles = meshcut::tessellate(polygon, bbox, 4);
    std::cout << "✅ Generated " << triangles.size() << " triangles\n";
    assert(triangles.size() > 0);
    
    // Test both orientations
    std::cout << "\n🔍 Testing triangle orientations...\n";
    auto triangles_tl_br = meshcut::tessellate(
        polygon, bbox, 4, meshcut::TriangleOrientation::TOPLEFT_BOTTOMRIGHT
    );
    auto triangles_tr_bl = meshcut::tessellate(
        polygon, bbox, 4, meshcut::TriangleOrientation::TOPRIGHT_BOTTOMLEFT  
    );
    
    std::cout << "✅ TL-BR orientation (\\): " << triangles_tl_br.size() << " triangles\n";
    std::cout << "✅ TR-BL orientation (/): " << triangles_tr_bl.size() << " triangles\n";
    assert(triangles_tl_br.size() == triangles_tr_bl.size());
    
    // Test tile tessellation 
    std::cout << "\n🔍 Testing tile tessellation...\n";
    meshcut::Polygon tile_polygon = {
        {0.25, 0.25},
        {0.75, 0.25}, 
        {0.75, 0.75},
        {0.25, 0.75}
    };
    
    auto tile_triangles = meshcut::tessellate_tile(tile_polygon, 8);
    std::cout << "✅ Tile tessellation: " << tile_triangles.size() << " triangles\n";
    assert(tile_triangles.size() > 0);
    
    std::cout << "\n🎉 All MeshCut tests passed!\n";
    std::cout << "📊 Library ready for MapLibre Native integration\n";
    
    return 0;
}