import math

def tile_to_latlon(z, x, y, center=True):
    """
    Convert slippy map tile (z/x/y) to latitude and longitude.
    If center=True returns the center of the tile, otherwise the top-left corner.
    """
    n = 2 ** z
    x_coord = x + 0.5 if center else x
    y_coord = y + 0.5 if center else y

    lon_deg = x_coord / n * 360.0 - 180.0
    lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * y_coord / n)))
    lat_deg = math.degrees(lat_rad)

    return lat_deg, lon_deg

# Example:
z, x, y = 14, 8710, 5744
lat, lon = tile_to_latlon(z, x, y)
print(f"Tile center at z={z}, x={x}, y={y} → lat={lat:.6f}, lon={lon:.6f}")
