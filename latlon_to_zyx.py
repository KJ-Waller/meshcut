import math

def latlon_to_tile(z, lat, lon, center=True):
    """
    Convert latitude and longitude to slippy map tile (z/x/y).
    If center=True returns the nearest tile center, otherwise the top-left tile indices.
    """
    n = 2 ** z
    x = (lon + 180.0) / 360.0 * n
    y = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n

    if center:
        x = int(x)
        y = int(y)
    else:
        x = math.floor(x)
        y = math.floor(y)

    return x, y

# Example:
z, lat, lon = 14, 47.31724, 11.43403
x, y = latlon_to_tile(z, lat, lon)
print(f"Lat={lat:.6f}, Lon={lon:.6f} → Tile z={z}, x={x}, y={y}")
