import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, DefaultDict
from collections import defaultdict
import numpy as np

def parse_svg_points(points_str: str) -> List[Tuple[float, float]]:
    """Parse SVG points string into list of (x,y) tuples."""
    points = []
    coords = points_str.split()
    for i in range(0, len(coords), 2):
        x = float(coords[i])
        y = float(coords[i+1])
        points.append((x, y))
    return points

def read_haken_svg(svg_path: str) -> List[Tuple[float, float]]:
    """Read Haken SVG file and return list of points."""
    tree = ET.parse(svg_path)
    root = tree.getroot()
    
    # Find the polyline with class "st1" which contains the main curve
    for polyline in root.findall('.//{http://www.w3.org/2000/svg}polyline'):
        if polyline.get('class') == 'st0':
            points_str = polyline.get('points')
            if points_str is None:
                raise ValueError("Points attribute is missing in the SVG polyline")
            return parse_svg_points(points_str)
    
    raise ValueError("Could not find the main curve in the SVG file")

def points_to_edges(points: List[Tuple[float, float]]) -> List[Tuple[int, int]]:
    """Convert list of points into list of edges (node indices)."""
    return [(i,(i+1)%len(points)) for i in range(len(points))]

def find_intersection(A, B, C, D):
    """Find the intersection point of line segments AB and CD.
    Returns the intersection point as a tuple (x,y) or None if lines are parallel."""
    
    # Handle vertical lines
    if A[0] == B[0] and C[0] == D[0]:  # Both lines are vertical
        return None
    
    if A[0] == B[0]:  # AB is vertical
        x = A[0]
        # Find y where CD crosses x
        m2 = (D[1] - C[1]) / (D[0] - C[0])
        y = m2 * (x - C[0]) + C[1]
    elif C[0] == D[0]:  # CD is vertical
        x = C[0]
        # Find y where AB crosses x
        m1 = (B[1] - A[1]) / (B[0] - A[0])
        y = m1 * (x - A[0]) + A[1]
    else:
        # Both lines are non-vertical
        m1 = (B[1] - A[1]) / (B[0] - A[0])
        m2 = (D[1] - C[1]) / (D[0] - C[0])
        
        if m1 == m2:  # Lines are parallel
            return None
            
        # Find intersection point
        x = (m1 * A[0] - m2 * C[0] + C[1] - A[1]) / (m1 - m2)
        y = m1 * (x - A[0]) + A[1]
    
    # Check if intersection point is within both line segments
    if (min(A[0], B[0]) <= x <= max(A[0], B[0]) and
        min(A[1], B[1]) <= y <= max(A[1], B[1]) and
        min(C[0], D[0]) <= x <= max(C[0], D[0]) and
        min(C[1], D[1]) <= y <= max(C[1], D[1])):
        return (x, y)
    return None

def sort_points_along_edge(first: Tuple[float, float], second: Tuple[float, float], points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """Sort points along the edge from first to second vertex.
    Returns points ordered by their distance from the first vertex."""
    # Calculate the direction vector of the edge
    dx = second[0] - first[0]
    dy = second[1] - first[1]
    
    # Calculate the parameter t for each point along the edge
    # t = 0 at first vertex, t = 1 at second vertex
    def get_t(point):
        if dx != 0:
            t = (point[0] - first[0]) / dx
        else:
            t = (point[1] - first[1]) / dy
        return t
    
    # Sort points by their t value
    sorted_points = sorted(points, key=get_t)
    return sorted_points

def find_self_crossings(points: List[Tuple[float, float]], edges: List[Tuple[int, int]]) -> Tuple[List[Tuple[Tuple[float, float], Tuple[int, int], Tuple[int, int]]], Dict[Tuple[int, int], List[Tuple[float, float]]]]:
    """Find all self-crossings in the curve.
    Returns:
        - List of tuples (intersection_point, edge1, edge2)
        - Dictionary mapping each edge to its sorted crossing points"""
    crossings_on_edge: Dict[Tuple[int, int], List[Tuple[float, float]]] = defaultdict(list)
    crossings = []
    n = len(edges)
    
    for i in range(n):
        for j in range(n):  # Check all pairs
            if i == j or i == (j+1)%n or i == (j-1)%n:
                continue
            edge1 = edges[i]
            edge2 = edges[j]
            
            # Get the points for each edge
            A = points[edge1[0]]
            B = points[edge1[1]]
            C = points[edge2[0]]
            D = points[edge2[1]]
            
            # Check if the edges intersect
            intersection = find_intersection(A, B, C, D)
            if intersection:
                crossings.append((intersection, edge1, edge2))
                crossings_on_edge[edge1].append(intersection)
    
    # Sort crossings along each edge
    for edge in crossings_on_edge:
        first = points[edge[0]]
        second = points[edge[1]]
        crossings_on_edge[edge] = sort_points_along_edge(first, second, crossings_on_edge[edge])
    
    return crossings, crossings_on_edge


def create_subdivided_edges(points: List[Tuple[float, float]], edges: List[Tuple[int, int]], crossings_on_edge: Dict[Tuple[int, int], List[Tuple[float, float]]]):
    """Create new edges that include all crossing points.
    Returns:
        - New list of points (original points + crossing points)
        - New list of edges connecting all points in order"""
    # Create a mapping from original points to their indices
    new_points = points.copy()
    new_edges = []
    is_crossing = [False]*len(points)
    
    # For each original edge, create new edges through its crossing points
    for edge in edges:
        first = points[edge[0]]
        second = points[edge[1]]
        crossings = crossings_on_edge.get(edge, [])
        
        # Get all points along this edge in order
        all_points = [first] + crossings + [second]

        if len(crossings) == 0:
            new_edges.append(edge)
            continue
    
        len_before = len(new_points)
        print(len_before)
        
        # Create edges between consecutive points
        for i in range(len(all_points) - 1):
            # Add points to new_points
            p1 = all_points[i]
            
            if i == 0:
                new_edges.append((edge[0], len(new_points)))
            elif i == len(all_points) - 2:
                new_points.append(p1)
                is_crossing.append(True)
                new_edges.append((len(new_points) - 1, edge[1]))
            else:
                new_points.append(p1)
                is_crossing.append(True)
                new_edges.append((len(new_points) - 1, len(new_points)))
            
        len_after = len(new_points)
        print(len_after)
        print(crossings)
        assert len_after == len_before + len(crossings)
    
    return new_points, new_edges, is_crossing

def assign_z_values(points: List[Tuple[float, float]], edges: List[Tuple[int, int]], is_crossing: List[bool], gauss_code: List[int]) -> Tuple[List[Tuple[float, float, float]], List[int]]:
    """Assign z-values to points based on Gauss code.
    Returns list of 3D points (x,y,z) where z is determined by Gauss code."""
    points_3d = [(x, y, 0.0) for x, y in points]
    count = 0
    
    # Start from edge 0
    current_edge_index = 0
    visited = set()
    is_crossing_array = np.array(is_crossing)
    num_crossings = np.sum(is_crossing_array)
    print(num_crossings)

    assigned_gauss_code = [0] * len(points)
    
    for edge in edges:
        prev_vertex = edge[0]
        next_vertex = edge[1]
        if is_crossing[prev_vertex]:
            print(count, gauss_code[count], prev_vertex)
            z_value = 100.0 if gauss_code[count] < 0 else -100.0
            points_3d[prev_vertex] = (points_3d[prev_vertex][0], points_3d[prev_vertex][1], z_value)
            assigned_gauss_code[prev_vertex] = gauss_code[count]
            count += 1

    return points_3d, assigned_gauss_code

if __name__ == "__main__":
    # Read the SVG and get points and edges
    points = read_haken_svg("ochiai2.svg")
    edges = points_to_edges(points)
    
    # Find self-crossings
    crossings, crossings_on_edge = find_self_crossings(points, edges)
    
    print(f"Number of nodes: {len(points)}")
    print(f"Number of edges: {len(edges)}")
    print(f"Number of self-crossings: {len(crossings)}")
    
    # Create new edges that include all crossing points
    new_points, new_edges, is_crossing = create_subdivided_edges(points, edges, crossings_on_edge)
    
    print(f"\nAfter subdivision:")
    print(f"Number of nodes: {len(new_points)}")
    print(f"Number of edges: {len(new_edges)}")
    
    # Assign z-values based on Gauss code
    gauss_code = "-1 2 -3 -4 5 -6 -7 8 -9 10 -11 -12 13 -14 15 16 4 -17 6 18 -19 11 20 21 -22 -23 24 25 -26 27 -28 29 14 30 -31 -13 -25 32 -33 -34 -27 28 35 36 -37 22 38 -39 12 -40 -41 42 -16 3 -2 1 -21 -38 23 37 45 -35 34 26 -29 -15 -42 43 -18 9 -8 7 17 -5 -43 41 -44 31 -30 44 40 19 -10 -20 39 -24 -32 33 -36 -45".split()
    gauss_code = [int(x) for x in gauss_code]
    print(gauss_code)
    
    points_3d, assigned_gauss_code = assign_z_values(new_points, new_edges, is_crossing, gauss_code)

    # detect consecutive [+, +, +, -, -, +]
    for i in range(len(gauss_code)):
        if gauss_code[i] > 0 and gauss_code[i+1] > 0 and gauss_code[i+2] > 0 and gauss_code[i+3] < 0 and gauss_code[i+4] < 0 and gauss_code[i+5] > 0:
            print(f"Three consecutive pluses at {gauss_code[i]}, {gauss_code[i+1]}, {gauss_code[i+2]}")
    
    # Visualize using polyscope
    points_array = np.array(points_3d)
    edges_array = np.array(new_edges)
    is_crossing_array = np.array(is_crossing)

    # shrink the points so that it fits in a -1 to 1 cube
    points_array = (points_array - np.min(points_array)) / (np.max(points_array) - np.min(points_array))
    points_array = points_array * 2 - 1
    gc_array = np.array(assigned_gauss_code)

    # save it to an obj file
    with open("./ochiai2.obj", "w") as f:
        for i in range(len(points_array)):
            f.write(f"v {points_array[i][0]} {points_array[i][1]} {points_array[i][2]}\n")
        for edge in new_edges:
            f.write(f"l {edge[0]+1} {edge[1]+1}\n")

    import polyscope as ps
    ps.init()
    curve = ps.register_curve_network("ochiai2_subdivided", points_array, edges_array)
    print(points_array.shape)
    print(is_crossing_array.shape)
    curve.add_scalar_quantity("is_crossing", is_crossing_array)
    curve.add_scalar_quantity("gauss_code", gc_array)
    
    ps.show()