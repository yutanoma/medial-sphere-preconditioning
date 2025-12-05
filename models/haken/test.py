import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

def gauss_code_to_planar(gauss_code):
    """
    Convert a Gauss code to a 2D planar layout.
    
    Args:
        gauss_code (str): Gauss code in the format "1+a,2-a,3+a,1-a,2+a,3-a"
        
    Returns:
        tuple: (points, crossings) where points is a list of (x,y) coordinates
               and crossings is a list of crossing information
    """
    # Parse the Gauss code
    elements = [elem.strip() for elem in gauss_code.split(',')]
    
    # Initialize lists to store points and crossings
    points = []
    crossings = {}
    
    # Generate points in a circular layout
    n = len(elements) // 2  # Number of crossings
    radius = 1.0
    for i in range(len(elements)):
        angle = 2 * np.pi * i / len(elements)
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points.append((x, y))
        
        # Store crossing information
        if '+' in elements[i] or '-' in elements[i]:
            crossing_num = int(elements[i].split('+')[0].split('-')[0])
            if crossing_num not in crossings:
                crossings[crossing_num] = {
                    'position': (x, y),
                    'type': '+' if '+' in elements[i] else '-'
                }
    
    return points, crossings

def plot_knot(points, crossings):
    """Plot the knot diagram using matplotlib."""
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot the curve
    x, y = zip(*points)
    ax.plot(x, y, 'b-', alpha=0.5)
    
    # Plot crossings
    for crossing_num, info in crossings.items():
        x, y = info['position']
        # Draw a circle at the crossing
        circle = Circle((x, y), 0.1, fill=False, color='black')
        ax.add_patch(circle)
        
        # Add crossing number
        ax.text(x, y, str(crossing_num), ha='center', va='center')
    
    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()

# Example usage
if __name__ == "__main__":
    # Example Gauss code for a trefoil knot
    gauss_code = "1+a,2-a,3+a,1-a,2+a,3-a"
    
    # Convert to planar layout
    points, crossings = gauss_code_to_planar(gauss_code)
    
    # Plot the result
    plot_knot(points, crossings)