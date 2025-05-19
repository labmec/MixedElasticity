import gmsh
import math

gmsh.initialize()

gmsh.model.add("honeycomb")

ndiv = 5

def generate_honeycomb_mesh(radius, rows, cols):
    hex_height = math.sqrt(3) * radius
    hex_width = 2 * radius
    y_offset = hex_height / 2

    for row in range(rows):
        for col in range(cols):
            x = col * hex_width * 0.75
            y = row * hex_height
            if col % 2 == 1:
                y += y_offset

            # Create hexagon points
            points = []
            for i in range(6):
                angle = math.pi / 3 * i
                px = x + radius * math.cos(angle)
                py = y + radius * math.sin(angle)
                points.append(gmsh.model.geo.addPoint(px, py, 0))

            # Create hexagon lines
            lines = []
            for i in range(6):
                lines.append(gmsh.model.geo.addLine(points[i], points[(i + 1) % 6]))

            # Create a closed curve and surface
            curve_loop = gmsh.model.geo.addCurveLoop(lines)
            gmsh.model.geo.addPlaneSurface([curve_loop])

def generate_honeycomb_mesh_hex(radius, divisions):

    cols = divisions*2-1
    hex_height = math.sqrt(3) * radius
    hex_width = 2 * radius
    y_offset = 0

    rows = divisions
    for col in range(cols):
        x = col * hex_width * 0.75
        for row in range(rows):
            y = row * hex_height + y_offset

            # Create hexagon points
            points = []
            for i in range(6):
                angle = math.pi / 3 * i
                px = x + radius * math.cos(angle)
                py = y + radius * math.sin(angle)
                points.append(gmsh.model.geo.addPoint(px, py, 0))

            # Create hexagon lines
            lines = []
            for i in range(6):
                lines.append(gmsh.model.geo.addLine(points[i], points[(i + 1) % 6]))

            # Create a closed curve and surface
            curve_loop = gmsh.model.geo.addCurveLoop(lines)
            gmsh.model.geo.addPlaneSurface([curve_loop])
        if (col >= divisions-1):
            y_offset += hex_height / 2
            rows -= 1
        else:
            y_offset -= hex_height / 2
            rows += 1

# Generate the honeycomb mesh
# generate_honeycomb_mesh(radius=1.0, rows=1, cols=2)
generate_honeycomb_mesh_hex(radius=1.0, divisions=ndiv)

gmsh.model.geo.synchronize()

# Add a physical group for the boundary (all lines)
boundary_entities = gmsh.model.getEntities(1)  # Get all 1D lines
gmsh.model.addPhysicalGroup(1, [line[1] for line in boundary_entities], tag=2)
gmsh.model.setPhysicalName(1, 2, "hexagon")

# Divide each hexagon into 6 equilateral triangles
entities = gmsh.model.getEntities(2)  # Get all 2D surfaces (hexagons)
for entity in entities:
    surface_tag = entity[1]
    boundary = gmsh.model.getBoundary([(2, surface_tag)], oriented=False)
    
    # Calculate the center of the hexagon
    x_sum, y_sum, z_sum = 0, 0, 0
    for b in boundary:
        point_coords = gmsh.model.getValue(0, b[1], [])
        x_sum += point_coords[0]
        y_sum += point_coords[1]
        z_sum += point_coords[2]
    num_points = len(boundary)
    center = gmsh.model.geo.addPoint(x_sum / num_points, y_sum / num_points, z_sum / num_points)
    
    # Create triangles
    for i in range(len(boundary)):
        p1 = boundary[i][1]
        p2 = boundary[(i + 1) % len(boundary)][1]
        gmsh.model.geo.addPlaneSurface([gmsh.model.geo.addCurveLoop([
            gmsh.model.geo.addLine(center, p1),
            gmsh.model.geo.addLine(p1, p2),
            gmsh.model.geo.addLine(p2, center)
        ])])

gmsh.model.geo.removeAllDuplicates()

gmsh.model.geo.synchronize()

# Create physical groups for the domain and boundary
surface_entities = gmsh.model.getEntities(2)  # Get all 2D surfaces

# Add a physical group for the domain (all triangle surfaces)
triangular_surfaces = []
for surface in surface_entities:
    boundary = gmsh.model.getBoundary([surface], oriented=False)
    if len(boundary) == 3:  # Check if the surface has 3 edges (triangle)
        triangular_surfaces.append(surface[1])

gmsh.model.addPhysicalGroup(2, triangular_surfaces, tag=1)
gmsh.model.setPhysicalName(2, 1, "domain")

# Remove all hexagon surfaces (non-triangular surfaces)
for surface in surface_entities:
    boundary = gmsh.model.getBoundary([surface], oriented=False)
    if len(boundary) != 3:  # Check if the surface is not a triangle
        gmsh.model.removeEntities([surface])

# Set transfinite line with 2 divisions for all lines
lines = gmsh.model.getEntities(1)  # Get all 1D entities (lines)
for line in lines:
    gmsh.model.mesh.setTransfiniteCurve(line[1], 2)

# Generate the 2D mesh
gmsh.model.mesh.generate(2)

gmsh.model.mesh.removeDuplicateElements()
gmsh.model.mesh.removeDuplicateNodes()

# Save the mesh to a file
gmsh.write("honeycomb.msh")

gmsh.fltk.run()
gmsh.finalize()