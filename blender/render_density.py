import struct
from pathlib import Path


# Assignment 3(b): read the binary VTK format written by our C++ exporter.
def read_particles_and_scalar(path):
    data = path.read_bytes()

    # Read particle positions stored as big-endian doubles.
    header_start = data.index(b"POINTS ")
    header_end = data.index(b"\n", header_start)
    header = data[header_start:header_end].split()

    count = int(header[1])
    if header[2] != b"double":
        raise ValueError("Expected double-precision particle positions.")

    point_start = header_end + 1
    point_end = point_start + count * 3 * 8

    coordinates = struct.unpack(
        f">{count * 3}d",
        data[point_start:point_end],
    )

    points = [
        coordinates[i:i + 3]
        for i in range(0, len(coordinates), 3)
    ]

    # Our exporter names the per-particle scalar field "scalar".
    # In the fluid files, this field contains density in kg/m^3.
    scalar_header = f"\nscalar 1 {count} double\n".encode("ascii")
    scalar_start = (
        data.index(scalar_header, point_end) + len(scalar_header)
    )

    scalars = struct.unpack(
        f">{count}d",
        data[scalar_start:scalar_start + count * 8],
    )

    return points, scalars


# Resolve input paths relative to the script, not the terminal location.
project_folder = Path(__file__).resolve().parent.parent
vtk_folder = project_folder / "output" / "assignment3"

for filename in (
    "fluid_density_without_boundary.vtk",
    "fluid_density_with_boundary.vtk",
):
    points, densities = read_particles_and_scalar(vtk_folder / filename)

    print(f"\n{filename}")
    print(f"Particle count: {len(points)}")
    print(f"Minimum density [kg/m^3]: {min(densities):.3f}")
    print(f"Maximum density [kg/m^3]: {max(densities):.3f}")
    
    import bpy
import bmesh
from mathutils import Vector


# Assignment 3(b): start a separate density visualization scene.
bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete(use_global=False)

output_folder = project_folder / "output" / "assignment3_blender"
output_folder.mkdir(parents=True, exist_ok=True)

# Use one fixed density scale for both cases.
density_min = 450.0
density_max = 1200.0
material_count = 64


def density_color(value):
    # Blue -> cyan -> green -> yellow -> red.
    t = max(0.0, min(1.0, (value - density_min)
                    / (density_max - density_min)))

    stops = [
        (0.05, 0.10, 0.85),
        (0.00, 0.75, 1.00),
        (0.10, 0.80, 0.20),
        (1.00, 0.85, 0.00),
        (0.90, 0.03, 0.02),
    ]

    position = t * (len(stops) - 1)
    index = min(int(position), len(stops) - 2)
    fraction = position - index

    return tuple(
        stops[index][channel] * (1.0 - fraction)
        + stops[index + 1][channel] * fraction
        for channel in range(3)
    ) + (1.0,)


materials = []

for index in range(material_count):
    value = density_min + (
        index / (material_count - 1)
    ) * (density_max - density_min)

    material = bpy.data.materials.new(f"Density_{value:.1f}")
    material.diffuse_color = density_color(value)
    materials.append(material)


# Reuse a small sphere template to avoid thousands of Blender objects.
template = bmesh.new()
bmesh.ops.create_icosphere(
    template,
    subdivisions=1,
    radius=0.0023,
)
template.verts.ensure_lookup_table()
template.verts.index_update()

sphere_vertices = [tuple(vertex.co) for vertex in template.verts]
sphere_faces = [
    tuple(vertex.index for vertex in face.verts)
    for face in template.faces
]
template.free()


def create_density_cutaway(filename, name, x_offset):
    points, densities = read_particles_and_scalar(vtk_folder / filename)

    vertices = []
    faces = []
    face_materials = []

    for point, density in zip(points, densities):
        # Hide the front half to expose the interior.
        # This affects visualization only.
        if point[1] < 0.05:
            continue

        base = len(vertices)

        vertices.extend(
            (
                point[0] + x_offset + vertex[0],
                point[1] + vertex[1],
                point[2] + vertex[2],
            )
            for vertex in sphere_vertices
        )

        faces.extend(
            tuple(base + index for index in face)
            for face in sphere_faces
        )

        normalized = max(
            0.0,
            min(1.0, (density - density_min)
                / (density_max - density_min)),
        )
        material_index = round(normalized * (material_count - 1))
        face_materials.extend([material_index] * len(sphere_faces))

    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(vertices, [], faces)
    mesh.update()

    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)

    for material in materials:
        mesh.materials.append(material)

    for polygon, material_index in zip(mesh.polygons, face_materials):
        polygon.material_index = material_index
        polygon.use_smooth = True

    return obj


create_density_cutaway(
    "fluid_density_without_boundary.vtk",
    "Fluid-only density",
    0.0,
)

create_density_cutaway(
    "fluid_density_with_boundary.vtk",
    "Density including floor",
    0.15,
)


# Show the ground plane only under the case with a boundary.
bpy.ops.mesh.primitive_plane_add(
    size=1.0,
    location=(0.20, 0.05, -0.0025),
)
floor = bpy.context.object
floor.name = "Boundary plane"
floor.scale = (0.124, 0.124, 1.0)

floor_material = bpy.data.materials.new("Boundary grey")
floor_material.diffuse_color = (0.35, 0.35, 0.35, 1.0)
floor.data.materials.append(floor_material)


# Workbench renders material colors without requiring lights or shaders.
scene = bpy.context.scene
scene.render.engine = "BLENDER_WORKBENCH"
scene.display.shading.light = "STUDIO"
scene.display.shading.color_type = "MATERIAL"
scene.display.shading.show_shadows = False
scene.display.shading.show_cavity = True
scene.display.shading.background_type = "WORLD"
scene.world.color = (0.06, 0.06, 0.06)

bpy.ops.object.camera_add(location=(0.125, -0.48, 0.25))
camera = bpy.context.object
target = Vector((0.125, 0.06, 0.045))
camera.rotation_euler = (
    target - camera.location
).to_track_quat("-Z", "Y").to_euler()
camera.data.type = "ORTHO"
camera.data.ortho_scale = 0.34
camera.data.clip_start = 0.001
scene.camera = camera

scene.render.resolution_x = 1600
scene.render.resolution_y = 900
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.filepath = str(output_folder / "density_comparison.png")
# Assignment 3(b): add labels aligned with the camera.
label_material = bpy.data.materials.new("Label white")
label_material.diffuse_color = (1.0, 1.0, 1.0, 1.0)


def add_camera_label(text, x, y, size):
    curve = bpy.data.curves.new("Figure label", type="FONT")
    curve.body = text
    curve.align_x = "CENTER"
    curve.size = size

    obj = bpy.data.objects.new(text, curve)
    bpy.context.collection.objects.link(obj)
    obj.data.materials.append(label_material)

    # Place labels in camera-local coordinates.
    obj.parent = camera
    obj.location = (x, y, -0.30)
    obj.rotation_euler = (0.0, 0.0, 0.0)


add_camera_label(
    "Fluid only", -0.075, 0.075, 0.007
)
add_camera_label(
    "Fluid + boundary", 0.075, 0.075, 0.007
)
add_camera_label(
    "Density [kg/m^3] | Shared scale: 450 to 1200",
    0.0, -0.077, 0.005
)
add_camera_label(
    "Blue: low | Yellow: near 1000 | Red: high",
    0.0, -0.087, 0.0045
)
# Save the editable scene and render the comparison.
bpy.ops.wm.save_as_mainfile(
    filepath=str(output_folder / "density_comparison.blend")
)
bpy.ops.render.render(write_still=True)

print("Left: fluid-only density")
print("Right: density including the floor contribution")
print("Shared density scale: 450–1200 kg/m^3")
print(f"Saved Blender scene and image to: {output_folder}")