import math
import struct
import sys
from pathlib import Path

import bpy
import mathutils
def read_particles(path):
    # Read positions from the binary VTK file produced by our exporter.
    data = path.read_bytes()
    header_start = data.index(b"POINTS ")
    header_end = data.index(b"\n", header_start)
    header = data[header_start:header_end].split()

    count = int(header[1])
    if header[2] != b"double":
        raise ValueError("Expected double-precision positions.")

    point_start = header_end + 1
    point_end = point_start + count * 3 * 8

    coordinates = struct.unpack(
        f">{count * 3}d", data[point_start:point_end]
    )
    points = [
        coordinates[i:i + 3]
        for i in range(0, len(coordinates), 3)
    ]

    # The exporter stores velocity in the FIELD array named "vector".
    vector_header = f"\nvector 3 {count} double\n".encode("ascii")
    velocity_start = (
        data.index(vector_header, point_end) + len(vector_header)
    )
    velocity_end = velocity_start + count * 3 * 8

    values = struct.unpack(
        f">{count * 3}d", data[velocity_start:velocity_end]
    )

    # Speed is the magnitude of the velocity vector, in m/s.
    speeds = [
        math.sqrt(
            values[i] ** 2
            + values[i + 1] ** 2
            + values[i + 2] ** 2
        )
        for i in range(0, len(values), 3)
    ]

    return points, speeds

def look_at(obj, point):
    direction = point - obj.location
    obj.rotation_euler = direction.to_track_quat("-Z", "Y").to_euler()


args = sys.argv[sys.argv.index("--") + 1:]
vtk_folder = Path(args[0]).resolve()
output_folder = Path(args[1]).resolve()
output_folder.mkdir(parents=True, exist_ok=True)
vtk_files = sorted(vtk_folder.glob("exercise1_*.vtk"))

if not vtk_files:
    raise RuntimeError("No VTK files found.")

# Start with an empty Blender scene.
bpy.ops.object.select_all(action="SELECT")
bpy.ops.object.delete(use_global=False)

# Make one mesh that stores all particle points.
particle_mesh = bpy.data.meshes.new("ParticlePoints")
initial_points, initial_speeds = read_particles(vtk_files[0])
particle_mesh.from_pydata(initial_points, [], [])

# Store one speed value per particle for subsequent color mapping.
speed_attribute = particle_mesh.attributes.new(
    name="speed",
    type="FLOAT",
    domain="POINT",
)
speed_attribute.data.foreach_set("value", initial_speeds)
particle_object = bpy.data.objects.new("Particles", particle_mesh)
bpy.context.collection.objects.link(particle_object)


# Turn every point into a small, smooth-looking sphere.
node_group = bpy.data.node_groups.new("ParticleSpheres", "GeometryNodeTree")
node_group.interface.new_socket(name="Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
node_group.interface.new_socket(name="Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
nodes = node_group.nodes
links = node_group.links
group_input = nodes.new("NodeGroupInput")
group_output = nodes.new("NodeGroupOutput")
mesh_to_points = nodes.new("GeometryNodeMeshToPoints")
mesh_to_points.inputs["Radius"].default_value = 0.0045
ico_sphere = nodes.new("GeometryNodeMeshIcoSphere")
ico_sphere.inputs["Radius"].default_value = 0.0045
ico_sphere.inputs["Subdivisions"].default_value = 2
instance = nodes.new("GeometryNodeInstanceOnPoints")
realize = nodes.new("GeometryNodeRealizeInstances")

links.new(group_input.outputs["Geometry"], mesh_to_points.inputs["Mesh"])
links.new(mesh_to_points.outputs["Points"], instance.inputs["Points"])
links.new(ico_sphere.outputs["Mesh"], instance.inputs["Instance"])
links.new(instance.outputs["Instances"], realize.inputs["Geometry"])

modifier = particle_object.modifiers.new("ParticleSpheres", "NODES")
modifier.node_group = node_group
# Color particles by actual speed using a fixed 0–10 m/s scale.
material = bpy.data.materials.new("ParticleSpeed")
material.use_nodes = True

shader_nodes = material.node_tree.nodes
shader_links = material.node_tree.links
shader_nodes.clear()

output = shader_nodes.new("ShaderNodeOutputMaterial")
emission = shader_nodes.new("ShaderNodeEmission")
emission.inputs["Strength"].default_value = 1.0

attribute = shader_nodes.new("ShaderNodeAttribute")
attribute.attribute_name = "speed"

map_range = shader_nodes.new("ShaderNodeMapRange")
map_range.inputs["From Min"].default_value = 0.0
# Use the final reference snapshot to set one shared scale for both images.
_, reference_speeds = read_particles(
    vtk_folder / "exercise1_0019.vtk"
)
speed_scale_max = max(max(reference_speeds), 1e-12)
map_range.inputs["From Max"].default_value = speed_scale_max
map_range.inputs["To Min"].default_value = 0.0
map_range.inputs["To Max"].default_value = 1.0
map_range.clamp = True

color_ramp = shader_nodes.new("ShaderNodeValToRGB")
color_ramp.color_ramp.interpolation = "LINEAR"
# Reference-style colors: zero speed is dark blue, maximum is white.
color_ramp.color_ramp.elements[0].color = (0.0, 0.01, 0.65, 1.0)
color_ramp.color_ramp.elements[1].color = (1.0, 1.0, 1.0, 1.0)

shader_links.new(attribute.outputs["Fac"], map_range.inputs["Value"])
shader_links.new(map_range.outputs["Result"], color_ramp.inputs["Fac"])
shader_links.new(color_ramp.outputs["Color"], emission.inputs["Color"])
shader_links.new(emission.outputs["Emission"], output.inputs["Surface"])

particle_object.data.materials.append(material)

set_material = nodes.new("GeometryNodeSetMaterial")
set_material.inputs["Material"].default_value = material
links.new(realize.outputs["Geometry"], set_material.inputs["Geometry"])
links.new(set_material.outputs["Geometry"], group_output.inputs["Geometry"])

# Show the simulated tank as an open-top cutaway.
# Hide the camera-facing x-min wall in the render only.
tank_vertices = [
    (-0.01, -0.01, -0.01),
    ( 0.16, -0.01, -0.01),
    ( 0.16,  0.81, -0.01),
    (-0.01,  0.81, -0.01),
    (-0.01, -0.01,  1.01),
    ( 0.16, -0.01,  1.01),
    ( 0.16,  0.81,  1.01),
    (-0.01,  0.81,  1.01),
]

tank_faces = [
    (0, 1, 2, 3),  # Bottom.
    (1, 5, 6, 2),  # Far side: maximum x.
    (0, 4, 5, 1),  # Minimum y.
    (3, 2, 6, 7),  # Maximum y.
]

tank_mesh = bpy.data.meshes.new("TankCutaway")
tank_mesh.from_pydata(tank_vertices, [], tank_faces)
tank_mesh.update()

tank = bpy.data.objects.new("TankCutaway", tank_mesh)
bpy.context.collection.objects.link(tank)

tank_material = bpy.data.materials.new("TankGrey")
tank_material.use_nodes = True
tank_shader = tank_material.node_tree.nodes.get("Principled BSDF")
tank_shader.inputs["Base Color"].default_value = (
    0.16, 0.18, 0.21, 1.0
)
tank_shader.inputs["Roughness"].default_value = 0.9
tank.data.materials.append(tank_material)

# Side view: y runs horizontally and z vertically.
# This exposes upward motion near the far end of the tank.
bpy.ops.object.camera_add(location=(-2.0, 0.40, 0.50))
camera = bpy.context.object
look_at(camera, mathutils.Vector((0.075, 0.40, 0.50)))
camera.data.type = "ORTHO"
camera.data.ortho_scale = 1.20
camera.data.clip_start = 0.001
bpy.context.scene.camera = camera

# Add two simple lights.
bpy.ops.object.light_add(type="AREA", location=(4.0, -3.0, 12.0))
bpy.context.object.data.energy = 1400
bpy.context.object.data.shape = "DISK"
bpy.context.object.data.size = 5.0
bpy.ops.object.light_add(type="AREA", location=(-4.0, 2.0, 5.0))
bpy.context.object.data.energy = 700
bpy.context.object.data.size = 4.0

scene = bpy.context.scene
scene.render.engine = "BLENDER_EEVEE"
scene.render.resolution_x = 1024
scene.render.resolution_y = 1024
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.fps = 50
scene.render.film_transparent = False
scene.world.color = (0.003, 0.008, 0.02)
# Inspect the initial fall, impact and spreading over one second.
# Inspect immediately before and after the velocity spike.
selected_frames = (0, 19)
for frame in selected_frames:
    vtk_file = vtk_folder / f"exercise1_{frame:04d}.vtk"

    if not vtk_file.is_file():
        raise FileNotFoundError(f"Missing snapshot: {vtk_file}")

    points, speeds = read_particles(vtk_file)
    if len(points) != len(particle_mesh.vertices):
        raise ValueError("Particle count changed between snapshots.")

    coordinates = [value for point in points for value in point]
    particle_mesh.vertices.foreach_set("co", coordinates)
    speed_attribute.data.foreach_set("value", speeds)
    particle_mesh.update()

    print(
        f"t={frame / 50.0:.6f} s, "
        f"particles={len(points)}, "
        f"maximum speed={max(speeds):.5f} m/s"
    )

    scene.frame_set(frame + 1)
    scene.render.filepath = str(
        output_folder / f"collision_{frame:04d}.png"
    )

    # Save each snapshot as an independently editable Blender scene.
    bpy.ops.wm.save_as_mainfile(
        filepath=str(output_folder / f"collision_{frame:04d}.blend")
    )

    bpy.ops.render.render(write_still=True)
    print(f"Rendered dam break at t = {frame / 50.0:.2f} s")