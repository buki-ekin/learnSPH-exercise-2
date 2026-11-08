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
# Scan all snapshots for a shared speed scale and fixed camera bounds.
speed_scale_max = 1e-12
bounds_min = [-0.01, -0.01, -0.01]
bounds_max = [0.16, 0.81, 1.01]

for vtk_file in vtk_files:
    scan_points, scan_speeds = read_particles(vtk_file)
    speed_scale_max = max(
        speed_scale_max, max(scan_speeds, default=0.0)
    )

    for point in scan_points:
        for axis in range(3):
            bounds_min[axis] = min(bounds_min[axis], point[axis])
            bounds_max[axis] = max(bounds_max[axis], point[axis])

print(f"Shared speed scale: 0 to {speed_scale_max:.3f} m/s")
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
# Fixed oblique camera framing the entire recorded trajectory.
view_min = mathutils.Vector(bounds_min)
view_max = mathutils.Vector(bounds_max)
view_center = (view_min + view_max) * 0.5
view_diagonal = (view_max - view_min).length

bpy.ops.object.camera_add(
    location=view_center + mathutils.Vector((-3.0, -0.6, 0.8))
)
camera = bpy.context.object
look_at(camera, view_center)
camera.data.type = "ORTHO"
camera.data.ortho_scale = max(1.2, view_diagonal * 1.15)
camera.data.clip_start = 0.001
camera.data.clip_end = 100.0
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
# Playback at 15 FPS gives 151 / 15 = 10.07 seconds.
# Simulation snapshots were exported at 50 FPS.
scene.render.fps = 15
scene.render.film_transparent = False
scene.world.color = (0.003, 0.008, 0.02)

for frame, vtk_file in enumerate(vtk_files):
    points, speeds = read_particles(vtk_file)

    # Rebuild only the point mesh to support particle deletion.
    # The object and its Geometry Nodes modifier remain unchanged.
    old_mesh = particle_object.data
    new_mesh = bpy.data.meshes.new(f"ParticlePoints_{frame:04d}")
    new_mesh.from_pydata(points, [], [])

    new_speed_attribute = new_mesh.attributes.new(
        name="speed",
        type="FLOAT",
        domain="POINT",
    )
    new_speed_attribute.data.foreach_set("value", speeds)
    new_mesh.materials.append(material)
    new_mesh.update()

    particle_object.data = new_mesh
    bpy.data.meshes.remove(old_mesh)

    scene.frame_set(frame + 1)
    scene.render.filepath = str(
        output_folder / f"custom_{frame:04d}.png"
    )

    print(
        f"Rendering {frame + 1}/{len(vtk_files)} | "
        f"simulation time={frame / 50.0:.2f} s | "
        f"particles={len(points)} | "
        f"max speed={max(speeds, default=0.0):.3f} m/s",
        flush=True,
    )

    bpy.ops.render.render(write_still=True)

# Save the final snapshot as an editable scene.
# This .blend contains the final state, not an animated mesh sequence.
bpy.ops.wm.save_as_mainfile(
    filepath=str(output_folder / "custom_final.blend")
)

print("All PNG frames rendered successfully.", flush=True)