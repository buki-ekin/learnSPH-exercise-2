import math
import struct
import sys
from pathlib import Path

import bpy
import mathutils


def read_points(path):
    data = path.read_bytes()
    header_start = data.index(b"POINTS ")
    header_end = data.index(b"\n", header_start)
    number_of_points = int(data[header_start:header_end].split()[1])
    point_start = header_end + 1
    point_end = point_start + number_of_points * 3 * 8
    values = struct.unpack(">" + "d" * number_of_points * 3, data[point_start:point_end])
    return [values[i:i + 3] for i in range(0, len(values), 3)]


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
particle_mesh.from_pydata(read_points(vtk_files[0]), [], [])
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
mesh_to_points.inputs["Radius"].default_value = 0.028
ico_sphere = nodes.new("GeometryNodeMeshIcoSphere")
ico_sphere.inputs["Radius"].default_value = 0.028
ico_sphere.inputs["Subdivisions"].default_value = 2
instance = nodes.new("GeometryNodeInstanceOnPoints")
realize = nodes.new("GeometryNodeRealizeInstances")

links.new(group_input.outputs["Geometry"], mesh_to_points.inputs["Mesh"])
links.new(mesh_to_points.outputs["Points"], instance.inputs["Points"])
links.new(ico_sphere.outputs["Mesh"], instance.inputs["Instance"])
links.new(instance.outputs["Instances"], realize.inputs["Geometry"])

modifier = particle_object.modifiers.new("ParticleSpheres", "NODES")
modifier.node_group = node_group

# A glossy blue-purple-pink-orange-yellow gradient like the reference video.
# The colour uses x+y, so it remains attached to a particle while it falls.
material = bpy.data.materials.new("RainbowPositionColor")
material.use_nodes = True
material.diffuse_color = (0.8, 0.08, 0.45, 1.0)
principled = material.node_tree.nodes.get("Principled BSDF")
principled.inputs["Roughness"].default_value = 0.2
principled.inputs["Metallic"].default_value = 0.12
shader_nodes = material.node_tree.nodes
shader_links = material.node_tree.links
geometry = shader_nodes.new("ShaderNodeNewGeometry")
separate = shader_nodes.new("ShaderNodeSeparateXYZ")
add_xy = shader_nodes.new("ShaderNodeMath")
add_xy.operation = "ADD"
map_range = shader_nodes.new("ShaderNodeMapRange")
color_ramp = shader_nodes.new("ShaderNodeValToRGB")
map_range.inputs["From Min"].default_value = -2.0
map_range.inputs["From Max"].default_value = 2.0
color_ramp.color_ramp.interpolation = "B_SPLINE"
color_ramp.color_ramp.elements[0].color = (0.005, 0.04, 0.8, 1.0)
purple = color_ramp.color_ramp.elements.new(0.27)
purple.color = (0.35, 0.02, 0.8, 1.0)
pink = color_ramp.color_ramp.elements.new(0.50)
pink.color = (1.0, 0.02, 0.35, 1.0)
orange = color_ramp.color_ramp.elements.new(0.73)
orange.color = (1.0, 0.25, 0.02, 1.0)
color_ramp.color_ramp.elements[-1].color = (1.0, 0.72, 0.08, 1.0)
shader_links.new(geometry.outputs["Position"], separate.inputs["Vector"])
shader_links.new(separate.outputs["X"], add_xy.inputs[0])
shader_links.new(separate.outputs["Y"], add_xy.inputs[1])
shader_links.new(add_xy.outputs[0], map_range.inputs["Value"])
shader_links.new(map_range.outputs["Result"], color_ramp.inputs["Fac"])
shader_links.new(color_ramp.outputs["Color"], principled.inputs["Base Color"])
particle_object.data.materials.append(material)

set_material = nodes.new("GeometryNodeSetMaterial")
set_material.inputs["Material"].default_value = material
links.new(realize.outputs["Geometry"], set_material.inputs["Geometry"])
links.new(set_material.outputs["Geometry"], group_output.inputs["Geometry"])


# Display the same finite floor used by the simulation.
bpy.ops.mesh.primitive_plane_add(
    size=5.0,
    location=(0.5, 0.5, -0.30),
)
floor = bpy.context.object
floor.name = "Simulation floor"

floor_material = bpy.data.materials.new("Floor grey")
floor_material.use_nodes = True
floor_shader = floor_material.node_tree.nodes.get("Principled BSDF")
floor_shader.inputs["Base Color"].default_value = (
    0.25, 0.28, 0.32, 1.0
)
floor_shader.inputs["Roughness"].default_value = 0.8
floor.data.materials.append(floor_material)
# Fixed oblique view: show approach along x and spreading along y and z.
bpy.ops.object.camera_add(location=(5.5, -8.0, 4.5))
camera = bpy.context.object
look_at(camera, mathutils.Vector((1.08, 0.48, 0.48)))
camera.data.type = "ORTHO"
camera.data.ortho_scale = 6.5
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
scene.render.fps = 30
scene.render.film_transparent = False
scene.world.color = (0.003, 0.008, 0.02)
# Inspect the initial fall, impact and spreading over one second.
selected_frames = (0, 6, 15, 30)

for frame in selected_frames:
    vtk_file = vtk_folder / f"exercise1_{frame:04d}.vtk"

    if not vtk_file.is_file():
        raise FileNotFoundError(f"Missing snapshot: {vtk_file}")

    points = read_points(vtk_file)

    if len(points) != len(particle_mesh.vertices):
        raise ValueError("Particle count changed between snapshots.")

    coordinates = [value for point in points for value in point]
    particle_mesh.vertices.foreach_set("co", coordinates)
    particle_mesh.update()

    scene.frame_set(frame + 1)
    scene.render.filepath = str(
        output_folder / f"collision_{frame:04d}.png"
    )

    # Save each snapshot as an independently editable Blender scene.
    bpy.ops.wm.save_as_mainfile(
        filepath=str(output_folder / f"collision_{frame:04d}.blend")
    )

    bpy.ops.render.render(write_still=True)
    print(f"Rendered collision at t = {frame / 30.0:.1f} s")