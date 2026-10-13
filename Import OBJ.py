import bpy
import json
from pathlib import Path


# ============================================================
# Configuration
# ============================================================

OBJ_PATH = r"C:\Users\Nikolas\Documents\Meshing of Botanical Trees\OBJs\test_tree\uv_test_6_metadata.obj"
JSON_PATH = r"C:\Users\Nikolas\Documents\Meshing of Botanical Trees\OBJs\test_tree\uv_test_6_metadata.json"

COLLECTION_NAME = "Imported OBJ"


# ============================================================
# Attribute definitions
# ============================================================

FLOAT_ATTRIBUTES = [
    "boundary_distance",
    "HC",
    "HL",
    "RW",
    "RB",
    "moisture",
    "root_distance",
    "uv_3",
]

VEC2_ATTRIBUTES = [
    "profile_position",
    "profile_polar_coordinate",
]

VEC3_ATTRIBUTES = [
    "position0",
    "direction0",
]

COLOR_ATTRIBUTES = [
    "color",
]

FACE_ATTRIBUTES = [
    "has_neighbor",
]


# ============================================================
# MTL parsing
# ============================================================

def parse_mtl(path):
    """
    Parse the useful subset of an MTL file.

    Returns:
        {
            material_name: {
                "Kd": [r, g, b],
                "Ks": [r, g, b],
                "Ns": value,
                "d": value,
                "map_Kd": Path or None,
            }
        }
    """

    materials = {}
    current = None

    path = Path(path)

    if not path.exists():
        print(f"WARNING: MTL file not found: {path}")
        return materials

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            # Strip full-line and trailing comments.
            line = line.split("#", 1)[0].strip()

            if not line:
                continue

            parts = line.split()
            keyword = parts[0]

            if keyword == "newmtl":
                name = " ".join(parts[1:])

                current = {
                    "Kd": [0.8, 0.8, 0.8],
                    "Ks": [0.0, 0.0, 0.0],
                    "Ns": 0.0,
                    "d": 1.0,
                    "map_Kd": None,
                }

                materials[name] = current

            elif current is None:
                continue

            elif keyword == "Kd" and len(parts) >= 4:
                current["Kd"] = [
                    float(parts[1]),
                    float(parts[2]),
                    float(parts[3]),
                ]

            elif keyword == "Ks" and len(parts) >= 4:
                current["Ks"] = [
                    float(parts[1]),
                    float(parts[2]),
                    float(parts[3]),
                ]

            elif keyword == "Ns" and len(parts) >= 2:
                current["Ns"] = float(parts[1])

            elif keyword == "d" and len(parts) >= 2:
                current["d"] = float(parts[1])

            elif keyword == "Tr" and len(parts) >= 2:
                # Tr is transparency, so d = 1 - Tr.
                current["d"] = 1.0 - float(parts[1])

            elif keyword == "map_Kd" and len(parts) >= 2:
                # This intentionally handles the common/simple case.
                texture_path = " ".join(parts[1:])
                current["map_Kd"] = path.parent / texture_path

    return materials


# ============================================================
# OBJ parsing
# ============================================================

def parse_obj(path):
    """
    Parse an OBJ file while preserving:

      * global vertex indices
      * texture coordinates
      * normals
      * object sections
      * face ordering
      * material assignment
      * smoothing groups

    Returns:

        vertices:
            [(x, y, z), ...]

        texcoords:
            [(u, v), ...]

        normals:
            [(x, y, z), ...]

        objects:
            [
                {
                    "name": str,
                    "faces": [
                        {
                            "corners": [
                                {
                                    "v": global vertex index,
                                    "vt": global texcoord index or None,
                                    "vn": global normal index or None,
                                },
                                ...
                            ],
                            "material": material name or None,
                            "smoothing": smoothing group,
                        }
                    ]
                }
            ]

        mtllibs:
            list of MTL paths
    """

    vertices = []
    texcoords = []
    normals = []

    objects = []
    mtllibs = []

    current_object = None
    current_material = None
    current_smoothing = None

    path = Path(path)

    def ensure_object():
        nonlocal current_object

        if current_object is None:
            current_object = {
                "name": "Object",
                "faces": [],
            }
            objects.append(current_object)

        return current_object

    def resolve_index(index, count):
        """
        Convert OBJ index to zero-based absolute index.
        """
        if index > 0:
            result = index - 1
        else:
            result = count + index

        if not (0 <= result < count):
            raise RuntimeError(
                f"Invalid OBJ index {index} "
                f"(count = {count})"
            )

        return result

    with open(path, "r", encoding="utf-8", errors="replace") as f:

        for line_number, line in enumerate(f, 1):

            # Strip full-line and trailing comments
            # (e.g. "f 1 2 3 # triangle").
            line = line.split("#", 1)[0].strip()

            if not line:
                continue

            # ------------------------------------------------
            # Vertex
            # ------------------------------------------------

            if line.startswith("v "):
                parts = line.split()

                if len(parts) < 4:
                    raise RuntimeError(
                        f"Invalid vertex at line {line_number}: {line}"
                    )

                vertices.append((
                    float(parts[1]),
                    float(parts[3]),
                    float(parts[2]),
                ))

            # ------------------------------------------------
            # Texture coordinate
            # ------------------------------------------------

            elif line.startswith("vt "):
                parts = line.split()

                if len(parts) < 2:
                    raise RuntimeError(
                        f"Invalid texture coordinate at "
                        f"line {line_number}: {line}"
                    )

                u = float(parts[1])
                v = float(parts[2]) if len(parts) >= 3 else 0.0

                texcoords.append((u, v))

            # ------------------------------------------------
            # Normal
            # ------------------------------------------------

            elif line.startswith("vn "):
                parts = line.split()

                if len(parts) < 4:
                    raise RuntimeError(
                        f"Invalid normal at line {line_number}: {line}"
                    )

                normals.append((
                    float(parts[1]),
                    float(parts[3]),
                    float(parts[2]),
                ))

            # ------------------------------------------------
            # Object
            # ------------------------------------------------

            elif line.startswith("o "):
                name = line[2:].strip()

                current_object = {
                    "name": name if name else f"Object_{len(objects):03d}",
                    "faces": [],
                }

                objects.append(current_object)

            # ------------------------------------------------
            # Material library
            # ------------------------------------------------

            elif line.startswith("mtllib "):
                mtl_name = line[len("mtllib "):].strip()
                # OBJ filenames may be quoted.
                if (
                    len(mtl_name) >= 2
                    and mtl_name[0] == '"'
                    and mtl_name[-1] == '"'
                ):
                    mtl_name = mtl_name[1:-1]
                mtllibs.append(path.parent / mtl_name)

            # ------------------------------------------------
            # Material
            # ------------------------------------------------

            elif line.startswith("usemtl "):
                current_material = line[len("usemtl "):].strip()

            # ------------------------------------------------
            # Smoothing
            # ------------------------------------------------

            elif line.startswith("s "):
                value = line[2:].strip()

                if value.lower() == "off":
                    current_smoothing = None
                else:
                    current_smoothing = value

            # ------------------------------------------------
            # Face
            # ------------------------------------------------

            elif line.startswith("f "):

                current_object = ensure_object()

                corners = []

                for vertex_spec in line.split()[1:]:

                    components = vertex_spec.split("/")

                    # v
                    v = resolve_index(
                        int(components[0]),
                        len(vertices),
                    )

                    # vt
                    vt = None

                    if len(components) >= 2 and components[1]:
                        vt = resolve_index(
                            int(components[1]),
                            len(texcoords),
                        )

                    # vn
                    vn = None

                    if len(components) >= 3 and components[2]:
                        vn = resolve_index(
                            int(components[2]),
                            len(normals),
                        )

                    corners.append({
                        "v": v,
                        "vt": vt,
                        "vn": vn,
                    })

                if len(corners) < 3:
                    raise RuntimeError(
                        f"Face with fewer than 3 vertices at "
                        f"line {line_number}"
                    )
                
                # Y/Z swap is a reflection and therefore reverses
                # the orientation of the triangles.
                corners.reverse()

                current_object["faces"].append({
                    "corners": corners,
                    "material": current_material,
                    "smoothing": current_smoothing,
                })

    return (
        vertices,
        texcoords,
        normals,
        objects,
        mtllibs,
    )


# ============================================================
# Blender material creation
# ============================================================

def create_blender_materials(mtl_data):
    """
    Create Blender materials from parsed MTL data.

    Returns:
        {
            material_name: bpy.types.Material
        }
    """

    result = {}

    for name, data in mtl_data.items():

        material = bpy.data.materials.get(name)

        if material is None:
            material = bpy.data.materials.new(name)

        material.use_nodes = True

        nodes = material.node_tree.nodes
        links = material.node_tree.links

        # Clear the default node setup.
        nodes.clear()

        output = nodes.new("ShaderNodeOutputMaterial")
        shader = nodes.new("ShaderNodeBsdfPrincipled")

        links.new(
            shader.outputs["BSDF"],
            output.inputs["Surface"],
        )

        # ----------------------------------------------------
        # Diffuse color
        # ----------------------------------------------------

        kd = data["Kd"]

        shader.inputs["Base Color"].default_value = (
            kd[0],
            kd[1],
            kd[2],
            1.0,
        )

        # ----------------------------------------------------
        # Roughness / specular
        # ----------------------------------------------------

        ns = data["Ns"]

        if "Roughness" in shader.inputs:
            if ns > 0:
                # Approximate OBJ Phong shininess with
                # Principled roughness.
                roughness = max(
                    0.0,
                    min(
                        1.0,
                        (2.0 / (ns + 2.0)) ** 0.5,
                    )
                )
            else:
                roughness = 1.0

            shader.inputs["Roughness"].default_value = roughness

        # ----------------------------------------------------
        # Transparency
        # ----------------------------------------------------

        alpha = data["d"]

        if "Alpha" in shader.inputs:
            shader.inputs["Alpha"].default_value = alpha

        # ----------------------------------------------------
        # Diffuse texture
        # ----------------------------------------------------

        texture_path = data["map_Kd"]

        if texture_path is not None:

            if texture_path.exists():

                try:
                    image = bpy.data.images.load(
                        str(texture_path),
                        check_existing=True,
                    )

                    texcoord = nodes.new(
                        "ShaderNodeTexCoord"
                    )

                    image_texture = nodes.new(
                        "ShaderNodeTexImage"
                    )

                    image_texture.image = image

                    links.new(
                        texcoord.outputs["UV"],
                        image_texture.inputs["Vector"],
                    )

                    links.new(
                        image_texture.outputs["Color"],
                        shader.inputs["Base Color"],
                    )

                except Exception as e:
                    print(
                        f"WARNING: Could not load texture "
                        f"{texture_path}: {e}"
                    )

            else:
                print(
                    f"WARNING: Texture not found: "
                    f"{texture_path}"
                )

        result[name] = material

    return result


# ============================================================
# Collection handling
# ============================================================

def create_collection(name):

    old_collection = bpy.data.collections.get(name)

    if old_collection is not None:

        for obj in list(old_collection.objects):
            bpy.data.objects.remove(
                obj,
                do_unlink=True,
            )

        bpy.data.collections.remove(old_collection)

    collection = bpy.data.collections.new(name)

    bpy.context.scene.collection.children.link(
        collection
    )

    return collection


# ============================================================
# Blender mesh creation
# ============================================================

def create_blender_objects(
    vertices,
    texcoords,
    normals,
    obj_data,
    collection,
    blender_materials,
):
    """
    Create Blender meshes.

    Returns:

        created_objects

        vertex_mapping:
            OBJ global vertex index ->
            list of (Blender object, local vertex index)

        face_mapping:
            OBJ global face index ->
            (Blender object, local polygon index)
    """

    created_objects = []

    vertex_mapping = [
        []
        for _ in vertices
    ]

    face_mapping = []

    global_face_index = 0

    for object_index, data in enumerate(obj_data):

        name = data["name"]
        faces = data["faces"]

        # ----------------------------------------------------
        # Determine local vertices.
        #
        # We preserve the order in which OBJ vertices are
        # first encountered by this object.
        # ----------------------------------------------------

        global_to_local = {}
        local_to_global = []

        for face in faces:
            for corner in face["corners"]:

                global_index = corner["v"]

                if global_index not in global_to_local:

                    local_index = len(local_to_global)

                    global_to_local[global_index] = local_index
                    local_to_global.append(global_index)

        local_vertices = [
            vertices[i]
            for i in local_to_global
        ]

        # ----------------------------------------------------
        # Create local face indices.
        # ----------------------------------------------------

        local_faces = []

        for face in faces:

            local_face = [
                global_to_local[corner["v"]]
                for corner in face["corners"]
            ]

            local_faces.append(local_face)

        # ----------------------------------------------------
        # Create mesh
        # ----------------------------------------------------

        mesh = bpy.data.meshes.new(
            name + "_Mesh"
        )

        mesh.from_pydata(
            local_vertices,
            [],
            local_faces,
        )

        mesh.update()

        obj = bpy.data.objects.new(
            name,
            mesh,
        )

        collection.objects.link(obj)

        created_objects.append(obj)

        # ----------------------------------------------------
        # Materials
        # ----------------------------------------------------

        material_slots = {}

        for face in faces:

            material_name = face["material"]

            if material_name is None:
                continue

            if material_name not in blender_materials:
                print(
                    f"WARNING: Material '{material_name}' "
                    f"was referenced by OBJ but not found "
                    f"in the MTL."
                )

                # Create a fallback material.
                material = bpy.data.materials.get(
                    material_name
                )

                if material is None:
                    material = bpy.data.materials.new(
                        material_name
                    )

                    material.use_nodes = True

                blender_materials[material_name] = material

            material = blender_materials[material_name]

            if material_name not in material_slots:

                slot_index = len(obj.data.materials)

                obj.data.materials.append(material)

                material_slots[material_name] = slot_index

        # Assign material indices.
        for local_face_index, face in enumerate(faces):

            material_name = face["material"]

            if material_name is not None:
                obj.data.polygons[
                    local_face_index
                ].material_index = material_slots[
                    material_name
                ]

        # ----------------------------------------------------
        # UV layer
        # ----------------------------------------------------

        has_uvs = any(
            corner["vt"] is not None
            for face in faces
            for corner in face["corners"]
        )

        if has_uvs:

            uv_layer = mesh.uv_layers.new(
                name="UVMap"
            )

            for polygon in mesh.polygons:

                face = faces[polygon.index]

                for loop_index, corner in zip(
                    polygon.loop_indices,
                    face["corners"],
                ):

                    if corner["vt"] is None:
                        continue

                    uv = texcoords[corner["vt"]]

                    uv_layer.data[
                        loop_index
                    ].uv = uv

        # ----------------------------------------------------
        # Custom normals
        # ----------------------------------------------------

        has_normals = any(
            corner["vn"] is not None
            for face in faces
            for corner in face["corners"]
        )

        if has_normals:

            loop_normals = []

            for polygon in mesh.polygons:

                face = faces[polygon.index]

                for corner in face["corners"]:

                    if corner["vn"] is None:
                        # Fallback to the polygon normal.
                        loop_normals.append(
                            polygon.normal[:]
                        )
                    else:
                        loop_normals.append(
                            normals[corner["vn"]]
                        )

            try:
                mesh.normals_split_custom_set(
                    loop_normals
                )
            except RuntimeError as e:
                print(
                    f"WARNING: Could not set custom normals "
                    f"for '{name}': {e}"
                )

        # ----------------------------------------------------
        # Smoothing
        # ----------------------------------------------------

        for polygon, face in zip(
            mesh.polygons,
            faces,
        ):
            polygon.use_smooth = (
                face["smoothing"] is not None
            )

        # ----------------------------------------------------
        # Vertex mapping
        # ----------------------------------------------------

        for local_index, global_index in enumerate(
            local_to_global
        ):
            vertex_mapping[
                global_index
            ].append(
                (obj, local_index)
            )

        # ----------------------------------------------------
        # Face mapping
        # ----------------------------------------------------

        for local_face_index in range(
            len(local_faces)
        ):
            face_mapping.append(
                (
                    obj,
                    local_face_index,
                )
            )

        global_face_index += len(local_faces)

    return (
        created_objects,
        vertex_mapping,
        face_mapping,
    )


# ============================================================
# JSON attribute import
# ============================================================

def import_attributes(
    objects,
    vertex_mapping,
    face_mapping,
    attr_path,
):

    with open(
        attr_path,
        "r",
        encoding="utf-8",
    ) as f:
        data = json.load(f)

    num_obj_vertices = len(vertex_mapping)
    num_obj_faces = len(face_mapping)

    print()
    print("JSON attribute statistics")
    print("=========================")

    for name in (
        FLOAT_ATTRIBUTES
        + VEC2_ATTRIBUTES
        + VEC3_ATTRIBUTES
        + COLOR_ATTRIBUTES
        + FACE_ATTRIBUTES
    ):
        if name in data:
            print(
                f"{name:30s}: {len(data[name])}"
            )

    print()
    print(f"OBJ vertices: {num_obj_vertices}")
    print(f"OBJ faces:    {num_obj_faces}")
    print()

    # --------------------------------------------------------
    # Generic vertex attribute writer
    # --------------------------------------------------------

    def write_vertex_attribute(
        name,
        data_type,
        property_name,
    ):

        if name not in data:
            return

        values = data[name]

        if len(values) != num_obj_vertices:
            raise RuntimeError(
                f"Vertex attribute '{name}' has "
                f"{len(values)} values, but OBJ contains "
                f"{num_obj_vertices} vertices."
            )

        attributes = {}

        for obj in objects:
            attributes[obj] = ensure_attribute(
                obj.data,
                name,
                data_type,
                "POINT",
            )

        for global_index, destinations in enumerate(
            vertex_mapping
        ):

            value = values[global_index]

            for obj, local_index in destinations:

                setattr(
                    attributes[obj].data[
                        local_index
                    ],
                    property_name,
                    value,
                )

    # --------------------------------------------------------
    # FLOAT
    # --------------------------------------------------------

    for name in FLOAT_ATTRIBUTES:
        write_vertex_attribute(
            name,
            "FLOAT",
            "value",
        )

    # --------------------------------------------------------
    # VEC2
    # --------------------------------------------------------

    for name in VEC2_ATTRIBUTES:
        write_vertex_attribute(
            name,
            "FLOAT2",
            "vector",
        )

    # --------------------------------------------------------
    # VEC3
    # --------------------------------------------------------

    for name in VEC3_ATTRIBUTES:
        write_vertex_attribute(
            name,
            "FLOAT_VECTOR",
            "vector",
        )

    # --------------------------------------------------------
    # COLOR
    # --------------------------------------------------------

    for name in COLOR_ATTRIBUTES:

        if name not in data:
            continue

        values = data[name]

        if len(values) != num_obj_vertices:
            raise RuntimeError(
                f"Vertex attribute '{name}' has "
                f"{len(values)} values, but OBJ contains "
                f"{num_obj_vertices} vertices."
            )

        attributes = {}

        for obj in objects:
            attributes[obj] = ensure_attribute(
                obj.data,
                name,
                "FLOAT_COLOR",
                "POINT",
            )

        for global_index, destinations in enumerate(
            vertex_mapping
        ):

            value = values[global_index]

            if len(value) == 3:
                value = [
                    value[0],
                    value[1],
                    value[2],
                    1.0,
                ]

            for obj, local_index in destinations:

                attributes[obj].data[
                    local_index
                ].color = value

    # --------------------------------------------------------
    # FACE
    # --------------------------------------------------------

    for name in FACE_ATTRIBUTES:

        if name not in data:
            continue

        values = data[name]

        if len(values) != num_obj_faces:
            raise RuntimeError(
                f"Face attribute '{name}' has "
                f"{len(values)} values, but OBJ contains "
                f"{num_obj_faces} faces."
            )

        attributes = {}

        for obj in objects:
            attributes[obj] = ensure_attribute(
                obj.data,
                name,
                "BOOLEAN",
                "FACE",
            )

        for global_index, (
            obj,
            local_index,
        ) in enumerate(face_mapping):

            attributes[obj].data[
                local_index
            ].value = values[global_index]

    print("Attributes imported successfully.")


# ============================================================
# Attribute helper
# ============================================================

def ensure_attribute(
    mesh,
    name,
    data_type,
    domain,
):

    attribute = mesh.attributes.get(name)

    if attribute is None:

        attribute = mesh.attributes.new(
            name=name,
            type=data_type,
            domain=domain,
        )

    else:

        if attribute.data_type != data_type:
            raise RuntimeError(
                f"Attribute '{name}' already exists with "
                f"type {attribute.data_type}, expected "
                f"{data_type}"
            )

        if attribute.domain != domain:
            raise RuntimeError(
                f"Attribute '{name}' already exists on "
                f"domain {attribute.domain}, expected "
                f"{domain}"
            )

    return attribute


# ============================================================
# Main
# ============================================================

def main():

    obj_path = Path(OBJ_PATH)
    json_path = Path(JSON_PATH)

    if not obj_path.exists():
        raise FileNotFoundError(
            f"OBJ file does not exist:\n{obj_path}"
        )

    if not json_path.exists():
        raise FileNotFoundError(
            f"JSON file does not exist:\n{json_path}"
        )

    print()
    print("========================================")
    print("Importing OBJ and JSON attributes")
    print("========================================")

    # --------------------------------------------------------
    # Parse OBJ
    # --------------------------------------------------------

    print(f"OBJ: {obj_path}")

    (
        vertices,
        texcoords,
        normals,
        obj_data,
        mtllibs,
    ) = parse_obj(obj_path)

    num_faces = sum(
        len(obj["faces"])
        for obj in obj_data
    )

    print(f"Vertices:    {len(vertices)}")
    print(f"Texcoords:   {len(texcoords)}")
    print(f"Normals:     {len(normals)}")
    print(f"Objects:     {len(obj_data)}")
    print(f"Faces:       {num_faces}")

    print()

    for i, obj in enumerate(obj_data):

        print(
            f"  [{i:3d}] "
            f"{obj['name']:30s} "
            f"{len(obj['faces'])} faces"
        )

    print()

    # --------------------------------------------------------
    # Parse MTL
    # --------------------------------------------------------

    mtl_data = {}

    for mtl_path in mtllibs:

        print(f"MTL: {mtl_path}")

        parsed = parse_mtl(mtl_path)

        mtl_data.update(parsed)

    print(
        f"Materials found: {len(mtl_data)}"
    )

    blender_materials = create_blender_materials(
        mtl_data
    )

    # --------------------------------------------------------
    # Create collection
    # --------------------------------------------------------

    collection = create_collection(
        COLLECTION_NAME
    )

    # --------------------------------------------------------
    # Create meshes
    # --------------------------------------------------------

    (
        objects,
        vertex_mapping,
        face_mapping,
    ) = create_blender_objects(
        vertices,
        texcoords,
        normals,
        obj_data,
        collection,
        blender_materials,
    )

    # --------------------------------------------------------
    # Import JSON attributes
    # --------------------------------------------------------

    import_attributes(
        objects,
        vertex_mapping,
        face_mapping,
        json_path,
    )

    # --------------------------------------------------------
    # Select imported objects
    # --------------------------------------------------------

    bpy.ops.object.select_all(
        action="DESELECT"
    )

    for obj in objects:
        obj.select_set(True)

    if objects:
        bpy.context.view_layer.objects.active = (
            objects[0]
        )

    print()
    print("========================================")
    print("Import complete")
    print("========================================")
    print(f"Collection: {COLLECTION_NAME}")
    print(f"Objects:    {len(objects)}")
    print()


main()