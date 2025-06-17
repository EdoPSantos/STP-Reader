# === Funções Globais (usadas por múltiplos ficheiros) ===

import os
import re
from collections import defaultdict

from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Plane
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.gp import gp_Pnt
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_FACE

# Utilitário para gerar chave de agrupamento de furos por eixo

def axis_key(adaptor):
    surface_type = adaptor.GetType()
    if surface_type == GeomAbs_Cylinder:
        axis = adaptor.Cylinder().Axis()
    elif surface_type == GeomAbs_Cone:
        axis = adaptor.Cone().Axis()
    else:
        raise ValueError("Tipo de superfície não suportado em axis_key")
    loc = axis.Location()
    dir = axis.Direction()
    return (
        round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2),
        round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)
    )

# Função para identificar se o furo atravessa a peça

def is_hole_through(faces, shape, tol=2.0):
    if not shape or shape.IsNull():
        raise ValueError("Shape inválido recebido em is_hole_through.")
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bounds = bbox.Get()
    hit_faces = set()
    for surface_type, face, adaptor in faces:
        vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while vertex_explorer.More():
            vertex = vertex_explorer.Current()
            pnt = BRep_Tool.Pnt(vertex)
            touched = get_bbox_faces_touched(pnt, bounds, tol)
            hit_faces.update(touched)
            vertex_explorer.Next()
    if faces:
        surface_type, face, adaptor = faces[0]
        if surface_type == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
        elif surface_type == GeomAbs_Cone:
            axis = adaptor.Cone().Axis()
        else:
            axis = None
        if axis:
            loc = axis.Location()
            dir = axis.Direction()
            for k in range(-3, 4):
                t = k * 0.25
                px = loc.X() + dir.X() * t * (bounds[3] - bounds[0])
                py = loc.Y() + dir.Y() * t * (bounds[4] - bounds[1])
                pz = loc.Z() + dir.Z() * t * (bounds[5] - bounds[2])
                pnt = gp_Pnt(px, py, pz)
                touched = get_bbox_faces_touched(pnt, bounds, tol)
                hit_faces.update(touched)
    axis_pairs = [('xmin','xmax'), ('ymin','ymax'), ('zmin','zmax')]
    for a, b in axis_pairs:
        if a in hit_faces and b in hit_faces:
            return True
    return False

# === utils para main_details.py ===

def extract_from_filename(filepath):
    base = os.path.basename(filepath)
    name, _ = os.path.splitext(base)
    parts = re.split(r"[_\-]", name)
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None, None

def extract_mold_and_part_from_step(filepath):
    file_mold, file_part = extract_from_filename(filepath)
    if file_mold and file_part:
        return file_mold, file_part

    mold, part = None, None
    product_mold, product_part = None, None
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "FILE_NAME" in line and (mold is None or part is None):
                match = re.search(r"FILE_NAME\('([^']+)'", line)
                if match:
                    filename = match.group(1)
                    base = os.path.basename(filename)
                    name, _ = os.path.splitext(base)
                    parts = re.split(r"[_\-]", name)
                    if len(parts) >= 2:
                        mold = parts[-2]
                        part = parts[-1]
            if "PRODUCT('" in line and (product_mold is None or product_part is None):
                match = re.search(r"PRODUCT\('([^']+)'", line)
                if match:
                    pname = match.group(1)
                    parts = re.split(r"[_\-]", pname)
                    if len(parts) >= 2:
                        product_mold = parts[-2]
                        product_part = parts[-1]

    if mold and part:
        return mold, part
    elif product_mold and product_part:
        return product_mold, product_part
    return None, None

def detect_positions_from_holes(axis_groups, shape):
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bounds = bbox.Get()
    zmin, zmax = bounds[2], bounds[5]
    thickness = zmax - zmin

    has_cone_up = False
    has_cone_down = False
    deep_blind_hole = False

    for faces in axis_groups:
        if is_hole_through(faces, shape):
            continue

        for stype, face, adaptor in faces:
            if stype == GeomAbs_Cone:
                dir_z = adaptor.Cone().Axis().Direction().Z()
                if dir_z > 0.5:
                    has_cone_up = True
                elif dir_z < -0.5:
                    has_cone_down = True
            elif stype == GeomAbs_Cylinder:
                face_bbox = Bnd_Box()
                brepbndlib.Add(face, face_bbox)
                _, _, fzmin, _, _, fzmax = face_bbox.Get()
                depth = abs(fzmax - fzmin)
                if depth > 0.9 * thickness:
                    deep_blind_hole = True

    if (has_cone_up and has_cone_down) or (deep_blind_hole and (has_cone_up or has_cone_down)):
        return 2
    if has_cone_up or has_cone_down or deep_blind_hole:
        return 2
    return 1

# === utils para cylindrical_details.py ===

def group_faces_by_axis(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    axis_groups = defaultdict(list)
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        surface_type = adaptor.GetType()
        if surface_type == GeomAbs_Cylinder:
            key = axis_key(adaptor)
            axis_groups[key].append((surface_type, face, adaptor))
        explorer.Next()
    return list(axis_groups.values())

def get_min_max_cylindrical_diameter(shape):
    min_d = None
    max_d = None
    axis_groups = group_faces_by_axis(shape)
    for faces in axis_groups:
        only_cylinders = all(stype == GeomAbs_Cylinder for stype, _, _ in faces)
        if only_cylinders:
            for stype, face, adaptor in faces:
                d = adaptor.Cylinder().Radius() * 2
                if (min_d is None) or (d < min_d):
                    min_d = d
                if (max_d is None) or (d > max_d):
                    max_d = d
    if min_d is not None and max_d is not None:
        return round(min_d, 2), round(max_d, 2)
    else:
        return None, None

def get_group_diameters(faces):
    diams = []
    for surface_type, face, adaptor in faces:
        if surface_type == GeomAbs_Cylinder:
            diams.append(round(adaptor.Cylinder().Radius() * 2, 2))
        elif surface_type == GeomAbs_Cone:
            base = round(adaptor.Cone().RefRadius() * 2, 2)
            diams.append(base)
    if not diams:
        return None, None
    return min(diams), max(diams)

def classify_hole_group_type(faces):
    types = set()
    for stype, _, _ in faces:
        if stype == GeomAbs_Cylinder:
            types.add("cilíndrico")
        elif stype == GeomAbs_Cone:
            types.add("cónico")
    if types == {"cilíndrico"}:
        return "cilíndrico"
    elif types == {"cónico"}:
        return "cónico"
    elif types == {"cilíndrico", "cónico"}:
        return "misto"
    else:
        return "desconhecido"

# === utils para rectangular_details.py ===

def group_rectangular_faces(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    rectangles = []
    while explorer.More():
        face = explorer.Current()
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        if adaptor.GetType() == GeomAbs_Plane:
            rectangles.append(face)
        explorer.Next()
    return rectangles

# === utils auxiliares (usados em várias partes) ===

def get_all_vertices_of_group(faces):
    vertices = []
    for face in faces:
        explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while explorer.More():
            v = explorer.Current()
            vertices.append(BRep_Tool.Pnt(v))
            explorer.Next()
    return vertices

def get_bbox_faces_touched(pnt, bbox, tol=0.5):
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    faces = []
    if abs(pnt.X() - xmin) < tol: faces.append("xmin")
    if abs(pnt.X() - xmax) < tol: faces.append("xmax")
    if abs(pnt.Y() - ymin) < tol: faces.append("ymin")
    if abs(pnt.Y() - ymax) < tol: faces.append("ymax")
    if abs(pnt.Z() - zmin) < tol: faces.append("zmin")
    if abs(pnt.Z() - zmax) < tol: faces.append("zmax")
    return faces

def project_point_to_bbox(loc, dir, bounds):
    px, py, pz = loc.X(), loc.Y(), loc.Z()
    dx, dy, dz = dir.X(), dir.Y(), dir.Z()
    if abs(dx) < 1e-8: dx = 1e-8
    if abs(dy) < 1e-8: dy = 1e-8
    if abs(dz) < 1e-8: dz = 1e-8

    t_values = []
    for i, (min_b, max_b, p, d) in enumerate(zip(bounds[:3], bounds[3:], [px, py, pz], [dx, dy, dz])):
        t1 = (min_b - p) / d
        t2 = (max_b - p) / d
        t_values.extend([t1, t2])
    t_values = [t for t in t_values if t > 0]
    if not t_values:
        return loc
    t = min(t_values)
    return gp_Pnt(px + t*dx, py + t*dy, pz + t*dz)
