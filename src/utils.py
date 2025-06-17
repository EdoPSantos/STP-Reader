# === Funções Globais (usadas por múltiplos ficheiros) ===

# === Módulos padrão ===
import os
import re
import math
from collections import defaultdict

# === OCC / pythonocc-core ===
from OCC.Core.gp import gp_Pnt
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBndLib import brepbndlib

from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Plane
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface

from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer

from OCC.Core.TopoDS import TopoDS_Shape, TopoDS_Face, topods
from OCC.Core.TopTools import TopTools_IndexedMapOfShape

# Utilitário para gerar chave de agrupamento de furos por eixo

def axis_key(adaptor, loc_tol=0.2, dir_tol=0.02):
    """
    Gera uma chave de agrupamento baseada no centro (location) e direção (direction) do eixo,
    ambos arredondados para permitir pequenas variações (tolerância).
    - loc_tol: tolerância para a localização do centro (em mm)
    - dir_tol: tolerância para a direção do eixo (unitário)
    """
    surface_type = adaptor.GetType()
    if surface_type == GeomAbs_Cylinder:
        axis = adaptor.Cylinder().Axis()
    elif surface_type == GeomAbs_Cone:
        axis = adaptor.Cone().Axis()
    else:
        raise ValueError("Tipo de superfície não suportado em axis_key")
    loc = axis.Location()
    dir = axis.Direction()
    # Arredonda para tolerância especificada
    return (
        round(loc.X() / loc_tol) * loc_tol,
        round(loc.Y() / loc_tol) * loc_tol,
        round(loc.Z() / loc_tol) * loc_tol,
        round(dir.X() / dir_tol) * dir_tol,
        round(dir.Y() / dir_tol) * dir_tol,
        round(dir.Z() / dir_tol) * dir_tol,
    )

# Função para identificar se o furo atravessa a peça
def is_hole_through(faces, shape, tol=2.0):
    """
    Verifica se pelo menos um vértice de qualquer face do grupo toca uma face extrema do bounding box,
    E pelo menos um vértice toca a outra face extrema (ou seja, o furo é passante).
    """
    if not shape or shape.IsNull():
        raise ValueError("Shape inválido recebido em is_hole_through.")
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bounds = bbox.Get()
    faces_touched = set()

    for _, face, _ in faces:
        vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while vertex_explorer.More():
            vertex = vertex_explorer.Current()
            pnt = BRep_Tool.Pnt(vertex)
            touched = get_bbox_faces_touched(pnt, bounds, tol)
            faces_touched.update(touched)
            vertex_explorer.Next()
    axis_pairs = [('xmin','xmax'), ('ymin','ymax'), ('zmin','zmax')]
    for a, b in axis_pairs:
        if a in faces_touched and b in faces_touched:
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

def group_connected_planar_faces(shape):
    """
    Retorna uma lista de todas as faces planas como grupos individuais.
    """
    groups = []

    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = topods.Face(explorer.Current())
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        if adaptor.GetType() == GeomAbs_Plane:
            groups.append(face)
        explorer.Next()
    return groups

def classify_rectangular_face(face, tol=1.0):
    if not isinstance(face, TopoDS_Face):
        try:
            face = topods.Face(face)
        except Exception as e:
            raise ValueError("Não foi possível converter a face para TopoDS_Face") from e

    bbox = Bnd_Box()
    brepbndlib.Add(face, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()

    dx = abs(xmax - xmin)
    dy = abs(ymax - ymin)
    dz = abs(zmax - zmin)

    sides = sorted([dx, dy, dz])
    largura, comprimento = round(sides[0], 2), round(sides[1], 2)

    tipo = "quadrado" if abs(largura - comprimento) <= tol else "retangular"
    return tipo, max(largura, comprimento), min(largura, comprimento)

def get_piece_dimensions(shape):
    """
    Calcula o comprimento e largura máximos da peça (bounding box global).
    Retorna (comprimento, largura).
    """
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    dx = abs(xmax - xmin)
    dy = abs(ymax - ymin)
    # No caso de peças planas, o maior valor é comprimento, o menor é largura
    comprimento = round(max(dx, dy), 2)
    largura = round(min(dx, dy), 2)
    return comprimento, largura

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

def group_faces_by_axis_and_proximity(shape, loc_tol=2.0, dir_tol=0.05):
    """
    Agrupa faces cilíndricas e cónicas por eixo (direção) e centro (location) com tolerância.
    Furos mistos (cilindro+cone) que estejam colineares e com centro próximo ficam no mesmo grupo.
    """
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    temp_groups = []

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        stype = adaptor.GetType()
        if stype not in (GeomAbs_Cylinder, GeomAbs_Cone):
            explorer.Next()
            continue

        # Obter centro e direção
        if stype == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
        elif stype == GeomAbs_Cone:
            axis = adaptor.Cone().Axis()
        loc = axis.Location()
        dir = axis.Direction()
        temp_groups.append({
            'loc': (loc.X(), loc.Y(), loc.Z()),
            'dir': (dir.X(), dir.Y(), dir.Z()),
            'stype': stype,
            'face': face,
            'adaptor': adaptor
        })
        explorer.Next()

    # Função de comparação de direção com tolerância (produto escalar)
    def is_dir_close(d1, d2, tol=dir_tol):
        dot = sum(a * b for a, b in zip(d1, d2))
        return abs(abs(dot) - 1.0) < tol

    # Agrupamento
    groups = []
    used = [False] * len(temp_groups)

    for i, g in enumerate(temp_groups):
        if used[i]:
            continue
        group = [g]
        used[i] = True
        for j, h in enumerate(temp_groups):
            if used[j]:
                continue
            # Mesma direção e centro suficientemente próximos
            if is_dir_close(g['dir'], h['dir'], dir_tol):
                dist = math.sqrt(sum((a - b) ** 2 for a, b in zip(g['loc'], h['loc'])))
                if dist < loc_tol:
                    group.append(h)
                    used[j] = True
        groups.append([(item['stype'], item['face'], item['adaptor']) for item in group])
    return groups

def get_auto_loc_tol(shape):
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    max_dim = max(abs(xmax - xmin), abs(ymax - ymin), abs(zmax - zmin))
    # Ajusta conforme os teus requisitos e experiência
    if max_dim < 200:
        return 2.0
    elif max_dim < 500:
        return 5.0
    else:
        return 10.0
