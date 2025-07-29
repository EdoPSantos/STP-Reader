from OCC.Core.TopTools import TopTools_IndexedMapOfShape, TopTools_IndexedDataMapOfShapeListOfShape, TopTools_ListIteratorOfListOfShape
from OCC.Core.TopoDS import TopoDS_Shape, topods, TopoDS_Face, TopoDS_Wire, TopoDS_Vertex
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Plane
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_FACE, TopAbs_EDGE
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.GeomAdaptor import GeomAdaptor_Curve
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Dir
from OCC.Core.GeomAbs import GeomAbs_Circle
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.BRepTools import breptools
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Bnd import Bnd_Box
from collections import Counter
import OCC.Core.TopExp
from collections import defaultdict
import numpy as np
import math
import os
import re

# ================================================== Code Recycled ==================================================
# ------------------------------------------------------ Main -------------------------------------------------------

def feat_key(feat, center_tol=1.0, d_tol=0.5):
    """
    Gera uma chave para agrupamento de features circulares/semicirculares.
    """
    return (round(feat['center'][0]/center_tol), round(feat['center'][1]/center_tol), round(feat['max_d']/d_tol))

def is_point_inside_bbox(point, bbox, tol=0.5):
    """
    Verifica se um ponto (x, y) está dentro do bounding box.
    """
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    x, y = point[:2]
    return (xmin - tol) <= x <= (xmax + tol) and (ymin - tol) <= y <= (ymax + tol)

def get_bbox(shape_or_face):
    """
    Retorna a bounding box de uma shape ou face: (xmin, ymin, zmin, xmax, ymax, zmax)
    """

    bbox = Bnd_Box()
    brepbndlib.Add(shape_or_face, bbox)
    return bbox.Get()

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

# -------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------- Name-File Extraction -----------------------------------------------

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

def extract_from_filename(filepath):
    base = os.path.basename(filepath)
    name, _ = os.path.splitext(base)
    parts = re.split(r"[_\-]", name)
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None, None

# -------------------------------------------------------------------------------------------------------------------

# ================================================== Code Details ===================================================
# ---------------------------------------------------- Circular -----------------------------------------------------

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

def get_all_radii_of_group(faces, min_allowed=3.0):
    """
    Retorna todos os raios (em mm) de todas as faces do grupo.
    Para cones, calcula sempre o raio nos dois extremos da face (base e topo).
    """

    radii = []
    has_cone = False
    for stype, face, adaptor in faces:
        if stype == GeomAbs_Cylinder:
            r = adaptor.Cylinder().Radius()
            if r >= min_allowed:
                radii.append(r)
        elif stype == GeomAbs_Cone:
            has_cone = True
            cone = adaptor.Cone()
            ref_radius = cone.RefRadius()
            semi_angle = cone.SemiAngle()
            axis = cone.Axis()
            axis_loc = axis.Location()
            axis_dir = axis.Direction()
            zs = []
            for pnt in get_all_vertices_of_group([face]):
                v = pnt.XYZ() - axis_loc.XYZ()
                z_rel = v.Dot(axis_dir.XYZ())
                zs.append(z_rel)
            if zs:
                z_min = min(zs)
                z_max = max(zs)
                r_min = abs(ref_radius + z_min * math.tan(semi_angle))
                r_max = abs(ref_radius + z_max * math.tan(semi_angle))
                if r_min >= min_allowed:
                    radii.append(r_min)
                if r_max >= min_allowed:
                    radii.append(r_max)
            else:
                for dz in [-10, 10]:
                    radius = abs(ref_radius + dz * math.tan(semi_angle))
                    if radius >= min_allowed:
                        radii.append(radius)
    radii = sorted(set(round(r, 2) for r in radii))
    return radii, has_cone

def get_all_vertices_of_group(faces):
    vertices = []
    for face in faces:
        explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while explorer.More():
            v = explorer.Current()
            vertices.append(BRep_Tool.Pnt(v))
            explorer.Next()
    return vertices

def merge_aligned_circular_features(features, axis_tol=30.0, center_tol=10.0, diameter_tol=2.0):
    """
    Agrupa features circulares com centros alinhados numa direção (por ex. Z), 
    com diâmetros semelhantes e próximos no espaço.
    """
    merged = []
    used = set()

    for i, f1 in enumerate(features):
        if i in used:
            continue
        group = [f1]
        for j, f2 in enumerate(features):
            if j <= i or j in used:
                continue
            # Verifica se os diâmetros são semelhantes
            if abs(f1['max_d'] - f2['max_d']) > diameter_tol:
                continue
            # Verifica alinhamento em X/Y (ou outra direção dominante)
            dx = abs(f1['center'][0] - f2['center'][0])
            dy = abs(f1['center'][1] - f2['center'][1])
            dz = abs(f1['center'][2] - f2['center'][2])
            # Aqui vamos supor alinhamento em Z
            if dx < axis_tol and dy < axis_tol and dz < center_tol:
                group.append(f2)
                used.add(j)
        used.add(i)

        # Agrupar os dados
        center_avg = tuple(round(sum(f['center'][k] for f in group)/len(group), 2) for k in range(3))
        all_min_d = [f['min_d'] for f in group]
        all_max_d = [f['max_d'] for f in group]
        merged.append({
            'center': center_avg,
            'min_d': round(min(all_min_d), 2),
            'max_d': round(max(all_max_d), 2)
        })
    return merged

def detect_positions_from_holes(axis_groups, shape):
    
    """
    Só considera grupos com pelo menos um cilindro (furos reais) para a lógica de posições.
    Ignora grupos só de cones (chanfros/cabeças de embutir).
    Retorna 2 posições se houver furos com fundo (não passantes) em direções opostas.
    Furos passantes são ignorados.
    """
    has_up = False
    has_down = False
    for i, faces in enumerate(axis_groups):
        if is_hole_through(faces, shape):
            # print(f"Grupo {i}: passante")
            continue
        # Só considera grupos com pelo menos um cilindro
        has_cylinder = any(stype == GeomAbs_Cylinder for stype, _, _ in faces)
        if not has_cylinder:
            # print(f"Grupo {i}: ignorado (só cone)")
            continue
        # Calcula direção e profundidade do cilindro principal
        for stype, face, adaptor in faces:
            if stype == GeomAbs_Cylinder:
                dir_z = adaptor.Cylinder().Axis().Direction().Z()
                # Calcula profundidade
                bbox = Bnd_Box()
                brepbndlib.Add(face, bbox)
                _, _, zmin, _, _, zmax = bbox.Get()
                depth = abs(zmax - zmin)
                # Considera só se profundidade for relevante (>2mm)
                if depth < 2.0:
                    continue
                if dir_z > 0.5:
                    has_up = True
                elif dir_z < -0.5:
                    has_down = True
                break
    if has_up and has_down:
        return 2
    return 1

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

def group_holes_by_center(features, center_tol=10.0):
    """
    Agrupa furos (circulares e semicirculares) com base na proximidade do centro (x, y).
    Para cada grupo, guarda todos os min_d e max_d.
    No resumo, mostra o menor min_d e o maior max_d de cada grupo.
    """
    grouped_by_center = []

    def dist2d(a, b):
        return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5

    for feat in features:
        cx, cy = feat['center'][:2]
        min_d = feat['min_d']
        max_d = feat['max_d']
        added = False
        for group in grouped_by_center:
            gx, gy = group['center']
            if dist2d((cx, cy), (gx, gy)) < center_tol:
                group['min_ds'].append(min_d)
                group['max_ds'].append(max_d)
                added = True
                break
        if not added:
            grouped_by_center.append({
                'center': (cx, cy),
                'min_ds': [min_d],
                'max_ds': [max_d]
            })

    all_counter = Counter()
    for group in grouped_by_center:
        min_d = min(group['min_ds'])
        max_d = max(group['max_ds'])
        key = (round(min_d, 1), round(max_d, 1))
        all_counter[key] += 1
    return all_counter

def get_face_area(face):
    """Retorna a área de uma face."""
    props = GProp_GProps()
    brepgprop.SurfaceProperties(face, props)
    return props.Mass()

# ----------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------- Semi-Cylindrical ---------------------------------------------------

def collect_semi_circular_arcs(shape, lateral_tol=3.0, angle_tol=30.0, radius_tol=1.0):
    """
    Recolhe todos os recortes semicirculares (arcos ~180º perto da lateral) e devolve lista de dicts:
    {'center': (x, y, z), 'radius': r, 'angle': angle}
    """

    results = []
    piece_bbox = get_bbox(shape)
    explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    while explorer.More():
        edge = explorer.Current()
        curve_data = BRep_Tool.Curve(edge)
        if not curve_data or len(curve_data) < 3:
            explorer.Next()
            continue
        curve_handle, first, last = curve_data
        adaptor = GeomAdaptor_Curve(curve_handle, first, last)
        if adaptor.GetType() == GeomAbs_Circle:
            circ = adaptor.Circle()
            center = circ.Location()
            r = circ.Radius()
            angle = abs(math.degrees(last - first))
            # Só arcos semicirculares (aprox 180º)
            if 180 - angle_tol <= angle <= 180 + angle_tol:
                x, y, z = center.X(), center.Y(), center.Z()
                xmin, ymin, zmin, xmax, ymax, zmax = piece_bbox
                near_lateral = (
                    abs(x - xmin) < lateral_tol or abs(x - xmax) < lateral_tol or
                    abs(y - ymin) < lateral_tol or abs(y - ymax) < lateral_tol
                )
                if near_lateral:
                    # Não duplica semicirculo já existente
                    found = False
                    for rec in results:
                        rx, ry, rz = rec['center']
                        rr = rec['radius']
                        if (abs(rx-x) < radius_tol and abs(ry-y) < radius_tol and abs(rz-z) < radius_tol and abs(rr-r) < radius_tol):
                            found = True
                            break
                    if not found:
                        results.append({'center': (x, y, z), 'radius': r, 'angle': angle})
        explorer.Next()
    return results

def group_semi_circular_arcs(arcs, group_tol=3.0):
    """
    Agrupa arcos semicirculares por centro X/Y (não por Z), tolerância em mm.
    """
    grouped = []
    for arc in arcs:
        x, y, z = arc['center']
        r = arc['radius']
        found = False
        for g in grouped:
            gx, gy = g['xy']
            if abs(x - gx) < group_tol and abs(y - gy) < group_tol:
                g['zs'].append(z)
                g['radii'].append(r)
                g['angles'].append(arc['angle'])
                found = True
                break
        if not found:
            grouped.append({
                'xy': (x, y),
                'zs': [z],
                'radii': [r],
                'angles': [arc['angle']]
            })
    return grouped

# ----------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------- Rectangular ------------------------------------------------------

def get_all_planar_faces_bbox(shape):
    """
    Retorna lista de dicts com centro, comprimento, largura, altura e bbox de todas as faces planas (inclusive circulares).
    """
    result = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        stype = adaptor.GetType()
        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        dx = abs(xmax - xmin)
        dy = abs(ymax - ymin)
        dz = abs(zmax - zmin)
        center = ((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
        tipo = None
        if stype == GeomAbs_Plane:
            tipo = 'plano'
        elif stype == GeomAbs_Cylinder:
            tipo = 'cilindro'
        elif stype == GeomAbs_Cone:
            tipo = 'cone'
        elif stype == GeomAbs_Sphere:
            tipo = 'esfera'
        elif stype == GeomAbs_Torus:
            tipo = 'toro'
        else:
            tipo = str(stype)
        result.append({
            'type': tipo,
            'center': tuple(round(c, 2) for c in center),
            'length': round(max(dx, dy), 2),
            'width': round(min(dx, dy), 2),
            'height': round(dz, 2),
            'bbox': (xmin, ymin, zmin, xmax, ymax, zmax)
        })
        explorer.Next()
    return result

# ----------------------------------------------------------------------------------------------------------------------