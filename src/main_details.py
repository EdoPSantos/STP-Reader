import os
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from collections import defaultdict
from .utils import axis_key, is_hole_through

def extract_from_filename(filepath):
    """
    Extrai o molde e a peça do nome do ficheiro físico.
    Exemplo: lixo_M251799_1-1.stp => molde=M251799, peça=1-1
    """
    base = os.path.basename(filepath)
    name, _ = os.path.splitext(base)
    import re
    parts = re.split(r"[_\-]", name)
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None, None

def extract_mold_and_part_from_step(filepath):
    """
    1. Tenta pelo nome do ficheiro físico (penúltima e última partes).
    2. Depois FILE_NAME no header.
    3. Depois PRODUCT.
    4. Se não encontrar, devolve (None, None)
    """
    import re

    # 1. Nome do ficheiro físico
    file_mold, file_part = extract_from_filename(filepath)
    if file_mold and file_part:
        return file_mold, file_part

    # 2. FILE_NAME no header
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

def group_faces_by_axis(shape):
    """Agrupa faces cilíndricas pelo eixo."""
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

    return axis_groups

def is_face_exterior(face, shape):
    """
    Verifica se uma face está na superfície exterior da peça (shell).
    Abordagem básica: compara bounding box dos vértices da face com bounding box global.
    """
    from OCC.Core.TopExp import TopExp_Explorer
    from OCC.Core.TopAbs import TopAbs_VERTEX
    from OCC.Core.BRep import BRep_Tool

    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    global_bounds = bbox.Get()  # (xmin, ymin, zmin, xmax, ymax, zmax)

    # Bounding box da face
    face_bbox = Bnd_Box()
    vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
    while vertex_explorer.More():
        vertex = vertex_explorer.Current()
        pnt = BRep_Tool.Pnt(vertex)
        face_bbox.Update(pnt.X(), pnt.Y(), pnt.Z())
        vertex_explorer.Next()
    face_bounds = face_bbox.Get()

    # Se algum vértice da face está quase colado ao bounding box global => é exterior
    tol = 1e-3
    for i in range(3):
        if abs(face_bounds[i] - global_bounds[i]) < tol or abs(face_bounds[i+3] - global_bounds[i+3]) < tol:
            return True
    return False

def classify_and_summarize_holes(shape, filepath=None):
    # --- Extrai nome de molde e peça ---
    mold, part = None, None
    if filepath:
        mold, part = extract_mold_and_part_from_step(filepath)
        mold_name = f"Molde: {mold}" if mold else "Molde: (não identificado)"
        part_name = f"Peça: {part}" if part else "Peça: (não identificada)"
    else:
        mold_name = "Molde: (desconhecido)"
        part_name = "Peça: (desconhecida)"

    # --- Bounding Box ---
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    length = abs(xmax - xmin)
    width = abs(ymax - ymin)

    # --- Lógica de furos ---
    axis_groups = group_faces_by_axis(shape)
    grouped_by_center = {}

    for faces in axis_groups.values():
        surface_type, _, adaptor = faces[0]
        if surface_type == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
            radius = adaptor.Cylinder().Radius()
        else:
            continue

        loc = axis.Location()
        center_key = (round(loc.X(), 2), round(loc.Y(), 2))

        if center_key not in grouped_by_center:
            grouped_by_center[center_key] = (faces, radius)
        else:
            if radius > grouped_by_center[center_key][1]:
                grouped_by_center[center_key] = (faces, radius)

    count_cylindrical_through = 0
    count_cylindrical_closed = 0
    count_rectangular_through = 0
    count_rectangular_closed = 0

    # Novo: verifica se as faces do grupo tocam superfícies externas
    for faces, _ in grouped_by_center.values():
    # Assume que todas as faces do grupo pertencem ao mesmo furo
        _, face, adaptor = faces[0]
        loc = adaptor.Cylinder().Axis().Location()
        dir = adaptor.Cylinder().Axis().Direction()
        if is_hole_through(loc, dir, shape):
            count_cylindrical_through += 1
        else:
            count_cylindrical_closed += 1

    # Furos retangulares: lógica ainda por implementar
    count_rectangular_through = 0
    count_rectangular_closed = 0

    # ------- PRINT UNIFICADO ---------
    print("="*26)
    print("      RESUMO GERAL")
    print("="*26)
    print(f"{mold_name}")
    print(f"{part_name}")
    print(f"Comprimento da peça: {length:.2f}")
    print(f"Largura da peça: {width:.2f}")
    print("-"*26)
    print(f"Furos redondos sem fundo (perfuram): {count_cylindrical_through}")
    print(f"Furos redondos com fundo (fechados): {count_cylindrical_closed}")
    print(f"Furos retangulares sem fundo (perfuram): {count_rectangular_through}")
    print(f"Furos retangulares com fundo (fechados): {count_rectangular_closed}")
    print("\n(Nota: Apenas furos com maior raio por centro são contabilizados)\n")
