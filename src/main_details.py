import os
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib_Add
from collections import defaultdict
from .utils import axis_key

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

def classify_and_summarize_holes(shape, filepath=None):
    # --- Info de Molde e Peça ---
    if filepath:
        basename = os.path.basename(filepath)
        name, _ = os.path.splitext(basename)
        mold_name = f"Molde: {name}"
        part_name = f"Peça: {name}"
    else:
        mold_name = "Molde: (desconhecido)"
        part_name = "Peça: (desconhecida)"

    # --- Bounding Box usando método recomendado ---
    bbox = Bnd_Box()
    brepbndlib_Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    length = abs(xmax - xmin)
    width = abs(ymax - ymin)

    # --- Lógica de furos (igual ao teu código) ---
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

    for faces, _ in grouped_by_center.values():
        if len(faces) >= 2:
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
