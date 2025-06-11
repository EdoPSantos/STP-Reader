from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Cone
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from collections import defaultdict
from .utils import axis_key, cone_axis_key

def count_cylindrical_holes_by_axis(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    axis_set = set()
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cylinder:
            axis_set.add(axis_key(adaptor))
        explorer.Next()
    return len(axis_set)

def count_conical_holes_by_axis(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    axis_set = set()
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cone:
            axis_set.add(cone_axis_key(adaptor))
        explorer.Next()
    return len(axis_set)

def analyze_general_summary(shape):
    from OCC.Core.gp import gp_Dir

    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    planar_faces = []

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)

        if adaptor.GetType() == GeomAbs_Plane:
            normal = adaptor.Plane().Axis().Direction()
            dir_vec = (round(abs(normal.X()), 1), round(abs(normal.Y()), 1), round(abs(normal.Z()), 1))
            if dir_vec in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]:
                planar_faces.append(face)
        explorer.Next()

    cylindrical_holes = count_cylindrical_holes_by_axis(shape)
    conical_holes = count_conical_holes_by_axis(shape)

    print("="*26)
    print("      RESUMO GERAL")
    print("="*26)
    print(f"Faces planas exteriores: {len(planar_faces)}")
    print(f"Furos cilíndricos distintos: {cylindrical_holes}")
    print(f"Furos cônicos distintos: {conical_holes}")
    print()

def analyze_cylindrical_holes(shape):
    print("="*26)
    print("    FUROS CILÍNDRICOS")
    print("="*26)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    holes = {}
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cylinder:
            key = axis_key(adaptor)
            cylinder = adaptor.Cylinder()
            radius = round(cylinder.Radius(), 2)
            diameter = round(2 * radius, 2)
            axis = cylinder.Axis()
            loc = axis.Location()
            dir = axis.Direction()
            if key not in holes:
                holes[key] = {
                    "radius": radius,
                    "diameter": diameter,
                    "center": (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2)),
                    "direction": (round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3))
                }
        explorer.Next()

    for i, (key, data) in enumerate(holes.items(), 1):
        print(f"Furo {i}:")
        print(f"  Centro: {data['center']}")
        print(f"  Direção: {data['direction']}")
        print(f"  Raio: {data['radius']:.2f}")
        print(f"  Diâmetro: {data['diameter']:.2f}")
        print()

def analyze_conical_holes(shape):
    print("="*26)
    print("     FUROS CÔNICOS")
    print("="*26)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    con_holes = {}
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cone:
            key = cone_axis_key(adaptor)
            cone = adaptor.Cone()
            radius = round(cone.RefRadius(), 2)
            diameter = round(2 * radius, 2)
            axis = cone.Axis()
            loc = axis.Location()
            dir = axis.Direction()
            if key not in con_holes:
                con_holes[key] = {
                    "center": (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2)),
                    "direction": (round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)),
                    "radius": radius,
                    "diameter": diameter
                }
        explorer.Next()

    for i, (key, data) in enumerate(con_holes.items(), 1):
        print(f"Furo cônico {i}:")
        print(f"  Centro: {data['center']}")
        print(f"  Direção: {data['direction']}")
        print(f"  Raio: {data['radius']:.2f}")
        print(f"  Diâmetro: {data['diameter']:.2f}")
        print()

    print(f"Total de furos cônicos: {len(con_holes)}\n")
