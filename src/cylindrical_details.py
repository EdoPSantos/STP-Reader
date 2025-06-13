import math
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_VERTEX
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.gp import gp_Pnt
from .utils import is_hole_through

def get_max_radius_of_cone_face(face, adaptor):
    # Calcula o maior raio real dos vértices de uma face cônica
    cone = adaptor.Cone()
    ref_radius = cone.RefRadius()
    semi_angle = cone.SemiAngle()
    z0 = cone.Location().Z()
    vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
    max_radius = 0
    while vertex_explorer.More():
        vertex = vertex_explorer.Current()
        pnt = BRep_Tool.Pnt(vertex)
        radius = abs(ref_radius + (pnt.Z() - z0) * math.tan(semi_angle))
        if radius > max_radius:
            max_radius = radius
        vertex_explorer.Next()
    return max_radius

def get_circular_hole_diameters_by_type(shape):
    """
    Agrupa faces cilíndricas e cônicas pelo centro e direção.
    Usa a função is_hole_through para distinguir furos passantes e fechados.
    """
    from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
    from OCC.Core.GeomAdaptor import GeomAdaptor_Surface

    center_to_data = dict()
    explorer = TopExp_Explorer(shape, TopAbs_FACE)

    while explorer.More():
        face = explorer.Current()
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        stype = adaptor.GetType()

        if stype == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
            loc = axis.Location()
            dir = axis.Direction()
            x, y, z = loc.X(), loc.Y(), loc.Z()
            radius = adaptor.Cylinder().Radius()
        elif stype == GeomAbs_Cone:
            axis = adaptor.Cone().Axis()
            loc = axis.Location()
            dir = axis.Direction()
            x, y, z = loc.X(), loc.Y(), loc.Z()
            radius = get_max_radius_of_cone_face(face, adaptor)
        else:
            explorer.Next()
            continue

        diameter = round(radius * 2, 2)
        key = (round(x, 2), round(y, 2), round(z, 2), round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3))
        if key not in center_to_data:
            center_to_data[key] = []
        center_to_data[key].append((face, diameter, loc, dir))
        explorer.Next()

    diameters_through = []
    diameters_closed = []

    for group in center_to_data.values():
        # Assume que todas as faces do grupo têm mesmo centro e direção
        face, diameter, loc, dir = group[0]
        if is_hole_through(loc, dir, shape):
            diameters_through.append(diameter)
        else:
            diameters_closed.append(diameter)
    return diameters_through, diameters_closed

def show_circular_hole_details(diameters_through, diameters_closed):
    from collections import Counter

    print("\n=== DETALHES DOS FUROS CIRCULARES SEM FUNDO ===")
    if not diameters_through:
        print("-- Nenhum furo circular sem fundo --\n")
    else:
        counter = Counter(diameters_through)
        for diameter, count in sorted(counter.items(), key=lambda x: x[0], reverse=True):
            print(f"-- {count} furo(s) com {diameter} de diâmetro --")

    print("\n=== DETALHES DOS FUROS CIRCULARES COM FUNDO ===")
    if not diameters_closed:
        print("-- Nenhum furo circular com fundo --\n")
    else:
        counter = Counter(diameters_closed)
        for diameter, count in sorted(counter.items(), key=lambda x: x[0], reverse=True):
            print(f"-- {count} furo(s) com {diameter} de diâmetro --")
