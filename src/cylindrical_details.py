import math
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_VERTEX
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone

def get_max_radius_of_cone_face(face, adaptor):
    """Calcula o maior raio real dos vértices de uma face cônica."""
    cone = adaptor.Cone()
    ref_radius = cone.RefRadius()
    semi_angle = cone.SemiAngle()
    z0 = cone.Location().Z()
    
    vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
    max_radius = 0
    while vertex_explorer.More():
        vertex = vertex_explorer.Current()
        pnt = BRep_Tool.Pnt(vertex)
        # Raio real naquele ponto (Z)
        radius = abs(ref_radius + (pnt.Z() - z0) * math.tan(semi_angle))
        if radius > max_radius:
            max_radius = radius
        vertex_explorer.Next()
    return max_radius

def get_circular_hole_diameters(shape):
    # Dicionário para guardar o maior diâmetro por centro (X, Y)
    center_to_max_diameter = dict()
    explorer = TopExp_Explorer(shape, TopAbs_FACE)

    while explorer.More():
        face = explorer.Current()
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        stype = adaptor.GetType()

        if stype == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
            loc = axis.Location()
            x, y = round(loc.X(), 2), round(loc.Y(), 2)
            radius = adaptor.Cylinder().Radius()
        elif stype == GeomAbs_Cone:
            axis = adaptor.Cone().Axis()
            loc = axis.Location()
            x, y = round(loc.X(), 2), round(loc.Y(), 2)
            radius = get_max_radius_of_cone_face(face, adaptor)
        else:
            explorer.Next()
            continue

        diameter = round(radius * 2, 2)

        # Guarda apenas o maior diâmetro para cada (x, y)
        if (x, y) not in center_to_max_diameter or diameter > center_to_max_diameter[(x, y)]:
            center_to_max_diameter[(x, y)] = diameter

        explorer.Next()

    # Retorna apenas os maiores diâmetros por centro
    return list(center_to_max_diameter.values())

def show_circular_hole_details(diameters):
    if not diameters:
        print("-- Nenhum furo circular --\n")
        return

    from collections import Counter
    counter = Counter(diameters)  # agrupa os diâmetros já arredondados

    print("\n=== DETALHES DOS FUROS CIRCULARES ===")
    for diameter, count in sorted(counter.items(), key=lambda x: x[0], reverse=True):
        print(f"-- {count} furo(s) com {diameter} de diâmetro --")
