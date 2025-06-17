import math
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_VERTEX
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from collections import defaultdict, Counter
from .utils import is_hole_through, classify_hole_group_type, axis_key, get_all_vertices_of_group


def group_faces_by_axis(shape):
    """
    Agrupa todas as faces cónicas e cilíndricas por centro e direção (eixo).
    Um grupo pode conter qualquer ordem/tipo (cilindro/cone/cone/cilindro/...)
    """
    hole_groups = defaultdict(list)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        stype = adaptor.GetType()
        if stype in (GeomAbs_Cylinder, GeomAbs_Cone):
            try:
                key = axis_key(adaptor)
                hole_groups[key].append((stype, face, adaptor))
            except:
                pass
        explorer.Next()
    return list(hole_groups.values())


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


def get_circular_hole_diameters_by_type(shape):
    """
    Devolve uma lista [(min_d, max_d, tipo), ...] para furos passantes e fechados.
    Cada grupo representa 1 furo real (pode ser misto cone/cilindro/cilindro+cone+cilindro).
    """
    hole_groups = group_faces_by_axis(shape)
    diameters_through = []
    diameters_closed = []

    for faces in hole_groups:
        radii, _ = get_all_radii_of_group(faces)
        if not radii:
            continue
        min_d = round(min(radii) * 2, 2)
        max_d = round(max(radii) * 2, 2)
        tipo = classify_hole_group_type(faces)
        result = is_hole_through(faces, shape)
        data = (min_d, max_d, tipo)
        if result:
            diameters_through.append(data)
        else:
            diameters_closed.append(data)
    return diameters_through, diameters_closed


def format_diameter_range(min_d, max_d, tipo):
    if tipo.lower() in ["cónico", "misto"] and min_d != max_d:
        return f"{min_d} a {max_d}"
    else:
        return f"{min_d}"


def show_circular_hole_details(diameters_through, diameters_closed):
    print("\n=== DETALHES DOS FUROS CIRCULARES SEM FUNDO ===")
    if not diameters_through:
        print("-- Nenhum furo circular sem fundo --\n")
    else:
        counter = Counter((tipo, format_diameter_range(min_d, max_d, tipo)) for min_d, max_d, tipo in diameters_through)
        def sort_key(tuplo): 
            val = tuplo[1]
            try:
                return float(val.split()[0])
            except:
                return 0
        for (tipo, diam_range), count in sorted(counter.items(), key=lambda x: sort_key(x[0]), reverse=True):
            print(f"-- {count} furo(s) [{tipo}] com {diam_range} mm de diâmetro --")

    print("\n=== DETALHES DOS FUROS CIRCULARES COM FUNDO ===")
    if not diameters_closed:
        print("-- Nenhum furo circular com fundo --\n")
    else:
        counter = Counter((tipo, format_diameter_range(min_d, max_d, tipo)) for min_d, max_d, tipo in diameters_closed)
        def sort_key(tuplo): 
            val = tuplo[1]
            try:
                return float(val.split()[0])
            except:
                return 0
        for (tipo, diam_range), count in sorted(counter.items(), key=lambda x: sort_key(x[0]), reverse=True):
            print(f"-- {count} furo(s) [{tipo}] com {diam_range} mm de diâmetro --")
