import math
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_VERTEX
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from collections import defaultdict, Counter
from .utils import is_hole_through, classify_hole_group_type, axis_key, get_all_vertices_of_group, group_faces_by_axis_and_proximity, get_auto_loc_tol, get_circular_hole_diameters_by_type

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
