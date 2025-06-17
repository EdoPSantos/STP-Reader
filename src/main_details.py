from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from collections import defaultdict
from .utils import extract_mold_and_part_from_step, detect_positions_from_holes, group_faces_by_axis_and_proximity, get_auto_loc_tol

import os
import re

def show_general_summary(shape, filepath=None):
    from collections import Counter

    mold, part = None, None
    if filepath:
        mold, part = extract_mold_and_part_from_step(filepath)
        mold_name = f"{mold}" if mold else "(não identificado)"
        part_name = f"{part}" if part else "(não identificada)"
    else:
        mold_name = "(desconhecido)"
        part_name = "(desconhecida)"

    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    length = abs(xmax - xmin)
    width = abs(ymax - ymin)

    loc_tol = get_auto_loc_tol(shape)
    axis_groups = group_faces_by_axis_and_proximity(shape, loc_tol=loc_tol)
    count_cylindrical = len(axis_groups)

    diameters = []  # Lista de tuplos: (min_diameter, max_diameter)
    for faces in axis_groups:
        radii = []
        for stype, face, adaptor in faces:
            if stype == GeomAbs_Cylinder:
                r = adaptor.Cylinder().Radius()
                if r > 0:
                    radii.append(round(r * 2, 2))
            elif stype == GeomAbs_Cone:
                base_r = adaptor.Cone().RefRadius()
                if base_r > 0:
                    radii.append(round(base_r * 2, 2))
        if radii:
            min_d = min(radii)
            max_d = max(radii)
            diameters.append((min_d, max_d))

    from collections import Counter
    counter = Counter(diameters)

    count_rectangular = 0

    positions = detect_positions_from_holes(axis_groups, shape)

    print("=" * 40)
    print("   RESUMO DA PEÇA PARA O GESTi")
    print("=" * 40)
    print(f"Molde: {mold_name}")
    print(f"Peça: {part_name}")
    print(f"Largura: {length:.0f} mm")
    print(f"Comprimento: {width:.0f} mm")
    print(f"Posições: {positions}")
    print(f"Qtd. furos cilíndricos: {count_cylindrical}")
    print(f"Qtd. furos quadrados: {count_rectangular}")
    print("-" * 40)

    if diameters:
        print("Diâmetro dos furos cilíndricos:")
        for (min_d, max_d), count in sorted(counter.items(), key=lambda x: -x[0][1]):  # Ordena pelo diâmetro máximo decrescente
            if min_d == max_d:
                print(f"  Tem {count} furo(s) com {min_d} mm de diâmetro")
            else:
                print(f"  Tem {count} furo(s) de {min_d} mm a {max_d} mm de diâmetro")
    else:
        print("Diâmetro dos furos cilíndricos: -")

    print("=" * 40)
