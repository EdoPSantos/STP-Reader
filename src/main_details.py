from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from collections import defaultdict
from .utils import extract_from_filename, extract_mold_and_part_from_step, detect_positions_from_holes, group_faces_by_axis

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

    axis_groups = group_faces_by_axis(shape)
    count_cylindrical = len(axis_groups)

    diameters = []
    for faces in axis_groups:
        radii = []
        for stype, face, adaptor in faces:
            if stype == GeomAbs_Cylinder:
                r = adaptor.Cylinder().Radius()
                if r > 0:
                    radii.append(round(r * 2, 2))
        if radii:
            diameters.append(max(radii))

    count_rectangular = 0

    positions = detect_positions_from_holes(axis_groups, shape)

    print("=" * 40)
    print("   RESUMO DA PEÇA PARA O GESTi")
    print("=" * 40)
    print(f"Molde: {mold_name}")
    print(f"Peça: {part_name}")
    print(f"Largura: {width:.0f} mm")
    print(f"Comprimento: {length:.0f} mm")
    print(f"Posições: {positions}")
    print(f"Qtd. furos cilíndricos: {count_cylindrical}")
    print(f"Qtd. furos quadrados: {count_rectangular}")
    print("-" * 40)

    if diameters:
        print("Diâmetro dos furos cilíndricos:")
        counter = Counter(diameters)
        for diameter, count in sorted(counter.items(), key=lambda x: -x[0]):
            print(f"  Tem {count} furo(s) com {diameter} mm")
    else:
        print("Diâmetro dos furos cilíndricos: -")

    print("=" * 40)
