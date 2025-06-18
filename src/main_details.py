from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from collections import Counter
from .utils import (
    extract_mold_and_part_from_step, detect_positions_from_holes,
    group_faces_by_axis_and_proximity, get_auto_loc_tol,
    group_non_circular_planar_faces, auto_filter_rectangular_faces,
    get_piece_dimensions, find_slots, get_bbox, is_open_edge
)

import os

def show_general_summary(shape, filepath=None):
    """
    Mostra um resumo geral da peça: dimensões, quantidade de furos circulares e retangulares/quadrados,
    e lista todos os diâmetros/medidas encontradas.
    """

    # 1. Identificação do molde e peça
    mold, part = None, None
    if filepath:
        mold, part = extract_mold_and_part_from_step(filepath)
        mold_name = f"{mold}" if mold else "(não identificado)"
        part_name = f"{part}" if part else "(não identificada)"
    else:
        mold_name = "(desconhecido)"
        part_name = "(desconhecida)"

    # 2. Dimensões globais da peça
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    length = abs(xmax - xmin)
    width = abs(ymax - ymin)

    # 3. Furos circulares (cilíndricos/conicos/mistos)
    loc_tol = get_auto_loc_tol(shape)
    axis_groups = group_faces_by_axis_and_proximity(shape, loc_tol=loc_tol)
    count_cylindrical = len(axis_groups)

    diameters = []
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
    counter_diam = Counter(diameters)

    # 4. Slots oblongos e furos retangulares/quadrados
    slots = find_slots(shape)
    slot_counter = Counter(slots)

    faces = group_non_circular_planar_faces(shape)
    filtered_faces = auto_filter_rectangular_faces(faces)
    piece_length, piece_width = get_piece_dimensions(shape)
    abs_tol = 2.0
    rel_tol = 0.97
    piece_bbox = get_bbox(shape)

    rectangular_counter = Counter()
    for face, shape_type, rlength, rwidth in filtered_faces:
        # Ignora slots abertos
        if is_open_edge(face, piece_bbox, tol=2.0):
            continue
        is_face_exterior = (
            (abs(rlength - piece_length) < abs_tol and abs(rwidth - piece_width) < abs_tol) or
            (abs(rlength - piece_width) < abs_tol and abs(rwidth - piece_length) < abs_tol)
        )
        if is_face_exterior and not is_open_edge(face, piece_bbox, tol=abs_tol):
            continue
        # Evita duplicados de slots
        if any(abs(rlength - slen) < 2.0 and abs(rwidth - swid) < 2.0 for (slen, swid) in slot_counter):
            continue
        # Elimina casos extremos (faces demasiado grandes)
        if (
            rlength > rel_tol * piece_length or
            rwidth > rel_tol * piece_width or
            rlength > rel_tol * piece_width or
            rwidth > rel_tol * piece_length
        ):
            continue
        key = (shape_type, rlength, rwidth)
        rectangular_counter[key] += 1

    # 5. Posições estimadas
    positions = detect_positions_from_holes(axis_groups, shape)

    # 6. Output formatado
    print("=" * 40)
    print("   RESUMO DA PEÇA PARA O GESTi")
    print("=" * 40)
    print(f"Molde: {mold_name}")
    print(f"Peça: {part_name}")
    print(f"Largura: {length:.0f} mm")
    print(f"Comprimento: {width:.0f} mm")
    print(f"Posições: {positions}")
    print(f"Qtd. furos circulares: {count_cylindrical}")
    print(f"Qtd. furos quadrados/retangulares: {sum(rectangular_counter.values())}")
    print("\n")
    print("-" * 40)

    # Furos circulares
    if diameters:
        print("Furos circulares:")
        for (min_d, max_d), count in sorted(counter_diam.items(), key=lambda x: -x[0][1]):  # Ordena pelo diâmetro máximo desc
            if min_d == max_d:
                print(f"  Tem {count} furo(s) com {min_d} mm de diâmetro")
            else:
                print(f"  Tem {count} furo(s) de {min_d} mm a {max_d} mm de diâmetro")
    else:
        print("Furos circulares: não extraídos")

    print("\n")
    print("-" * 40)
    # Furos retangulares/quadrados internos
    if rectangular_counter:
        print("Furos retangulares/quadrados internos:")
        for (shape_type, rlength, rwidth), count in sorted(rectangular_counter.items(), key=lambda x: (-x[0][1], -x[0][2])):
            if shape_type == "square":
                print(f"  Tem {count} furo(s) [quadrado] com {rwidth} mm de lado")
            else:
                print(f"  Tem {count} furo(s) [retangular] com {rlength} mm de comprimento e {rwidth} mm de largura")
    else:
        print("Furos retangulares/quadrados internos: não extraídos")
    print("\n")
    print("=" * 40)
