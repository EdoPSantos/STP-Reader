from OCC.Core.BRep import BRep_Tool
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepTools import breptools
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from collections import Counter
from .utils import (
    extract_mold_and_part_from_step, detect_positions_from_holes,
    group_faces_by_axis_and_proximity, get_auto_loc_tol,
    group_non_circular_planar_faces, auto_filter_rectangular_faces,
    get_piece_dimensions, find_slots, get_bbox, is_open_edge, 
    get_circular_hole_diameters_by_type, is_face_circular_complete,
    get_all_radii_of_group, is_face_half_circle, get_face_area
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
    bbox_coords = bbox.Get()  # <- IMPORTANTE: obter coordenadas!
    xmin, ymin, zmin, xmax, ymax, zmax = bbox_coords
    dim1 = abs(xmax - xmin)
    dim2 = abs(ymax - ymin)
    # O maior vai para "Largura", o menor para "Comprimento"
    largura = max(dim1, dim2)
    comprimento = min(dim1, dim2)

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

    # Agrupa furos circulares por eixo e proximidade (como nos detalhes cilíndricos)
    axis_groups = group_faces_by_axis_and_proximity(shape, loc_tol=get_auto_loc_tol(shape))
    positions = detect_positions_from_holes(axis_groups, shape)

    diameters = []

    # Processamento dos furos circulares/meia-lua
    for faces in axis_groups:
        radii, _ = get_all_radii_of_group(faces)
        if not radii:
            continue
        # Seleciona a face de maior área do grupo para testar o tipo
        main_face = max((f[1] for f in faces), key=get_face_area)
        outer_wire = breptools.OuterWire(main_face)
        edge_explorer = TopExp_Explorer(outer_wire, TopAbs_EDGE)
        edge_count = 0
        while edge_explorer.More():
            edge_count += 1
            edge_explorer.Next()
        print(f"Face: {main_face}, área: {get_face_area(main_face):.2f}, nº de arestas no wire: {edge_count}, raio mínimo: {min(radii)}")

        # AQUI: Passar as coordenadas, não o objeto!
        if is_face_half_circle(main_face, piece_bbox=bbox_coords):
            min_r = round(min(radii), 2)
            max_r = round(max(radii), 2)
            print(f"DEBUG: ENCONTREI UM MEIA-LUA com raio(s) {min_r} a {max_r}, face={main_face}")
            diameters.append((min_r, max_r))
        else:
            min_d = round(min(radii) * 2, 2)
            max_d = round(max(radii) * 2, 2)
            diameters.append((min_d, max_d))

    count_cylindrical = len(diameters)
    diam_counter = Counter()
    for min_d, max_d in diameters:
        diam_counter[(min_d, max_d)] += 1

    # 6. Output formatado (mantido igual)
    print("=" * 40)
    print("   RESUMO DA PEÇA PARA O GESTi")
    print("=" * 40)
    print(f"Molde: {mold_name}")
    print(f"Peça: {part_name}")
    print(f"Largura: {largura:.0f} mm")
    print(f"Comprimento: {comprimento:.0f} mm")
    print(f"Posições: {positions}")
    print(f"Qtd. furos circulares: {count_cylindrical}")
    print(f"Qtd. furos quadrados/retangulares: {sum(rectangular_counter.values())}")
    print("\n")
    print("-" * 40)

    # Furos circulares
    if diameters:
        print("Furos circulares:")
        for (min_d, max_d), count in sorted(diam_counter.items(), key=lambda x: -x[0][1]):  # Ordena pelo diâmetro máximo desc
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
