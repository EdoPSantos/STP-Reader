from OCC.Core.BRep import BRep_Tool
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_EDGE
from collections import Counter
from .utils import (
    extract_mold_and_part_from_step, detect_positions_from_holes,
    group_faces_by_axis_and_proximity, get_auto_loc_tol,
    group_non_circular_planar_faces, auto_filter_rectangular_faces,
    get_piece_dimensions, find_slots, get_bbox, is_open_edge,
    get_all_radii_of_group, get_face_area,
    collect_semi_circular_arcs, group_semi_circular_arcs
)
import os

# Função utilitária para comparar centros e diâmetros com tolerância
def is_duplicate_circle(center1, d1, list_centers, list_ds, center_tol=3.0, d_tol=2.0):
    for i, c2 in enumerate(list_centers):
        if abs(center1[0] - c2[0]) < center_tol and abs(center1[1] - c2[1]) < center_tol:
            if abs(d1 - list_ds[i]) < d_tol:
                return True
    return False

# 1. Furos retangulares/quadrados
def get_rectangular_holes(shape):
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
        if is_open_edge(face, piece_bbox, tol=2.0):
            continue
        is_face_exterior = (
            (abs(rlength - piece_length) < abs_tol and abs(rwidth - piece_width) < abs_tol) or
            (abs(rlength - piece_width) < abs_tol and abs(rwidth - piece_length) < abs_tol)
        )
        if is_face_exterior and not is_open_edge(face, piece_bbox, tol=abs_tol):
            continue
        if any(abs(rlength - slen) < 2.0 and abs(rwidth - swid) < 2.0 for (slen, swid) in slot_counter):
            continue
        if (
            rlength > rel_tol * piece_length or
            rwidth > rel_tol * piece_width or
            rlength > rel_tol * piece_width or
            rwidth > rel_tol * piece_length
        ):
            continue
        key = (shape_type, rlength, rwidth)
        rectangular_counter[key] += 1
    return rectangular_counter

# 2. Furos circulares (mas só únicos)
def get_circular_holes(shape, semi_features=None):
    """
    Retorna todos os furos circulares completos (ignorando duplicados semicirculares, se fornecido).
    """
    axis_groups = group_faces_by_axis_and_proximity(shape, loc_tol=get_auto_loc_tol(shape))
    features = []

    for faces in axis_groups:
        main_face = max((f[1] for f in faces), key=get_face_area)
        radii, _ = get_all_radii_of_group(faces)
        if not radii:
            continue

        min_d = round(min(radii) * 2, 2)
        max_d = round(max(radii) * 2, 2)

        from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
        from OCC.Core.GeomAbs import GeomAbs_Cylinder
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(main_face))

        if adaptor.GetType() == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
            loc = axis.Location()
            center_tuple = (round(loc.X(), 2), round(loc.Y(), 2))
        else:
            bbox = get_bbox(main_face)
            cx = (bbox[0] + bbox[3]) / 2
            cy = (bbox[1] + bbox[4]) / 2
            center_tuple = (round(cx, 2), round(cy, 2))

        # Verificar duplicado com semicirculares (comparar centro e max_d)
        if semi_features and is_grouped_with_semi_arc(center_tuple, max_d, semi_features):
            continue

        features.append({
            'center': center_tuple,
            'min_d': min_d,
            'max_d': max_d
        })

    return features

# 3. Recortes semicirculares (guardar centros e raios)
def get_semi_circular_holes(shape, group_tol=3.0, lateral_tol=3.0):
    raw_semicircles = collect_semi_circular_arcs(shape, lateral_tol=lateral_tol)
    grouped_semicircles = group_semi_circular_arcs(raw_semicircles, group_tol=group_tol)
    semi_features = []
    for g in grouped_semicircles:
        rmin, rmax = round(min(g['radii']), 2), round(max(g['radii']), 2)
        center = (round(g['xy'][0],2), round(g['xy'][1],2))
        semi_features.append({
            'center': center,
            'min_d': rmin,
            'max_d': rmax,
            'group': g  # Para debug/impressão
        })
    return semi_features

def is_near(center1, center2, tol=3.0):
    return abs(center1[0] - center2[0]) < tol and abs(center1[1] - center2[1]) < tol

def is_same_diameter(d1, d2, tol=2.0):
    return abs(d1 - d2) < tol

# Função para verificar se o círculo é duplicado por algum semicircular
def is_grouped_with_semi_arc(center, d, semi_features, center_tol=3.0, d_tol=2.0):
    """
    Verifica se um círculo completo está dentro de um grupo de recortes semicirculares.
    """
    for sf in semi_features:
        group = sf.get('group')
        if not group:
            continue
        gx, gy = group['xy']
        if abs(center[0] - gx) < center_tol and abs(center[1] - gy) < center_tol:
            for radius in group['radii']:
                if abs(d - 2 * radius) < d_tol:
                    return True
    return False


# Função principal para o sumário geral
def show_general_summary(shape, filepath=None):
    semi_features = get_semi_circular_holes(shape)
    circular_features = get_circular_holes(shape, semi_features=semi_features)
    rectangular_counter = get_rectangular_holes(shape)

    # Detectar número de posições (Top/Bottom) com base nos furos
    hole_groups = group_faces_by_axis_and_proximity(shape, loc_tol=get_auto_loc_tol(shape))
    positions = detect_positions_from_holes(hole_groups, shape)

    # Contadores por faixa de diâmetro para circulares e semicirculares juntos
    all_counter = Counter()
    # Circulares (ignora duplicados com semicirculares)
    for feat in circular_features:
        all_counter[(feat['min_d'], feat['max_d'])] += 1
    # Semicirculares (cada um como se fosse um furo circular extra, como pediste)
    for semi in semi_features:
        all_counter[(semi['min_d'], semi['max_d'])] += 1

    # --- Mold e peça
    if filepath:
        mold, part = extract_mold_and_part_from_step(filepath)
        mold_name = f"{mold}" if mold else "(não identificado)"
        part_name = f"{part}" if part else "(não identificada)"
    else:
        mold_name = "(desconhecido)"
        part_name = "(desconhecida)"

    # --- Dimensões globais
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bbox_coords = bbox.Get()
    xmin, ymin, zmin, xmax, ymax, zmax = bbox_coords
    dim1 = abs(xmax - xmin)
    dim2 = abs(ymax - ymin)
    largura = max(dim1, dim2)
    comprimento = min(dim1, dim2)

    # --- OUTPUT
    print("=" * 40)
    print("   RESUMO DA PEÇA PARA O GESTi")
    print("=" * 40)
    print(f"Molde: {mold_name}")
    print(f"Peça: {part_name}")
    print(f"Largura: {largura:.1f} mm")
    print(f"Comprimento: {comprimento:.1f} mm")
    print(f"Posições: {positions}")
    print(f"Qtd. furos circulares (inclui semicirculares): {sum(all_counter.values())}")
    print(f"Qtd. furos quadrados/retangulares: {sum(rectangular_counter.values())}")
    print("\n" + "-" * 40)

    print("Furos circulares (inclui semicirculares):")
    if all_counter:
        for (min_d, max_d), count in sorted(all_counter.items(), key=lambda x: -x[0][1]):
            if min_d == max_d:
                print(f"  Tem {count} furo(s) com {min_d} mm de diâmetro")
            else:
                print(f"  Tem {count} furo(s) de {min_d} mm a {max_d} mm de diâmetro")
    else:
        print("  Nenhum furo circular encontrado.")
        
    print("\n" + "-" * 40)

    if rectangular_counter:
        print("Furos retangulares/quadrados internos:")
        for (shape_type, rlength, rwidth), count in sorted(rectangular_counter.items(), key=lambda x: (-x[0][1], -x[0][2])):
            if shape_type == "square":
                print(f"  Tem {count} furo(s) [quadrado] com {rwidth} mm de lado")
            else:
                print(f"  Tem {count} furo(s) [retangular] com {rlength} mm de comprimento e {rwidth} mm de largura")
    else:
        print("Furos retangulares/quadrados internos: não extraídos")
    print("\n" + "=" * 40)