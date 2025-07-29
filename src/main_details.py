from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRep import BRep_Tool
from collections import defaultdict
from OCC.Core.Bnd import Bnd_Box
import OCC.Core.TopExp
from .utils import (
    get_bbox, get_face_area,
    group_faces_by_axis_and_proximity, get_auto_loc_tol,
    get_all_radii_of_group, collect_semi_circular_arcs,
    group_semi_circular_arcs, group_holes_by_center,
    extract_mold_and_part_from_step, detect_positions_from_holes,
    merge_aligned_circular_features,
    get_all_planar_faces_bbox,
    feat_key, is_point_inside_bbox
)

#-------------------------------------------------------------------------------------------------
# 1. Furos retangulares/quadrados
def get_rectangular_features(shape):
    all_faces = get_all_planar_faces_bbox(shape)
    chapa_bbox = get_bbox(shape)
    xmin_c, ymin_c, zmin_c, xmax_c, ymax_c, zmax_c = chapa_bbox

    def face_touches_border(face_bbox, plate_bbox, tol=0.5):
        xmin, ymin, _, xmax, ymax, _ = face_bbox
        xmin_c, ymin_c, _, xmax_c, ymax_c, _ = plate_bbox
        if abs(xmin - xmax) < tol:
            if abs(xmin - xmin_c) < tol or abs(xmin - xmax_c) < tol:
                return True
        if abs(ymin - ymax) < tol:
            if abs(ymin - ymin_c) < tol or abs(ymin - ymax_c) < tol:
                return True
        if (abs(xmin - xmin_c) < tol and abs(xmax - xmax_c) < tol and
                abs(ymin - ymin_c) < tol and abs(ymax - ymax_c) < tol):
            return True
        return False

    faces_xy = [f for f in all_faces if zmin_c <= f['center'][2] <= zmax_c and not face_touches_border(f['bbox'], chapa_bbox)]
    tol = 1e-3
    n = len(faces_xy)
    conexoes = [[] for _ in range(n)]
    for i, f1 in enumerate(faces_xy):
        bbox1 = f1['bbox']
        x1min, y1min, _, x1max, y1max, _ = bbox1
        for j, f2 in enumerate(faces_xy):
            if j == i:
                continue
            bbox2 = f2['bbox']
            x2min, y2min, _, x2max, y2max, _ = bbox2
            x_conectado = (
                abs(x1min - x2max) < tol or abs(x1max - x2min) < tol or
                abs(x1min - x2min) < tol or abs(x1max - x2max) < tol
            )
            y_conectado = (
                abs(y1min - y2max) < tol or abs(y1max - y2min) < tol or
                abs(y1min - y2min) < tol or abs(y1max - y2max) < tol
            )
            if x_conectado and y_conectado:
                conexoes[i].append(j)

    usados = set()
    grupos = []

    def dfs(idx, grupo):
        usados.add(idx)
        grupo.append(faces_xy[idx])
        for viz in conexoes[idx]:
            if viz not in usados:
                dfs(viz, grupo)

    for i in range(n):
        if i not in usados:
            grupo = []
            dfs(i, grupo)
            grupos.append(grupo)

    grupos_filtrados = []
    for grupo in grupos:
        grupo_filtrado = [f for f in grupo if abs(f['bbox'][5] - f['bbox'][2]) > 1.0]
        if len(grupo_filtrado) >= 2:
            grupos_filtrados.append(grupo_filtrado)

    grupos_planos = []
    grupos_cilindros = []
    for grupo in grupos_filtrados:
        tipos = set(f.get('type', '') for f in grupo)
        if tipos <= {'plano'}:
            grupos_planos.append(grupo)
        elif tipos <= {'cilindro', 'cone'} or tipos == {'cilindro'} or tipos == {'cone'}:
            grupos_cilindros.append(grupo)

    grupos_excluidos_centros = [set(tuple(f['center']) for f in gr) for gr in grupos_planos + grupos_cilindros]

    rectangular_counter = []
    grupo_id = 1
    tol_bbox = 1.0  # tolerância para comparar bbox do grupo com bbox da chapa
    for grupo in grupos:
        grupo_filtrado = [f for f in grupo if abs(f['bbox'][5] - f['bbox'][2]) > 1.0]
        if len(grupo_filtrado) < 2:
            continue
        centros_grupo = set(tuple(f['center']) for f in grupo_filtrado)
        if any(centros_grupo == centros_excl for centros_excl in grupos_excluidos_centros):
            continue
        xmins = [f['bbox'][0] for f in grupo_filtrado]
        ymins = [f['bbox'][1] for f in grupo_filtrado]
        xmaxs = [f['bbox'][3] for f in grupo_filtrado]
        ymaxs = [f['bbox'][4] for f in grupo_filtrado]
        xmin_g = min(xmins)
        xmax_g = max(xmaxs)
        ymin_g = min(ymins)
        ymax_g = max(ymaxs)
        # Verifica se bbox do grupo coincide com bbox da chapa
        if (abs(xmin_g - xmin_c) < tol_bbox and abs(xmax_g - xmax_c) < tol_bbox and
            abs(ymin_g - ymin_c) < tol_bbox and abs(ymax_g - ymax_c) < tol_bbox):
            continue  # ignora grupo igual à chapa
        comprimento_g = round(xmax_g - xmin_g, 3)
        largura_g = round(ymax_g - ymin_g, 3)
        rectangular_counter.append({
            'grupo_id': grupo_id,
            'num_faces': len(grupo_filtrado),
            'comprimento': comprimento_g,
            'largura': largura_g,
            'faces': grupo_filtrado
        })
        grupo_id += 1
    return rectangular_counter


# 2. Furos circulares (mas só únicos)
def get_circular_features(shape, semi_features=None):
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
            center_tuple = (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2))
        else:
            bbox = get_bbox(main_face)
            cx = (bbox[0] + bbox[3]) / 2
            cy = (bbox[1] + bbox[4]) / 2
            center_tuple = (round(cx, 2), round(cy, 2), 0.0)

        # Verificar duplicado com semicirculares (comparar centro e max_d)
        if semi_features and is_grouped_with_semi_arc(center_tuple, max_d, semi_features):
            continue

        # (Removido filtro de duplicados por centro/diâmetro)

        features.append({
            'center': center_tuple,
            'min_d': min_d,
            'max_d': max_d
        })

    return features

# 3. Recortes semicirculares (guardar centros e raios)
def get_semi_circular_features(shape, group_tol=3.0, lateral_tol=3.0):
    """
    Retorna features semicirculares agrupadas.
    """
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

# =====================
# FUNÇÃO PRINCIPAL DE SUMÁRIO
# =====================

def summarize_piece(shape, filepath=None):
    """
    Processa e retorna todos os dados sumarizados da peça para o show_general_summary.
    """
    semi_features = get_semi_circular_features(shape)
    circular_features = get_circular_features(shape, semi_features=semi_features)
    circular_features = merge_aligned_circular_features(circular_features)
    rectangular_counter = get_rectangular_features(shape)

    rectangular_bboxes = []
    for grupo in rectangular_counter:
        xmins = [f['bbox'][0] for f in grupo['faces']]
        ymins = [f['bbox'][1] for f in grupo['faces']]
        xmaxs = [f['bbox'][3] for f in grupo['faces']]
        ymaxs = [f['bbox'][4] for f in grupo['faces']]
        xmin_g = min(xmins)
        ymin_g = min(ymins)
        xmax_g = max(xmaxs)
        ymax_g = max(ymaxs)
        rectangular_bboxes.append((xmin_g, ymin_g, 0, xmax_g, ymax_g, 0))

    circular_inside_rect = []
    circular_outside_rect = []
    for feat in circular_features:
        center = feat['center']
        inside_any_rect = any(is_point_inside_bbox(center, bbox) for bbox in rectangular_bboxes)
        if inside_any_rect:
            circular_inside_rect.append(feat)
        else:
            circular_outside_rect.append(feat)

    hole_groups = group_faces_by_axis_and_proximity(shape, loc_tol=get_auto_loc_tol(shape))
    positions = detect_positions_from_holes(hole_groups, shape)

    # Agrupamento antigo: deduplicar por centro e diâmetro
    unique_circ = {}
    for feat in circular_outside_rect:
        key = feat_key(feat)
        if key not in unique_circ:
            unique_circ[key] = feat
    circ_counter = group_holes_by_center(list(unique_circ.values()))
    unique_circ_list = list(unique_circ.values())

    unique_semi = {}
    for feat in semi_features:
        key = feat_key(feat)
        if key not in unique_semi:
            unique_semi[key] = feat
    semi_counter = group_holes_by_center(list(unique_semi.values()))

    if filepath:
        mold, part = extract_mold_and_part_from_step(filepath)
        mold_name = f"{mold}" if mold else "(não identificado)"
        part_name = f"{part}" if part else "(não identificada)"
    else:
        mold_name = "(desconhecido)"
        part_name = "(desconhecida)"

    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bbox_coords = bbox.Get()
    xmin, ymin, zmin, xmax, ymax, zmax = bbox_coords
    dim1 = abs(xmax - xmin)
    dim2 = abs(ymax - ymin)
    largura = max(dim1, dim2)
    comprimento = min(dim1, dim2)
    grupo_id = len(rectangular_counter)

    return {
        'mold_name': mold_name,
        'part_name': part_name,
        'largura': largura,
        'comprimento': comprimento,
        'positions': positions,
        'circ_counter': circ_counter,
        'semi_counter': semi_counter,
        'rectangular_counter': rectangular_counter,
        'grupo_id': grupo_id,
        'unique_circ_list': unique_circ_list
    }

def show_general_summary(shape, filepath=None):
    """
    Mostra um sumário geral da peça, incluindo furos, dimensões e agrupamentos.
    """
    summary = summarize_piece(shape, filepath)
    print("=" * 40)
    print("      RESUMO DA PEÇA PARA O GESTi")
    print("=" * 40)
    print(f"Molde: {summary['mold_name']}")
    print(f"Peça: {summary['part_name']}")
    print(f"Largura: {summary['largura']:.1f} mm")
    print(f"Comprimento: {summary['comprimento']:.1f} mm")
    print(f"Posições: {summary['positions']}")
    print(f"Qtd. furos circulares: {sum(summary['circ_counter'].values())}")
    print(f"Qtd. furos semicirculares: {sum(summary['semi_counter'].values())}")
    print(f"Qtd. furos quadrados/retangulares: {summary['grupo_id']}")
    print("\n" + "-" * 40)

    print("Furos circulares:")
    circ_counter = summary['circ_counter']
    if circ_counter:
        for (min_d, max_d), count in sorted(circ_counter.items(), key=lambda x: -x[0][1]):
            if min_d == max_d:
                print(f"  Tem {count} furo(s) com {min_d:.1f} mm de diâmetro")
            else:
                print(f"  Tem {count} furo(s) de {min_d:.1f} mm a {max_d:.1f} mm de diâmetro")
    else:
        print("  Nenhum furo circular encontrado.")
    print("\n" + "-" * 40)

    print("Furos semicirculares:")
    semi_counter = summary['semi_counter']
    if semi_counter:
        for (min_d, max_d), count in sorted(semi_counter.items(), key=lambda x: -x[0][1]):
            if min_d == max_d:
                print(f"  Tem {count} furo(s) com {min_d:.1f} mm de diâmetro")
            else:
                print(f"  Tem {count} furo(s) de {min_d:.1f} mm a {max_d:.1f} mm de diâmetro")
    else:
        print("  Nenhum furo semicircular encontrado.")
    print("\n" + "-" * 40)

    print("Furos Retangulares:")
    rectangular_counter = summary['rectangular_counter']
    if rectangular_counter:
        agrupados = defaultdict(int)
        for grupo in rectangular_counter:
            c = grupo['comprimento']
            l = grupo['largura']
            maior = max(c, l)
            menor = min(c, l)
            key = (maior, menor)
            agrupados[key] += 1  # Conta grupos, não faces
        for (comprimento, largura), total in sorted(agrupados.items(), key=lambda x: (-x[1], -x[0][0], -x[0][1])):
            print(f"  Tem {total} furo(s) de {comprimento:.1f} mm de comprimento e {largura:.1f} mm de largura")                      
    else:
        print("  Nenhum furo retangular encontrado.")

    print("\n" + "=" * 40)

    print("\n" + "-" * 16 + "Testes" + "-" * 16)

    # DEBUG: Mostrar detalhes de todos os furos circulares do resumo
    print("\n[DEBUG] Detalhes completos dos furos circulares do resumo:")
    circ_counter = summary['circ_counter']
    unique_circ_list = summary.get('unique_circ_list', [])
    if circ_counter:
        for (min_d, max_d), count in sorted(circ_counter.items(), key=lambda x: -x[0][1]):
            print(f"\nGrupo: min_d={min_d}, max_d={max_d}, quantidade={count}")
            # Listar todos os furos deste grupo
            for feat in unique_circ_list:
                if round(feat['min_d'],1) == min_d and round(feat['max_d'],1) == max_d:
                    print(f"  - centro={feat['center']}, min_d={feat['min_d']}, max_d={feat['max_d']}")
    else:
        print("Nenhum furo circular encontrado.")
    
    print("\n" + "-" * 40)
