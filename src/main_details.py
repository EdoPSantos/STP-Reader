from OCC.Core.BRep import BRep_Tool
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.TopExp import TopExp_Explorer
import OCC.Core.TopExp
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_Plane
from collections import Counter, defaultdict
from .utils import (
    find_slots, group_non_circular_planar_faces,
    auto_filter_rectangular_faces, get_piece_dimensions,
    get_bbox, is_open_edge, get_face_area,
    group_faces_by_axis_and_proximity, get_auto_loc_tol,
    get_all_radii_of_group, collect_semi_circular_arcs,
    group_semi_circular_arcs, group_holes_by_center,
    extract_mold_and_part_from_step, detect_positions_from_holes,
    merge_aligned_circular_features,
    get_all_planar_faces_bbox,
)

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
            center_tuple = (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2))
        else:
            bbox = get_bbox(main_face)
            cx = (bbox[0] + bbox[3]) / 2
            cy = (bbox[1] + bbox[4]) / 2
            center_tuple = (round(cx, 2), round(cy, 2), 0.0)

        # Verificar duplicado com semicirculares (comparar centro e max_d)
        if semi_features and is_grouped_with_semi_arc(center_tuple, max_d, semi_features):
            continue

        # Ignora duplicados entre furos circulares (comparando centro e diâmetro)
        if any(
            abs(center_tuple[0] - f['center'][0]) < 3.0 and
            abs(center_tuple[1] - f['center'][1]) < 3.0 and
            abs(max_d - f['max_d']) < 2.0
            for f in features
        ):
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

def is_point_inside_bbox(point, bbox, tol=0.5):
    """
    Verifica se um ponto (x, y) está dentro do bounding box.
    """
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    x, y = point[:2]  # Corrigido para aceitar tupla de 3 elementos
    return (xmin - tol) <= x <= (xmax + tol) and (ymin - tol) <= y <= (ymax + tol)

# Função principal para o sumário geral
def show_general_summary(shape, filepath=None):
    semi_features = get_semi_circular_holes(shape)
    circular_features = get_circular_holes(shape, semi_features=semi_features)
    # Aplica agrupamento/filtragem para evitar duplicidades
    circular_features = merge_aligned_circular_features(circular_features)
    rectangular_counter = get_rectangular_holes(shape)

    # --- NOVO: Extrair bounding box de cada retângulo realmente contado
    rectangular_bboxes = []
    faces = group_non_circular_planar_faces(shape)
    filtered = auto_filter_rectangular_faces(faces)
    for (shape_type, rlength, rwidth), count in rectangular_counter.items():
        for face, stype, length, width in filtered:
            if stype == shape_type and abs(length - rlength) < 0.01 and abs(width - rwidth) < 0.01:
                bbox = get_bbox(face)
                rectangular_bboxes.append(bbox)

    # --- Separar os furos circulares
    circular_inside_rect = []
    circular_outside_rect = []

    for feat in circular_features:
        center = feat['center']
        inside_any_rect = any(is_point_inside_bbox(center, bbox) for bbox in rectangular_bboxes)
        if inside_any_rect:
            circular_inside_rect.append(feat)
        else:
            circular_outside_rect.append(feat)

    # Detectar número de posições (Top/Bottom) com base nos furos
    hole_groups = group_faces_by_axis_and_proximity(shape, loc_tol=get_auto_loc_tol(shape))
    positions = detect_positions_from_holes(hole_groups, shape)

    # --- Remover duplicados entre circulares e semicirculares (por centro e diâmetro)
    def feat_key(feat, center_tol=1.0, d_tol=0.5):
        return (round(feat['center'][0]/center_tol), round(feat['center'][1]/center_tol), round(feat['max_d']/d_tol))

    unique_feats = {}
    for feat in circular_outside_rect:
        key = feat_key(feat)
        if key not in unique_feats:
            unique_feats[key] = feat
    for feat in semi_features:
        key = feat_key(feat)
        if key not in unique_feats:
            unique_feats[key] = feat
    all_counter = group_holes_by_center(list(unique_feats.values()))

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
# --------------------------------------- Teste --------------------------------------- 
    """
    print(f"=== DEBUG : Lista de furos circulares encontrados ===")
    for f in circular_features:
        print(f"Centro: {f['center']}, min_d: {f['min_d']}, max_d: {f['max_d']}")
    print(f"=== FIM DEBUG ===")

    print(f"=== DEBUG : Lista de furos semicirculares encontrados ===")
    for f in semi_features:
        print(f"Centro: {f['center']}, min_d: {f['min_d']}, max_d: {f['max_d']}")
    print(f"=== FIM DEBUG ===")

    print(f"=== DEBUG : Retângulos identificados ===")
    for face, shape_type, length, width in filtered:
        print(f"Tipo: {shape_type}, comprimento: {length}, largura: {width}")
    print(f"=== FIM DEBUG ===")

    print(f"=== DEBUG : Retângulos contados ===")
    for (shape_type, rlength, rwidth), count in rectangular_counter.items():
        print(f"Tipo: {shape_type}, comprimento: {rlength}, largura: {rwidth}, count: {count}")
    print(f"=== FIM DEBUG ===")
    """
# -------------------------------------------------------------------------------------
    
    print("Furos circulares (inclui semicirculares):")
    if all_counter:
        for (min_d, max_d), count in sorted(all_counter.items(), key=lambda x: -x[0][1]):
            if min_d == max_d:
                print(f"  Tem {count} furo(s) com {min_d:.1f} mm de diâmetro")
            else:
                print(f"  Tem {count} furo(s) de {min_d:.1f} mm a {max_d:.1f} mm de diâmetro")
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
    
    print("\n" + "-" * 40)

    all_faces = get_all_planar_faces_bbox(shape)
    """
    print("\n=== Todas as formas geométricas (planas, curvas, circulares, etc.) ===")
    if all_faces:
        for idx, f in enumerate(all_faces, 1):
            bbox = f['bbox']
            xmin, ymin, zmin, xmax, ymax, zmax = bbox
            tipo = f.get('type', '')
            center = f.get('center', (None, None, None))
            # Corrige apenas para cilindros: garante que ymin < ymax e o centro Y está entre eles
            if tipo == 'cilindro' and ymin > ymax:
                ymin, ymax = ymax, ymin
            # Opcional: garantir que o centro Y está entre ymin e ymax (para casos raros)
            if tipo == 'cilindro' and not (ymin <= center[1] <= ymax):
                ymin, ymax = min(ymin, ymax), max(ymin, ymax)
            bbox_fmt = (round(xmin, 3), round(ymin, 3), round(zmin, 3), round(xmax, 3), round(ymax, 3), round(zmax, 3))
            print(f"\n{idx}) Tipo: {tipo}, Centro: {center}, \n  1- Comprimento: {f['length']} mm, Largura: {f['width']} mm, Altura: {f['height']} mm, \n  2- BBox: {bbox_fmt}")
    else:
        print("Nenhuma face geométrica encontrada.")"""

# -------------------------------------------------------------------------------------
    # === Grupos de furos retangulares (apenas faces de all_faces, agrupando por extremos x/y) ===
    print("\n=== Grupos de furos retangulares (por conectividade XY das faces exibidas) ===")

    chapa_bbox = get_bbox(shape)
    xmin_c, ymin_c, zmin_c, xmax_c, ymax_c, zmax_c = chapa_bbox

    # Filtra faces para garantir que estão dentro do z da chapa
    faces_xy = [f for f in all_faces if zmin_c <= f['center'][2] <= zmax_c]

    # Agrupa faces por conectividade XY (se algum extremo x/y de uma face "toca" outro de outra face)
    # --- Agrupamento por conectividade transitiva (DFS) ---
    tol = 1e-3
    n = len(faces_xy)
    conexoes = [[] for _ in range(n)]
    # Monta lista de conexões
    for i, f1 in enumerate(faces_xy):
        bbox1 = f1['bbox']
        x1min, y1min, _, x1max, y1max, _ = bbox1
        for j, f2 in enumerate(faces_xy):
            if j == i:
                continue
            bbox2 = f2['bbox']
            x2min, y2min, _, x2max, y2max, _ = bbox2
            x_conectado = (
                abs(x1min - x2max) < tol or abs(x1max - x2min) < tol or abs(x1min - x2min) < tol or abs(x1max - x2max) < tol
            )
            y_conectado = (
                abs(y1min - y2max) < tol or abs(y1max - y2min) < tol or abs(y1min - y2min) < tol or abs(y1max - y2max) < tol
            )
            if x_conectado and y_conectado:
                conexoes[i].append(j)

    # Busca em profundidade para encontrar componentes conexos
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

    print(f"DEBUG: {len(grupos)} grupos encontrados por conectividade XY")

    # --- Resumo dos grupos ---
    print("\nResumo dos grupos:")
    resumo = []
    for grupo in grupos:
        if len(grupo) < 2:
            continue
        xmins = [f['bbox'][0] for f in grupo]
        ymins = [f['bbox'][1] for f in grupo]
        xmaxs = [f['bbox'][3] for f in grupo]
        ymaxs = [f['bbox'][4] for f in grupo]
        xmin_g = min(xmins)
        xmax_g = max(xmaxs)
        ymin_g = min(ymins)
        ymax_g = max(ymaxs)
        comprimento_g = round(xmax_g - xmin_g, 3)
        largura_g = round(ymax_g - ymin_g, 3)
        resumo.append((len(grupo), comprimento_g, largura_g))
    for idx, (nfaces, comp, larg) in enumerate(resumo, 1):
        print(f"  Grupo {idx}: {nfaces} faces | Comprimento: {comp} mm | Largura: {larg} mm")

    # --- Detalhamento dos grupos ---
    grupo_id = 1
    for grupo in grupos:
        if len(grupo) < 2:
            continue
        xmins = [f['bbox'][0] for f in grupo]
        ymins = [f['bbox'][1] for f in grupo]
        xmaxs = [f['bbox'][3] for f in grupo]
        ymaxs = [f['bbox'][4] for f in grupo]
        xmin_g = min(xmins)
        xmax_g = max(xmaxs)
        ymin_g = min(ymins)
        ymax_g = max(ymaxs)
        comprimento_g = round(xmax_g - xmin_g, 3)
        largura_g = round(ymax_g - ymin_g, 3)
        print(f"\nGrupo {grupo_id} (({len(grupo)}) faces conectadas XY): Comprimento grupo: {comprimento_g} mm, Largura grupo: {largura_g} mm\n")
        grupo_id += 1
        for idx_grupo, info in enumerate(grupo, 1):
            tipo = info.get('type', '')
            center = info.get('center', (None, None, None))
            bbox = info['bbox']
            xmin, ymin, zmin, xmax, ymax, zmax = bbox
            # Corrige apenas para cilindros: garante que ymin < ymax e o centro Y está entre eles
            if tipo == 'cilindro' and ymin > ymax:
                ymin, ymax = ymax, ymin
            if tipo == 'cilindro' and not (ymin <= center[1] <= ymax):
                ymin, ymax = min(ymin, ymax), max(ymin, ymax)
            bbox_fmt = (round(xmin, 3), round(ymin, 3), round(zmin, 3), round(xmax, 3), round(ymax, 3), round(zmax, 3))
            comprimento_bbox = round(xmax - xmin, 3)
            largura_bbox = round(ymax - ymin, 3)
            print(f"  Face (idx {idx_grupo} de {len(grupo)}): Tipo: {tipo}, Centro: {center}, "
                  f"\n    Comprimento: {info['length']} mm, Largura: {info['width']} mm, Altura: {info['height']} mm"
                  f"\n    BBox: {bbox_fmt}\n")
    
    print("\n" + "=" * 40)