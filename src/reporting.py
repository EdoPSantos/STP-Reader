"""
Módulo de formatação e relatórios para análise de peças.
Separa a lógica de formatação da análise geométrica.
"""

import math
from typing import Dict, List, Any, Tuple
from collections import defaultdict, Counter
from .data_structures import AnalysisWarning
from .utils import (
    angle_to_word, calculate_face_angle_from_dimensions,
    get_bbox, get_face_area, is_hole_through
)
import logging

logger = logging.getLogger(__name__)

def format_piece_header(summary: Dict[str, Any]) -> str:
    """Formata o cabeçalho do relatório da peça."""
    bbox = summary.get('bbox')
    if not bbox:
        return "ERRO: Informações de bbox não disponíveis"
    
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    comprimento = abs(xmax - xmin)
    largura = abs(ymax - ymin)
    altura = abs(zmax - zmin)
    
    header = []
    header.append("=" * 40)
    header.append("      RESUMO DA PEÇA PARA O GESTi")
    header.append("=" * 40)
    header.append(f"Molde: {summary['mold_name']}")
    header.append(f"Peça: {summary['part_name']}")
    header.append("\n" + "-" * 40)
    header.append(f"Comprimento: {comprimento:.1f} mm")
    header.append(f"Largura: {largura:.1f} mm")
    header.append(f"Altura: {altura:.1f} mm")
    header.append(f"Posições: {summary['positions']}")
    header.append("\n" + "-" * 40)
    
    # Contar furos
    valores_circ = list(summary['circ_counter'].values())
    valores_filtrados = [v for v in valores_circ if v <= 20]  # Remover valores fantasma
    total_grouped_circular = sum(valores_filtrados)
    
    header.append(f"Qtd. furos circulares: {total_grouped_circular}")
    header.append(f"Qtd. furos semicirculares: {sum(summary['semi_counter'].values())}")
    header.append(f"Qtd. furos quadrados/retangulares: {summary['grupo_id']}")
    header.append("\n" + "-" * 40)
    
    return "\n".join(header)

def format_warnings(warnings: List[AnalysisWarning]) -> str:
    """Formata warnings se existirem."""
    output = []
    
    # Agrupar warnings por tipo
    warnings_by_type = defaultdict(list)
    for warning in warnings:
        warnings_by_type[warning.warning_type].append(warning)
    
    # Formatar cada tipo
    for warning_type, type_warnings in warnings_by_type.items():
        if warning_type == 'circular':
            output.append("WARNING -> só aparece se for do tipo Circulo")
        elif warning_type == 'semicircular':
            output.append("WARNING -> só aparece se for do tipo semicirculares")
        elif warning_type == 'rectangular':
            output.append("WARNING -> só aparece se for do tipo Retangulares")
        
        for warning in type_warnings:
            output.append(f"  {warning.message}")
        output.append("-----")
    
    return "\n".join(output) if output else ""

def format_circular_holes(summary: Dict[str, Any], shape, axis_groups: List[Any]) -> Tuple[str, List[str]]:
    """
    Formata a seção de furos circulares.
    Retorna (texto_formatado, lista_de_direcoes_globais)
    """
    output = []
    todas_direcoes_globais = []
    
    # Warnings primeiro
    warnings_text = format_warnings([w for w in summary.get('warnings', []) if w.warning_type == 'circular'])
    if warnings_text:
        output.append(warnings_text)
    
    output.append("Furos circulares:")
    
    circ_counter = summary['circ_counter']
    unique_circ_list = summary.get('unique_circ_list', [])
    
    if circ_counter:
        bbox = summary.get('bbox')
        altura_chapa = abs(bbox[5] - bbox[2]) if bbox else 0
        
        for (min_d, max_d), count_from_counter in sorted(circ_counter.items(), key=lambda x: -x[0][1]):
            grupo_feats = [feat for feat in unique_circ_list 
                          if round(feat['min_d'], 1) == min_d and round(feat['max_d'], 1) == max_d]
            
            if not grupo_feats:
                continue
            
            # Processar cada feature
            features_detalhadas = []
            
            for feat in grupo_feats:
                feature_info = _process_circular_feature(feat, axis_groups, shape, altura_chapa)
                if feature_info:  # Filtrar furos com profundidade < 1mm
                    features_detalhadas.append(feature_info)
                    if feature_info.get('direcao'):
                        todas_direcoes_globais.append(feature_info['direcao'])
            
            # Agrupar e formatar
            if features_detalhadas:
                subgrupos_text = _format_circular_subgroups(features_detalhadas)
                output.extend(subgrupos_text)
    else:
        output.append("  Nenhum furo circular encontrado.")
    
    return "\n".join(output), todas_direcoes_globais

def _process_circular_feature(feat: Dict[str, Any], axis_groups: List[Any], 
                             shape, altura_chapa: float) -> Dict[str, Any]:
    """Processa uma feature circular individual."""
    center = feat['center']
    
    # Encontrar grupo correspondente
    matching_group = None
    for faces in axis_groups:
        main_face = max((f[1] for f in faces), key=get_face_area)
        from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
        from OCC.Core.BRep import BRep_Tool
        from OCC.Core.GeomAbs import GeomAbs_Cylinder
        
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(main_face))
        
        if adaptor.GetType() == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
            loc = axis.Location()
            face_center = (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2))
        else:
            bbox = get_bbox(main_face)
            cx = (bbox[0] + bbox[3]) / 2
            cy = (bbox[1] + bbox[4]) / 2
            face_center = (round(cx, 2), round(cy, 2), 0.0)
        
        if (abs(center[0] - face_center[0]) < 0.1 and 
            abs(center[1] - face_center[1]) < 0.1):
            matching_group = faces
            break
    
    # Calcular profundidade e tipo
    profundidade, tipo = _calculate_hole_depth_and_type(feat, matching_group, shape, altura_chapa)
    
    # Filtrar furos muito pequenos
    if isinstance(profundidade, (int, float)) and profundidade < 1.0:
        return None
    
    # Calcular diâmetros corretos
    actual_min_d, actual_max_d, eh_furo_escalonado = _calculate_feature_diameters(feat)
    
    # Calcular área e perímetro
    if 'total_area' in feat and feat['total_area'] > 0:
        area_real = feat['total_area']
    else:
        raio = actual_max_d / 2
        area_real = math.pi * raio ** 2
    
    raio = actual_max_d / 2
    perimetro = 2 * math.pi * raio
    area_display = area_real if tipo == "cego" else None
    
    # Analisar componentes geométricos
    components_info, direcao = _analyze_geometric_components(feat, eh_furo_escalonado, tipo, altura_chapa)
    
    return {
        'actual_min_d': actual_min_d,
        'actual_max_d': actual_max_d,
        'profundidade': profundidade,
        'tipo': tipo,
        'area_display': area_display,
        'perimetro': perimetro,
        'area_real': area_real,
        'components_info': components_info,
        'direcao': direcao
    }

def _calculate_hole_depth_and_type(feat: Dict[str, Any], matching_group, shape, altura_chapa: float) -> Tuple[Any, str]:
    """Calcula profundidade e tipo do furo."""
    profundidade = 'n/d'
    tipo = "passante"
    
    if matching_group:
        is_through = is_hole_through(matching_group, shape)
        
        if 'coordinates' in feat and feat['coordinates'] and 'z_range' in feat['coordinates']:
            profundidade = feat['coordinates']['z_range'].get('height', 'n/d')
            if isinstance(profundidade, (int, float)):
                profundidade = round(profundidade, 2)
        else:
            # Fallback: calcular baseado no bbox
            all_z_values = []
            for face_tuple in matching_group:
                face = face_tuple[1]
                face_bbox = get_bbox(face)
                all_z_values.extend([face_bbox[2], face_bbox[5]])
            if all_z_values:
                profundidade = round(abs(max(all_z_values) - min(all_z_values)), 2)
        
        # Determinar tipo baseado na profundidade
        if isinstance(profundidade, (int, float)):
            tolerancia = 0.1
            if profundidade < (altura_chapa - tolerancia):
                tipo = "cego"
            else:
                tipo = "passante"
        else:
            tipo = "passante" if is_through else "cego"
    
    return profundidade, tipo

def _calculate_feature_diameters(feat: Dict[str, Any]) -> Tuple[float, float, bool]:
    """Calcula diâmetros corretos da feature."""
    actual_min_d = feat['min_d']
    actual_max_d = feat['max_d']
    eh_furo_escalonado = False
    
    if 'geometric_components' in feat and feat['geometric_components']:
        cones = [c for c in feat['geometric_components'] if c['geometric_type'] == 'conica']
        cilindros = [c for c in feat['geometric_components'] if c['geometric_type'] == 'cilindrica']
        
        # Lógica complexa de determinação de diâmetros (simplificada aqui)
        if len(cones) >= 2 and len(cones) > len(cilindros):
            # Verificar cones simétricos
            alturas_cones = [c['altura'] for c in cones]
            altura_referencia = alturas_cones[0]
            cones_mesma_altura = sum(1 for altura in alturas_cones if abs(altura - altura_referencia) < 0.1)
            
            if cones_mesma_altura >= 2 and len(cilindros) >= 1:
                # Mostrar apenas diâmetro do cilindro
                diametros_cilindros = [c['diameter'] for c in cilindros]
                if diametros_cilindros:
                    diametro_cilindro_menor = min(diametros_cilindros)
                    actual_min_d = diametro_cilindro_menor
                    actual_max_d = diametro_cilindro_menor
        
        eh_furo_escalonado = (abs(actual_max_d - actual_min_d) > 0.1)
    
    return actual_min_d, actual_max_d, eh_furo_escalonado

def _analyze_geometric_components(feat: Dict[str, Any], eh_furo_escalonado: bool, 
                                 tipo: str, altura_chapa: float) -> Tuple[List[str], str]:
    """Analisa componentes geométricos e determina direção."""
    components_info = []
    direcao = None
    
    if 'geometric_components' in feat and feat['geometric_components']:
        for component in feat['geometric_components']:
            if component['geometric_type'] == 'cilindrica':
                components_info.append(f"Cilíndrico ⌀{component['diameter']:.1f}mm h{component['altura']:.1f}mm")
            elif component['geometric_type'] == 'conica':
                components_info.append(f"Cônico ⌀{component['diameter_min']:.1f}-{component['diameter_max']:.1f}mm h{component['altura']:.1f}mm")
                
                if eh_furo_escalonado or tipo == "cego":
                    if 'direction' in component and direcao is None:
                        if component['direction'] == 'expanding':
                            direcao = "↓"
                        elif component['direction'] == 'contracting':
                            direcao = "↑"
    
    # Determinar direção para furos cegos baseado na posição Z
    if tipo == "cego" and 'coordinates' in feat and feat['coordinates']:
        coords = feat['coordinates']
        if 'z_range' in coords:
            z_min = coords['z_range'].get('min', 0)
            z_max = coords['z_range'].get('max', altura_chapa)
            if abs(z_max - altura_chapa) < abs(z_min - 0):
                direcao = "↓"
            else:
                direcao = "↑"
    
    return components_info, direcao

def _format_circular_subgroups(features_detalhadas: List[Dict[str, Any]]) -> List[str]:
    """Formata subgrupos de features circulares."""
    output = []
    
    # Agrupar features idênticas
    subgrupos = {}
    for feat_info in features_detalhadas:
        components_str = " + ".join(feat_info['components_info']) if feat_info['components_info'] else "Simples"
        key = (feat_info['actual_min_d'], feat_info['actual_max_d'], feat_info['profundidade'], 
               feat_info['tipo'], components_str)
        
        if key not in subgrupos:
            subgrupos[key] = {'count': 0, 'direcoes': [], 'areas': [], 'perimetros': []}
        subgrupos[key]['count'] += 1
        subgrupos[key]['areas'].append(feat_info['area_real'])
        subgrupos[key]['perimetros'].append(feat_info['perimetro'])
        if feat_info['direcao']:
            subgrupos[key]['direcoes'].append(feat_info['direcao'])
    
    # Formatar cada subgrupo
    for (d_min, d_max, profundidade, tipo, components_str), info in sorted(subgrupos.items(), key=lambda x: (-x[1]['count'], x[0])):
        if d_min == d_max:
            diametro_str = f"{d_min:.1f} mm"
        else:
            diametro_str = f"{d_min:.1f} a {d_max:.1f} mm"
        
        line = f"  Quantidade: {info['count']}:\n     Diâmetro: {diametro_str} | Profundidade: {profundidade if profundidade != 'n/d' else 'n/d'} mm | Tipo: {tipo}"
        
        if components_str != "Simples":
            line += f"\n       Geometria: {components_str}"
        
        # Formatar direções
        if info['direcoes']:
            direcoes_count = {}
            for dir in info['direcoes']:
                direcoes_count[dir] = direcoes_count.get(dir, 0) + 1
            
            if len(direcoes_count) == 1:
                dir_symbol = list(direcoes_count.keys())[0]
                line += f"\n       Direção: {dir_symbol}"
            else:
                dir_info = []
                total_down = direcoes_count.get("↓", 0)
                total_up = direcoes_count.get("↑", 0)
                if total_down >= total_up:
                    if total_down > 0:
                        dir_info.append(f"{total_down}×↓")
                    if total_up > 0:
                        dir_info.append(f"{total_up}×↑")
                else:
                    if total_up > 0:
                        dir_info.append(f"{total_up}×↑")
                    if total_down > 0:
                        dir_info.append(f"{total_down}×↓")
                line += f"\n       Direção: {' + '.join(dir_info)}"
        
        output.append(line)
    
    return output

def format_semicircular_holes(summary: Dict[str, Any], shape) -> str:
    """Formata a seção de furos semicirculares."""
    from .main_details import get_semi_circular_features  # Import local para evitar circular
    
    output = []
    
    # Warnings primeiro
    warnings_text = format_warnings([w for w in summary.get('warnings', []) if w.warning_type == 'semicircular'])
    if warnings_text:
        output.append(warnings_text)
    
    output.append("Furos semicirculares:")
    
    semi_counter = summary['semi_counter']
    if semi_counter:
        bbox = summary.get('bbox')
        altura_chapa = abs(bbox[5] - bbox[2]) if bbox else 0
        semi_features = get_semi_circular_features(shape)
        
        for (min_d, max_d), _ in sorted(semi_counter.items(), key=lambda x: -x[0][1]):
            grupo_semi_feats = [feat for feat in semi_features 
                               if round(feat['min_d'], 1) == min_d and round(feat['max_d'], 1) == max_d]
            
            if not grupo_semi_feats:
                continue
            
            # Processar features semicirculares (similar aos circulares mas mais simples)
            subgrupos_text = _format_semicircular_subgroups(grupo_semi_feats, altura_chapa)
            output.extend(subgrupos_text)
    else:
        output.append("  Nenhum furo semicircular encontrado.")
    
    return "\n".join(output)

def _format_semicircular_subgroups(grupo_semi_feats: List[Dict[str, Any]], altura_chapa: float) -> List[str]:
    """Formata subgrupos de features semicirculares."""
    output = []
    
    # Processar cada feature
    features_detalhadas = []
    for feat in grupo_semi_feats:
        actual_min_d = feat['min_d']
        profundidade = altura_chapa
        tipo = "passante"
        
        # Calcular área e perímetro semicirculares
        raio = actual_min_d / 2
        area_semicircular = (math.pi * raio ** 2) / 2
        perimetro_semicircular = (math.pi * raio) + actual_min_d
        
        # Analisar componentes geométricos
        components_info = []
        if 'geometric_components' in feat and feat['geometric_components']:
            for component in feat['geometric_components']:
                if component['geometric_type'] == 'cilindrica':
                    components_info.append(f"Cilíndrico ⌀{component['diameter']:.1f}mm h{component['altura']:.1f}mm")
                elif component['geometric_type'] == 'conica':
                    components_info.append(f"Cônico ⌀{component['diameter_min']:.1f}-{component['diameter_max']:.1f}mm h{component['altura']:.1f}mm")
        
        features_detalhadas.append({
            'actual_min_d': actual_min_d,
            'profundidade': profundidade,
            'tipo': tipo,
            'perimetro': perimetro_semicircular,
            'area_real': area_semicircular,
            'components_info': components_info
        })
    
    # Agrupar e formatar
    subgrupos = {}
    for feat_info in features_detalhadas:
        components_str = " + ".join(feat_info['components_info']) if feat_info['components_info'] else "Semicircular simples"
        key = (feat_info['actual_min_d'], feat_info['profundidade'], 
               feat_info['tipo'], feat_info['perimetro'], feat_info['area_real'], components_str)
        
        if key not in subgrupos:
            subgrupos[key] = {'count': 0}
        subgrupos[key]['count'] += 1
    
    # Formatar cada subgrupo
    for (d_min, profundidade, tipo, perimetro, area_real, components_str), info in sorted(subgrupos.items(), key=lambda x: (-x[1]['count'], x[0])):
        diametro_str = f"{d_min:.1f} mm"
        
        line = f"  Quantidade: {info['count']}:"
        line += f"\n     Diâmetro: {diametro_str} | Profundidade: {profundidade:.1f} mm | Tipo: {tipo}"
        line += f"\n       Geometria: {components_str}"
        
        output.append(line)
    
    return output

def format_rectangular_holes(summary: Dict[str, Any], shape) -> str:
    """Formata a seção de furos retangulares."""
    output = []
    
    # Warnings primeiro
    warnings_text = format_warnings([w for w in summary.get('warnings', []) if w.warning_type == 'rectangular'])
    if warnings_text:
        output.append(warnings_text)
    
    output.append("Furos Retangulares:")
    
    rectangular_counter = summary['rectangular_counter']
    if rectangular_counter:
        bbox = summary.get('bbox')
        altura = abs(bbox[5] - bbox[2]) if bbox else 0
        
        for i, grupo in enumerate(rectangular_counter, 1):
            rect_text = _format_single_rectangular_hole(grupo, i, shape, altura)
            output.append(rect_text)
    else:
        output.append("  Nenhum furo retangular encontrado.")
    
    return "\n".join(output)

def _format_single_rectangular_hole(grupo: Dict[str, Any], index: int, shape, altura_chapa: float) -> str:
    """Formata um único furo retangular."""
    comprimento = grupo['comprimento']
    largura = grupo['largura']
    
    # Calcular profundidade e tipo (lógica simplificada)
    profundidade = altura_chapa
    tipo = "passante"
    direcao = ""
    
    # Análise básica de tipo baseado nas faces
    if 'faces' in grupo and grupo['faces']:
        # Lógica simplificada - pode ser expandida conforme necessário
        faces_group = grupo['faces']
        if isinstance(faces_group, list) and len(faces_group) > 0:
            primeiro_item = faces_group[0]
            if isinstance(primeiro_item, dict) and 'height' in primeiro_item:
                # Calcular com base nos metadados
                alturas_elementos = [item['height'] for item in faces_group if item.get('height', 0) > 0.1]
                if alturas_elementos:
                    altura_maxima = max(alturas_elementos)
                    if altura_maxima < altura_chapa * 0.95:
                        tipo = "cego"
                        profundidade = round(altura_maxima, 1)
                        direcao = " ↓"  # Assumir direção padrão
    
    # Formatar saída principal
    tipo_completo = f"{tipo}{direcao}" if direcao else tipo
    output = f"  Furo {index}: {comprimento:.1f} x {largura:.1f} mm - Profundidade: {profundidade:.1f} mm - Tipo: {tipo_completo}"
    
    # Adicionar detalhes das faces se disponível
    faces_details = _format_rectangular_face_details(grupo)
    if faces_details:
        output += "\n" + faces_details
    
    return output

def _format_rectangular_face_details(grupo: Dict[str, Any]) -> str:
    """Formata detalhes das faces de um furo retangular."""
    if 'faces' not in grupo or not grupo['faces']:
        return ""
    
    # Obter limites do furo retangular
    hole_bbox = grupo.get('bbox', None)
    
    output = []
    output.append("       Detalhes das faces:")
    
    faces_group = grupo['faces']
    
    if isinstance(faces_group[0], dict):
        # Metadados das faces
        faces_agrupadas = {}
        for idx, face_data in enumerate(faces_group):
            length = face_data.get('length', 0)
            width = face_data.get('width', 0)
            height = face_data.get('height', 0)
            face_type = face_data.get('type', 'desconhecida')
            
            # Calcular ângulo usando novas regras
            angle = calculate_face_angle_from_dimensions(length, width, height, face_type)
            
            # Normalizar tipo para exibição
            if face_type in ['cone', 'cilindro']:
                display_type = face_type
            elif face_type == 'plano':
                display_type = "plana"
            else:
                display_type = face_type
            
            # Ordenar dimensões
            dim1 = max(length, width)
            dim2 = min(length, width)
            
            # Chave de agrupamento
            key = (round(dim1, 1), round(dim2, 1), round(height, 1), round(angle, 1), display_type)
            
            if key not in faces_agrupadas:
                faces_agrupadas[key] = []
            
            # Criar entrada da face com todos os campos necessários
            face_entry = {
                'index': idx + 1,
                'length': length,
                'width': width,
                'height': height,
                'angle': angle,
                'type': display_type
            }
            
            # Copiar campos de diâmetro se existirem
            if 'diameter' in face_data:
                face_entry['diameter'] = face_data['diameter']
            if 'ref_diameter' in face_data:
                face_entry['ref_diameter'] = face_data['ref_diameter']
            
            # Copiar centro se existir
            if 'center' in face_data:
                face_entry['center'] = face_data['center']
                
            faces_agrupadas[key].append(face_entry)
        
        # Organizar por orientação
        faces_por_orientacao = _group_faces_by_orientation(faces_agrupadas, hole_bbox)
        
        # Formatar saída
        ordem_orientacao = ['horizontal', 'diagonal', 'vertical']
        for orientacao in ordem_orientacao:
            if orientacao in faces_por_orientacao:
                output.append(f"       Faces {orientacao}s:")
                for face_info in faces_por_orientacao[orientacao]:
                    output.append(f"         {face_info['descricao']}")
                output.append("")  # linha em branco
    
    return "\n".join(output) if len(output) > 1 else ""

def _group_faces_by_orientation(faces_agrupadas: Dict, hole_bbox=None) -> Dict[str, List[Dict[str, Any]]]:
    """Agrupa faces por orientação (horizontal, diagonal, vertical)."""
    
    def get_position_label(face_center, hole_bbox):
        """Determina se o centro da face está 'In' (dentro) ou 'Out' (fora) do furo."""
        if not hole_bbox or not face_center:
            return ""
        
        xmin, ymin, zmin, xmax, ymax, zmax = hole_bbox
        fx, fy, fz = face_center
        
        # Verificar se o centro está dentro dos limites do furo (com uma pequena tolerância)
        tol = 1.0  # tolerância em mm
        if (xmin - tol <= fx <= xmax + tol and 
            ymin - tol <= fy <= ymax + tol):
            return " (In)"
        else:
            return " (Out)"
    
    faces_por_orientacao = {}
    
    for (dim1, dim2, altura_face, angulo, tipo), faces_do_grupo in faces_agrupadas.items():
        angle_word = angle_to_word(angulo)
        
        if angle_word not in faces_por_orientacao:
            faces_por_orientacao[angle_word] = []
        
        if len(faces_do_grupo) == 1:
            face = faces_do_grupo[0]
            descricao = f"Face {face['index']}: {face['length']:.1f} x {face['width']:.1f} x {face['height']:.1f} mm ({face['type']})"
            # Adicionar informações de diâmetro para faces cilíndricas
            if face['type'] == 'cilindro' and 'diameter' in face:
                descricao += f" - Diâmetro: {face['diameter']:.2f} mm"
            # Adicionar informações de diâmetro para faces cônicas
            if face['type'] == 'cone' and 'ref_diameter' in face:
                descricao += f" - Diâmetro base: {face['ref_diameter']:.2f} mm"
            # Adicionar label de posição para cilindros e cones
            if face['type'] in ['cilindro', 'cone'] and 'center' in face:
                # Lógica estrita: centro deve estar estritamente dentro dos limites
                xmin, ymin, zmin, xmax, ymax, zmax = hole_bbox if hole_bbox else (None, None, None, None, None, None)
                fx, fy, fz = face['center']
                if hole_bbox and (xmin < fx < xmax and ymin < fy < ymax):
                    descricao += " (In)"
                else:
                    descricao += " (Out)"
            faces_por_orientacao[angle_word].append({
                'descricao': descricao,
                'count': 1
            })
        else:
            indices = [str(face['index']) for face in faces_do_grupo]
            descricao = f"Faces {', '.join(indices)}: {dim1:.1f} x {dim2:.1f} x {altura_face:.1f} mm ({tipo}) (x{len(faces_do_grupo)})"
            # Para múltiplas faces, verificar se são cilíndricas e adicionar info de diâmetro
            primeira_face = faces_do_grupo[0]
            if tipo == 'cilindro' and 'diameter' in primeira_face:
                descricao += f" - Diâmetro: {primeira_face['diameter']:.2f} mm"
            if tipo == 'cone' and 'ref_diameter' in primeira_face:
                descricao += f" - Diâmetro base: {primeira_face['ref_diameter']:.2f} mm"
            # Adicionar label de posição se disponível
            if tipo in ['cilindro', 'cone'] and 'center' in primeira_face:
                xmin, ymin, zmin, xmax, ymax, zmax = hole_bbox if hole_bbox else (None, None, None, None, None, None)
                fx, fy, fz = primeira_face['center']
                if hole_bbox and (xmin < fx < xmax and ymin < fy < ymax):
                    descricao += " (In)"
                else:
                    descricao += " (Out)"
            faces_por_orientacao[angle_word].append({
                'descricao': descricao,
                'count': len(faces_do_grupo)
            })
    
    return faces_por_orientacao

def format_piece_position(todas_direcoes_globais: List[str]) -> str:
    """Calcula e formata a posição da chapa baseado nas direções dos furos."""
    direcoes_unicas = set(todas_direcoes_globais)
    if len(direcoes_unicas) <= 1:
        posicao_chapa = 1  # Mesma direção ou sem direção
    else:
        posicao_chapa = 2  # Direções diferentes
    
    return f"\nPosição da chapa: {posicao_chapa}"

def generate_full_report(summary: Dict[str, Any], shape) -> str:
    """
    Gera o relatório completo da peça.
    Esta é a função principal que substitui a lógica de formatação do show_general_summary.
    """
    output_parts = []
    
    # Cabeçalho
    header = format_piece_header(summary)
    output_parts.append(header)
    
    # Furos circulares
    circular_text, todas_direcoes_globais = format_circular_holes(summary, shape, summary.get('hole_groups', []))
    output_parts.append(circular_text)
    
    # Posição da chapa baseada nas direções
    position_text = format_piece_position(todas_direcoes_globais)
    output_parts.append(position_text)
    output_parts.append("\n" + "-" * 40)
    
    # Furos semicirculares
    semicircular_text = format_semicircular_holes(summary, shape)
    output_parts.append(semicircular_text)
    output_parts.append("\n" + "-" * 40)
    
    # Furos retangulares
    rectangular_text = format_rectangular_holes(summary, shape)
    output_parts.append(rectangular_text)
    
    # Rodapé
    output_parts.append("\n" + "=" * 40)
    
    return "\n".join(output_parts)
