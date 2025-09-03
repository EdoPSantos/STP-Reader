"""
Utilitários geométricos consolidados - VERSÃO UNIFICADA
Consolida utils.py + geometry_utils.py para eliminar duplicações.

Este ficheiro contém TODAS as funções geométricas numa só localização,
eliminando duplicações de constantes, importações e funções.
"""

import math
import os
import re
import logging
from typing import List, Dict, Optional, Tuple, Any
from collections import Counter, defaultdict

# Importações OpenCASCADE unificadas
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone, GeomAbs_Plane, GeomAbs_Sphere, GeomAbs_Torus, GeomAbs_Circle
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface, GeomAdaptor_Curve
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_FACE, TopAbs_EDGE
from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.Bnd import Bnd_Box
import OCC.Core.TopExp

logger = logging.getLogger(__name__)
logger.disabled = True

# ============================================================================
# CONSTANTES UNIFICADAS - Eliminação total de duplicações
# ============================================================================
# Tolerâncias para análise geométrica
TOL_BBOX = 1.0               # Tolerância para bounding box (main_details.py tinha 1.0)
TOL_IS_HOLE = 2.0            # Tolerância para is_hole_through
TOL_CENTER = 1.0             # Tolerância para centro (main_details.py tinha 1.0)
TOL_D = 0.5                  # Tolerância para diâmetro (agrupamento de features)
TOL_FACE_BORDER = 0.5        # Tolerância para faces na borda
TOL_PLANAR_HEIGHT = 1.0      # Tolerância para altura mínima de face plana
TOL_CONEXAO = 0.5            # Tolerância para ligação de faces

# Tolerâncias específicas para semicirculares
GROUP_TOL = 3.0              # Tolerância para agrupamento de semicirculares
LATERAL_TOL = 150.0          # Tolerância para lateralidade de semicirculares
ANGLE_TOL = 30.0             # Tolerância angular para arcos semicirculares
RADIUS_TOL = 1.0             # Tolerância para raio em arcos semicirculares

# Strings de tipo geométrico unificadas
TIPO_PLANO = 'plano'
TIPO_CILINDRO = 'cilindro'
TIPO_CONE = 'cone'
TIPO_ESFERA = 'esfera'
TIPO_TORO = 'toro'

# ============================================================================
# FUNÇÕES DE EXTRAÇÃO DE NOMES
# ============================================================================

def extract_mold_and_part_from_step(filename):
    base = os.path.splitext(os.path.basename(filename))[0]
    parts = base.split("_")
    if len(parts) < 2:
        raise ValueError("Nome do ficheiro não contém '_' suficiente para extrair Molde e Peça")
    molde = parts[-2]
    peca = parts[-1]
    return molde, peca

def extract_from_filename(filepath):
    """Extrai nome do molde e peça do nome do ficheiro."""
    base = os.path.basename(filepath)
    name, _ = os.path.splitext(base)
    parts = re.split(r"[_\-]", name)
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None, None

# ============================================================================
# FUNÇÕES DE ÂNGULOS E ORIENTAÇÃO
# ============================================================================

def calculate_face_angle_from_dimensions(length: float, width: float, height: float, face_type: str) -> float:
    """
    Calcula o ângulo da face baseado nas dimensões e tipo conforme as regras corretas:
    - 0°: Para qualquer face com altura (Z) = 0
    - 45°: Para faces cônicas (sempre) ou planas com X e Y > 0 e altura > 0  
    - 90°: Para cilindros com altura > 0 ou planas com apenas uma dimensão (X ou Y) e altura > 0
    """
    TOL_ZERO = 0.01
    
    length = length if length > TOL_ZERO else 0
    width = width if width > TOL_ZERO else 0  
    height = height if height > TOL_ZERO else 0
    
    if height == 0:
        return 0.0
    
    if face_type == 'cone':
        return 45.0
    elif face_type == 'cilindro':
        return 90.0
    elif face_type in ['plana', 'plano']:
        has_length = length > 0
        has_width = width > 0
        
        if has_length and has_width:
            return 45.0
        elif has_length or has_width:
            return 90.0
        else:
            return 90.0
    
    return 0.0

def angle_to_word(angle: float) -> str:
    """Converte ângulo numérico para palavra correspondente."""
    if abs(angle - 90.0) < 0.1:
        return "vertical"
    elif abs(angle - 45.0) < 0.1:
        return "diagonal"
    elif abs(angle - 0.0) < 0.1:
        return "horizontal"
    else:
        return f"{angle:.1f}°"

def get_real_face_angle(face, face_type: str) -> float:
    """
    FUNÇÃO NÃO UTILIZADA - Mantida para compatibilidade.
    Obter ângulo real baseado no tipo de face.
    """
    try:
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        surface_type = adaptor.GetType()
        
        if surface_type == GeomAbs_Cone:
            cone = adaptor.Cone()
            semi_angle_rad = cone.SemiAngle()
            semi_angle_deg = math.degrees(semi_angle_rad)
            inclination_angle = 90.0 - semi_angle_deg
            return max(0.0, min(90.0, abs(inclination_angle)))
            
        elif surface_type == GeomAbs_Cylinder:
            cylinder = adaptor.Cylinder()
            axis = cylinder.Axis().Direction()
            
            vertical_vec = [0, 0, 1]
            axis_vec = [axis.X(), axis.Y(), axis.Z()]
            
            dot_product = sum(a * b for a, b in zip(axis_vec, vertical_vec))
            angle_rad = math.acos(abs(dot_product))
            angle_deg = math.degrees(angle_rad)
            
            return round(angle_deg, 1)
            
        elif surface_type == GeomAbs_Plane:
            plane = adaptor.Plane()
            normal = plane.Axis().Direction()
            normal_vec = [normal.X(), normal.Y(), normal.Z()]
            
            vertical_vec = [0, 0, 1]
            dot_product = sum(a * b for a, b in zip(normal_vec, vertical_vec))
            
            angle_rad = math.acos(abs(dot_product))
            angle_deg = math.degrees(angle_rad)
            
            return round(angle_deg, 1)
            
        else:
            return calculate_face_angle_and_type(face)['angle']
            
    except Exception as e:
        return 0.0

def calculate_face_angle_and_type(face) -> Dict[str, Any]:
    """
    FUNÇÃO NÃO UTILIZADA - Mantida para compatibilidade.
    Calcula o ângulo de inclinação e tipo de uma face.
    """
    try:
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        surface_type = adaptor.GetType()
        
        is_curved = surface_type not in [GeomAbs_Plane]
        face_type = "curva" if is_curved else "plana"
        
        if surface_type == GeomAbs_Plane:
            plane = adaptor.Plane()
            normal = plane.Axis().Direction()
            normal_vec = [normal.X(), normal.Y(), normal.Z()]
            
        else:
            try:
                u_min, u_max, v_min, v_max = surface.Bounds()
                u_mid = (u_min + u_max) / 2
                v_mid = (v_min + v_max) / 2
                
                point = gp_Pnt()
                du = gp_Vec()
                dv = gp_Vec()
                surface.D1(u_mid, v_mid, point, du, dv)
                
                normal_vec = [
                    du.Y() * dv.Z() - du.Z() * dv.Y(),
                    du.Z() * dv.X() - du.X() * dv.Z(), 
                    du.X() * dv.Y() - du.Y() * dv.X()
                ]
                
                mag = (sum(x**2 for x in normal_vec))**0.5
                if mag > 0:
                    normal_vec = [x/mag for x in normal_vec]
                else:
                    normal_vec = [0, 0, 1]
            except Exception as e:
                normal_vec = [0, 0, 1]
        
        vertical_vec = [0, 0, 1]
        dot_product = sum(a * b for a, b in zip(normal_vec, vertical_vec))
        
        dot_product = max(-1.0, min(1.0, dot_product))
        
        angle_rad = math.acos(abs(dot_product))
        angle_deg = math.degrees(angle_rad)
        
        return {
            "angle": round(angle_deg, 1),
            "type": face_type
        }
        
    except Exception as e:
        return {"angle": 0.0, "type": "desconhecida"}

# ============================================================================
# FUNÇÕES DE BOUNDING BOX
# ============================================================================

def get_bbox(shape_or_face):
    """Retorna a bounding box de uma shape ou face: (xmin, ymin, zmin, xmax, ymax, zmax)"""
    bbox = Bnd_Box()
    brepbndlib.Add(shape_or_face, bbox)
    return bbox.Get()

def is_point_inside_bbox(point, bbox, tol=0.5):
    """Verifica se um ponto (x, y) está dentro do bounding box."""
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    x, y = point[:2]
    return (xmin - tol) <= x <= (xmax + tol) and (ymin - tol) <= y <= (ymax + tol)

def get_bbox_faces_touched(pnt, bbox, tol=0.5):
    """Identifica quais faces do bounding box um ponto toca."""
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    faces = []
    if abs(pnt.X() - xmin) < TOL_BBOX: faces.append("xmin")
    if abs(pnt.X() - xmax) < TOL_BBOX: faces.append("xmax")
    if abs(pnt.Y() - ymin) < TOL_BBOX: faces.append("ymin")
    if abs(pnt.Y() - ymax) < TOL_BBOX: faces.append("ymax")
    if abs(pnt.Z() - zmin) < TOL_BBOX: faces.append("zmin")
    if abs(pnt.Z() - zmax) < TOL_BBOX: faces.append("zmax")
    return faces

# ============================================================================
# FUNÇÕES DE AGRUPAMENTO E CHAVES
# ============================================================================

def feat_key(feat, center_tol=1.0, d_tol=0.5):
    """Gera uma chave para agrupamento de features circulares/semicirculares."""
    return (round(feat['center'][0]/TOL_CENTER), round(feat['center'][1]/TOL_CENTER), round(feat['max_d']/TOL_D))

def normalize_grouping_values(center: Tuple[float, float, float], diameter: float, 
                            center_tol: float = 15.0, d_tol: float = 0.5) -> Tuple[int, int, int]:
    """Normaliza valores para agrupamento consistente."""
    return (
        round(center[0] / center_tol),
        round(center[1] / center_tol),
        round(diameter / d_tol)
    )

# ============================================================================
# FUNÇÕES DE ÁREA E GEOMETRIA
# ============================================================================

def get_face_area(face):
    """Retorna a área de uma face."""
    props = GProp_GProps()
    brepgprop.SurfaceProperties(face, props)
    return props.Mass()

def has_real_void_in_center(faces, center_point: Tuple[float, float, float], 
                           min_radius: float, depth_info: Optional[float] = None) -> bool:
    """Verifica se há um vazio real no centro do suposto furo."""
    try:
        if depth_info and isinstance(depth_info, (int, float)):
            if depth_info < 1.0:
                return False
        
        if len(faces) < 3:
            if (len(faces) == 1 or len(faces) == 2) and min_radius >= 2.0:
                if min_radius <= 10.0:
                    return True
                else:
                    return True
            return False
            
        face_z_ranges = []
        face_types = []
        
        for face_info in faces:
            if len(face_info) >= 3:
                face_data = face_info[2]
                if 'z_range' in face_data:
                    face_z_ranges.append((face_data['z_range']['min'], face_data['z_range']['max']))
                    face_types.append(face_data.get('surface_type', 'unknown'))
        
        if len(face_z_ranges) >= 2:
            face_z_ranges.sort()
            
            for i in range(len(face_z_ranges) - 1):
                z_max_atual = face_z_ranges[i][1]
                z_min_proximo = face_z_ranges[i + 1][0]
                gap = z_min_proximo - z_max_atual
                
                if gap > 0.5:
                    return False
            
            return True
        
        return len(faces) >= 3
        
    except Exception as e:
        return True

# ============================================================================
# FUNÇÕES DE ORIENTAÇÃO E VALIDAÇÃO
# ============================================================================

def check_face_orientation_for_hole(face, adaptor, group_center):
    """
    Verifica se uma face curva está orientada corretamente para ser parte de um furo.
    Versão mais conservativa - só rejeita casos muito claros de faces "de costas".
    """
    try:
        orientation = face.Orientation()
        
        bbox = get_bbox(face)
        face_center_x = (bbox[0] + bbox[3]) / 2
        face_center_y = (bbox[1] + bbox[4]) / 2
        face_center_z = (bbox[2] + bbox[5]) / 2
        
        dist_to_group = (
            (face_center_x - group_center[0])**2 + 
            (face_center_y - group_center[1])**2 + 
            (face_center_z - group_center[2])**2
        )**0.5
        
        if dist_to_group > 50.0:
            return False
            
        return True
            
    except Exception as e:
        return True
    
    return True

# ============================================================================
# FUNÇÕES PRINCIPAIS DE ANÁLISE CIRCULAR
# ============================================================================

def group_faces_by_axis_and_proximity(shape, loc_tol=2.0, dir_tol=0.05):
    """
    NOVA ABORDAGEM: Recolhe todas as faces curvas, verifica limites da chapa,
    agrupa similar aos retângulos, depois analisa geometria para identificar cilindros/cones.
    """
    # 1. RECOLHER todas as faces curvas da shape
    all_curved_faces = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        stype = adaptor.GetType()
        
        if stype in [GeomAbs_Cylinder, GeomAbs_Cone]:
            bbox = get_bbox(face)
            xmin, ymin, zmin, xmax, ymax, zmax = bbox
            center_x = (xmin + xmax) / 2
            center_y = (ymin + ymax) / 2
            center_z = (zmin + zmax) / 2
            
            face_data = {
                'face': face,
                'type': stype,
                'adaptor': adaptor,
                'bbox': bbox,
                'center': (round(center_x, 2), round(center_y, 2), round(center_z, 2)),
                'area': get_face_area(face)
            }
            
            if stype == GeomAbs_Cylinder:
                cylinder = adaptor.Cylinder()
                axis = cylinder.Axis()
                face_data.update({
                    'radius': cylinder.Radius(),
                    'axis_location': (axis.Location().X(), axis.Location().Y(), axis.Location().Z()),
                    'axis_direction': (axis.Direction().X(), axis.Direction().Y(), axis.Direction().Z())
                })
            elif stype == GeomAbs_Cone:
                cone = adaptor.Cone()
                axis = cone.Axis()
                face_data.update({
                    'ref_radius': cone.RefRadius(),
                    'semi_angle': cone.SemiAngle(),
                    'axis_location': (axis.Location().X(), axis.Location().Y(), axis.Location().Z()),
                    'axis_direction': (axis.Direction().X(), axis.Direction().Y(), axis.Direction().Z())
                })
            
            all_curved_faces.append(face_data)
        
        explorer.Next()
    
    # 2. VERIFICAR LIMITES DA CHAPA
    chapa_bbox = get_bbox(shape)
    xmin_c, ymin_c, zmin_c, xmax_c, ymax_c, zmax_c = chapa_bbox
    
    largura_chapa = abs(xmax_c - xmin_c)
    altura_chapa = abs(ymax_c - ymin_c)
    espessura_chapa = abs(zmax_c - zmin_c)
    
    tol_limite_x = min(10.0, largura_chapa * 0.02)
    tol_limite_y = min(10.0, altura_chapa * 0.02)
    tol_limite_z = min(5.0, espessura_chapa * 0.5)
    
    faces_dentro_limites = []
    for face_data in all_curved_faces:
        cx, cy, cz = face_data['center']
        face_bbox = face_data['bbox']
        
        face_xmin, face_ymin, face_zmin, face_xmax, face_ymax, face_zmax = face_bbox
        
        if (face_xmin >= xmin_c - tol_limite_x and face_xmax <= xmax_c + tol_limite_x and
            face_ymin >= ymin_c - tol_limite_y and face_ymax <= ymax_c + tol_limite_y and
            face_zmin >= zmin_c - tol_limite_z and face_zmax <= zmax_c + tol_limite_z):
            faces_dentro_limites.append(face_data)
    
    # 3. AGRUPAR faces
    n = len(faces_dentro_limites)
    conexoes = [[] for _ in range(n)]
    
    tol_conexao = 3.0
    
    for i in range(n):
        face1 = faces_dentro_limites[i]
        bbox1 = face1['bbox']
        
        for j in range(i + 1, n):
            face2 = faces_dentro_limites[j]
            bbox2 = face2['bbox']
            
            dx = abs(face1['center'][0] - face2['center'][0])
            dy = abs(face1['center'][1] - face2['center'][1])
            
            if dx > loc_tol * 2 or dy > loc_tol * 2:
                continue
            
            x_overlap = not (bbox1[3] < bbox2[0] - tol_conexao or bbox2[3] < bbox1[0] - tol_conexao)
            y_overlap = not (bbox1[4] < bbox2[1] - tol_conexao or bbox2[4] < bbox1[1] - tol_conexao)
            z_overlap = not (bbox1[5] < bbox2[2] - tol_conexao or bbox2[5] < bbox1[2] - tol_conexao)
            
            overlaps = sum([x_overlap, y_overlap, z_overlap])
            if overlaps >= 2:
                conexoes[i].append(j)
                conexoes[j].append(i)
    
    # 4. DFS para agrupar
    groups = []
    used = [False] * n

    def dfs(idx, grupo):
        used[idx] = True
        grupo.append(faces_dentro_limites[idx])
        for viz in conexoes[idx]:
            if not used[viz]:
                dfs(viz, grupo)

    for i in range(n):
        if not used[i]:
            grupo = []
            dfs(i, grupo)
            if grupo:
                groups.append(grupo)
    
    # 5. ANALISAR GEOMETRIA
    analyzed_groups = []
    
    for i, grupo in enumerate(groups):
        grupo.sort(key=lambda x: x['area'], reverse=True)
        
        if grupo:
            avg_x = sum(f['center'][0] for f in grupo) / len(grupo)
            avg_y = sum(f['center'][1] for f in grupo) / len(grupo)
            avg_z = sum(f['center'][2] for f in grupo) / len(grupo)
            group_center = (avg_x, avg_y, avg_z)
            
            faces_corretas = []
            for face_data in grupo:
                if check_face_orientation_for_hole(face_data['face'], face_data['adaptor'], group_center):
                    faces_corretas.append(face_data)
            
            if faces_corretas:
                geometric_analysis = analyze_group_geometry(faces_corretas)
                
                group_with_analysis = []
                for face_data in faces_corretas:
                    group_with_analysis.append((
                        face_data['type'],
                        face_data['face'], 
                        face_data['adaptor'],
                        face_data,
                        geometric_analysis
                    ))
                
                analyzed_groups.append(group_with_analysis)
    
    return analyzed_groups

def analyze_group_geometry(grupo):
    """
    Analisa um grupo de faces para determinar a geometria real:
    - Se o raio é constante ao longo de Z = CILINDRO
    - Se o raio varia ao longo de Z = CONE
    """
    if not grupo:
        return {'type': 'unknown', 'components': []}
    
    face_info = []
    for face_data in grupo:
        bbox = face_data['bbox']
        zmin, zmax = bbox[2], bbox[5]
        z_center = (zmin + zmax) / 2
        
        if face_data['type'] == GeomAbs_Cylinder:
            radius = face_data['radius']
            face_info.append({
                'z_min': zmin,
                'z_max': zmax,
                'z_center': z_center,
                'radius': radius,
                'type': 'cylinder',
                'face_data': face_data
            })
        elif face_data['type'] == GeomAbs_Cone:
            cone = face_data['adaptor'].Cone()
            ref_radius = cone.RefRadius()
            semi_angle = cone.SemiAngle()
            apex = cone.Apex()
            axis = cone.Axis()
            
            height = abs(zmax - zmin)
            
            apex_z = apex.Z()
            
            dist_to_zmin = abs(zmin - apex_z)
            dist_to_zmax = abs(zmax - apex_z)
            
            r_at_zmin = dist_to_zmin * math.tan(semi_angle)
            r_at_zmax = dist_to_zmax * math.tan(semi_angle)
            
            if r_at_zmin < 0:
                r_at_zmin = abs(r_at_zmin)
            if r_at_zmax < 0:
                r_at_zmax = abs(r_at_zmax)
            
            apex_z = face_data.get('apex_z', z_center)
            z_mid = (zmin + zmax) / 2
            
            if abs(r_at_zmax - r_at_zmin) < 0.01:
                if apex_z < z_mid:
                    direction = 'expanding'
                else:
                    direction = 'contracting'
            else:
                direction = 'expanding' if r_at_zmax > r_at_zmin else 'contracting'
            
            face_info.append({
                'z_min': zmin,
                'z_max': zmax,
                'z_center': z_center,
                'radius_min': min(r_at_zmin, r_at_zmax),
                'radius_max': max(r_at_zmin, r_at_zmax),
                'radius_at_zmin': r_at_zmin,
                'radius_at_zmax': r_at_zmax,
                'type': 'cone',
                'direction': direction,
                'face_data': face_data
            })
    
    face_info.sort(key=lambda x: x['z_center'])
    
    components = []
    
    for info in face_info:
        if info['type'] == 'cylinder':
            components.append({
                'geometric_type': 'cilindrica',
                'z_min': round(info['z_min'], 1),
                'z_max': round(info['z_max'], 1),
                'radius': round(info['radius'], 2),
                'diameter': round(info['radius'] * 2, 2),
                'altura': round(abs(info['z_max'] - info['z_min']), 2),
                'area_lateral': round(2 * math.pi * info['radius'] * abs(info['z_max'] - info['z_min']), 2),
                'volume': round(math.pi * info['radius']**2 * abs(info['z_max'] - info['z_min']), 2),
                'face_data': info['face_data']
            })
        elif info['type'] == 'cone':
            components.append({
                'geometric_type': 'conica',
                'z_min': round(info['z_min'], 1),
                'z_max': round(info['z_max'], 1),
                'radius_min': round(info['radius_min'], 2),
                'radius_max': round(info['radius_max'], 2),
                'diameter_min': round(info['radius_min'] * 2, 2),
                'diameter_max': round(info['radius_max'] * 2, 2),
                'altura': round(abs(info['z_max'] - info['z_min']), 2),
                'direction': info['direction'],
                'apex_angle': round(math.degrees(info['face_data']['semi_angle'] * 2), 1),
                'conicidade': round(math.degrees(info['face_data']['semi_angle']), 1),
                'face_data': info['face_data']
            })
    
    cylinder_count = sum(1 for c in components if c['geometric_type'] == 'cilindrica')
    cone_count = sum(1 for c in components if c['geometric_type'] == 'conica')
    
    if cylinder_count > 0 and cone_count == 0:
        hole_type = 'cylindrical'
    elif cone_count > 0 and cylinder_count == 0:
        hole_type = 'conical'
    elif cylinder_count > 0 and cone_count > 0:
        hole_type = 'mixed'
    else:
        hole_type = 'unknown'
    
    total_area = calculate_hole_total_area(components)
    
    if components:
        avg_x = sum(c.get('face_data', {}).get('center', [0, 0, 0])[0] for c in components) / len(components)
        avg_y = sum(c.get('face_data', {}).get('center', [0, 0, 0])[1] for c in components) / len(components)
        avg_z = sum(c.get('face_data', {}).get('center', [0, 0, 0])[2] for c in components) / len(components)
        group_center = (round(avg_x, 2), round(avg_y, 2), round(avg_z, 2))
    else:
        group_center = (0.0, 0.0, 0.0)
    
    coordinates = extract_hole_coordinates(components, group_center)
    
    return {
        'type': hole_type,
        'components': components,
        'total_components': len(components),
        'cylinder_count': cylinder_count,
        'cone_count': cone_count,
        'total_area': total_area,
        'coordinates': coordinates
    }

def calculate_hole_total_area(geometric_components):
    """Calcula a área total de um furo baseado nos seus componentes geométricos."""
    if not geometric_components:
        return 0.0
    
    total_area = 0.0
    
    for component in geometric_components:
        if component['geometric_type'] == 'cilindrica':
            r = component['radius']
            area = math.pi * r * r
            total_area += area
            
        elif component['geometric_type'] == 'conica':
            r_max = component['radius_max']
            area = math.pi * r_max * r_max
            total_area += area
    
    return round(total_area, 2)

def extract_hole_coordinates(geometric_components, center):
    """Extrai todas as coordenadas relevantes de um furo."""
    coordinates = {
        'center_principal': center,
        'components_coords': [],
        'z_range': {'min': None, 'max': None},
        'bounding_coords': {}
    }
    
    if not geometric_components:
        return coordinates
    
    z_values = []
    
    for i, component in enumerate(geometric_components):
        comp_coords = {
            'component_id': i + 1,
            'geometric_type': component['geometric_type'],
            'z_min': component.get('z_min', 0),
            'z_max': component.get('z_max', 0)
        }
        
        if component['geometric_type'] == 'cilindrica':
            comp_coords.update({
                'radius': component['radius'],
                'center_z': (component.get('z_min', 0) + component.get('z_max', 0)) / 2
            })
        elif component['geometric_type'] == 'conica':
            comp_coords.update({
                'radius_min': component['radius_min'],
                'radius_max': component['radius_max'],
                'center_z': (component.get('z_min', 0) + component.get('z_max', 0)) / 2
            })
        
        coordinates['components_coords'].append(comp_coords)
        z_values.extend([component.get('z_min', 0), component.get('z_max', 0)])
    
    if z_values:
        coordinates['z_range']['min'] = round(min(z_values), 2)
        coordinates['z_range']['max'] = round(max(z_values), 2)
        coordinates['z_range']['height'] = round(abs(max(z_values) - min(z_values)), 2)
    
    if geometric_components:
        max_radius = 0
        for component in geometric_components:
            if component['geometric_type'] == 'cilindrica':
                max_radius = max(max_radius, component['radius'])
            elif component['geometric_type'] == 'conica':
                max_radius = max(max_radius, component['radius_max'])
        
        if max_radius > 0:
            cx, cy, cz = center
            coordinates['bounding_coords'] = {
                'x_min': round(cx - max_radius, 2),
                'x_max': round(cx + max_radius, 2),
                'y_min': round(cy - max_radius, 2),
                'y_max': round(cy + max_radius, 2),
                'z_min': coordinates['z_range']['min'],
                'z_max': coordinates['z_range']['max'],
                'diameter': round(max_radius * 2, 2)
            }
    
    return coordinates

def get_all_radii_of_group(faces, min_allowed=0.0):
    """
    Retorna todos os raios (em mm) e detalhes geométricos de todas as faces do grupo.
    NOVA VERSÃO: Trabalha com análise geométrica completa.
    """
    radii = []
    has_cone = False
    
    if len(faces) > 0 and len(faces[0]) == 5:
        geometric_analysis = faces[0][4]
        
        geometric_components = geometric_analysis.get('components', [])
        
        for component in geometric_components:
            if component['geometric_type'] == 'cilindrica':
                r = component['radius']
                if r >= min_allowed:
                    radii.append(r)
            elif component['geometric_type'] == 'conica':
                has_cone = True
                r_min = component['radius_min']
                r_max = component['radius_max']
                if r_min >= min_allowed:
                    radii.append(r_min)
                if r_max >= min_allowed:
                    radii.append(r_max)
        
        radii = sorted(set(round(r, 2) for r in radii))
        return radii, has_cone, geometric_components
    
    else:
        geometric_components = []
        
        for item in faces:
            if len(item) >= 3:
                stype, face, adaptor = item[:3]
                
                component_info = {
                    'type': stype,
                    'face': face,
                    'adaptor': adaptor
                }

                if stype == GeomAbs_Cylinder:
                    r = adaptor.Cylinder().Radius()
                    if r >= min_allowed:
                        radii.append(r)
                    
                    bbox = get_bbox(face)
                    altura = abs(bbox[5] - bbox[2])
                    
                    component_info.update({
                        'geometric_type': 'cilindrica',
                        'radius': round(r, 2),
                        'diameter': round(r * 2, 2),
                        'altura': round(altura, 2),
                        'area_lateral': round(2 * 3.14159 * r * altura, 2),
                        'volume': round(3.14159 * r * r * altura, 2)
                    })
                    
                elif stype == GeomAbs_Cone:
                    has_cone = True
                    cone = adaptor.Cone()
                    ref_radius = cone.RefRadius()
                    semi_angle = cone.SemiAngle()
                    
                    bbox = get_bbox(face)
                    altura_cone = abs(bbox[5] - bbox[2])
                    
                    r_min = abs(ref_radius - altura_cone/2 * math.tan(semi_angle))
                    r_max = abs(ref_radius + altura_cone/2 * math.tan(semi_angle))
                    
                    if r_min >= min_allowed:
                        radii.append(r_min)
                    if r_max >= min_allowed:
                        radii.append(r_max)
                    
                    component_info.update({
                        'geometric_type': 'conica',
                        'ref_radius': round(ref_radius, 2),
                        'radius_min': round(r_min, 2),
                        'radius_max': round(r_max, 2),
                        'diameter_min': round(r_min * 2, 2),
                        'diameter_max': round(r_max * 2, 2),
                        'apex_angle': round(math.degrees(semi_angle * 2), 1),
                        'altura': round(altura_cone, 2),
                        'conicidade': round(math.degrees(semi_angle), 1)
                    })
                
                geometric_components.append(component_info)
        
        radii = sorted(set(round(r, 2) for r in radii))
        return radii, has_cone, geometric_components

# ============================================================================
# FUNÇÕES DE ANÁLISE DE FUROS
# ============================================================================

def is_hole_through(faces, shape, tol=2.0):
    """
    Verifica se pelo menos um vértice de qualquer face do grupo toca uma face extrema do bounding box.
    """
    if not shape or shape.IsNull():
        raise ValueError("Shape inválido recebido em is_hole_through.")
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bounds = bbox.Get()
    faces_touched = set()

    for item in faces:
        if len(item) == 5:
            _, face, _, _, _ = item
        elif len(item) == 4:
            _, face, _, _ = item
        else:
            _, face, _ = item
            
        vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while vertex_explorer.More():
            vertex = vertex_explorer.Current()
            pnt = BRep_Tool.Pnt(vertex)
            touched = get_bbox_faces_touched(pnt, bounds, TOL_IS_HOLE)
            faces_touched.update(touched)
            vertex_explorer.Next()
    axis_pairs = [('xmin','xmax'), ('ymin','ymax'), ('zmin','zmax')]
    for a, b in axis_pairs:
        if a in faces_touched and b in faces_touched:
            return True
    return False

def detect_positions_from_holes(axis_groups, shape):
    """Determina o número de posições baseado na direção dos furos."""
    has_up = False
    has_down = False
    for i, faces in enumerate(axis_groups):
        if is_hole_through(faces, shape):
            continue
        
        has_cylinder = False
        for item in faces:
            if len(item) == 5:
                stype, _, _, _, _ = item
            elif len(item) == 4:
                stype, _, _, _ = item
            else:
                stype, _, _ = item
            if stype == GeomAbs_Cylinder:
                has_cylinder = True
                break
                
        if not has_cylinder:
            continue
            
        for item in faces:
            if len(item) == 5:
                stype, face, adaptor, dados, analise = item
            elif len(item) == 4:
                stype, face, adaptor, dados = item
            else:
                stype, face, adaptor = item
                
            if stype == GeomAbs_Cylinder:
                dir_z = adaptor.Cylinder().Axis().Direction().Z()
                bbox = Bnd_Box()
                brepbndlib.Add(face, bbox)
                _, _, zmin, _, _, zmax = bbox.Get()
                depth = abs(zmax - zmin)
                if depth < 2.0:
                    continue
                if dir_z > 0.5:
                    has_up = True
                elif dir_z < -0.5:
                    has_down = True
                break
    if has_up and has_down:
        return 2
    return 1

def group_holes_by_center(features, center_tol=10.0):
    """Agrupa furos (circulares e semicirculares) com base na proximidade do centro."""
    grouped_by_center = []

    def dist2d(a, b):
        return ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) ** 0.5

    for feat in features:
        cx, cy = feat['center'][:2]
        min_d = feat['min_d']
        max_d = feat['max_d']
        added = False
        for group in grouped_by_center:
            gx, gy = group['center']
            if dist2d((cx, cy), (gx, gy)) < center_tol:
                group['min_ds'].append(min_d)
                group['max_ds'].append(max_d)
                added = True
                break
        if not added:
            grouped_by_center.append({
                'center': (cx, cy),
                'min_ds': [min_d],
                'max_ds': [max_d]
            })

    all_counter = Counter()
    for group in grouped_by_center:
        min_d = min(group['min_ds'])
        max_d = max(group['max_ds'])
        key = (round(min_d, 1), round(max_d, 1))
        all_counter[key] += 1
    return all_counter

# ============================================================================
# FUNÇÕES SEMICIRCULARES
# ============================================================================

def collect_semi_circular_arcs(shape, lateral_tol=3.0, angle_tol=30.0, radius_tol=1.0):
    """
    Recolhe todos os recortes semicirculares (arcos ~180º perto da lateral).
    """
    results = []
    piece_bbox = get_bbox(shape)
    xmin, ymin, zmin, xmax, ymax, zmax = piece_bbox
    
    explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    while explorer.More():
        edge = explorer.Current()
        curve_data = BRep_Tool.Curve(edge)
        if not curve_data or len(curve_data) < 3:
            explorer.Next()
            continue
        curve_handle, first, last = curve_data
        adaptor = GeomAdaptor_Curve(curve_handle, first, last)
        if adaptor.GetType() == GeomAbs_Circle:
            circ = adaptor.Circle()
            center = circ.Location()
            r = circ.Radius()
            angle = abs(math.degrees(last - first))
            x, y, z = center.X(), center.Y(), center.Z()
            
            is_large_semicircle = r >= 100.0 and angle >= 90.0
            is_normal_semicircle = 130.0 <= angle <= 180.0
            
            if (is_normal_semicircle or is_large_semicircle) and r >= 10.0:
                tolerance_alignment = 2.0
                
                dist_to_left = abs(x - xmin)
                dist_to_right = abs(x - xmax)
                dist_to_bottom = abs(y - ymin)
                dist_to_top = abs(y - ymax)
                
                aligned_with_left = dist_to_left <= tolerance_alignment
                aligned_with_right = dist_to_right <= tolerance_alignment
                aligned_with_bottom = dist_to_bottom <= tolerance_alignment
                aligned_with_top = dist_to_top <= tolerance_alignment
                
                is_aligned_with_edge = aligned_with_left or aligned_with_right or aligned_with_bottom or aligned_with_top
                
                if is_aligned_with_edge:
                    found = False
                    for rec in results:
                        rx, ry, rz = rec['center']
                        rr = rec['radius']
                        if (abs(rx-x) < RADIUS_TOL and abs(ry-y) < RADIUS_TOL and 
                            abs(rz-z) < RADIUS_TOL and abs(rr-r) < RADIUS_TOL):
                            found = True
                            break
                    
                    if not found:
                        results.append({'center': (x, y, z), 'radius': r, 'angle': angle})
                else:
                    min_edge_distance = min(dist_to_left, dist_to_right, dist_to_bottom, dist_to_top)
                    edge_tolerance = 10.0
                    
                    is_large_center_semicircle = (r >= 100.0 and 
                                                min_edge_distance > 100.0 and
                                                angle >= 170.0)
                    
                    if min_edge_distance <= edge_tolerance:
                        found = False
                        for rec in results:
                            rx, ry, rz = rec['center']
                            rr = rec['radius']
                            if (abs(rx-x) < RADIUS_TOL and abs(ry-y) < RADIUS_TOL and 
                                abs(rz-z) < RADIUS_TOL and abs(rr-r) < RADIUS_TOL):
                                found = True
                                break
                        
                        if not found:
                            results.append({'center': (x, y, z), 'radius': r, 'angle': angle})
                    elif is_large_center_semicircle:
                        found = False
                        for rec in results:
                            rx, ry, rz = rec['center']
                            rr = rec['radius']
                            if (abs(rx-x) < RADIUS_TOL and abs(ry-y) < RADIUS_TOL and 
                                abs(rz-z) < RADIUS_TOL and abs(rr-r) < RADIUS_TOL):
                                found = True
                                break
                        
                        if not found:
                            results.append({'center': (x, y, z), 'radius': r, 'angle': angle})
        
        explorer.Next()
    
    # Filtrar duplicados e erros
    unique_results = []
    for i, result in enumerate(results):
        rx, ry, rz = result['center']
        rr = result['radius']
        
        is_duplicate = False
        for j, existing in enumerate(unique_results):
            ex, ey, ez = existing['center']
            er = existing['radius']
            
            diff_x = abs(rx - ex)
            diff_y = abs(ry - ey) 
            diff_z = abs(rz - ez)
            diff_r = abs(rr - er)
            
            if (diff_x < 0.5 and diff_y < 0.5 and diff_r < 0.5):
                is_duplicate = True
                break
        
        if not is_duplicate:
            unique_results.append(result)
    
    filtered_results = []
    for i, result in enumerate(unique_results):
        is_manufacturing_error = False
        rx, ry, rz = result['center']
        rr = result['radius']
        
        for j, other in enumerate(unique_results):
            if i == j:
                continue
            ox, oy, oz = other['center']
            or_radius = other['radius']
            
            distance = math.sqrt((rx - ox)**2 + (ry - oy)**2)
            
            if (rr < or_radius * 0.3 and
                distance < 20.0 and
                rr < 20.0):
                is_manufacturing_error = True
                break
        
        if not is_manufacturing_error:
            filtered_results.append(result)
    
    return filtered_results

def group_semi_circular_arcs(arcs, group_tol=3.0):
    """Agrupa arcos semicirculares por centro X/Y."""
    grouped = []
    for arc in arcs:
        x, y, z = arc['center']
        r = arc['radius']
        found = False
        for g in grouped:
            gx, gy = g['xy']
            if abs(x - gx) < GROUP_TOL and abs(y - gy) < GROUP_TOL:
                g['zs'].append(z)
                g['radii'].append(r)
                g['angles'].append(arc['angle'])
                found = True
                break
        if not found:
            grouped.append({
                'xy': (x, y),
                'zs': [z],
                'radii': [r],
                'angles': [arc['angle']]
            })
    return grouped

# ============================================================================
# FUNÇÕES RETANGULARES/PLANARES
# ============================================================================

def get_all_planar_faces_bbox(shape):
    """
    Retorna lista de dicts com centro, comprimento, largura, altura e bbox 
    de todas as faces planas (inclusive circulares).
    """
    result = []
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        adaptor = GeomAdaptor_Surface(BRep_Tool.Surface(face))
        stype = adaptor.GetType()
        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        dx = abs(xmax - xmin)
        dy = abs(ymax - ymin)
        dz = abs(zmax - zmin)
        center = ((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2)
        
        face_info = {
            'type': None,
            'center': tuple(round(c, 2) for c in center),
            'length': round(max(dx, dy), 2),
            'width': round(min(dx, dy), 2),
            'height': round(dz, 2),
            'bbox': (xmin, ymin, zmin, xmax, ymax, zmax)
        }
        
        if stype == GeomAbs_Plane:
            face_info['type'] = TIPO_PLANO
        elif stype == GeomAbs_Cylinder:
            face_info['type'] = TIPO_CILINDRO
            try:
                cylinder = adaptor.Cylinder()
                radius = cylinder.Radius()
                face_info['radius'] = round(radius, 2)
                face_info['diameter'] = round(radius * 2, 2)
            except Exception as e:
                pass  # Se não conseguir extrair, continua sem diâmetro
        elif stype == GeomAbs_Cone:
            face_info['type'] = TIPO_CONE
            try:
                cone = adaptor.Cone()
                ref_radius = cone.RefRadius()
                semi_angle = cone.SemiAngle()
                face_info['ref_radius'] = round(ref_radius, 2)
                face_info['semi_angle'] = round(semi_angle, 4)
                face_info['ref_diameter'] = round(ref_radius * 2, 2)
            except Exception as e:
                pass  # Se não conseguir extrair, continua sem diâmetro
        elif stype == GeomAbs_Sphere:
            face_info['type'] = TIPO_ESFERA
        elif stype == GeomAbs_Torus:
            face_info['type'] = TIPO_TORO
        else:
            face_info['type'] = str(stype)
            
        result.append(face_info)
        explorer.Next()
    return result

# ============================================================================
# FUNÇÕES NÃO UTILIZADAS - Mantidas para compatibilidade
# ============================================================================

def merge_aligned_circular_features(features, axis_tol=30.0, center_tol=10.0, diameter_tol=0.5):
    """
    FUNÇÃO NÃO UTILIZADA - Mantida para compatibilidade.
    Agrupa features circulares com centros alinhados numa direção.
    """
    merged = []
    used = set()

    for i, f1 in enumerate(features):
        if i in used:
            continue
        group = [f1]
        for j, f2 in enumerate(features):
            if j <= i or j in used:
                continue
            if abs(f1['max_d'] - f2['max_d']) > diameter_tol:
                continue
            dx = abs(f1['center'][0] - f2['center'][0])
            dy = abs(f1['center'][1] - f2['center'][1])
            dz = abs(f1['center'][2] - f2['center'][2])
            if dx < axis_tol and dy < axis_tol and dz < center_tol:
                group.append(f2)
                used.add(j)
        used.add(i)

        center_avg = tuple(round(sum(f['center'][k] for f in group)/len(group), 2) for k in range(3))
        all_min_d = [f['min_d'] for f in group]
        all_max_d = [f['max_d'] for f in group]
        
        all_geometric_components = []
        has_cone_combined = False
        total_area_combined = 0.0
        all_coordinates = []
        
        for f in group:
            if 'geometric_components' in f:
                all_geometric_components.extend(f['geometric_components'])
            if f.get('has_cone', False):
                has_cone_combined = True
            if 'total_area' in f:
                total_area_combined += f['total_area']
            if 'coordinates' in f:
                all_coordinates.append(f['coordinates'])
        
        combined_coordinates = {}
        if all_coordinates:
            combined_coordinates = all_coordinates[0].copy() if all_coordinates else {}
            
            if len(all_coordinates) > 1:
                all_z_mins = [coord.get('z_range', {}).get('min', 0) for coord in all_coordinates if coord.get('z_range', {}).get('min') is not None]
                all_z_maxs = [coord.get('z_range', {}).get('max', 0) for coord in all_coordinates if coord.get('z_range', {}).get('max') is not None]
                
                if all_z_mins and all_z_maxs:
                    combined_coordinates['z_range'] = {
                        'min': round(min(all_z_mins), 2),
                        'max': round(max(all_z_maxs), 2),
                        'height': round(abs(max(all_z_maxs) - min(all_z_mins)), 2)
                    }
        
        merged_feature = {
            'center': center_avg,
            'min_d': round(min(all_min_d), 2),
            'max_d': round(max(all_max_d), 2),
            'geometric_components': all_geometric_components,
            'num_components': len(all_geometric_components),
            'has_cone': has_cone_combined,
            'total_area': round(total_area_combined, 2),
            'coordinates': combined_coordinates
        }
        
        merged.append(merged_feature)
        
    return merged
