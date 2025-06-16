from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.gp import gp_Pnt
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_VERTEX, TopAbs_FACE
from collections import defaultdict

def axis_key(adaptor):
    surface_type = adaptor.GetType()
    if surface_type == GeomAbs_Cylinder:
        axis = adaptor.Cylinder().Axis()
    elif surface_type == GeomAbs_Cone:
        axis = adaptor.Cone().Axis()
    else:
        raise ValueError("Tipo de superfície não suportado em axis_key")
    loc = axis.Location()
    dir = axis.Direction()
    return (
        round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2),
        round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)
    )

def get_all_vertices_of_group(faces):
    vertices = []
    for face in faces:
        explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while explorer.More():
            v = explorer.Current()
            vertices.append(BRep_Tool.Pnt(v))
            explorer.Next()
    return vertices

def get_bbox_faces_touched(pnt, bbox, tol=0.5):
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    faces = []
    if abs(pnt.X() - xmin) < tol: faces.append("xmin")
    if abs(pnt.X() - xmax) < tol: faces.append("xmax")
    if abs(pnt.Y() - ymin) < tol: faces.append("ymin")
    if abs(pnt.Y() - ymax) < tol: faces.append("ymax")
    if abs(pnt.Z() - zmin) < tol: faces.append("zmin")
    if abs(pnt.Z() - zmax) < tol: faces.append("zmax")
    return faces

def project_point_to_bbox(loc, dir, bounds):
    """
    Projeta o ponto central 'loc' ao longo da direção 'dir' até tocar no bounding box 'bounds'.
    Retorna o ponto na face mais próxima.
    """
    px, py, pz = loc.X(), loc.Y(), loc.Z()
    dx, dy, dz = dir.X(), dir.Y(), dir.Z()
    # Protege contra direção nula
    if abs(dx) < 1e-8: dx = 1e-8
    if abs(dy) < 1e-8: dy = 1e-8
    if abs(dz) < 1e-8: dz = 1e-8

    t_values = []
    for i, (min_b, max_b, p, d) in enumerate(zip(bounds[:3], bounds[3:], [px, py, pz], [dx, dy, dz])):
        t1 = (min_b - p) / d
        t2 = (max_b - p) / d
        t_values.extend([t1, t2])
    t_values = [t for t in t_values if t > 0]
    if not t_values:
        return loc
    t = min(t_values)
    return gp_Pnt(px + t*dx, py + t*dy, pz + t*dz)

def is_hole_through(faces, shape, tol=2.0):
    if not shape or shape.IsNull():
        raise ValueError("Shape inválido recebido em is_hole_through.")
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bounds = bbox.Get()
    hit_faces = set()
    for surface_type, face, adaptor in faces:
        vertex_explorer = TopExp_Explorer(face, TopAbs_VERTEX)
        while vertex_explorer.More():
            vertex = vertex_explorer.Current()
            pnt = BRep_Tool.Pnt(vertex)
            touched = get_bbox_faces_touched(pnt, bounds, tol)
            hit_faces.update(touched)
            vertex_explorer.Next()
    if faces:
        surface_type, face, adaptor = faces[0]
        if surface_type == GeomAbs_Cylinder:
            axis = adaptor.Cylinder().Axis()
        elif surface_type == GeomAbs_Cone:
            axis = adaptor.Cone().Axis()
        else:
            axis = None
        if axis:
            loc = axis.Location()
            dir = axis.Direction()
            for k in range(-3, 4):
                t = k * 0.25
                px = loc.X() + dir.X() * t * (bounds[3] - bounds[0])
                py = loc.Y() + dir.Y() * t * (bounds[4] - bounds[1])
                pz = loc.Z() + dir.Z() * t * (bounds[5] - bounds[2])
                pnt = gp_Pnt(px, py, pz)
                touched = get_bbox_faces_touched(pnt, bounds, tol)
                hit_faces.update(touched)
    axis_pairs = [('xmin','xmax'), ('ymin','ymax'), ('zmin','zmax')]
    for a, b in axis_pairs:
        if a in hit_faces and b in hit_faces:
            return True
    return False

def get_group_diameters(faces):
    diams = []
    for surface_type, face, adaptor in faces:
        if surface_type == GeomAbs_Cylinder:
            diams.append(round(adaptor.Cylinder().Radius() * 2, 2))
        elif surface_type == GeomAbs_Cone:
            # Podes usar o raio de referência ou calcular ambos os extremos
            base = round(adaptor.Cone().RefRadius() * 2, 2)
            semi_angle = adaptor.Cone().SemiAngle()
            # Se quiseres, calcula o maior diâmetro do cone pelas extremidades (opcional)
            diams.append(base)
        # outros tipos ignorados
    if not diams:
        return None, None
    return min(diams), max(diams)

def classify_hole_group_type(faces):
    types = set()
    for stype, _, _ in faces:
        if stype == GeomAbs_Cylinder:
            types.add("cilíndrico")
        elif stype == GeomAbs_Cone:
            types.add("cónico")
    if types == {"cilíndrico"}:
        return "cilíndrico"
    elif types == {"cónico"}:
        return "cónico"
    elif types == {"cilíndrico", "cónico"}:
        return "misto"
    else:
        return "desconhecido"
    
def detect_positions(shape):
    """
    Deteta se é preciso virar a peça (1 ou 2 posições) com base nos furos cónicos não passantes.
    """
    from .cylindrical_details import group_faces_by_axis
    from OCC.Core.GeomAbs import GeomAbs_Cone, GeomAbs_Cylinder

    axis_groups = group_faces_by_axis(shape)
    orientations = []

    for faces in axis_groups:
        has_cone = any(stype == GeomAbs_Cone for stype, _, _ in faces)
        # Só analisar furos que têm cone (cabeça de embutir)
        if has_cone:
            for stype, face, adaptor in faces:
                if stype == GeomAbs_Cone:
                    axis = adaptor.Cone().Axis()
                    dir_z = axis.Direction().Z()
                    # Marca como para cima (+Z) ou para baixo (-Z)
                    if dir_z > 0.5:
                        orientations.append('up')
                    elif dir_z < -0.5:
                        orientations.append('down')
    # Se existir 'up' e 'down', são necessárias 2 posições
    if 'up' in orientations and 'down' in orientations:
        return 2
    # Só uma orientação -> 1 posição
    return 1

def group_faces_by_axis(shape):
    """
    Agrupa todas as faces cilíndricas pelo eixo.
    """
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    axis_groups = defaultdict(list)
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        surface_type = adaptor.GetType()
        if surface_type == GeomAbs_Cylinder:
            key = axis_key(adaptor)
            axis_groups[key].append((surface_type, face, adaptor))
        explorer.Next()
    return list(axis_groups.values())

def get_min_max_cylindrical_diameter(shape):
    """
    Calcula o diâmetro mínimo e máximo entre todos os furos cilíndricos da peça.
    """
    from OCC.Core.GeomAbs import GeomAbs_Cylinder
    min_d = None
    max_d = None
    axis_groups = group_faces_by_axis(shape)
    for faces in axis_groups:
        only_cylinders = all(stype == GeomAbs_Cylinder for stype, _, _ in faces)
        if only_cylinders:
            # Considera todos os cilindros no grupo
            for stype, face, adaptor in faces:
                d = adaptor.Cylinder().Radius() * 2
                if (min_d is None) or (d < min_d):
                    min_d = d
                if (max_d is None) or (d > max_d):
                    max_d = d
    if min_d is not None and max_d is not None:
        return round(min_d, 2), round(max_d, 2)
    else:
        return None, None
