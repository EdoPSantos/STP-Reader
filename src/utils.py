from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.gp import gp_Pnt
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib

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

def get_bbox_faces_touched(pnt, bbox, tol=0.5):
    """
    Retorna as faces do bounding box tocadas por um ponto.
    Possíveis valores: 'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax'
    """
    xmin, ymin, zmin, xmax, ymax, zmax = bbox
    faces = []
    if abs(pnt.X() - xmin) < tol:
        faces.append("xmin")
    if abs(pnt.X() - xmax) < tol:
        faces.append("xmax")
    if abs(pnt.Y() - ymin) < tol:
        faces.append("ymin")
    if abs(pnt.Y() - ymax) < tol:
        faces.append("ymax")
    if abs(pnt.Z() - zmin) < tol:
        faces.append("zmin")
    if abs(pnt.Z() - zmax) < tol:
        faces.append("zmax")
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

def is_hole_through(loc, dir, shape, tol=0.5):
    """
    Devolve True se o furo for passante (sem fundo), False se for fechado (com fundo).
    A verificação é feita apenas quando os dois extremos tocam faces opostas do mesmo eixo (ex: zmin/zmax).
    """
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    bounds = bbox.Get()

    pnt1 = project_point_to_bbox(loc, dir, bounds)
    pnt2 = project_point_to_bbox(loc, dir.Reversed(), bounds)

    faces1 = get_bbox_faces_touched(pnt1, bounds, tol)
    faces2 = get_bbox_faces_touched(pnt2, bounds, tol)

    # Para cada eixo, se um extremo toca uma face e o outro a face oposta, é passante
    axis_faces = [('xmin', 'xmax'), ('ymin', 'ymax'), ('zmin', 'zmax')]
    for face1, face2 in axis_faces:
        if (face1 in faces1 and face2 in faces2) or (face2 in faces1 and face1 in faces2):
            return True  # FURO PASSANTE (sem fundo)
    return False  # FURO FECHADO (com fundo)
