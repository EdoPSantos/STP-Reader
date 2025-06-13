from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone

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
