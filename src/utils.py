# Funções utilitárias para gerar chaves únicas com base nos eixos dos furos

def axis_key(adaptor):
    axis = adaptor.Cylinder().Axis()
    loc = axis.Location()
    dir = axis.Direction()
    return (
        round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2),
        round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)
    )

def cone_axis_key(adaptor):
    axis = adaptor.Cone().Axis()
    loc = axis.Location()
    dir = axis.Direction()
    return (
        round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2),
        round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)
    )
