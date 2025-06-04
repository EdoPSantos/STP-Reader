import os
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.gp import gp_Ax1


def load_step_file(filename):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise Exception("Erro ao ler ficheiro STEP")
    step_reader.TransferRoots()
    shape = step_reader.OneShape()
    return shape

def axis_key(adaptor):
    """Gera uma chave única simplificada com centro e direção."""
    axis = adaptor.Cylinder().Axis()
    loc = axis.Location()
    dir = axis.Direction()
    return (
        round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2),
        round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)
    )

def count_cylindrical_holes_by_axis(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    axis_set = set()

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cylinder:
            key = axis_key(adaptor)
            axis_set.add(key)
        explorer.Next()

    return len(axis_set)

def analyze_faces(shape):
    from collections import defaultdict

    radius_count = defaultdict(int)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        geom_type = adaptor.GetType()

        if geom_type == GeomAbs_Cylinder:
            cylinder = adaptor.Cylinder()
            radius = round(cylinder.Radius(), 2)
            radius_count[radius] += 1

        explorer.Next()

    if radius_count:
        print("Resumo dos furos cilíndricos por raio:")
        for radius, count in sorted(radius_count.items()):
            print(f"Existem {count} circunferências com {radius:.2f} de raio")
    else:
        print("Nenhuma face cilíndrica encontrada.")

def choose_file_from_folder(folder):
    files = [f for f in os.listdir(folder) if f.lower().endswith('.stp')]
    if not files:
        print("Nenhum ficheiro .stp encontrado na pasta.")
        return None

    print("Escolha um ficheiro .stp:")
    for i, f in enumerate(files):
        print(f"{i+1}: {f}")

    while True:
        escolha = input("Número do ficheiro: ")
        if escolha.isdigit():
            escolha = int(escolha)
            if 1 <= escolha <= len(files):
                return os.path.join(folder, files[escolha - 1])
        print("Escolha inválida. Tente novamente.")

if __name__ == "__main__":
    folder = "./viewers"
    filepath = choose_file_from_folder(folder)
    if filepath:
        shape = load_step_file(filepath)
        num_furos_cil = count_cylindrical_holes_by_axis(shape)
        print(f"Número estimado de furos cilíndricos: {num_furos_cil}")
        analyze_faces(shape)
