import os
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface

def load_step_file(filename):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise Exception("Erro ao ler ficheiro STEP")

    step_reader.TransferRoots()
    shape = step_reader.OneShape()
    return shape

def count_cylindrical_holes(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    cylinder_faces = 0
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        geom_type = adaptor.GetType()
        if geom_type == GeomAbs_Cylinder:
            cylinder_faces += 1
        explorer.Next()
    return cylinder_faces // 2  # Assume 2 faces por furo atravessante

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
        num_furos = count_cylindrical_holes(shape)
        print(f"Número de furos (cilíndricos) no ficheiro '{os.path.basename(filepath)}': {num_furos}")
