import os
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone

# Função para carregar ficheiro STEP
def load_step_file(filename):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise Exception("Erro ao ler ficheiro STEP")
    step_reader.TransferRoots()
    shape = step_reader.OneShape()
    return shape

# Função para permitir seleção do ficheiro no terminal
def choose_file_from_folder(folder):
    files = [f for f in os.listdir(folder) if f.lower().endswith('.stp')]
    if not files:
        print("Nenhum ficheiro .stp encontrado na pasta.")
        return None

    print("Escolha um ficheiro .stp:")
    for i, f in enumerate(files):
        print(f"{i+1}: {f}")

    while True:
        choice = input("Número do ficheiro: ")
        if choice.isdigit():
            choice = int(choice)
            if 1 <= choice <= len(files):
                return os.path.join(folder, files[choice - 1])
        print("Escolha inválida. Tente novamente.")
