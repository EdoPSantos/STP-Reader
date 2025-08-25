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
    if shape.IsNull():
        raise Exception("Shape inválido ou vazio")
    return shape

# Função para permitir seleção do ficheiro no terminal
def choose_file_from_folder(folder):
    files = [f for f in os.listdir(folder) if f.lower().endswith(('.stp', '.step', '.igs'))]
    if not files:
        print("Nenhum ficheiro .stp/.step/.igs encontrado na pasta.")
        return None

    print("Escolha um ficheiro .stp/.step/.igs:")
    for i, f in enumerate(files):
        print(f"{i+1}: {f}")
    print("a: Processar TODOS os ficheiros")

    while True:
        choice = input("Número do ficheiro ou 'a' para todos: ").strip().lower()
        if choice == 'a':
            return 'ALL'
        elif choice.isdigit():
            choice = int(choice)
            if 1 <= choice <= len(files):
                return os.path.join(folder, files[choice - 1])
        print("Escolha inválida. Tente novamente.")

def get_all_files_from_folder(folder):
    """Retorna lista de todos os ficheiros .stp/.step/.igs na pasta"""
    files = [f for f in os.listdir(folder) if f.lower().endswith(('.stp', '.step', '.igs'))]
    return [os.path.join(folder, f) for f in files]
