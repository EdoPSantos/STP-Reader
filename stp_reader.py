import os
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Cone
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer, topexp_MapShapesAndAncestors
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape


def load_step_file(filename):
    step_reader = STEPControl_Reader()
    status = step_reader.ReadFile(filename)
    if status != IFSelect_RetDone:
        raise Exception("Erro ao ler ficheiro STEP")
    step_reader.TransferRoots()
    shape = step_reader.OneShape()
    return shape

def count_holes_by_type(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    cylinder_faces = 0
    cone_faces = 0

    faces = []
    while explorer.More():
        face = explorer.Current()
        faces.append(face)
        explorer.Next()

    # Mapa de todas as arestas e as faces associadas
    map_faces = TopTools_IndexedDataMapOfShapeListOfShape()
    topexp_MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, map_faces)

    for face in faces:
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        geom_type = adaptor.GetType()

        if geom_type not in (GeomAbs_Cylinder, GeomAbs_Cone):
            continue

        has_cone_neighbor = False
        edge_explorer = TopExp_Explorer(face, TopAbs_EDGE)
        while edge_explorer.More():
            edge = edge_explorer.Current()

            try:
                neighbor_faces = map_faces.FindFromKey(edge)
            except:
                edge_explorer.Next()
                continue

            neighbors = []
            i = 1
            while True:
                try:
                    neighbors.append(neighbor_faces.Value(i))
                    i += 1
                except:
                    break

            for neighbor in neighbors:
                if neighbor.IsSame(face):
                    continue
                neighbor_surface = BRep_Tool.Surface(neighbor)
                neighbor_adaptor = GeomAdaptor_Surface(neighbor_surface)
                neighbor_type = neighbor_adaptor.GetType()
                if neighbor_type == GeomAbs_Cone:
                    has_cone_neighbor = True
            edge_explorer.Next()

        if geom_type == GeomAbs_Cylinder and not has_cone_neighbor:
            cylinder_faces += 1
        elif geom_type == GeomAbs_Cone:
            cone_faces += 1

    return cylinder_faces // 2, cone_faces // 2


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
        num_furos_cil, num_furos_cone = count_holes_by_type(shape)
        print(f"Número de furos cilíndricos: {num_furos_cil}")
        print(f"Número de furos cônicos: {num_furos_cone}")
        analyze_faces(shape)
