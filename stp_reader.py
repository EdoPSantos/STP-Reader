import os
from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Core.IFSelect import IFSelect_RetDone
from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder, GeomAbs_Plane, GeomAbs_Cone
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

def cone_axis_key(adaptor):
    """Gera uma chave única simplificada para cones, usando o eixo."""
    axis = adaptor.Cone().Axis()
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

def count_conical_holes_by_axis(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    axis_set = set()

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cone:
            key = cone_axis_key(adaptor)
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

def analyze_general_summary(shape):
    from OCC.Core.gp import gp_Dir

    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    planar_faces = []

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)

        if adaptor.GetType() == GeomAbs_Plane:
            normal = adaptor.Plane().Axis().Direction()
            # Consideramos superfícies planas exteriores se a normal for ±Z ou ±X ou ±Y
            dir_vec = (round(abs(normal.X()), 1), round(abs(normal.Y()), 1), round(abs(normal.Z()), 1))
            if dir_vec in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]:
                planar_faces.append(face)

        explorer.Next()

    # Furos reais
    cylindrical_holes = count_cylindrical_holes_by_axis(shape)
    conical_holes = count_conical_holes_by_axis(shape)

    print("="*26)
    print("      RESUMO GERAL")
    print("="*26)
    print(f"Faces planas exteriores: {len(planar_faces)}")
    print(f"Furos cilíndricos distintos: {cylindrical_holes}")
    print(f"Furos cônicos distintos: {conical_holes}")
    print()


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

def analyze_cylindrical_holes(shape):
    print("="*26)
    print("    FUROS CILÍNDRICOS")
    print("="*26)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    holes = {}
    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cylinder:
            key = axis_key(adaptor)
            cylinder = adaptor.Cylinder()
            radius = round(cylinder.Radius(), 2)
            diameter = round(2 * radius, 2)
            axis = cylinder.Axis()
            loc = axis.Location()
            dir = axis.Direction()
            # Agrupa por eixo (um furo por eixo)
            if key not in holes:
                holes[key] = {
                    "radius": radius,
                    "diameter": diameter,
                    "center": (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2)),
                    "direction": (round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)),
                    "faces": []
                }
            holes[key]["faces"].append(face)
        explorer.Next()

    for i, (key, data) in enumerate(holes.items(), 1):
        print(f"Furo {i}:")
        print(f"  Centro: {data['center']}")
        print(f"  Direção: {data['direction']}")
        print(f"  Raio: {data['radius']:.2f}")
        print(f"  Diâmetro: {data['diameter']:.2f}")
        # Comprimento e tipo (passante/cego) podem ser estimados com mais geometria
        print()

        def analyze_conical_holes(shape):
    print("="*26)
    print("     FUROS CÔNICOS")
    print("="*26)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    con_holes = {}

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cone:
            key = cone_axis_key(adaptor)
            cone = adaptor.Cone()
            radius = round(cone.RefRadius(), 2)
            diameter = round(2 * radius, 2)
            axis = cone.Axis()
            loc = axis.Location()
            dir = axis.Direction()
            if key not in con_holes:
                con_holes[key] = {
                    "center": (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2)),
                    "direction": (round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)),
                    "radius": radius,
                    "diameter": diameter
                }
        explorer.Next()

    for i, (key, data) in enumerate(con_holes.items(), 1):
        print(f"Furo cônico {i}:")
        print(f"  Centro: {data['center']}")
        print(f"  Direção: {data['direction']}")
        print(f"  Raio: {data['radius']:.2f}")
        print(f"  Diâmetro: {data['diameter']:.2f}")
        print()

    print(f"Total de furos cônicos: {len(con_holes)}\n")


def analyze_cylindrical_and_conical_holes(shape):
    print("="*26)
    print("    FUROS CILÍNDRICOS")
    print("="*26)
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    cyl_holes = {}
    con_holes = {}

    while explorer.More():
        face = explorer.Current()
        surface = BRep_Tool.Surface(face)
        adaptor = GeomAdaptor_Surface(surface)
        if adaptor.GetType() == GeomAbs_Cylinder:
            key = axis_key(adaptor)
            cylinder = adaptor.Cylinder()
            radius = round(cylinder.Radius(), 2)
            diameter = round(2 * radius, 2)
            axis = cylinder.Axis()
            loc = axis.Location()
            dir = axis.Direction()
            if key not in cyl_holes:
                cyl_holes[key] = {
                    "center": (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2)),
                    "direction": (round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)),
                    "radius": radius,
                    "diameter": diameter
                }
        elif adaptor.GetType() == GeomAbs_Cone:
            key = cone_axis_key(adaptor)
            cone = adaptor.Cone()
            radius = round(cone.RefRadius(), 2)
            diameter = round(2 * radius, 2)
            axis = cone.Axis()
            loc = axis.Location()
            dir = axis.Direction()
            if key not in con_holes:
                con_holes[key] = {
                    "center": (round(loc.X(), 2), round(loc.Y(), 2), round(loc.Z(), 2)),
                    "direction": (round(dir.X(), 3), round(dir.Y(), 3), round(dir.Z(), 3)),
                    "radius": radius,
                    "diameter": diameter
                }
        explorer.Next()

if __name__ == "__main__":
    folder = "viewers"
    filepath = choose_file_from_folder(folder)
    if filepath:
        shape = load_step_file(filepath)
        while True:
            print("\nO que deseja visualizar?")
            print("1 - Resumo geral")
            print("2 - Detalhes dos furos cilíndricos")
            print("3 - Detalhes dos furos cônicos")
            print("0 - Sair")
            opcao = input("Escolha uma opção: ").strip()
            if opcao == "1":
                analyze_general_summary(shape)
            elif opcao == "2":
                analyze_cylindrical_holes(shape)
            elif opcao == "3":
                analyze_conical_holes(shape)
            elif opcao == "0":
                print("Saindo...")
                break
            else:
                print("Opção inválida. Tente novamente.")

            # Submenu após mostrar detalhes
            while True:
                print("1 - Voltar ao menu principal")
                print("0 - Sair")
                subop = input("Escolha uma opção: ").strip().lower()
                if subop == "1":
                    break
                elif subop == "0":
                    print("Saindo...")
                    exit()
                else:
                    print("Opção inválida. Tente novamente.")