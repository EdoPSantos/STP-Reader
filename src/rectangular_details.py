from .utils import group_rectangular_faces

def show_rectangular_details(shape):
    holes = group_rectangular_faces(shape)
    if not holes:
        print("-- Nenhum furo retangular --\\n")
        return

    print("-- Furos retangulares --")
    for i, face in enumerate(holes, 1):
        print(f"Furo {i}: 1 face (detalhes ainda não implementados)")
    print()
