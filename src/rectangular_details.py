from collections import Counter
from .utils import group_connected_planar_faces, classify_rectangular_face

def show_rectangular_details(shape):
    holes = group_connected_planar_faces(shape)
    if not holes:
        print("\n=== DETALHES DOS FUROS RETANGULARES ===")
        print("-- Nenhum furo retangular --\n")
        return

    # Agrupa por tipo e medidas
    counter = Counter()
    for face in holes:
        tipo, comprimento, largura = classify_rectangular_face(face)
        key = (tipo, comprimento, largura)
        counter[key] += 1

    print("=== DETALHES DOS FUROS RETANGULARES ===")
    for (tipo, comprimento, largura), count in sorted(counter.items(), key=lambda x: (-x[0][1], -x[0][2])):
        if tipo == "quadrado":
            print(f"-- {count} furo(s) [quadrado] com {largura} mm de lado --")
        else:
            print(f"-- {count} furo(s) [retangular] com {comprimento} mm de comprimento e {largura} mm de largura --")
    print()
