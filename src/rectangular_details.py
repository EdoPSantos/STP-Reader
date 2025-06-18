from collections import Counter
from .utils import group_non_circular_planar_faces, classify_rectangular_face, auto_filter_rectangular_faces, get_piece_dimensions

def show_rectangular_details(shape):
    faces = group_non_circular_planar_faces(shape)
    filtered_faces = auto_filter_rectangular_faces(faces)
    if not filtered_faces:
        print("\n=== DETALHES DOS FUROS RETANGULARES ===")
        print("-- Nenhum furo retangular --\n")
        return

    piece_length, piece_width = get_piece_dimensions(shape)
    # Parâmetros universais de filtragem para ignorar bordas/chanfros/exteriores
    abs_tol = 2.0  # tolerância absoluta em mm
    rel_tol = 0.97  # tolerância relativa, 97% do tamanho total

    counter = Counter()
    for face, shape_type, length, width in filtered_faces:
        # Ignorar qualquer face em que pelo menos um lado é igual ou quase igual ao exterior da peça
        if (
            abs(length - piece_length) < abs_tol or
            abs(length - piece_width) < abs_tol or
            abs(width - piece_length) < abs_tol or
            abs(width - piece_width) < abs_tol or
            length > rel_tol * piece_length or
            width > rel_tol * piece_width or
            length > rel_tol * piece_width or
            width > rel_tol * piece_length
        ):
            continue
        key = (shape_type, length, width)
        counter[key] += 1

    print("=== DETALHES DOS FUROS RETANGULARES ===")
    if not counter:
        print("-- Nenhum furo retangular --\n")
    else:
        for (shape_type, length, width), count in sorted(counter.items(), key=lambda x: (-x[0][1], -x[0][2])):
            if shape_type == "square":
                print(f"-- {count} furo(s) [quadrado] com {width} mm de lado --")
            else:
                print(f"-- {count} furo(s) [retangular] com {length} mm de comprimento e {width} mm de largura --")
    print()
