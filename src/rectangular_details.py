from collections import Counter

def show_rectangular_details(shape):
    """
    Mostra furos retangulares (incluindo slots com extremos arredondados) e quadrados,
    diferenciando entre furos internos e slots abertos à borda.
    """
    # Detetar slots oblongos (retangulares arredondados)
    slots = find_slots(shape)
    slot_counter = Counter(slots)

    if slots:
        print("=== DETALHES DOS SLOTS OBLONGOS (RETANGULARES ARREDONDADOS) ===")
        for (length, width), count in sorted(slot_counter.items(), key=lambda x: (-x[0][0], -x[0][1])):
            print(f"-- {count} slot(s) oblongos com {length} mm de comprimento e {width} mm de largura --")
        print()

    # Faces planas não circulares (retangulares e quadradas)
    faces = group_non_circular_planar_faces(shape)
    filtered_faces = auto_filter_rectangular_faces(faces)
    piece_length, piece_width = get_piece_dimensions(shape)
    abs_tol = 2.0
    rel_tol = 0.97

    counter = Counter()
    piece_bbox = get_bbox(shape)
    for face, shape_type, length, width in filtered_faces:
        # Ignora faces que são slots abertos (à borda)
        if is_open_edge(face, piece_bbox, tol=2.0):
            continue
        # Só elimina faces totalmente exteriores (face igual à da chapa, chanfro)
        is_face_exterior = (
            (abs(length - piece_length) < abs_tol and abs(width - piece_width) < abs_tol) or
            (abs(length - piece_width) < abs_tol and abs(width - piece_length) < abs_tol)
        )
        # Corrigido: Mantém slots abertos à borda, mas agora a função espera o bbox!
        if is_face_exterior and not is_open_edge(face, piece_bbox, tol=abs_tol):
            continue

        # Ignora duplicados que já são slots
        if any(abs(length - slen) < 2.0 and abs(width - swid) < 2.0 for (slen, swid) in slot_counter):
            continue

        # Elimina casos extremos (faces quase do tamanho da chapa)
        if (
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
