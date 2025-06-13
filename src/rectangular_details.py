# codigo parecido com o cylindrical details.py mas para retangulos ou formatos que não sejão cilindricos

def show_rectangular_details(rectangular_holes):
    if not rectangular_holes:
        print("-- Nenhum furo retangular --\n")
        return

    print("-- Furos retangulares --")
    for i, faces in enumerate(rectangular_holes, 1):
        print(f"Furo {i}: {len(faces)} face(s) (detalhes não implementados)")
    print()
