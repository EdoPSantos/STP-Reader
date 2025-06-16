from OCC.Core.BRep import BRep_Tool
from OCC.Core.GeomAbs import GeomAbs_Cylinder
from OCC.Core.GeomAdaptor import GeomAdaptor_Surface
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from collections import defaultdict
from .utils import get_min_max_cylindrical_diameter, group_faces_by_axis

import os
import re

def extract_from_filename(filepath):
    """
    Extrai o nome do molde e da peça a partir do nome do ficheiro.
    """
    base = os.path.basename(filepath)
    name, _ = os.path.splitext(base)
    parts = re.split(r"[_\-]", name)
    if len(parts) >= 2:
        return parts[-2], parts[-1]
    return None, None

def extract_mold_and_part_from_step(filepath):
    """
    Tenta extrair molde e peça a partir do ficheiro STEP.
    Procura por nome do ficheiro e por entradas PRODUCT no texto.
    """
    file_mold, file_part = extract_from_filename(filepath)
    if file_mold and file_part:
        return file_mold, file_part

    mold, part = None, None
    product_mold, product_part = None, None
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "FILE_NAME" in line and (mold is None or part is None):
                match = re.search(r"FILE_NAME\('([^']+)'", line)
                if match:
                    filename = match.group(1)
                    base = os.path.basename(filename)
                    name, _ = os.path.splitext(base)
                    parts = re.split(r"[_\-]", name)
                    if len(parts) >= 2:
                        mold = parts[-2]
                        part = parts[-1]
            if "PRODUCT('" in line and (product_mold is None or product_part is None):
                match = re.search(r"PRODUCT\('([^']+)'", line)
                if match:
                    pname = match.group(1)
                    parts = re.split(r"[_\-]", pname)
                    if len(parts) >= 2:
                        product_mold = parts[-2]
                        product_part = parts[-1]

    if mold and part:
        return mold, part
    elif product_mold and product_part:
        return product_mold, product_part
    return None, None

def show_general_summary(shape, filepath=None):
    """
    Mostra apenas os campos principais para preenchimento do GESTi:
    Molde, Peça, Largura, Comprimento, Posições, Qtd. furos cilíndricos, Qtd. furos retangulares.
    """
    # Extrai o molde e a peça
    mold, part = None, None
    if filepath:
        mold, part = extract_mold_and_part_from_step(filepath)
        mold_name = f"{mold}" if mold else "(não identificado)"
        part_name = f"{part}" if part else "(não identificada)"
    else:
        mold_name = "(desconhecido)"
        part_name = "(desconhecida)"

    # Calcula bounding box
    bbox = Bnd_Box()
    brepbndlib.Add(shape, bbox)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    length = abs(xmax - xmin)
    width = abs(ymax - ymin)

    # Conta furos cilíndricos (por eixo)
    axis_groups = group_faces_by_axis(shape)
    count_cylindrical = len(axis_groups)

    # Lista dos diâmetros de cada furo
    diameters = []
    from OCC.Core.GeomAbs import GeomAbs_Cylinder
    for faces in axis_groups:
        only_cylinders = all(stype == GeomAbs_Cylinder for stype, _, _ in faces)
        if only_cylinders:
            for stype, face, adaptor in faces:
                diameters.append(round(adaptor.Cylinder().Radius() * 2, 2))

    # Conta furos retangulares - placeholder para implementar
    count_rectangular = 0

    # Detectar automaticamente número de posições (exemplo: sempre 1)
    positions = 1  # Podes substituir por uma função tipo detect_positions(shape)

    print("="*40)
    print("   RESUMO DA PEÇA PARA O GESTi")
    print("="*40)
    print(f"Molde: {mold_name}")
    print(f"Peça: {part_name}")
    print(f"Largura: {width:.0f} mm")
    print(f"Comprimento: {length:.0f} mm")
    print(f"Posições: {positions}")
    print(f"Qtd. furos cilíndricos: {count_cylindrical}")
    print(f"Qtd. furos quadrados: {count_rectangular}")
    print("-"*40)
    if diameters:
        print("Diâmetro dos furos cilíndricos:", ', '.join(str(d) + " mm" for d in sorted(set(diameters))))
    else:
        print("Diâmetro  dos furos cilíndricos: -")
    print("="*40)
