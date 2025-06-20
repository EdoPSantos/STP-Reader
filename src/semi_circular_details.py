from .utils import get_bbox, collect_semi_circular_arcs, group_semi_circular_arcs

def show_semi_circular_details(shape):
    """
    Mostra todos os recortes semicirculares agrupados (evita duplicados por Z).
    """
    lateral_tol = 3.0  # tolerância em mm para considerar na lateral
    group_tol = 3.0    # tolerância para agrupar por centro X/Y (mm)

    # 1. Recolhe todos os arcos semicirculares candidatos
    raw_results = collect_semi_circular_arcs(shape, lateral_tol=lateral_tol)
    # 2. Agrupa por X/Y
    grouped = group_semi_circular_arcs(raw_results, group_tol=group_tol)

    print("\n=== DETALHES DOS RECORTES SEMICIRCULARES ===")
    if not grouped:
        print("Nenhum recorte semicircular identificado automaticamente.")
    else:
        for idx, g in enumerate(grouped, 1):
            x, y = g['xy']
            zmin, zmax = min(g['zs']), max(g['zs'])
            rmin, rmax = min(g['radii']), max(g['radii'])
            angle_avg = sum(g['angles']) / len(g['angles'])
            area_avg = sum(g['areas']) / len(g['areas'])
            print(f"\nSemicírculo {idx}:")
            print(f"  Centro: X={x:.2f}  Y={y:.2f}")
            print(f"  Z (mín/máx): {zmin:.2f} / {zmax:.2f}")
            print(f"  Raio (mín/máx): {rmin:.2f} / {rmax:.2f} mm")
            print(f"  Ângulo médio do arco: {angle_avg:.1f}º")
            print(f"  Área média (aprox): {area_avg:.2f} mm²")
            print(f"  Ocorrências agrupadas: {len(g['zs'])}")
        print(f"\nTotal: {len(grouped)} recorte(s) semicircular(es) agrupado(s).")
    print("=" * 50)
