from .loader import load_step_file, choose_file_from_folder 
from .main_details import classify_and_summarize_holes
from .cylindrical_details import get_circular_hole_diameters, show_circular_hole_details
from .rectangular_details import show_rectangular_details


def main():
    folder = "stp_igs_files"
    filepath = choose_file_from_folder(folder)
    if not filepath:
        return

    shape = load_step_file(filepath)

    while True:
        print("\n=== MENU PRINCIPAL ===")
        print("1 - Resumo geral")
        print("2 - Detalhes dos furos redondos")
        print("3 - Detalhes dos furos retangulares")
        print("4 - Escolher outro ficheiro")
        print("0 - Sair")

        opcao = input("Escolha uma opção: ").strip()

        if opcao == "1":
            classify_and_summarize_holes(shape, filepath)

        elif opcao == "2":
            circular = get_circular_hole_diameters(shape)
            show_circular_hole_details(circular)

        elif opcao == "3":
            rectangular = []  # A ser implementado futuramente com deteção geométrica de furos retangulares
            show_rectangular_details(rectangular)

        elif opcao == "4":
            new_path = choose_file_from_folder(folder)
            if new_path:
                filepath = new_path
                shape = load_step_file(filepath)
                print(f"Novo ficheiro carregado: {filepath}")

        elif opcao == "0":
            print("Saindo...")
            break

        else:
            print("Opção inválida. Tente novamente.")


if __name__ == "__main__":
    main()