from .loader import load_step_file, choose_file_from_folder 
from .main_details import show_general_summary
from .cylindrical_details import get_circular_hole_diameters_by_type, show_circular_hole_details
from .rectangular_details import show_rectangular_details
from .operations_details import show_operations_details  # NOVO!

def main():
    folder = "stp_igs_files"
    filepath = choose_file_from_folder(folder)
    if not filepath:
        return

    shape = load_step_file(filepath)

    if not shape or shape.IsNull():
        print("Erro: shape inválido ou ficheiro STEP corrompido.")
        return

    while True:
        print("\n=== MENU PRINCIPAL ===")
        print("1 - Resumo da peça para o GESTi")
        print("2 - Detalhes dos furos redondos")
        print("3 - Detalhes dos furos retangulares")
        print("4 - Detalhes de posição e operações")
        print("5 - Escolher outro ficheiro")
        print("0 - Sair")

        opcao = input("Escolha uma opção: ").strip()

        if opcao == "1":
            show_general_summary(shape, filepath)

        elif opcao == "2":
            diameters_through, diameters_closed = get_circular_hole_diameters_by_type(shape)
            show_circular_hole_details(diameters_through, diameters_closed)

        elif opcao == "3":
            show_rectangular_details(shape)
        
        elif opcao == "4":
            show_operations_details(shape, filepath)  # Nova função, ver abaixo
        
        elif opcao == "5":
            new_path = choose_file_from_folder(folder)
            if new_path:
                new_shape = load_step_file(new_path)
                if not new_shape or new_shape.IsNull():
                    print("Erro ao carregar novo ficheiro STEP.")
                else:
                    filepath = new_path
                    shape = new_shape
                    print(f"Novo ficheiro carregado: {filepath}")
        
        elif opcao == "0":
            print("Saindo...")
            break

        else:
            print("Opção inválida. Tente novamente.")

if __name__ == "__main__":
    main()
