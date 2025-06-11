from .loader import load_step_file, choose_file_from_folder
from .analyzer import analyze_general_summary, analyze_cylindrical_holes, analyze_conical_holes

def main():
    folder = "stp_igs_files"
    filepath = choose_file_from_folder(folder)
    if not filepath:
        return

    shape = load_step_file(filepath)

    while True:
        print("\n=== MENU PRINCIPAL ===")
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

if __name__ == "__main__":
    main()
