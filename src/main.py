from .loader import load_step_file, choose_file_from_folder, get_all_files_from_folder
from .main_details import show_general_summary

def process_all_files(folder):
    """Processa todos os ficheiros na pasta"""
    all_files = get_all_files_from_folder(folder)
    
    if not all_files:
        print("Nenhum ficheiro encontrado na pasta.")
        return
    
    # Perguntar formato preferido
    print("\nEscolha o formato de exportação:")
    print("1 - TXT (padrão)")
    print("2 - Excel")
    format_choice = input("Formato (1 ou 2): ").strip()
    
    export_format = 'excel' if format_choice == '2' else 'txt'
    print(f"Formato escolhido: {export_format.upper()}")
    
    print(f"\n=== PROCESSANDO {len(all_files)} FICHEIROS ===")
    
    for i, filepath in enumerate(all_files, 1):
        filename = filepath.split('\\')[-1]  # Extrair apenas o nome do ficheiro
        print(f"\n[{i}/{len(all_files)}] Processando: {filename}")
        print("=" * 60)
        
        try:
            shape = load_step_file(filepath)
            
            if not shape or shape.IsNull():
                print(f"ERRO: Ficheiro {filename} inválido ou corrompido.")
                continue
            
            # Processar com formato escolhido
            show_general_summary(shape, filepath, save_to_file=True, export_format=export_format)
            
        except Exception as e:
            print(f"ERRO ao processar {filename}: {str(e)}")
            continue
    
    print(f"\n=== PROCESSAMENTO CONCLUÍDO ({len(all_files)} ficheiros) ===")

def main():
    folder = "stp_igs_files"
    filepath = choose_file_from_folder(folder)
    if not filepath:
        return
    
    # Se escolheu processar todos os ficheiros
    if filepath == 'ALL':
        process_all_files(folder)
        return

    shape = load_step_file(filepath)

    if not shape or shape.IsNull():
        print("Erro: shape inválido ou ficheiro STEP corrompido.")
        return
    
    # Processar automaticamente o resumo da peça (gera TXT automaticamente)
    show_general_summary(shape, filepath, save_to_file=True, export_format='txt')
    
    while True:
        print("\n=== MENU PRINCIPAL ===")
        print("1 - Exportar para Excel")
        print("2 - Escolher outro ficheiro")
        print("0 - Sair")

        opcao = input("Escolha uma opção: ").strip().lower()

        if opcao == "1":
            print("Exportando para Excel...")
            show_general_summary(shape, filepath, save_to_file=True, export_format='excel')

        elif opcao == "2":
            new_path = choose_file_from_folder(folder)
            if new_path:
                if new_path == 'ALL':
                    process_all_files(folder)
                    break
                else:
                    new_shape = load_step_file(new_path)
                    if not new_shape or new_shape.IsNull():
                        print("Erro ao carregar novo ficheiro STEP.")
                    else:
                        filepath = new_path
                        shape = new_shape
                        print(f"Novo ficheiro carregado: {filepath}")
                        # Processar automaticamente o novo ficheiro (gera TXT automaticamente)
                        show_general_summary(shape, filepath, save_to_file=True, export_format='txt')
        
        elif opcao == "0":
            print("Saindo...")
            break

        else:
            print("Opção inválida. Tente novamente.")

if __name__ == "__main__":
    main()
