import file_paths as FILE_PATHS

def run_analysis(KOID: str, target_organism: str):
    if target_organism.endswith('.csv'):
        target_organism = target_organism.replace('.csv', "")
    data = FILE_PATHS.crate_target_analysis_file(target_organism)
