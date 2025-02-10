import file_paths as FILE_PATHS
import os
import shutil

def remove_empty_KOID(targetOrganism: str):
    def get_K_folders(directory: str):
        return [f for f in os.listdir(directory) if f.startswith("K") and os.path.isdir(os.path.join(directory, f))]
    
    if targetOrganism.endswith('.csv'):
        targetOrganism = targetOrganism.replace('.csv', "")

    KOIDs = get_K_folders(FILE_PATHS.GET_ORGANISM_FOLDER_PATH(targetOrganism))

    for KOID in KOIDs:
        folderPath = os.path.join(FILE_PATHS.GET_KOID_FOLDER_PATH(KOID, targetOrganism))
        MSAResultsPath = FILE_PATHS.GET_KOID_MSA_ANALYSIS_DF_FILE_PATH(KOID, targetOrganism)  # Ensure this gets the correct path

        if not os.path.exists(MSAResultsPath):  # Check if the file does not exist
            try:
                shutil.rmtree(folderPath)  # Delete the folder and its contents
                print(f"Deleted folder: {folderPath}")
            except Exception as e:
                print(f"Error deleting folder {folderPath}: {e}")


if __name__ == "__main__":
    remove_empty_KOID("K25706", "test_target_bacteria")