import os
import pandas as pd

HOME_DIRECTORY = os.getcwd()
DATAFILES = os.path.join(HOME_DIRECTORY, "DataFiles")
APIS = os.path.join(HOME_DIRECTORY, "API")
PERSISTANCE_FILE_LOCATION = os.path.join(HOME_DIRECTORY, "Persistance_Files")

def GET_RELATIVE_PATH(FILE_PATH: str, start = os.getcwd()) -> str:
    return os.path.relpath(FILE_PATH, start)

def CHECK_FOLDER_EXISTS(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"New Folder '{folder_path}' created.")
        return True
    else:
        return False
    
def GET_ORGANISM_FOLDER_PATH(targetOrganism: str) -> str:
    folder_path = os.path.join(DATAFILES, targetOrganism)
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_ORGANISM_ANALYSIS_FOLDER_PATH(targetOrganism: str) -> str:
    folder_path = os.path.join(GET_ORGANISM_FOLDER_PATH(targetOrganism), "Analysis")
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_KOID_FOLDER_PATH(KOID: str, targetOrganism: str) -> str:
    folder_path = os.path.join(GET_ORGANISM_FOLDER_PATH(targetOrganism), KOID)
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_KOID_MSA_FOLDER_PATH(KOID: str, targetOrganism: str) -> str:
    folder_path = os.path.join(GET_KOID_FOLDER_PATH(KOID, targetOrganism), "MSA")
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_KOID_MSA_ANALYSIS_FILE_PATH(KOID:str, targetOrganism: str) -> str:
    return os.path.join(GET_KOID_FOLDER_PATH(KOID, targetOrganism), "MSA", f"{KOID}_MSA_RESULTS.aln")

def GET_KOID_MSA_ANALYSIS_DF_FILE_PATH(KOID:str, targetOrganism: str) -> str:
    return os.path.join(GET_KOID_FOLDER_PATH(KOID, targetOrganism), "MSA", f"Conservation_DF.csv")
    
def GET_KOID_UNIPROT_ENTRIES_PATH(KOID: str, targetOrganism: str) -> str:
    folder_path = os.path.join(GET_KOID_FOLDER_PATH(KOID, targetOrganism), "UniProt_Entries")
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_UNPROT_ENTRY_FILE_PATH(KOID: str, UniprotID: str, type: str, targetOrganism: str) -> str:
    return os.path.join(GET_KOID_UNIPROT_ENTRIES_PATH(KOID, targetOrganism), UniprotID) + type

def CHECK_PATH_EXISTS(file_path: str) -> bool:
    return os.path.exists(file_path)

def crate_target_analysis_file(targetOrganism: str):
    if targetOrganism.endswith('.csv'):
        targetOrganism = targetOrganism.replace('.csv', "")
    KOIDs = os.listdir(GET_ORGANISM_FOLDER_PATH(targetOrganism))
    all_data = []
    for KOID in KOIDs:
        if os.path.exists(GET_KOID_MSA_ANALYSIS_DF_FILE_PATH(KOID, targetOrganism)):
            data = pd.read_csv(GET_KOID_MSA_ANALYSIS_DF_FILE_PATH(KOID, targetOrganism))
            first_uniprotID = data.at[0, 'UniProtID']
            uniprot_to_koid_mapping = {first_uniprotID: KOID}
            filtered_data = data[data['UniProtID'] == first_uniprotID].copy()
            filtered_data['UniProtID'] = filtered_data['UniProtID'].map(uniprot_to_koid_mapping)
            filtered_data.rename(columns={'UniProtID': 'KOID'}, inplace=True)
            filtered_data.drop(columns=["Unnamed: 0", "MSA_Binding_Location", "Binding_Location"], inplace = True)
            all_data.append(filtered_data)
    if all_data:
        # Concatenate all DataFrames in the list into a single DataFrame
        analysisDF = pd.concat(all_data, ignore_index=True)
        
        # Ensure you're passing a full path with filename to save the CSV
        output_file_path = GET_ORGANISM_ANALYSIS_FOLDER_PATH(targetOrganism) + "/analysis_result.csv"
        analysisDF.to_csv(output_file_path, index=False)
        print(f"Analysis saved to {output_file_path}")
    else:
        print("No valid data found to analyze.")

#Concatenates all files from a specified folder
def readDirectoryContents(folder_path):
    fileList = os.listdir(folder_path)
    data = ""
    for fileName in fileList:
        filePath = os.path.join(folder_path, fileName)
        with open(filePath, 'r', encoding="utf-8") as inF:
            data += inF.read()
            data += "\n"
    return data

if __name__ == "__main__":
    crate_target_analysis_file("target_bacteria.csv")