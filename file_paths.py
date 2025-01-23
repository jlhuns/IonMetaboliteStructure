import os

HOME_DIRECTORY = os.getcwd()
DATAFILES = os.path.join(HOME_DIRECTORY, "DataFiles")
APIS = os.path.join(HOME_DIRECTORY, "API")

def check_Folder_Exsits(folder_path):
    # Check if the folder exists
    if not os.path.exists(folder_path):
        # Create the folder if it doesn't exist
        os.makedirs(folder_path)
        print(f"New Folder '{folder_path}' created.")
        return True
    else:
        return False

def GET_KOID_FOLDER_PATH(KOID: str) -> str:
    folder_path = os.path.join(DATAFILES, KOID)
    check_Folder_Exsits(folder_path)
    return folder_path

def GET_KOID_MSA_PATH(KOID: str) -> str:
    folder_path = os.path.join(GET_KOID_FOLDER_PATH(KOID), "MSA")
    check_Folder_Exsits(folder_path)
    return folder_path

def GET_KOID_UNIPROT_ENTRYS_PATH(KOID: str) -> str:
    folder_path = os.path.join(GET_KOID_FOLDER_PATH(KOID), "UniProt_Entries")
    check_Folder_Exsits(folder_path)
    return folder_path




if __name__ == "__main__":
    print(HOME_DIRECTORY)
    print(DATAFILES)
    print(getKOID("K00937"))
