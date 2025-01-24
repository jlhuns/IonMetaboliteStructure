import os

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
    
def GET_ORAGANISM_FOLDER_PATH(targetOrganism: str) -> str:
    folder_path = os.path.join(DATAFILES, targetOrganism)
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_KOID_FOLDER_PATH(KOID: str, targetOrganism: str) -> str:
    folder_path = os.path.join(GET_ORAGANISM_FOLDER_PATH(targetOrganism), KOID)
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_KOID_MSA_PATH(KOID: str, targetOrganism: str) -> str:
    folder_path = os.path.join(GET_KOID_FOLDER_PATH(KOID, targetOrganism), "MSA")
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_KOID_UNIPROT_ENTRYS_PATH(KOID: str, targetOrganism: str) -> str:
    folder_path = os.path.join(GET_KOID_FOLDER_PATH(KOID, targetOrganism), "UniProt_Entries")
    CHECK_FOLDER_EXISTS(folder_path)
    return folder_path

def GET_UNPROT_ENTRY_FILE_PATH(KOID: str, UniprotID: str, type: str, targetOrganism: str) -> str:
    return os.path.join(GET_KOID_UNIPROT_ENTRYS_PATH(KOID, targetOrganism), UniprotID) + type

def CHECK_PATH_EXISTS(file_path: str) -> bool:
    return os.path.exists(file_path)