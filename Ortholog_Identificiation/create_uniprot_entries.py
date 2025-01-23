from typing import List
from API import Uniprot_API
import file_paths as FILEPATH
import os
#To test run in root directory using: python -m Ortholog_Identificiation.create_uniprot_entries

def create_uniprot_entires(uniprotIDs: List[str], KOID: str):
    for uniprotID in uniprotIDs:
        filePath = FILEPATH.GET_UNPROT_ENTRY_FILE_PATH(KOID, uniprotID, ".txt")
        if(FILEPATH.CHECK_FILE_EXISTS(filePath)):
            print("File exsits: ", filePath)
        else:
            Uniprot_API.fetch_uniprot_file(uniprotID, FILEPATH.GET_KOID_UNIPROT_ENTRYS_PATH(KOID), "txt")
