from typing import List
from API import Uniprot_API
import file_paths as FILEPATH
import os
import pandas as pd
import shutil

PERSISTENCE_FILE = FILEPATH.PERSISTANCE_FILE_LOCATION + "/uniprot_entry_file_locations.csv"

def create_uniprot_entires(uniprotIDs: List[str], KOID: str, targetOrganism: str) -> None:
    if targetOrganism.endswith('.csv'):
        targetOrganism = targetOrganism.replace('.csv', "")

    if os.path.exists(PERSISTENCE_FILE):
        uniprot_entry_df = pd.read_csv(PERSISTENCE_FILE, index_col="UniProtID")
    else:
        # DataFrame with UniProtID as the index for fast lookups
        uniprot_entry_df = pd.DataFrame(columns=["UniProtID", "File_Location"]).set_index("UniProtID")

    for uniprotID in uniprotIDs:
        if uniprotID in uniprot_entry_df.index:
            # If UniProt ID already exists, print its location
            sourcePath = uniprot_entry_df.loc[uniprotID, "File_Location"]
            destinationPath = FILEPATH.GET_KOID_UNIPROT_ENTRIES_PATH(KOID, targetOrganism)

            #Check if the same file exsists in the current folder
            try:
                os.path.samefile(sourcePath, destinationPath + "\\" + uniprotID + ".txt")
                # print(f"Source and destination are the same, skipping copy: {sourcePath}")
                continue
            except Exception as e:
                print("")
            try:
                #if the file exsists elsewhere, copy that file instead of doing an API call
                shutil.copy(sourcePath, destinationPath)
                print(f"File copied successfully from {sourcePath} to {destinationPath}")
            except Exception as e:
                print(f"An error occurred: {e}")

        else:
            filePath = FILEPATH.GET_UNPROT_ENTRY_FILE_PATH(KOID, uniprotID, ".txt", targetOrganism)

            if(FILEPATH.CHECK_PATH_EXISTS(filePath)):
                print("File exsits: ", filePath)
            else:
                Uniprot_API.fetch_uniprot_file(uniprotID, FILEPATH.GET_KOID_UNIPROT_ENTRIES_PATH(KOID, targetOrganism), "txt")

            # Add the new entry to the DataFrame
            relative_file_path = FILEPATH.GET_RELATIVE_PATH(filePath)
            uniprot_entry_df.loc[uniprotID] = relative_file_path

    uniprot_entry_df.to_csv(PERSISTENCE_FILE)
