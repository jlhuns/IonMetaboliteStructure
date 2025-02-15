import requests
import os

BASE_DOWNLOAD_URL = "https://alphafold.ebi.ac.uk/files/AF-{UNIPROTID}-F1-model_v4.pdb"

def get_pdb_file(uniprotID):
    url = BASE_DOWNLOAD_URL.replace("{UNIPROTID}", uniprotID)
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        print(f"Error getting PDB file, code: {response.status_code}")
        return
