from typing import List
import requests

#URL FORM = https://rest.kegg.jp/<operation>/<argument>[/<argument2[/<argument3> ...]]
BASE_URL = "https://rest.kegg.jp"

"""
Output format
The output of all operations is in a text format:
tab-delimited text returned from list, find, conv and link
flat file database format returned from get
text message returned from info
Status code
The HTTP status code can be used to check if the operation was successful.

Code  	Meaning
200	Success
400	Bad request (syntax error, wrong database name, etc.)
404	Not found (e.g., requesting amino acid sequence for RNA)
"""


def info(database):
    databases = [
        "kegg", "pathway", "brite", "module",  "ko",  "genes", "vg", "vp", "ag",
        "genome", "ligand", "compound", "glycan", "reaction", "rclass", "enzyme",
        "network", "variant", "disease", "drug", "dgroup"
    ]
    if database in databases:
        url = f"{BASE_URL}/info/{database}"
        response = requests.get(url)
        if response.status_code == 200:
            return response.text
        else:
            return f'Error: {response.status_code}'
    else:
        return f'Database: {database} not found in {databases}'

def list_entries(database, *args):
    """
    Obtain a list of entries from a specific KEGG database.
    
    :param database: KEGG database to query (e.g., 'pathway', 'ko', 'gene', etc.)
    :param args: Additional optional arguments like organism code or database-specific options
    :return: List of entries in the specified database
    """
    valid_databases = [
        "pathway", "brite", "module", "ko", "vg", "vp", "ag", "genome", "compound",
        "glycan", "reaction", "rclass", "enzyme", "network", "variant", "disease",
        "drug", "dgroup", "organism"
    ]
    
    if database not in valid_databases:
        return f"Invalid database '{database}'"
    
    # Form URL based on database and optional arguments
    url = f"{BASE_URL}/list/{database}"
    
    if args:
        url += f"/{'/'.join(args)}"
    
    response = requests.get(url)
    
    if response.status_code == 200:
        return response.text
    else:
        return f"Error: {response.status_code}"
    
def find_genes(KOID:str):
    url = f'{BASE_URL}/find/genes/{KOID}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return f'Error: {response.status_code}'

    

def convert_KeggID_to_UniprotID(GeneIdentifier: List[str], DataBase = "uniprot"):
    genes = "+".join(GeneIdentifier)
    url = f'{BASE_URL}/conv/{DataBase}/{genes}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return f'Error: {response.status_code}'

if __name__ == "__main__":
    # Test examples
    genes = ["tma:TM0862", "aeo:O23A_p0654", "dra:DR_A0031"]
    print(find_genes("k00937"))
