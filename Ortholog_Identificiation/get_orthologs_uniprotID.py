#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#This will just be an array of the KO's
#May need faster search function.
from io import StringIO
from typing import List
import requests
import re
import pandas as pd
from Bio.KEGG import REST
import os
import file_paths as FILEPATH
from API import Kegg_API

# KEGG API base URL
KEGG_API_BASE_URL = "http://rest.kegg.jp"
PERSISTENCE_FILE = FILEPATH.PERSISTANCE_FILE_LOCATION + "/orthologs_uniprotID.csv"

def get_genes(KOID: str) -> List[str]:
    print("Getting KOID Genes")
    data = Kegg_API.find_genes(KOID)
    pattern = r"^\S+"
    genes = re.findall(pattern, data, flags=re.MULTILINE) 
    #about 20 seconds
    return genes

def find_common_genes(genes, TARGET_ORGANISMS_FILEPATH):
    print("Getting common Genes")
    organismsDf = pd.read_csv(TARGET_ORGANISMS_FILEPATH) 
    keggCodes = organismsDf['kegg_organism_code'].tolist()
    common_genes = [gene for gene in genes if any(gene.split(':')[0] == kegg_code for kegg_code in keggCodes)]
    return common_genes

def find_uniprot_ids(genes: List[str]):
    print("Finding Uniprot Ids for common genes")
    data = Kegg_API.convert_KeggID_to_UniprotID(genes)
    pattern = r'up:(\S+)'
    uniprot_ids = re.findall(pattern, data, flags=re.MULTILINE) 
    return uniprot_ids

def get_uniprot_ids(KOID, target_organism):
    genes = get_genes(KOID)
    common_genes = find_common_genes(genes, target_organism)
    print(f"Number of common genes found: {len(common_genes)}")
    uniprot_ids = find_uniprot_ids(common_genes)
    print(f"For {KOID}, returning uniprot IDs : {uniprot_ids}")
    return uniprot_ids