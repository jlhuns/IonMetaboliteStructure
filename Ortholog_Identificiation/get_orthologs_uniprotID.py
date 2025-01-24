#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#This will just be an array of the KO's
#May need faster search function.
import requests
import re
import pandas as pd
from Bio.KEGG import REST
import os
import file_paths as FILEPATH

# KEGG API base URL
KEGG_API_BASE_URL = "http://rest.kegg.jp"
PERSISTENCE_FILE = FILEPATH.PERSISTANCE_FILE_LOCATION + "/orthologs_uniprotID.csv"

# returns an DICTIONARY of the the genes in the ortholog (K0) and their identifier number thing
def get_genes_for_kolid(kolid):
    """Fetch gene information for a given KO ID from KEGG API."""
    url = f"{KEGG_API_BASE_URL}/get/{kolid}"
    try:
        # Make a GET request to the KEGG API
        response = requests.get(url)
        response.raise_for_status()  # Raises an HTTPError if the response code isn't 2xx
        
        data = response.text  # Get the response text

        # Check if the response contains the entry for KO ID (i.e., d ata contains "ENTRY")
        if "ENTRY" not in data:
            print(f"KO ID {kolid} not found.")
            return []

        # Parse the genes section from the response text
        genesSection = {}
        isGenesSection = False
        for line in data.splitlines():
            if line.startswith("GENES"):
                isGenesSection = True
            if isGenesSection:
                pattern = r':\s*([^(\s]+)'
                if line.strip() == "///":  # End of GENES section
                    break
                # Extract the part before the colon in each line
                partBeforeColon = line.split(":")[0].strip()  # Split by ":" and get the first part
                partAfterCollon = ""
                match = re.search(pattern, line)
                if match:
                    partAfterCollon = match.group(1)
                genesSection[partBeforeColon] = partAfterCollon
    
        return genesSection  # Return the genes section

    except requests.exceptions.RequestException as e:
        print(f"Error fetching KO ID {kolid}: {e}")
        return []  # Return an empty list on error

def get_target_gene_codes(TARGET_ORGANISMS_FILEPATH):
    # filter the list of genes for the appropriate organisms (if applicable)
    organismsDf = pd.read_csv(TARGET_ORGANISMS_FILEPATH) 
    keggCodes = organismsDf['kegg_organism_code'].tolist()

    return(keggCodes)

def get_similar_gene_ids(kolid, TARGET_ORGANISMS_FILEPATH):
    #get the target genes dict
    targetGenes = get_target_gene_codes(TARGET_ORGANISMS_FILEPATH)

    #get genes IDS from K0 id
    koGenes = get_genes_for_kolid(kolid)

    koGeneIDs = koGenes.keys()

    #make sure all lowercase for comparison and change ko IDs to a set
    koSet = set(gene.lower() for gene in koGeneIDs)  
    targetGenesLower = [gene.lower() for gene in targetGenes]

    return [gene for gene in targetGenesLower if gene in koSet], koGenes

def get_uniprot_from_kegg(gene,kegg_id):
    # KEGG REST API endpoint
    url = f"https://rest.kegg.jp/get/{gene}:{kegg_id}"
    # Send a GET request to the endpoint
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Read the response data
        data = response.text
        match = re.search(r"UniProt:\s*(\S+)", data)
        if match:
            return match.group(1)
        else:
            return f"Error: No {gene}:{kegg_id} match to Uniprot ID"
    else:
        return f"Error: {response.status_code} - Unable to retrieve data"

#returns an ARRAY of unipro IDs
def get_uniprot_ids(koID, TARGET_FILEPATH):
    if os.path.exists(PERSISTENCE_FILE):
        kegg_entry_df = pd.read_csv(PERSISTENCE_FILE, index_col="gene")
    else:
        # DataFrame with UniProtID as the index for fast lookups
        kegg_entry_df = pd.DataFrame(columns=["gene","gene_identifier", "UniprotID"]).set_index("gene")

    similarGenes, dictIdAndIdentifier = get_similar_gene_ids(koID,TARGET_FILEPATH)


    orthologKOIDs = []

    for gene in similarGenes:
        if gene in kegg_entry_df.index:
            uniprotID = kegg_entry_df.loc[gene, "UniprotID"]
            if not uniprotID.startswith("Error"):
                print(uniprotID)
                orthologKOIDs.append(uniprotID)
        else:
            gene_identifier = dictIdAndIdentifier[gene.upper()]
            id = get_uniprot_from_kegg(gene, gene_identifier)
            kegg_entry_df.loc[gene] = {"gene_identifier": gene_identifier, "UniprotID": id}
        
            if not id.startswith("Error"):
                orthologKOIDs.append(id)
    
    kegg_entry_df.to_csv(PERSISTENCE_FILE)        
    return orthologKOIDs
