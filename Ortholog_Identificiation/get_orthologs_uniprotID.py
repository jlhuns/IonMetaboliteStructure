#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#This will just be an array of the KO's
#May need faster search function.
import requests
import re
import pandas as pd
from Bio.KEGG import REST

# KEGG API base URL
KEGG_API_BASE_URL = "http://rest.kegg.jp"

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

# def get_uniprot_ids_by_gene(gene_name):
#     gene_response = REST.kegg_get(gene_name).read()
#     # Regular expression pattern to match UniProt ID
#     pattern = r'UniProt:\s+(\S+)'

#     # Extracting UniProt ID using regular expression
#     match = re.search(pattern, gene_response)
#     if match:
#         uniprot_id = match.group(1)
#     else:
#         uniprot_id = None

#     # print(f'gene: {gene_name} --> uniprot_id: {uniprot_id}') # TODO: remove after testing
#     return uniprot_id

def get_uniprot_from_kegg(gene,kegg_id):
    # KEGG REST API endpoint
    #url = f"http://rest.kegg.jp/get/{kegg_id}"
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
def get_unipro_ids(koID, TARGET_FILEPATH):
    similarGenes, dictIdAndIdentifier = get_similar_gene_ids(koID,TARGET_FILEPATH)

    orthologKOIDs = []

    for gene in similarGenes:
        #gene = gene.upper()
        id = get_uniprot_from_kegg(gene, dictIdAndIdentifier[gene.upper()])
       
        if not id.startswith("Error"):
            orthologKOIDs.append(id)
        
    return orthologKOIDs

def main():
    """Main function to test fetching genes for a specified KO ID."""
    ko_id = 'K00973'  # Example KO ID (can change to any KO ID you want to query)
    
    #This must be the relative file path
    targetFilePath = "/Users/kristinbillings/Our Structure data/IonMetaboliteStructure/target_prokaryotes.csv"
    
    # Get the genes for the specified KO ID
    genes = get_genes_for_kolid(ko_id)
    
    # Print the result (the genes associated with the KO ID)
    # if genes:
    #     print(f"Genes for KO ID {ko_id}:")
    #     for gene in genes:
    #         print(gene)
    # else:
    #     print(f"No genes found for KO ID {ko_id} or there was an error.")

    #mylist, yes = get_similar_gene_ids(ko_id, targetFilePath)
    #print(mylist)
    yeet = get_unipro_ids(ko_id, targetFilePath)
    print(yeet)

if __name__ == "__main__":
    main()

