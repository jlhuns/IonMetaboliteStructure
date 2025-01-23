#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#This will just be an array of the KO's
#May need faster search function.
import requests
import pandas as pd

# KEGG API base URL
KEGG_API_BASE_URL = "http://rest.kegg.jp"

# returns an ARRAY of the the genes in the ortholog (K0)
def get_genes_for_kolid(kolid):
    """Fetch gene information for a given KO ID from KEGG API."""
    url = f"{KEGG_API_BASE_URL}/get/{kolid}"
    try:
        # Make a GET request to the KEGG API
        response = requests.get(url)
        response.raise_for_status()  # Raises an HTTPError if the response code isn't 2xx
        
        data = response.text  # Get the response text

        # Check if the response contains the entry for KO ID (i.e., data contains "ENTRY")
        if "ENTRY" not in data:
            print(f"KO ID {kolid} not found.")
            return []

        # Parse the genes section from the response text
        genesSection = []
        isGenesSection = False
        for line in data.splitlines():
            if line.startswith("GENES"):
                isGenesSection = True
            if isGenesSection:
                if line.strip() == "///":  # End of GENES section
                    break
                # Extract the part before the colon in each line
                partBeforeColon = line.split(":")[0].strip()  # Split by ":" and get the first part
                genesSection.append(partBeforeColon)

        return genesSection  # Return the genes section

    except requests.exceptions.RequestException as e:
        print(f"Error fetching KO ID {kolid}: {e}")
        return []  # Return an empty list on error

def get_target_gene_codes(TARGET_ORGANISMS_FILEPATH):
    # filter the list of genes for the appropriate organisms (if applicable)
    organismsDf = pd.read_csv(TARGET_ORGANISMS_FILEPATH) 
    keggCodes = organismsDf['kegg_organism_code'].tolist()
    
    return(keggCodes)

def get_orthologs(kolid, TARGET_ORGANISMS_FILEPATH):
    #get the target genes
    targetGenes = get_target_gene_codes(TARGET_ORGANISMS_FILEPATH)

    #get genes from K0 id
    get_genes_for_kolid(kolid)

    #get IDs for the K0 ids



    print(targetGenes)
    return "not done yet"


def main():
    """Main function to test fetching genes for a specified KO ID."""
    ko_id = 'K00973'  # Example KO ID (can change to any KO ID you want to query)
    
    # Get the genes for the specified KO ID
    genes = get_genes_for_kolid(ko_id)
    
    # Print the result (the genes associated with the KO ID)
    if genes:
        print(f"Genes for KO ID {ko_id}:")
        for gene in genes:
            print(gene)
    else:
        print(f"No genes found for KO ID {ko_id} or there was an error.")

if __name__ == "__main__":
    main()

