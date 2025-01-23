#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#This will just be an array of the KO's
#May need faster search function.
import requests

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
        genes_section = []
        is_genes_section = False
        for line in data.splitlines():
            if line.startswith("GENES"):
                is_genes_section = True
            if is_genes_section:
                if line.strip() == "///":  # End of GENES section
                    break
                genes_section.append(line.strip())

        return genes_section  # Return the genes section

    except requests.exceptions.RequestException as e:
        print(f"Error fetching KO ID {kolid}: {e}")
        return []  # Return an empty list on error

def get_target_gene_codes():
    # filter the list of genes for the appropriate organisms (if applicable)
    organisms_df = pd.read_csv(target_organisms_filepath) 
    ortholog_genes_df['kegg_organism_code'] = ortholog_genes_df.apply(get_organism_code_from_row, axis=1)
    filtered_genes = ortholog_genes_df.loc[ortholog_genes_df['kegg_organism_code'].isin(organisms_df['kegg_organism_code'])]

def get_orthologs(kolid):
    #get the target genes



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

