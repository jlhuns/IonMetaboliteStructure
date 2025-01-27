#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#Use UNIPROT ID to get uniprot Entries (CREATE FILES UNDER DATAFILES MSA UNIPROT)
#Make a fasta file using UNIPROT ENTRIES (ADD FILE UNDER KOID MSA)
#RUN MSA (ADD file to MSA under correct KOID)
#Analyze MSA (binding or not, return DF )

from Ortholog_Identificiation import create_uniprot_entries, get_orthologs_uniprotID
from MSA import analyze_active_sites, run_msa
import time



if __name__ == "__main__":
    start_time = time.time()
    
    uniprot_ids = get_orthologs_uniprotID.get_uniprot_ids("K00937", "target_prokaryotes.csv")
    create_uniprot_entries.create_uniprot_entires(uniprot_ids, "K00937", "target_prokaryotes.csv")
    run_msa.create_msa_file("K00937", "target_prokaryotes.csv")
    analyze_active_sites.Create_Analysis_DF("K00937", "target_prokaryotes.csv")
    
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"The program ran for {elapsed_time:.6f} seconds.")

