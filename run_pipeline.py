#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#Use UNIPROT ID to get uniprot Entries (CREATE FILES UNDER DATAFILES MSA UNIPROT)
#Make a fasta file using UNIPROT ENTRIES (ADD FILE UNDER KOID MSA)
#RUN MSA (ADD file to MSA under correct KOID)
#Analyze MSA (binding or not, return DF )

from Ortholog_Identificiation import create_uniprot_entries, get_orthologs_uniprotID
from MSA import analyze_active_sites, run_msa
import time
import sys

def main(KOID_file, target_organism_file):
    start_time = time.time()

    with open(KOID_file, 'r') as inF:
        lines = inF.readlines()

        for KOID_line in lines:
            KOID = KOID_line.replace("\n", "")
            # uniprot_ids = get_orthologs_uniprotID.get_uniprot_ids(KOID, target_organism_file)
            # create_uniprot_entries.create_uniprot_entires(uniprot_ids, KOID, target_organism_file)
            # run_msa.create_msa_file(KOID, target_organism_file)
            analyze_active_sites.analyze_MSA(KOID, target_organism_file) 

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"The program ran for {elapsed_time:.6f} seconds.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Please provide at least two arguments.")

