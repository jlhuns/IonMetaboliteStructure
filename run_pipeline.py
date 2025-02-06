#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#Use UNIPROT ID to get uniprot Entries (CREATE FILES UNDER DATAFILES MSA UNIPROT)
#Make a fasta file using UNIPROT ENTRIES (ADD FILE UNDER KOID MSA)
#RUN MSA (ADD file to MSA under correct KOID)
#Analyze MSA (binding or not, return DF )
import tempfile
from Ortholog_Identificiation import create_uniprot_entries, get_orthologs_uniprotID
from MSA import analyze_active_sites, run_msa
from Analysis import conservation_analysis
import time
import sys
import file_paths as FILE_PATHS

def main(KOID_file, target_organism_file):
    start_time = time.time()

        # Create a temporary file for logging errors
    with tempfile.NamedTemporaryFile(mode='w+', delete=False) as temp_log:
        log_file_path = temp_log.name  # Store the path of the log file
        temp_log.write("Error Log Started\n")  # Optional header for the log file

    with open(KOID_file, 'r') as inF:
        lines = inF.readlines()

        for KOID_line in lines:
            KOID = KOID_line.replace("\n", "")
            try:
                # Attempt each step, and skip to the next KOID if any error occurs
                uniprot_ids = get_orthologs_uniprotID.get_uniprot_ids(KOID, target_organism_file)
                create_uniprot_entries.create_uniprot_entires(uniprot_ids, KOID, target_organism_file)
                run_msa.create_msa_file(KOID, target_organism_file)
                analyze_active_sites.analyze_MSA(KOID, target_organism_file)
            except Exception as e:
                # Log the error and skip to the next KOID
                temp_log.write(f"Error processing KOID {KOID}: {e}\n")
                continue  # Continue to the next KOID if an error occurs
        FILE_PATHS.crate_target_analysis_file(target_organism_file)


    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"The program ran for {elapsed_time:.6f} seconds.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Please provide at least two arguments.")

