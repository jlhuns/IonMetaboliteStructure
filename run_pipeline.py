from Ortholog_Identificiation import create_uniprot_entries, get_orthologs_uniprotID
from MSA import analyze_active_sites, run_msa
from Analysis import conservation_analysis
from FileCare import CleanKOIDs

import time
import sys
import file_paths as FILE_PATHS

from API import Kegg_API
import re
import os
def load_seen_koids(filename="seenKOIDS.txt"):
    """Load previously seen KOIDs from a file."""
    if not os.path.exists(filename):
        return set()
    with open(filename, "r") as f:
        return set(line.strip() for line in f)

def save_seen_koid(koid, filename="seenKOIDS.txt"):
    """Append a new KOID to the file."""
    with open(filename, "a") as f:
        f.write(koid + "\n")

def remove_seen_koids_from_file(KOID_FILE, seenKOIDS):
    """Remove all seen KOIDs from the KOID_FILE."""
    with open(KOID_FILE, 'r') as f:
        lines = f.readlines()

    with open(KOID_FILE, 'w') as f:
        for line in lines:
            # Filter out KOIDs that are in seenKOIDS
            new_line = re.sub(r"'(K\d+)'", lambda m: '' if m.group(1) in seenKOIDS else m.group(0), line)
            if new_line.strip():  # Don't write empty lines back to the file
                f.write(new_line)

def main(KOID_FILE, target_organism_file):
    start_time = time.time()
    seenKOIDS = load_seen_koids()
    all_koids = []

    # Regular expression to match KOIDs like 'K12345'
    koid_pattern = re.compile(r"'(K\d+)'")
    with open(KOID_FILE, 'r') as inF:
        for line in inF:
            # Find all KOIDs in the current line
            koids_in_line = koid_pattern.findall(line)
            all_koids.extend(koids_in_line)  # Add found KOIDs to the list

    # all_koids will now contain all the KOIDs from the entire file
    seenKOIDS = []
    for KOID in all_koids:
        try:
            seenKOIDS = load_seen_koids()

            folder_path = os.path.join(FILE_PATHS.GET_ORGANISM_FOLDER_PATH(target_organism_file).replace('.csv', ""), KOID)
            if KOID not in seenKOIDS:
                seenKOIDS.add(KOID)
                save_seen_koid(KOID)
            else:
                print(f"Skipping {KOID}: Already Seen KOID.")
                analyze_active_sites.analyze_MSA(KOID, target_organism_file) #Temporary Line
                continue
            
            if os.path.exists(folder_path) and os.path.isdir(folder_path):
                print(f"Skipping {KOID}: Folder already exists.")
                continue
            # Attempt each step, and skip to the next KOID if any error occurs
            uniprot_ids = get_orthologs_uniprotID.get_uniprot_ids(KOID, target_organism_file)
            if len(uniprot_ids) < 20:
                continue
            create_uniprot_entries.create_uniprot_entires(uniprot_ids, KOID, target_organism_file)
            run_msa.create_msa_file(KOID, target_organism_file)
            analyze_active_sites.analyze_MSA(KOID, target_organism_file)
        except:
            print(f"{KOID} ran into an error")
            continue
    #After each KOID runs
    CleanKOIDs.remove_empty_KOID(target_organism_file)
    FILE_PATHS.crate_target_analysis_file(target_organism_file)

        

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"The program ran for {elapsed_time:.6f} seconds.")


if __name__ == "__main__":
    if len(sys.argv) > 2:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Please provide at least two argument.")
