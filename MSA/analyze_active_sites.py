import os
import re
from typing import Counter
import pandas as pd
import file_paths as FILE_PATH
from Bio import AlignIO

CWD = os.getcwd()
PARENT_DIR = os.path.dirname(CWD)

def create_analysis_df(KOID: str, target_organism: str):
    uniprotEntryFolder = FILE_PATH.GET_KOID_UNIPROT_ENTRIES_PATH(KOID, target_organism)
    data = []

    for file_name in os.listdir(uniprotEntryFolder):
        file_path = os.path.join(uniprotEntryFolder, file_name)
        proteinKOID = ""
        siteType = ""
        location = ""
        conservationScore = "TBD"
        locationDescription = ""

        if os.path.isfile(file_path):
            with open(file_path, 'r') as inF:
                lines = inF.readlines()
                for i, line in enumerate(lines[:-1]):
                    nextLine = lines[i + 1] if i + 1 < len(lines) else ""
                    if(line.startswith("ID")):
                        proteinKOID = line.split()[1]
                    elif(line.startswith("FT")):
                        if(line.split()[1] == "BINDING" or line.split()[1] == "ACT_SITE"):
                            siteType = (line.split()[1])
                            location = (line.split()[2])
                            quotes_match = re.findall(r'"([^"]+)"', nextLine)
                            locationDescription = quotes_match[0] if quotes_match else ""

                            data.append({
                                "KOID": proteinKOID,
                                "Type": siteType,
                                "Binding_Location": location,
                                "conservationScore": conservationScore,
                                "Description": locationDescription
                            })
    analysisDF = pd.DataFrame(data)
    return analysisDF

def get_conservation_score(analysisDF, KOID: str, target_organism: str):
    # Function to determine Clustal-style conservation symbols
    def clustal_symbol(column):
        counts = Counter(column)
        most_common_residue, max_count = counts.most_common(1)[0]
        
        if max_count == len(column):
            return '*'  # Fully conserved
        elif max_count >= len(column) * 0.8:  # Threshold for strong similarity
            return ':'
        elif max_count >= len(column) * 0.6:  # Threshold for weak similarity
            return '.'
        else:
            return ' '  # No significant conservation
        
    def get_positions(seq_position, seq):
        residue_index = 0
        for MSA_index, value in enumerate(seq):
            if value not in ("-", "\n", " "):
                residue_index += 1
            if(residue_index == seq_position):
                return MSA_index, value     
            
    MSA_RESULTS_FILE_PATH = FILE_PATH.GET_KOID_MSA_ANALYSIS_FILE_PATH(KOID, target_organism)
    alignment = AlignIO.read(MSA_RESULTS_FILE_PATH, "clustal")

    for index, row in analysisDF.iterrows():


        seq_id = row["KOID"]  # Replace with the actual sequence ID
        seq_record = next((record for record in alignment if record.id == seq_id), None)
        binding_location = row["Binding_Location"]

        try:
            # If it converts successfully, treat it as a single value
            binding_location_int = int(binding_location)
            # print(f"Processing single binding location: {binding_location_int}")

            # Get position and conservation score for this single position
            MSA_index, value = get_positions(binding_location_int, seq_record.seq)
            column = [record.seq[MSA_index] for record in alignment]
            conservation_score = clustal_symbol(column)

            # Update the conservation score in the DataFrame
            analysisDF.at[index, 'conservationScore'] = conservation_score

        except ValueError:
            # If conversion fails, it's likely a range (e.g., "15..23")
            print(f"Binding location is a range or non-numeric: {binding_location}")
    
    print(analysisDF.to_string())





    

                



def analyze_MSA(KOID: str, target_organism: str):
    if target_organism.endswith('.csv'):
        target_organism = target_organism.replace('.csv', "")

    analysisDF = create_analysis_df(KOID, target_organism)
    get_conservation_score(analysisDF, KOID, target_organism)
