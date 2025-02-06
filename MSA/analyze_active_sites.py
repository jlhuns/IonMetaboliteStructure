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
                        if(line.split()[1] == "BINDING" or line.split()[1] == "ACT_SITE" or line.split()[1] == "SITE"):
                            siteType = (line.split()[1])
                            location = (line.split()[2])
                            quotes_match = re.findall(r'"([^"]+)"', nextLine)

                            if(nextLine.count('"') == 1):
                                quotes_match = re.findall(r'(?<=\").*', nextLine)
                                locationDescription = quotes_match[0] + "..." if quotes_match else ""
                                
                            else:
                                locationDescription = quotes_match[0] if quotes_match else ""

                            data.append({
                                "UniProtID": proteinKOID,
                                "Type": siteType,
                                "MSA_Binding_Location": "",
                                "Binding_Location": location,
                                "conservationScore": conservationScore,
                                "Value": "",
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
            return '-'  # No significant conservation
        
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


        seq_id = row["UniProtID"]
        seq_record = next((record for record in alignment if record.id == seq_id), None)
        binding_location = row["Binding_Location"]

        pattern = r'(\d+)\.\.(\d+)'
        if re.match(pattern, binding_location):
            conservation_score_data = ""
            values = ""
            binding_location_start = int(re.match(pattern, binding_location).group(1))
            binding_location_end = int(re.match(pattern, binding_location).group(2))

            MSA_index_start, value = get_positions(binding_location_start, seq_record.seq)
            MSA_index_end, value = get_positions(binding_location_end, seq_record.seq)
            residue_counts = []
            for MSA_index in range(MSA_index_start, MSA_index_end + 1):
                column = [record.seq[MSA_index] for record in alignment]
                values += seq_record.seq[MSA_index]
                conservation_score_data += (clustal_symbol(column))
                residue_counts.append(dict(Counter(column)))

            

            analysisDF.at[index, 'conservationScore'] = conservation_score_data
            analysisDF.at[index, 'Value'] = values
            analysisDF.at[index, 'Residue_Counts'] = str(residue_counts)
            analysisDF.at[index, 'MSA_Binding_Location'] = f'{MSA_index_start}..{MSA_index_end}'
            
        else:
            binding_location_int = int(binding_location)

            # Get position and conservation score for this single position
            MSA_index, value = get_positions(binding_location_int, seq_record.seq)
            column = [record.seq[MSA_index] for record in alignment]
            conservation_score = clustal_symbol(column)

            residue_counts = dict(Counter(column))
            # Update the conservation score in the DataFrame
            analysisDF.at[index, 'conservationScore'] = conservation_score
            analysisDF.at[index, 'Value'] = value
            analysisDF.at[index, 'Residue_Counts'] = str(residue_counts)
            analysisDF.at[index, 'MSA_Binding_Location'] = f'{MSA_index}'
    return analysisDF


def analyze_MSA(KOID: str, target_organism: str):
    if target_organism.endswith('.csv'):
        target_organism = target_organism.replace('.csv', "")
    if not (os.path.exists(os.path.join(FILE_PATH.GET_KOID_MSA_FOLDER_PATH(KOID, target_organism), f"{KOID}_MSA_Results.aln"))):
       print(f"MSA file does not exist, skipping MSA analysis")
       return
    if target_organism.endswith('.csv'):
        target_organism = target_organism.replace('.csv', "")
    
    analysisDF = create_analysis_df(KOID, target_organism)
    resultsDF = get_conservation_score(analysisDF, KOID, target_organism)    
    if(resultsDF.empty):
        print(f"analysisDF is empty for K0: {KOID}. No Binding or Active Sites Found")
        with open(os.path.join(FILE_PATH.GET_KOID_MSA_FOLDER_PATH(KOID, target_organism), "EmptyDF.txt"), 'w') as file:
            file.write(f"analysisDF is empty for K0: {KOID}. No Binding or Active Sites Found")
        return resultsDF
    clear_results = resultsDF.drop_duplicates(subset="MSA_Binding_Location")
    clear_results = clear_results.drop("Binding_Location", axis=1)
    clear_results.to_csv(os.path.join(FILE_PATH.GET_KOID_MSA_FOLDER_PATH(KOID, target_organism), "Simple_Conservation_DF.csv"))
    resultsDF.to_csv(os.path.join(FILE_PATH.GET_KOID_MSA_FOLDER_PATH(KOID, target_organism), "Conservation_DF.csv"))


if __name__ == "__main__":
    analyze_MSA("K00031", "target_bacteria.csv")