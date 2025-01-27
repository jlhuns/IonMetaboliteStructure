import os
import re
import pandas as pd
import file_paths as FILE_PATH

CWD = os.getcwd()
PARENT_DIR = os.path.dirname(CWD)

# PATH_TO_DATAFILES = os.path.join(PARENT_DIR, "datafiles")
# PATH_TO_UNIPROT_ENTRIES = os.path.join(PATH_TO_DATAFILES, "uniprot_entries")
# PATH_TO_EVAL_FILES = os.path.join(PATH_TO_DATAFILES, "eval_files")
# PATH_TO_MSA_FILE = os.path.join(PATH_TO_DATAFILES, "muscle_data", "alignment.aln")

def convert_uniprot_data_to_position(entry):
    value_str = entry[2]
    # Check if the value string contains periods
    if re.search(r'\.+', value_str):
        # Split the string using regular expression to handle any number of periods
        start_str, end_str = re.split(r'\.+', value_str)

        # Convert the start and end strings to integers
        start = int(start_str)
        end = int(end_str)

        # Generate the range of values
        return list(range(start, end + 1))
    else:
        # Convert the string to an integer
        return int(value_str)


def read_uniprot_file_to_analyze_active_sites(directory, filename):
    data = []
    act_str = 'ACT_SITE'
    bind_str = 'BINDING'
    for file in os.listdir(directory):
        if file == filename:
            filepath = os.path.join(directory, filename)
            with open(filepath, 'r') as inF:
                lines = inF.readlines()
                for i in range(len(lines) - 1):  # Iterate through lines except the last one
                    line = lines[i].strip()
                    next_line = lines[i + 1].strip()

                    if line.startswith('ID'):
                        uniProt_ID = line.split()[1].replace(';', '')
                        # print(uniProt_ID)

                    if line.startswith('FT'):
                        entry = line.split()               
                        # Process the current FT line (active site or binding site)
                        if entry[1] == act_str:
                            position = convert_uniprot_data_to_position(entry)
                            # if isinstance(position, list):
                            #     print(position[0], position[-1])
                            data.append({
                                "UniProt_ID": uniProt_ID,
                                "Type": act_str,
                                "Position": position, #if there are two numbers like [450, 455] then the active sites goes from 450 to 455
                                "Conservation_Score": "TBD",
                                "Description": re.search(r'="([^"]+)"', ' '.join(next_line.split())).group(1)
                            })                         
                        if entry[1] == bind_str:
                            position = convert_uniprot_data_to_position(entry)
                            data.append({
                                "UniProt_ID": uniProt_ID,
                                "Type": bind_str,
                                "Position": position,
                                "Conservation_Score": "TBD",
                                "Description": re.search(r'="([^"]+)"', next_line).group(1)
                            })
    activeSitesDF = pd.DataFrame(data)
    # print(activeSitesDF)
    return activeSitesDF

def FindActiveSitesInMSA(activeSitesDF, MSAFile):
    # print(activeSitesDF)
    def find_char_position(full_string: str, non_hyphen_index) -> int:
        non_hyphen_count = 0
        for pos, char in enumerate(full_string):
            if char != '-' and char != "\n":
                if non_hyphen_count == non_hyphen_index:
                    return pos - 1 #base 0 indexing
                non_hyphen_count += 1
        # print(sequence)
        raise ValueError("Index out of range for non-hyphenated string.")
    
    with open(MSAFile, 'r') as inF:
        lines = inF.readlines()
        uniProtIDs = set(activeSitesDF.loc[:, "UniProt_ID"])
        uniProtIDs = list(uniProtIDs)
        results = []
        
        for uniProtID in uniProtIDs:
            conservationScoreString = ""
            sequence = ""
            charactersToRemove = 0
            

            for i in range(len(lines)):
                line = lines[i]
                if(line.startswith(uniProtID)):
                    newLineCount = line.count("\n")
                    spaces = len([char for char in line if char.isspace()])
                    match = re.search(r'\S+', line)
                    charactersToRemove = spaces + match.end() - newLineCount
                    line = line[charactersToRemove:]
                    line = line.replace(" ", "")
                    sequence += line

                if(line.startswith(" ")):
                    line = line[charactersToRemove:]
                    conservationScoreString += (line)

            filtered_df = activeSitesDF[activeSitesDF['UniProt_ID'].str.startswith(uniProtID)]
            for index, row in filtered_df.iterrows():
                if row['Type'] == "BINDING":
                    position = find_char_position(sequence, row['Position'])  # Define or import this function
                    #handles edge case where if conservation score is the last char then it returns "\n" but going back one will give the correct conservation score.
                    if conservationScoreString[position] == "\n":
                        filtered_df.loc[index, 'Conservation_Score'] = conservationScoreString[position-1]
                    else:
                        filtered_df.loc[index, 'Conservation_Score'] = conservationScoreString[position]



            results.append(filtered_df)
            # print(results)
        final_results = pd.concat(results, ignore_index=True)
        return final_results

def Create_Analysis_DF(KOID: str, targetOrganism: str):
    if targetOrganism.endswith('.csv'):
        targetOrganism = targetOrganism.replace('.csv', "")
    
    PATH_TO_UNIPROT_ENTRIES = FILE_PATH.GET_KOID_UNIPROT_ENTRIES_PATH(KOID, targetOrganism)
    PATH_TO_MSA_FILE = FILE_PATH.GET_KOID_MSA_ANALYSIS_FILE_PATH(KOID, targetOrganism)

    file_names = [file for file in os.listdir(PATH_TO_UNIPROT_ENTRIES) if os.path.isfile(os.path.join(PATH_TO_UNIPROT_ENTRIES, file))]
    dataFrames = []
    for file in file_names:
        dataFrames.append(read_uniprot_file_to_analyze_active_sites(PATH_TO_UNIPROT_ENTRIES, file))
    final_dataframe = pd.concat(dataFrames, ignore_index=True)

    ActiveSitesDF = FindActiveSitesInMSA(final_dataframe, PATH_TO_MSA_FILE)
    ActiveSitesDF.to_csv(os.path.join(FILE_PATH.GET_KOID_MSA_FOLDER_PATH(KOID, targetOrganism), "ActiveSitesDF"))