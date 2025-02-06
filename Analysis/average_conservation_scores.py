import file_paths as FILE_PATH
import pandas as pd
import os

def average_conservation_score(target_organism):
    target_organism_path = FILE_PATH.GET_ANALYSIS_RESULT_FILE(target_organism)
    conservation_df = pd.read_csv(target_organism_path)
    conservation_df = conservation_df.drop_duplicates(subset = "MSA_Binding_Location")
    descriptionAndScore = {}

    for index, row in conservation_df.iterrows():
        if row["Description"] not in descriptionAndScore:
            scoreAndCount = {}
            scoreAndCount[row["conservationScore"]] = 1
            descriptionAndScore[row["Description"]] = scoreAndCount
        else:
            if row["conservationScore"] in descriptionAndScore[row["Description"]]:
                descriptionAndScore[row["Description"]][row["conservationScore"]] = descriptionAndScore[row["Description"]][row["conservationScore"]]+1
            else:
                descriptionAndScore[row["Description"]][row["conservationScore"]] = 1
        
    rows = []

    # Loop through dictionary to create rows
    for substance, sub_data in descriptionAndScore.items():
        for key, count in sub_data.items():
            rows.append([substance, key, count])

    # Create pandas DataFrame
    df = pd.DataFrame(rows, columns=['Description', 'conservationScore', 'Count'])

    df.to_csv(os.path.join(FILE_PATH.GET_ORGANISM_ANALYSIS_FOLDER_PATH(target_organism), "All_Average_Scores.csv"))
                

if __name__ == "__main__":
    average_conservation_score("target_bacteria")

    
