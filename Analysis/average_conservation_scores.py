import file_paths as FILE_PATH
import pandas as pd

def average_conservation_score(file):
    conservation_df = pd.read_csv(file)
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
        
    print(descriptionAndScore)  


                


if __name__ == "__main__":
    test_file = FILE_PATH.GET_ANALYSIS_RESULT_FILE("target_bacteria")
    print(average_conservation_score(test_file))

    
