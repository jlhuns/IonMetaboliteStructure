import file_paths as FILEPATH
from API import Muscle_API
import os

def create_msa_file(KOID: str, targetOrganism: str):
    if targetOrganism.endswith('.csv'):
        targetOrganism = targetOrganism.replace('.csv', "")

    folderPath = FILEPATH.GET_KOID_UNIPROT_ENTRIES_PATH(KOID, targetOrganism)
    sequence_data = FILEPATH.readDirectoryContents(folderPath)

    inputData = {
        "email": "joshua.l.hunsaker@gmail.com",
        "title": f"{targetOrganism}: {KOID}",
        "data": sequence_data
    }

    filePath = os.path.join(FILEPATH.GET_KOID_MSA_FOLDER_PATH(KOID, targetOrganism), f"{KOID}_MSA_Results.aln")
    if FILEPATH.CHECK_PATH_EXISTS(filePath):
        print(f"MSA File aready exists at: {filePath}")
    else:
        jobID = Muscle_API.submit_job(inputData.get("email"), inputData.get("title"), inputData.get("data"))
        Muscle_API.checkStatus(jobID)
        response = Muscle_API.get_results(jobID, "aln-clustalw")

        #create file with result
        with open(filePath, 'w') as file:
            file.write(response)
        file.close()