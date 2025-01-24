import time
import requests
import xml.etree.ElementTree as ET
import os

# Base URL for the Clustal Omega Job Dispatcher API
BASE_URL = 'https://www.ebi.ac.uk/Tools/services/rest/muscle'
DEFAULT_HEADERS = {
    'Accept': 'application/json'
}

# Function to list available parameters
def get_parameters():
    url = f'{BASE_URL}/parameters'
    response = requests.get(url, headers=DEFAULT_HEADERS)
    if response.status_code == 200:
        return response.json()
    else:
        return f'Error: {response.status_code}'

# Function to get detailed information about a specific parameter
def get_parameter_details(parameter):
    url = f'{BASE_URL}/parameterdetails/{parameter}'
    response = requests.get(url, headers=DEFAULT_HEADERS)
    if response.status_code == 200:
        return response.json()
    else:
        return f'Error: {response.status_code}'

#Returns Job_ID
def submit_job(email, title, sequence, parameters = None)-> str:
    url = f'{BASE_URL}/run'
    data = {
        'email': email,
        'title': title,
        'sequence': sequence,
        # 'parameters': parameters  # Replace with actual parameters
    }
    response = requests.post(url, data=data) 
    if response.status_code == 200:
        return response.text  # Return job ID
    else:
        return f'Error: {response.status_code}'

# Function to check the status of a job
def check_job_status(job_id):
    url = f'{BASE_URL}/status/{job_id}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        return f'Error: {response.status_code}'

# Function to get result types for a job
def get_result_types(job_id):
    url = f'{BASE_URL}/resulttypes/{job_id}'
    response = requests.get(url, headers=DEFAULT_HEADERS)
    if response.status_code == 200:
        return response.json()
    else:
        return f'Error: {response.status_code}'

# Function to get results from a job
def get_results(job_id, result_type):
    url = f'{BASE_URL}/result/{job_id}/{result_type}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text  # Returns the result data
    else:
        return f'Error: {response.status_code}'
 
    #MOVE LATER
def readDirectoryContents(folder_path):
    fileList = os.listdir(folder_path)
    sequences = ""
    for fileName in fileList:
        filePath = os.path.join(folder_path, fileName)
        with open(filePath, 'r', encoding="utf-8") as inF:
            sequences += inF.read()
            sequences += "\n"
    return sequences

def checkStatus(jobID):
    keepChecking = True
    while keepChecking:
        time.sleep(5)
        status = check_job_status(jobID)
        if status == "FINISHED":
            keepChecking = False
        else:
            print("Waiting for results...")


# Example usage
if __name__ == '__main__':
    # Step 1: Get parameters
    print('Getting available parameters...')
    parameters = get_parameters()
    print(parameters)

    print("checking a parameter description")
    description = get_parameter_details('sequence')
    print(description)
    
    folderPath = r"C:\Users\joshu\Desktop\Structure Research\IonMetaboliteStructure\DataFiles\target_prokaryotes\K00973\UniProt_Entries"

    sequence = readDirectoryContents(folderPath)

    # Step 2: Submit a job (replace with actual input data)
    print('Submitting job...')
    job_id = submit_job("joshua.l.hunsaker@gmail.com", "test", sequence)
    print(job_id)
    
    # Step 3: Check job status (replace with actual job ID)
    print('Checking job status...')
    job_status = checkStatus(job_id)
    print(job_status)
    
    # Step 4: Get result types (replace with actual job ID)
    print('Getting result types...')
    result_types = get_result_types(job_id)
    print(result_types)
    
    # Step 5: Get results (replace with actual result type)
    result_type = 'aln-clustalw'  # Replace with actual result type
    print('Getting results...')
    results = get_results(job_id, result_type)
    print(results)
