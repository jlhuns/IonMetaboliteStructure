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

#USE INSTEAD OF CHECK_JOB_STATUS to prevent ERRORS!
def checkStatus(jobID):
    keepChecking = True
    while keepChecking:
        time.sleep(2)
        status = check_job_status(jobID)
        if status == "FINISHED":
            keepChecking = False
        else:
            print("Waiting for results...")