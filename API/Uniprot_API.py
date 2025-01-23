import requests
import os

def fetch_uniprot_file(uniprot_id: str, output_dir: str = ".", file_format: str = "txt"):
    """
    Fetches a UniProt file for a given UniProt ID and saves it to a specified directory.

    Parameters:
        uniprot_id (str): The UniProt ID of the protein (e.g., "P12345").
        output_dir (str): Directory where the file will be saved. Default is the current directory.
        file_format (str): Format of the file to retrieve (e.g., "txt", "fasta", "xml"). Default is "txt".

    Returns:
        str: Path to the saved file.

    Raises:
        ValueError: If the response from UniProt is not successful.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/"
    url = f"{base_url}{uniprot_id}.{file_format}"

    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise ValueError(f"Failed to fetch UniProt file for ID {uniprot_id}: {e}")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save the file
    file_path = os.path.join(output_dir, f"{uniprot_id}.{file_format}")
    with open(file_path, "wb") as file:
        file.write(response.content)

    print(f"File saved to: {file_path}")
    return file_path

def fetch_and_extract_uniprot_sections(uniprot_id: str, sections: list[str], output_dir: str = "."):
    """
    Fetches a UniProt .txt file for a given UniProt ID and extracts specific sections.

    Parameters:
        uniprot_id (str): The UniProt ID of the protein (e.g., "P12345").
        sections (list[str]): The sections to extract (e.g., ["ID", "SQ"]).
        output_dir (str): Directory where the output file will be saved. Default is the current directory.

    Returns:
        str: Path to the saved file with the extracted sections.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"

    try:
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        raise ValueError(f"Failed to fetch UniProt file for ID {uniprot_id}: {e}")

    # Extract the sections
    lines = response.text.splitlines()
    extracted_lines = []
    for section in sections:
        extracted_lines.extend([line for line in lines if line.startswith(section)])
    
    if not extracted_lines:
        raise ValueError(f"None of the specified sections ({', '.join(sections)}) were found in the file for ID {uniprot_id}.")

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save the extracted sections
    output_path = os.path.join(output_dir, f"{uniprot_id}.txt")
    with open(output_path, "w") as file:
        file.write("\n".join(extracted_lines))

    print(f"Extracted sections saved to: {output_path}")
    return output_path