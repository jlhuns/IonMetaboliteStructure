import file_paths as FILEPATH
import re
import pandas as pd


def get_file_lines(file_name: str) -> list[str]:
    data = []
    with open(file_name, 'r') as inF:
        data = inF.readlines()
    return data

def get_section_data(lines: list[str], section_header: str) -> list[str]:
    section_header = section_header.upper()
    data = []
    for line in lines:
        if line.startswith(section_header):
            line = line.replace(section_header, "")
            line = line.strip()
            data.append(line)
    return data

def check_for_potential_allosteric(KOID, target_organism, line: str, uniprotID, file_count) -> None:
    df = FILEPATH.LOAD_PERSISTANCE_FILE_ALLOSTERIC()  # Ensure this function returns a DataFrame

    key_words = ["allosteric", "Allosteric enzyme", "Activity regulation"]
    
    # Convert line to lowercase and check which keywords are present
    line_lower = line.lower()
    matched_keywords = [keyword for keyword in key_words if keyword.lower() in line_lower]

    if matched_keywords:
        if KOID in df['KOID'].values:
            # Append UniProt ID only if it's not already in the list
            existing_uniprot_ids = set(df.loc[df['KOID'] == KOID, 'UniProtID'].values[0].split(', '))
            if uniprotID not in existing_uniprot_ids:
                df.loc[df['KOID'] == KOID, 'UniProtID'] = ', '.join(existing_uniprot_ids | {uniprotID})
            
            # Count unique UniProt IDs for this KOID only
            unique_uniprot_count = len(df.loc[df['KOID'] == KOID, 'UniProtID'].values[0].split(', '))
            print(f"KOID {KOID}: Unique UniProt IDs: {unique_uniprot_count}")
            df.loc[df['KOID'] == KOID, 'mentioned_files'] = unique_uniprot_count
            df.loc[df['KOID'] == KOID, 'total_files'] = file_count
        else:
            # Create new entry
            new_data = pd.DataFrame({
                'KOID': [KOID], 
                "target_organism": [target_organism], 
                "mentions": [", ".join(matched_keywords)],
                "UniProtID": [uniprotID]
            })
            df = pd.concat([df, new_data], ignore_index=True)

        FILEPATH.SAVE_PERSISTANCE_FILE_ALLOSTERIC(df)  # Ensure this function properly saves

    return  # Explicit return

def get_cofactors(lines: list[str]) -> list[str]:
    lines = get_section_data(lines, "CC")
    cofactors = []
    pattern = r'Name=([^;]+)'

    for line in lines:
        test = re.findall(pattern, line)
        if test:
            cofactors.extend(test)
    return cofactors

def check_ion_metabolite(df):
    metals = [
    "Al", "Sb", "As", "Ba", "Be", "Bi", "Ca", "Cd", "Cr", "Co",
    "Cu", "Dy", "Er", "Eu", "Ga", "Au", "Hf", "Ho", "Fe", "Pb",
    "Li", "Mn", "Hg", "Mo", "Nd", "Ni", "Os", "Pd", "Pt", "Pu",
    "Po", "K", "Rh", "Rb", "Ru", "Sm", "Sc", "Si", "Ag", "Na",
    "Ta", "Tl", "Th", "Sn", "Ti", "U", "V", "Yb", "Zn", "Zr", "Mg",
    "Se", "W"
    ]
    pattern = r'\b(?:' + '|'.join(metals) + r')\b'
    
    for index, row in df.iterrows():
        description = row["Description"]
        if re.search(pattern, description):
            df.at[index, "Ligand_Type"] = "Ion"
        else:
            df.at[index, "Ligand_Type"] = "Metabolite"










if __name__ == "__main__":
    # koid = "K01679"
    # lines = get_file_lines(FILEPATH.GET_UNPROT_ENTRY_FILE_PATH(koid, "A0R2U8", ".txt", "target_bacteria"))
    # section_cc = get_section_data(lines, 'cC')


    # #A Test
    # df = pd.read_csv(FILEPATH.GET_KOID_MSA_ANALYSIS_DF_FILE_PATH(koid, "target_bacteria"))
    # ions = [
    # "Na(+)", "K(+)", "Ca(2+)", "Mg(2+)", "Fe(3+)", "fe-su", "Mn(2+) 2"
    # ]

    # # List of random metabolites
    # metabolites = [
    #     "Glucose", "Lactate", "Pyruvate", "Citrate", "ATP", "NADH", "FADH2", "Urea", "Creatinine", "Acetyl-CoA"
    # ]

    # # Combine ions and metabolites
    # data = ions + metabolites

    # # Create DataFrame dictionary
    # dataDF = {
    #     "Description": data
    # }

    # # Convert to DataFrame
    # test_df = pd.DataFrame(dataDF)
    # check_ion_metabolite(test_df)

    file_path = FILEPATH.GET_ANALYSIS_RESULT_FILE("test_target_bacteria")
    df = pd.read_csv(file_path)
    check_ion_metabolite(df)
    df.to_csv(file_path)

    
