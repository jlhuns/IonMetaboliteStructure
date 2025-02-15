from API import PDB_API
import pandas as pd
from Bio import PDB
import Structure.graphStrucutre as graphStrucutre

import file_paths as FILEPATH

def get_average_coordinate(coordinates):
    # Calculate average of X, Y, Z coordinates
    total_x, total_y, total_z = 0, 0, 0
    count = 0
    
    for x, y, z in coordinates:
        total_x += x
        total_y += y
        total_z += z
        count += 1

    if count == 0:
        return (0, 0, 0)  # Return (0, 0, 0) if no coordinates are provided

    return (total_x / count, total_y / count, total_z / count)




def create_structure_atoms_df(KOID):
    pdb_data = PDB_API.get_pdb_file(KOID)
    parser = PDB.PDBParser(QUIET=True)

    with open("temp.pdb", "w") as f:
        f.write(pdb_data)

    structure = parser.get_structure(KOID, "temp.pdb")
    
    atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append({
                        "Model": model.id,
                        "Chain": chain.id,
                        "Residue": residue.resname,
                        "Residue_ID": residue.id[1],
                        "Atom": atom.name,
                        "X": atom.coord[0],
                        "Y": atom.coord[1],
                        "Z": atom.coord[2]
                    })

    df = pd.DataFrame(atoms)
    return df

def find_binding_atoms_X_Y_Z(KOID, atomsDF, conservationDF):
    filtered_df = conservationDF[conservationDF["UniProtID"].str.startswith(KOID)]
    results = []

    for index, row in filtered_df.iterrows():
        Type = row["Type"]
        Binding_location = row["Binding_Location"]
        coordinates = []

        print(f"Binding Location: {Binding_location}")

        if ".." in Binding_location:
            # If Binding_Location contains a range (e.g., 59..62), split and check for multiple residues
            start, end = Binding_location.split("..")
            start = int(start)
            end = int(end)
            
            # Filter for residues within the range
            matching_atoms = atomsDF[(atomsDF["Residue_ID"] >= start) & (atomsDF["Residue_ID"] <= end)]
            if not matching_atoms.empty:
                coordinates.extend(matching_atoms[["X", "Y", "Z"]].values)
                print(f"Found atoms in range {start}..{end}:")
                print(matching_atoms[["Residue_ID", "Atom", "X", "Y", "Z"]])
        else:
            # If Binding_Location is a single residue, just filter for it
            matching_atoms = atomsDF[atomsDF["Residue_ID"].astype(str).str.strip() == str(Binding_location).strip()]
            if not matching_atoms.empty:
                coordinates.extend(matching_atoms[["X", "Y", "Z"]].values)
                print(f"Found atom: {matching_atoms[["Residue_ID", "Atom", "X", "Y", "Z"]]}")

        if coordinates:
            avg_coordinates = get_average_coordinate(coordinates)
            results.append({
                "KOID": KOID,
                "Type": Type,
                "Binding_Location": Binding_location,
                "Coordinates": avg_coordinates
            })
        else:
            results.append({
                "KOID": KOID,
                "Type": Type,
                "Binding_Location": Binding_location,
                "Coordinates": None
            })
    
    # Create DataFrame from results
    results_df = pd.DataFrame(results)
    return results_df


        



if __name__ == "__main__":
    KOID = "A0A0H4VDQ2"
    atomsdf = create_structure_atoms_df(KOID)
    conservation_df_file_path = FILEPATH.GET_KOID_MSA_ANALYSIS_DF_FILE_PATH("K00611", "target_bacteria")
    conservation_df = pd.read_csv(conservation_df_file_path)
    df = find_binding_atoms_X_Y_Z(KOID, atomsdf, conservation_df)
    graphStrucutre.plot_binding_locations(df)
    print(df)
    
