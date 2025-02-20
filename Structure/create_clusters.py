from API import PDB_API
import pandas as pd
from Bio import PDB
import numpy as np
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import pdist

import file_paths as FILEPATH

def get_average_coordinate(coordinates):
    """Calculate average X, Y, Z coordinates for binding residues."""
    if not coordinates:
        return None
    return tuple(np.mean(np.array(coordinates), axis=0))

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
                # print(f"Found atoms in range {start}..{end}:")
                # print(matching_atoms[["Residue_ID", "Atom", "X", "Y", "Z"]])
        else:
            # If Binding_Location is a single residue, just filter for it
            matching_atoms = atomsDF[atomsDF["Residue_ID"].astype(str).str.strip() == str(Binding_location).strip()]
            if not matching_atoms.empty:
                coordinates.extend(matching_atoms[["X", "Y", "Z"]].values)
                # print(f"Found atom: {matching_atoms[["Residue_ID", "Atom", "X", "Y", "Z"]]}")

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

    results_df = pd.DataFrame(results)
    return results_df

def compute_ligand_distances(structure, binding_sites):
    """Computes minimum ligand distance to binding sites."""
    ligand_distances = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                ligand_atoms = [atom for atom in residue if atom.element != "H"]
                
                for _, row in binding_sites.iterrows():
                    if row["Coordinates"] is None:
                        continue
                    
                    residue_coords = np.array(row["Coordinates"])
                    min_distance = min(np.linalg.norm(atom.coord - residue_coords) for atom in ligand_atoms)
                    
                    ligand_distances.append({
                        "KOID": row["KOID"],
                        "Ligand_Residue": residue.resname,
                        "Ligand_Chain": chain.id,
                        "Binding_Type": row["Type"],
                        "Binding_Location": row["Binding_Location"],
                        "Ligand_Distance": min_distance
                    })

    return pd.DataFrame(ligand_distances)

def hierarchical_clustering(binding_sites_df, threshold):
    # Extract coordinates
    coordinates = np.array(binding_sites_df['Coordinates'].dropna().values.tolist())
    
    # Calculate pairwise Euclidean distances between residues
    distance_matrix = pdist(coordinates, metric='euclidean')
    
    # Perform hierarchical clustering (Ward's method is commonly used)
    linkage_matrix = sch.linkage(distance_matrix, method='ward')

    clusters = sch.fcluster(linkage_matrix, threshold, criterion='distance')

    # Add cluster labels to your DataFrame
    binding_sites_df['Cluster'] = pd.Series(clusters, index=binding_sites_df.index)
    
    return binding_sites_df, linkage_matrix

import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt

def plot_dendrogram(linkage_matrix, threshold):
    """
    Visualizes the hierarchical clustering dendrogram.
    """
    plt.figure(figsize=(10, 8))
    
    # Plot dendrogram with the linkage matrix
    sch.dendrogram(linkage_matrix)
    
    # Add horizontal line to represent the threshold (for visualizing cluster separation)
    plt.axhline(y=threshold, color='r', linestyle='--')
    
    # Customize labels and title
    plt.title(f"Hierarchical Clustering Dendrogram: {KOID}")
    plt.xlabel("Binding Sites")
    plt.ylabel("Distance (angstrom)")
    
    # Show plot
    plt.show()

# After performing hierarchical clustering

if __name__ == "__main__":
    KOID = "A0A3R9E617"
    try:
        atomsdf = create_structure_atoms_df(KOID)
    except Exception as e:
        print("Error creating structure df: ", e)
        exit()
    conservation_df_file_path = FILEPATH.GET_KOID_MSA_ANALYSIS_DF_FILE_PATH("K01903", "target_bacteria")
    conservation_df = pd.read_csv(conservation_df_file_path)
    binding_sites_df = find_binding_atoms_X_Y_Z(KOID, atomsdf, conservation_df)
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(KOID, "temp.pdb")
    ligand_distances_df = compute_ligand_distances(structure, binding_sites_df)

    print("Binding Sites DataFrame:")
    print(binding_sites_df)

    # Perform Hierarchical Clustering
    threshold = 20
    clustered_binding_sites_df, linkage_matrix = hierarchical_clustering(binding_sites_df, threshold)
    print("Clustered Binding Sites DataFrame:")
    print(linkage_matrix)
    print(clustered_binding_sites_df)
    print(linkage_matrix)

    plot_dendrogram(linkage_matrix, threshold)

