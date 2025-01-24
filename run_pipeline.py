#Get Orthologs from KEGG using KOID and target csv. (RETURN VALUE: ARRAY of uniprot ID)
#Use UNIPROT ID to get uniprot Entries (CREATE FILES UNDER DATAFILES MSA UNIPROT)
#Make a fasta file using UNIPROT ENTRIES (ADD FILE UNDER KOID MSA)
#RUN MSA (ADD file to MSA under correct KOID)
#Analyze MSA (binding or not, return DF )

from Ortholog_Identificiation import create_uniprot_entries, get_orthologs_uniprotID



if __name__ == "__main__":
    uniprot_ids = get_orthologs_uniprotID.get_uniprot_ids("K00973", "target_prokaryotes.csv")
    create_uniprot_entries.create_uniprot_entires(uniprot_ids, "K00973", "target_prokaryotes.csv")
