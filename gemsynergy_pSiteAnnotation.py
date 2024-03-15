import psite_annotation as pa
import pandas as pd

path = "p10_evidence.txt"
evidence_df = pd.read_csv(path, sep="\t")

pa.addPeptideAndPsitePositions(evidence_df, "Phosphosite_seq.fasta", pspInput=True)

#evidence_df[['Proteins', 'Modified sequence', 'Start positions', 'End positions', 'Site sequence context', 'Site positions']]

evidence_df.to_csv("evidence_annotated.txt", sep="\t", index=None)
