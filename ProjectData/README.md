# Torok-Storb 2018 - SHIP program

The Torok-Storb lab project involved comparing co-cultures of a canine DS1 cell line of stromal cells with normal or AML hematopoietic bone marrow cells in PBS. The goal of the project was to determine the influence of hematopoietic cells on the non-hematopoietic environment.

Normal cells were CD45+ and CD34+. Stromal cells are connective tissue cells of any organ the most common types of such being fibroblasts and pericytes. Hematopoietic cells are blood forming cells. PBS stands for Phosphate-Buffered Saline.

In this experiment DS1 cells co-cultured with normal bone marrow are the control (labeled NormBM) versus those co-cultured with the AML bone marrow cell line ML3.  The data from the two sets of three samples was generated on the same day. All the RNA extracted via column extractions was high quality (based on RIN). The Illumina TruSeq RNA library prep v2 was used to generate libraries from the RNA samples for 50 base, paired end, unstranded sequencing.

## Metadata:
The Translational Genomics Repository (TGR) provides a data dictionary for metadata associated with genomic data sets. The specific scientific metadata for this study can be found in the file: `torokstorbMetadata.csv`
The definitions of the annotations (the column names) and the values (the values in the columns) can be found here: https://ontology.fredhutch.org/ if you click on the "Translational Genomics Repository" tab.


## RNA sequencing data:
Each folder contains files associated with a single dataset the metadata for which is in the above csv. Each dataset in the TGR receives a unique molecular_id, as listed in the metadata csv file. This id is also the name of each of the sub-folders here that contain the results from each individual dataset's processing via the STAR workflow.
