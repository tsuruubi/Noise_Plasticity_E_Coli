# Coupling between noise and plasticity in protein expression level in E. coli

## Descriptions of folders
Zenodo: Read the "Readme.txt" file in the Zenodo folder and follow the instruction. Download three datasets (YFP_Fused_DS01–03) from the Zenodo repository (https://doi.org/10.5281/zenodo.13318136; https://doi.org/10.5281/zenodo.13318501; https://doi.org/10.5281/zenodo.13319277). Unzip and place the folders (ex.T02191218) into the Zenodo folder.  
FlowData_Output: output files of the "FlowData_Analysis.Rmd" script  
ENV_DataSet: the Env dataset comprising of transcriptome profiles (PRECISE 2.0) obtained from iModulonDB (http://imodulondb.org/)  
RegulonDB: known regulatory relationship between transcriptional regulators and target genes obtained from RegulonDB (https://regulondb.ccg.unam.mx)  
W3110_genes: essential genes and gene names of E. coli W3110 used in the Mut dataset  
MG1655_genes: codes for extracting gene names and locus tags of E. coli MG1655  
Meta: metadata of fcs file, medium information  
Sigmulon: known ptomoters of genes in E. coli MG1655 obtained from RegulonDB (https://regulondb.ccg.unam.mx)  
Growth_Curve: kime course data of cell densities under different nutrient conditions  
Metabolite_TF_Interactions: kelationship between metabolic effectors and transcirptional regulators, obtained from Lempp et al (https://doi.org/10.1038/s41467-019-12474-1)  
ECOCYC: known regulatory relationship between transcriptional regulators and target genes obtained from EcoCyc (https://ecocyc.org)  
NOISE_CH_DataSet: protein noise profiles obtained from Taniguchi et al (https://doi.org/10.1126/science.1188308)  
NOISE_PL_Dataset: promoter-mediated noise profiles obtained from Urchueguia et al (https://doi.org/10.1371/journal.pbio.3001491)  
Output: output files of the "Main.Rmd" script  

## Descriptions of scripts
FlowData_Analysis.Rmd: codes for calculating the means and sds of the GFP distributions  
Main.Rmd: codes for analysis and figures  

## Requirements:
R 4.3.2 or greator  