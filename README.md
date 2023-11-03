OCT 10 
CPTAC3 - several diseases in it, has other data related to it, proteomics project in addition to genomics = would be best for the project 
Proteomic data commons where these cases also exist and be cool to later tie in that information
Overall, around 3000 cases/chunks of data that he will process through the workflow; will take around 2 months to get totally done (continuously integrate data as they are moving) and call an endpoint 
Use to start building the parts to match records, organizing, based on what is in the directory so we can iteratively add on to it 
Workflow is not running as fast as we hoped - only 1 node running 
Tgz file has 130 read groups/results, 
Merge tgz file with query data 
In tgz = results of workflow (2 files for each bit that goes through; 1 is counts file and the other is multiqc results file)
Main concern is counts file, but would be cool to open multiqc file to assess quality of the counts data (data generated)

Query (2 tables) = raw output from query, will have meta data about file itself and about the patient (case) 
In GDC, use graph database - has diff entities that relate to each other 
Case (patient), sample (biopsy), alliquot (subsample of analyte, which is RNA isolated from the biopsy and used to build library for sequencing), read group
Two query files because we have 2 separate but related databases 
Case ID and read group and file info = read group level (one case could appear multiple times) 
case UID, clinical variables (age, sex, disease state, etc) = this will be at the case level only
Tgz file = output of the workflow; workflow will take read group ID and UIDs for the actual fastq files 
Push files through fastqc and salmon
Multiqc will take outputs and wrap into data file 
What is included in tqz file = quantification file from salmon, data file from multiqc 
For 130 runs (read groups) 
Naming convention is read group EUID = should be able to use those IDs to match it back to those tables 
And can correlate all the way back to the case 

What to do now: all in Python jupyter notebook (can use pandas because there is tabular data) 
Start by exploring tsv files, figuring how to match everything together 
Once we have feel for that, open up tgz file and extract all the content 
Become familiar with the content (because these have transcript expression data in it) 
Start playing around with clustering those 
How many results do i have from these types of clusters
OCT 17 
Read counts data itself
Dont worry about multiqc files, just the quant.sf (output from salmon)
Open quant.sf and will be tsv file
Look at documentation for salmon for format 
Will give counts info (raw counts, tpm, num reads, etc) 
Use TPM for our purposes (tags per million) - like normalization 
For each of the files, open and extract TPM and transcript name/identifier 
Concatenate all of them into 1 big dataframe for now 
Can then work on clustering them
Maybe look into doing PCA (such as sklearn) 
Visualization of the data (heatmaps, 2D or 3D cluster visualizations) 
Metadata analysis 
1st way: visualizations, use metadata as a visual representation; something that would color all of the dots by age, or race, or disease type etc 
Heatmap or cluster map: programs will accept another column/row based on categorical column you give it (like race) 


Git - just use code, don’t keep the data there 
Use box for the data 

Waqaas Notes:

Read the counts data itself
open the sf file
extract gene name, and then the TPM
concatenate it all into one dataframe
PCA on the counts data, sklearn > visualization of the data, heatmaps, 2d-3d cluster visualizations
visualizations would be like age/disease type, or there could be clusters for categorical variables, such as race.
then work on clustering them
future: slice out particular disease data from all the available data
query data has all the clinical variables
if you wanted to get out some disease symptoms, then you could make a selection out of that query data
make functions that can make it an interactive process

OCT 22 
TPM = transcripts per million 
For every case, there is SF file containing all of the genes and a TPM for every gene 
Filtering out the genes in question 
May have to make function
PCA

What we did : 
Create code that reads the tsv files in the data directory and joins the files ending with sur.tsv into one data frame and those not ending with sur.tsv into another 
Create a dataframe from the files not ending in sur.tsv 
Read the quant.sf file 
In the quant.sf file: 
Transcript name, length, effective length, TPM, NumReads 
of all sf files: create dictionary with name and TPM 
Append this to the dataframe 
** note: changed “name” to “transcript_id” for clarity 

Questions to ask:
In the combined dataframe, we have 2 entries with the same read_group_id (one being R1, the other R2), but there is only one corresponding sf file, how to interpret?
The sf file is stupid large, so I think we must prefilter for relevant genes before any merging. How do we decide? Do we use biopython or other packages to enable searching the transcript IDs using gene IDs?
Is there an alternative to a combined dataframe? We are considering turning the sf files into wide format dataframes for easier integration, but we keep them separate until the scope of the analysis is determined.

Oct 25
How/where do we process SF file
1 SF file per read group; 1 row per read file 
Separate record for R1 and R2, but rest of variables should be the same 
Can remove unaligned read IDs and R1 R2 
Don’t bring all of the counts to the same table 
Prefiltering? - don’t filter the SF file at all
Keeps the SF/counts separate from metadata, what we will want to do is make a selection based off metadata 
Then we can choose the counts for a set of read groups 
Work on that filtering mechanism first 
The readgroup ID is the key to link those 
Seurat to make clusters? Can we do something similar?
Sort of; seurat specifically for single cell RNAseq
Figure out what diseases are represented in the metadata
Pick the most populous one
experiment?
Sf files into 1 large DF
PCA
Prior to PCA - matrix of cases/samples -> samples and new principle components 
Carried out across transcript IDs; would want to end up with 
heat map, k means clustering
Filtering (due to large file) 
Variants filtration; 1 row corresponds to 1 transcript, calc variants and pick arbitrary cutoff that leaves a manageable amt of rows (ask matt to elaborate)

Oct 30 

read in the sf files one by one
use the Name column in the files to generate columns in the dataframe, and the content of the cell will be the TPM value
for the next sf file, if the Name column already exists, then add the TPM value to the existing column, otherwise create a new column
Z test, t test, anova 

Oct 31
Q’s for Garcia
taking data that we have then using that to do the metadata analysis 
What other variables of interest should we look at
Need to find which variables have enough data to be able to statistically investigate 
Do column by column missingness - do we have at least 80-90% data there? 
That will refine about 10-20 columns that actually have majority of the data
What type of test for gender and race? 
Still trying to turn all the sf files into a dataframe for the data we have 
Next part of the code was meant to match the name and TPM and populate all of the columns for the 133 cases
Transcript IDs as column names - we don’t need to do that. We want the transcript IDs as the rows for the orientation 
Take the tpm column and keep that and get rid of everything else. So we just need the transcript ID and TPM column for 1 file, set the index to be the transcript IDs, then maybe replace the column name ‘TPM’ with the read group ID. And have that be the dataframe
Then open the next file and basically do the same except the name column is the index, change the name of the TPM to the new read group ID then concatenate that column to the last one - then build it up column by column
No need to use numpy array or anything, just let it do it automatically 

Nov 2
Separated code into a table generator and a metadata analysis file.
Cleaned code to introduce 'already-executed' contingencies and standardized the variable names for dataframes.
Successfully compiled the TPM dataframe
