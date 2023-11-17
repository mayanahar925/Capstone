**October 10:** Discussed the project, which involves CPTAC3 data containing several diseases and both genomics and proteomics information. The plan is to process around 3,000 cases' data using Python and integrate it as it arrives. This will take approximately two months. Explored the data and planned to match records and organize based on directory structure. Noted that the workflow was running slowly due to only one node.

**October 17:** Discussed processing the raw counts data. Focused on TPM (Transcripts Per Million) data in quant.sf files generated by Salmon. Planned to open and extract TPM values and transcript identifiers. The goal was to concatenate all this data into one large dataframe for further analysis, potentially including clustering.

**October 22:** Discussed the specifics of TPM data, focusing on filtering relevant genes. Planned to extract TPM data and concatenate it into a single dataframe. Mentioned the possibility of using PCA (Principal Component Analysis) for data visualization.

**October 25:** Addressed how to process Salmon (SF) files. Each SF file corresponds to a read group and contains gene data. Discussed prefiltering for relevant genes and deciding on scope. Mentioned the option of turning SF files into wide-format dataframes for easier integration.

**October 30:** Discussed code to read SF files one by one, extract TPM values, and generate a single dataframe. The idea was to create columns with gene names and populate them with TPM values from different SF files. Mentioned statistical tests like Z-tests, t-tests, and ANOVA.

**October 31:** Planned to use metadata for analysis. Discussed selecting variables of interest based on available data and assessing missing data. Considered which statistical tests to use for variables like gender and race. Discussed matching gene names and TPM values for 133 cases. Emphasized that TPM values should be stored in a dataframe with transcript IDs as rows and read group IDs as columns.

**November 2:** Separated the code into a table generator and a metadata analysis file. Cleaned the code, introduced 'already-executed' contingencies, and standardized variable names. Successfully compiled the TPM dataframe.

**November 4:** Managed to perform a PCA and visualize the samples using 2d and 3d scatterplots. Also processed the dataframes for further clustering analysis.

**November 11:** Imported pyensembl to map transcript IDs to gene IDs. Defined a query maker function, which outputs two dataframes that captures the cases as rows along with the transcripts as columns with the TPM as the value. Used T test to compare group A and group B to identify significant and nonsignificant transcript IDs based on the given query. Defined a second function that interprets the variance at the gene level between group A and group B.

**November 11** Imported pyensembl to map transcript IDs to gene IDs. Defined a query maker function, which outputs two dataframes that captures the cases as rows along with the transcripts as columns with the TPM as the value. Used T test to compare group A and group B to identify significant and nonsignificant transcript IDs based on the given query. Defined a second function that interprets the variance at the gene level between group A and group B. 