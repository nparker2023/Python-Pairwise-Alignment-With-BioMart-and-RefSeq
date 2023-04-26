# Python Pipeline Tutorial to Query and Manipulate Data from BioMart and RefSeq in Order to Perform a Pairwise Alignment Between Two Species 
### Overview
The following tutorial gives a step by step guide on how to successfully use the pipeline in order to get the desired results. The Python script without the explanations has also been provided and can be found above.
```Python
!pip install pybiomart --quiet
!pip install biopython --quiet
```
The following lines for this block of code install packages that aren't already pre-installed in the Google Collab Python notebook. This is done by using the code "!pip install", which is then followed by the specified package. The use of "--quiet" prevents the code block from outputting any information related to the package installation. However, if some other type of IDE is being used, this code does not need to be used. Due to this, the process of installing packages may vary depending on the IDE being used.

```Python
from pybiomart import Dataset, Server
import pandas as pd
from Bio import Entrez, Align, SeqIO
from Bio.Align import substitution_matrices
```
The following lines for this block of code will load the packages that are required for coding. The use of import allows for an installed package to be accessed and used within the code. The use of "as" allows for the package to be given a nickname or shortened. In this example, the package pandas is shortened to pd. The use of "from" followed by import allows for specified subpackages to be accessed from a specified package. If numerous subpackages are to be accessed, then each one can be listed after import and separated with a comma.

```Python
def mart_finder(file_name_1):
    server = Server(host='http://www.ensembl.org')
    list_1 = server.list_marts()
    list_1.to_csv(file_name_1, index=False)
```
The following lines create a function that will return a csv that contains a list of all the possible BioMart Ensembl marts. The function has one argument. The argument file_name_1 refers to the name of the csv file for all possible marts. The list of marts was saved as list_1 and converted to a csv file. The list of all possible marts can be used to determine what mart will be entered for the mart_name argument of the database_finder function.

```Python
def database_finder(mart_name, file_name_2):
    server = Server(host='http://www.ensembl.org')
    mart = server[mart_name]
    list_2 = mart.list_datasets()
    list_2.to_csv(file_name_2, index=False)
```
The following lines create a function that will return a csv that contains a list of all the possible BioMart Ensembl databases. The function has two arguments. The first argument mart_name refers to the name of the selected mart. This mart can be found from the list that was created from the mart_finder function. The second argument file_name_2 refers to the name of the csv file for possible datasets. The list of datasets is saved as list_2 and converted to a csv file. The list of all possible databases can be used to determine what species will be entered for the species argument of the dataset_retrieve function.

```Python
def filters_attributes(species, file_1, file_2):
    species_dataset = Dataset(name=species, host='http://www.ensembl.org')
    list_1 = species_dataset.filters
    with open(file_1, 'w') as f:
        for item in list_1.keys():
            f.write("%s,%s\n" % (item, list_1[item]))

    list_2 = species_dataset.attributes
    with open(file_2, 'w') as f:
        for item in list_2.keys():
            f.write("%s,%s\n" % (item, list_2[item]))
```  
The following lines create a function that will return csvs containing all the possible attributes and filters for a specified dataset. The argument species refers to the selected database and the arguments of file_1 and file_2 refer to what the csvs will be called. The species_dataset variable is set equal to the Dataset subpackage, which requires a name and host to access a specified database. The name variable is set equal to the argument of species and the host variable is set equal to the Ensembl website. The attributes and filters for the species_dataset variable can be found by adding .filters or .attributes after it. The filters and attributes are then saved as list_1 and list_2. Since these two variables are dictionaries, they are then written to csv files through the use of for loops and .key() functions. It should be noted however, that some listed filters and attributes are incomplete and missing information due to the way that the package was created. These attributes and filters can be used in the dataset_retrieve function. While there are some attributes and filters that are the same from dataset to dataset, there are some attributes and filters that are unique to each dataset. For the purposes of this tutorial, all the attributes and filters used here are present for both of the selected datasets.

```Python
def dataset_retrieve(species, chrom, file_name):
    species_dataset = Dataset(name=species, host='http://www.ensembl.org')
    species_query = species_dataset.query(
        attributes=['refseq_mrna', 'refseq_peptide', 'ensembl_gene_id', 'external_gene_name', 'description',
                    'start_position', 'end_position', 'strand', 'chromosome_name', 'name_1006'],
        filters={'chromosome_name': [chrom]})
    filtered_set = species_query.dropna()
    filtered_set.columns = filtered_set.columns.str.replace(' ', '_')
    filtered_set.to_csv(file_name, index=False)
```
The following lines create a function that will return a csv file containing queried data from BioMart. In order to get this file, three arguments must be entered. These arguments are species, chrom, and file_name. The species argument refers to the database from which data will be queried. Furthermore, this database contains genetic information related to a particular species. A list of possible species can be found using the database_finder function mentioned above. The chrom argument refers to the chromosomal location in which the genes are located. The file_name argument refers to what the output csv file will be named. The species_dataset refers to the variable that is set equal to the dataset that will be selected from BioMart. Using the Dataset subpackage, which requires name and host, it is able to access a specified database. The name variable is set equal to the argument of species and the host variable is set equal to the Ensembl website. The species_query refers to the variable that is set equal to the queries that will be retrieved based on the specified attributes and filters. The queries are filtered by chromosomal location, which is set equal to the argument of chrom. All the attributes and filters can be found using the filters_attributes function mentioned above. The variable filtered_set is equal to the species_query dataframe that has been filtered to remove queries that are missing information. The filtered_set is then updated so that all spaces in column names are replaced with underscores. This allows for easier manipulation of the data. The filtered_set is then converted to a csv file and given a name based on the file_name argument.

```Python
def gene_list(species_1, chrom, species_2_id, species_2_gene_name, file_name):
    species_dataset = Dataset(name=species_1, host='http://www.ensembl.org')
    gene_list_query = species_dataset.query(
        attributes=['ensembl_gene_id', 'external_gene_name', species_2_id, species_2_gene_name],
        filters={'chromosome_name': [chrom]})
    filtered_set = gene_list_query.dropna()
    filtered_set.columns = filtered_set.columns.str.replace(' ', '_')
    filtered_set.to_csv(file_name, index=False)
```
The following lines create a function that lists the genes of a particular species and its corresponding homologs to another specified species. This information is then saved to a csv file. This function has five arguments that must be entered. The first argument is the name of the database from which the information will be queried. The second argument is the chromosomal location in which the genes are located. The third argument is the gene ids for the second specified species, which is used as an attribute when querying the data. The fourth is the gene names for the second specified species, which is also used as an attribute when querying the data. The fifth argument is file_name, which is used to name the newly created csv. To query them, it accesses a database as previously mentioned in the dataset_retrieve function. While this function is like the dataset_retrieve function, the attributes for this function are different. In particular, the attributes list the gene id and gene name for the specified species and homologs. All the queries are filtered by chromosomal location, which is set equal to the argument of chrom. All of this information is set equal to the variable of gene_list_query. The gene_list_query is then filtered to drop any queries that have blank/missing information. This filtered list is called filtered_set. The filter_set then has its column names renamed to convert any spaces to underscores. This filtered set is then saved to a csv and the name is specified by the argument file_name.

```Python
def gene_list_dataset_1_filter(species, gene_list, species_filter, filter_gene):
    dataset = pd.read_csv(species)
    genes = pd.read_csv(gene_list)
    list_1 = genes['Gene_name'].unique()
    dataset_query = dataset.query("Gene_name in @list_1")
    dataset_query.to_csv(species_filter, index=False)
    list_2 = dataset_query['Gene_name'].unique()
    genes_filter = genes.query("Gene_name in @list_2")
    genes_filter.to_csv(filter_gene, index=False)
```

The following lines create a function that filters out queries that don't appear on the gene list for the first species dataset. This function has four arguments, which are species, gene_list, species_filter, and species_gene. The first argument refers to the filename of the species dataset that is read in and the second argument refers to the filename of the gene list that is read in. The argument species_filter refers to the file that is created after absent genes are filtered out of the dataset. The argument filter_gene refers to the file that is created after genes from the gene_list are filtered out if they are absent from the dataset. The gene list and the dataset for the first species are read in and set equal to the variables of dataset and genes. To get all the unique genes for the first species, a variable called list_1 is created. This variable is a dataframe that is made up of genes from the gene variable. This removes any duplicate gene names from the dataframe and helps to make the dataframe as small as possible. The variable dataset_query is created and returns queries from the dataset variable that have a gene name that is present in list_1. These queries are then saved to a csv file and given a name with the species_filter argument. The gene list then must be updated to filter out genes that were not present in the dataset. A variable called list_2 is created and makes a dataframe of all the gene names from the filtered dataset. Like list_1, this dataframe removes duplicate names that might be present in the list. The variable gene_filter is created and returns queries from the genes variable that have a gene name that is present in list_2. These queries are then saved to a csv file and given a name through the use of the filter_gene argument.

```Python
def gene_list_dataset_2_filter(species, gene_list, column_name, file_name_1, file_name_2):
    dataset = pd.read_csv(species)
    genes = pd.read_csv(gene_list)
    list_1 = genes[column_name].unique()
    dataset_query = dataset.query("Gene_stable_ID in @list_1")
    dataset_query.to_csv(file_name_1, index=False)
    list_2 = dataset_query['Gene_stable_ID'].unique()
    dataset_query_2 = genes[genes[column_name].isin(list_2)]
    dataset_query_2.to_csv(file_name_2, index=False)
```

The following lines create a function that filters out queries that don't appear on the gene list for the second species dataset. This function has five arguments, which are species, gene_list, column_name, file_name_1, and file_name_2. Like the arguments in the previously mentioned function, the first argument refers to the filename of the species dataset that is read in. The second argument refers to the filename of the filtered gene list that is read in. This filtered gene list was created from the gene_filter variable in the previous function. The argument column_name refers to the column name that corresponds to the homologs. The variable dataset_query is created and returns queries from the dataset that have a gene name that is present in list_1. These queries are then saved to a csv file and given a name with the file_name_1 argument. For the filtered gene list to reflect the genes that are present for species_2, it must be queried against the queried dataset. To do this, all the unique gene ids for the second species must be retrieved. This is done through the creation of a variable called list_2 which creates a dataframe of unique gene ids from the dataset_query variable. This removes any duplicate gene ids from the dataframe and helps to make the dataframe as small as possible. Since the name of the homolog can vary depending on the species that is chosen for analysis, the argument of column_name allows for different parameters to be entered for the function. The variable query_2 is created and returns queries from the genes variable that have a gene id that is present in list_2 based of off the parameter entered for column_name. These queries are then saved to a csv file and given a name with the file_name_2 argument.

```Python
def dataset_1_final_filter(species, gene_list, file_name):
    dataset = pd.read_csv(species)
    genes = pd.read_csv(gene_list)
    list_1 = list(genes['Gene_name'].unique())
    dataset_query = dataset.query("Gene_name in @list_1")
    dataset_query.to_csv(file_name, index=False)
```
The following lines create a function that updates the first dataset for the first species. As shown in the gene_list_dataset_2_filter function, the gene list was updated to remove any genes that were not present in the dataset for the second species. This process also removed the corresponding homologs for the first species. Due to the removal of the corresponding homologs, the dataset for the first species is updated for a final time to reflect the filtering that occurred for the dataset of the second species. This function has three arguments, which are species, gene_list, and file_name. The first argument refers to the filename of the filtered dataset that is read in for the first species. The second argument refers to the filename of the filtered gene list that is read in. This filtered gene list was created in the previous function and shows all the genes present for both datasets. The third argument refers to what the filtered dataset will be called. After reading in the csvs for the genes and dataset variables, a variable called list_1 is created. The variable creates a dataframe of unique gene names from the dataset variable. The variable dataset_query is created and returns queries from the dataset that have a gene name that is present in list_1. These queries are then saved to a csv file and given a name with the file_name argument.

```Python
def gene_ontology_filter(file, go_term, go_name_filter):
    filtered_species = pd.read_csv(file)
    query = filtered_species[filtered_species['GO_term_name'].isin([go_term])]
    query.to_csv(go_name_filter, index=False)
```
The following lines create a function that returns queries that have a particular gene ontology term. The function has three arguments, which are file, go_term, and go_name filter. The first argument refers to the filtered species file that is read in and the second parameter refers to the selected gene ontology term. The third argument refers to the name that the output csv file is given. The filtered species file is read in and saved as filtered_species. The query variable is then created, which returns queries if the specified go_term is present in the GO_term_name column. These queries are then saved to a csv file and given a name with the go_name_filter argument.

```Python
def ref_seq_list(file_name, gene, column_name, name):
    file = pd.read_csv(file_name)
    filtered = file.loc[:, [column_name, "Gene_name"]]
    filtered.drop_duplicates(keep='first', inplace=True)
    list = [gene]
    desired_gene = filtered.query("Gene_name in @list")
    desired_gene.to_csv(name, index=False)
```
The following line creates a function that returns Refseqs for a particular gene. The function has four arguments, which are file_name, column_name, gene, and name. The first argument refers to the selected file that contains queries filtered by gene ontology name for a specified species. The second argument column_name refers to the column in which the desired Refseqs are. The Refseqs can be selected from either the RefSeq_mRNA_ID or RefSeq_peptide_ID column. These columns contain the sequences for a specified gene in mRNA and protein form. The third argument gene refers to the gene of interest, and the fourth argument refers to what the output csv will be called. The filtered variable is created by selecting one of the two previously mentioned columns and the column Gene_name. The filtered variable is then manipulated so that duplicated Refseqs are dropped. A variable called list is created and creates a list which contains whatever was entered for the gene argument. The variable desired gene is created, which returns queries that have a gene name that is present in the list variable. The desired gene variable is then saved to a csv file and given a name with the name argument.

```Python
def ref_seq_sequence(email, db_type, id, file_name):
    Entrez.email = email
    net_handle = Entrez.efetch(db=db_type, id=id, rettype='fasta', retmode='text')
    out_handle = open(file_name, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
```
The following lines create a function that outputs a RefSeq based on its id and database type. The function has four arguments, which are email, db_type, id, and file. The argument email refers to the email address that is entered. To access the package for that can retrieve a RefSeq, an email address must be entered each time. If no email address is entered, a warning will come up that one needs to be entered. The second parameter refers to the database type, which can either be nucleotide or protein. The third parameter id refers to the RefSeq id that is entered. The fourth parameter file_name refers to what the output file will be named. After entering a valid email address, a variable called net_handle is created, which retrieves information for a specified RefSeq depending on the type of database and RefSeq id that is entered. A variable called out_handle is created, which opens the file that the RefSeq information will be written to. The name of this file is given by the file_name argument. The information from net_handle is then read and written into the out_handle variable. The out_handle and net_handle are then closed for the information to be saved to the file. It should be noted, however, that while rettype is set to fasta, it could be changed to 'gb', which is for GenBank.

```Python
def matrix():
    matrix_list = substitution_matrices.load()
    print('The following pre-defined matrices of', ', '.join(matrix_list), 'are available.')
```
The following lines of code return a function that prints out all the available pre-defined substitution matrices.

```Python
def possible_pairwise_alignment(open_gap, extend_gap, matrix, file_1, file_2):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    sequence_1 = SeqIO.read(file_1, "fasta")
    sequence_2 = SeqIO.read(file_2, "fasta")
    alignments = aligner.align(sequence_1, sequence_2)
    print("There are", len(alignments), "possible alignments.")
```
The following lines create a function that will calculate all the possible pairwise alignments for two sequences. The function has five arguments, which are open_gap, extend_gap, matrix, file_1, and file_2. The first and second arguments refer to the penalties for the open and extended gaps. It is important to note that all values for these penalties must be entered as negative numbers. The third argument is the type of matrix that will be used for the alignment. The fourth and fifth arguments refer to the sequence files for the first and second species. The aligner variable is created and uses PairWiseAligner(), which automatically creates an object with all the needed information to perform a pairwise alignment. In addition to this, the open gap and extended gap scores are adjusted for the aligner object. An open gap score refers to the penalty that occurs for the first residue in a gap. An extended gap score refers to the penalty that occurs for any additional residues that occur in a gap. The type of matrix that will be used is entered through the use of the matrix argument. The variable sequence_1 is set equal to reading in file_1 in fasta format. The variable sequence_2 is also set equal to reading in file_2 in fasta format. Since this pairwise alignment is going to align the homologs of different species, file_1 should be one species and file_2 should be another. The variable alignments are set equal to creating all the possible alignments for sequences 1 and 2 for a specific alignment type. The number of possible alignments is then printed out to show how many exist.

```Python
def pairwise_alignment(open_gap, extend_gap, matrix, file_1, file_2, file_name, alignment):
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    sequence_1 = SeqIO.read(file_1, "fasta")
    sequence_2 = SeqIO.read(file_2, "fasta")
    score = aligner.score(sequence_1.seq, sequence_2.seq)
    alignments = aligner.align(sequence_1, sequence_2)
    with open(file_name, 'w') as file:
        file.writelines(["Matrix:", ' ', str(matrix), '\n'])
        file.writelines(["Gap penalty:", ' ', str(abs(open_gap)), '\n'])
        file.writelines(["Extended penalty:", ' ', str(abs(extend_gap)), '\n'])
        file.writelines(["Score:", ' ', str(score), '\n'])
        file.writelines(['\n', str(alignments[alignment])])
```
The following lines create a function that will calculate all the possible pairwise alignments for two sequences and save them to a text file. This function is very similar to the one above, but this one allows for an alignment to be written in a file along with its score. Since there can be hundreds of alignments, it is not reasonable to write them all out at once, since there would be an overload of work for whatever device this code is run on. Due to this, this function allows for a specific alignment to be selected. For this example, there are 480 alignments to choose from. These alignments are selected by using their index. Since there are 480 choices, the index would be 0 to 479. It is important to remember that an index always starts with 0 and is one less than the actual amount. The function has seven arguments, which are open_gap, extend_gap, matrix, file_1, file_2, alignment, file_name. The first and second arguments refer to the penalties for the open and extended gaps. It is important to note that all values for these penalties must be entered as negative numbers. The third argument is the type of matrix that will be used for the alignment. The fourth and fifth arguments refer to the sequence files for the first and second species. The sixth argument refers to the index value entered to get the selected alignment. The seventh argument refers to the name given to the output txt file. The aligner variable is created and uses PairWiseAligner(), which automatically creates an object with all the needed information to perform a pairwise alignment. In addition to this, the open gap and extended gap scores are adjusted for the aligner object. An open gap score refers to the penalty that occurs for the first residue in a gap. An extended gap score refers to the penalty that occurs for any additional residues that occur in a gap. The type of matrix that will be used is entered with the matrix argument. The variable sequence_1 is set equal to reading in file_1 in fasta format. The variable sequence_2 is also set equal to reading in file_2 in fasta format. The variable score is set equal to calculating the alignment score between the specified sequences. Since this pairwise alignment is going to align the homologs of different species, file_1 should be one species and file_2 should be another. The variable alignments are set equal to creating all the possible alignments for sequences 1 and 2 for a specific alignment type. A text file is then opened and given a name through the file_name argument and the specified alignment based on index is written to the file as a string along with its score, gap penalty, extended penalty, and matrix type.

### Function Arguments

```Python
if __name__ == '__main__':
    mart_finder('mart_list.csv')
    database_finder('ENSEMBL_MART_ENSEMBL', 'database_list.csv')
    filters_attributes('hsapiens_gene_ensembl', 'h_filter.csv', 'h_attrib.csv')
    filters_attributes('mmusculus_gene_ensembl', 'm_filter.csv', 'm_attrib.csv')
    dataset_retrieve('hsapiens_gene_ensembl', '5', 'species_1.csv')
    dataset_retrieve('mmusculus_gene_ensembl', '18', 'species_2.csv')
    gene_list('hsapiens_gene_ensembl', '5', 'mmusculus_homolog_ensembl_gene', 'mmusculus_homolog_associated_gene_name',
              'genes.csv')
    gene_list_dataset_1_filter('species_1.csv', 'genes.csv', 'species_1_filter.csv', 'filtered_gene.csv')
    gene_list_dataset_2_filter('species_2.csv', 'filtered_gene.csv', 'Mouse_gene_stable_ID',
                               'species_2_filter_final.csv', 'filtered_gene_final.csv')
    dataset_1_final_filter('species_1_filter.csv', 'filtered_gene_final.csv', 'species_1_filter_final.csv')
    gene_ontology_filter('species_1_filter_final.csv', 'plasma membrane', 'species_1_go.csv')
    gene_ontology_filter('species_2_filter_final.csv', 'plasma membrane', 'species_2_go.csv')
    ref_seq_list('species_1_go.csv', 'APC', 'RefSeq_peptide_ID', 'species_1_ref.csv')
    ref_seq_list('species_2_go.csv', 'Apc', 'RefSeq_peptide_ID', 'species_2_ref.csv')
    ref_seq_sequence('nickolas.parker@maine.edu', 'protein', 'NP_001394379', 'H_APC_ref_seq.fasta')
    ref_seq_sequence('nickolas.parker@maine.edu', 'protein', 'NP_001347909', 'M_Apc_ref_seq.fasta')
    matrix()
    possible_pairwise_alignment(-10, -0.5, 'BLOSUM62', 'H_APC_ref_seq.fasta', 'M_Apc_ref_seq.fasta')
    pairwise_alignment(-10, -0.5, 'BLOSUM62', 'H_APC_ref_seq.fasta', 'M_Apc_ref_seq.fasta', 'alignment.txt', 0)
```    
The following code calls all the previously mentioned functions. The order in which the functions are called and the order in which the arguments are entered into the functions is important. If the functions are called in the wrong order, this could cause the code to fail and possibly give inaccurate results as well as extra or unwanted files. The entered information for each function corresponds to a specified argument within the function. Entering wrong or extra information into the arguments will interfere with the pipeline process, as it will cause the functions to return an error message. Some of the functions are called more than once, since they deal with both species. Other functions deal with one species, so their functionality is more limited in comparison to the ones that can deal with both species. It should be noted, however, that some functions have the ability to be run multiple times without interfering with anything. The best example is the function pairwise_alignment, which can be called 480 times to get all the possible alignments.
