# Python Pipeline Tutorial to Query and Manipulate Data from BioMart and RefSeq in Order to Perform a Pairwise Alignment Between Two Species 

### Overview

The following tutorial gives a step by step guide on how to successfully use the pipeline in order to get the desired results. The Python script without the explanations has also been provided and can be found above.

### Requirements

The following packages are required for this pipeline.

```Python
# Install packages
!pip install pybiomart --quiet
!pip install biopython --quiet
```
```Python
# Load packages
from pybiomart import Dataset, Server
import pandas as pd
from Bio import Entrez, Align, SeqIO
from Bio.Align import substitution_matrices
```
### Find All Marts

```Python
def mart_finder(file_name_1):
    server = Server(host='http://www.ensembl.org')
    list_1 = server.list_marts()
    list_1.to_csv(file_name_1, index=False)
```
The following lines create a function that will return a csv that contains a list of all the possible BioMart Ensembl marts. The function has one argument. The argument file_name_1 refers to the name of the csv file for all possible marts. The list of marts was saved as list_1 and converted to a csv file. The list of all possible marts can be used to determine what mart will be entered for the mart_name argument of the database_finder function.

### Find All BioMart Ensembl Databases 
```Python
def database_finder(mart_name, file_name_2):
    server = Server(host='http://www.ensembl.org')
    mart = server[mart_name]
    list_2 = mart.list_datasets()
    list_2.to_csv(file_name_2, index=False)
```
The following lines create a function that will return a csv that contains a list of all the possible BioMart Ensembl databases. The function has two arguments. The first argument mart_name refers to the name of the selected mart. This mart can be found from the list that was created from the mart_finder function. The second argument file_name_2 refers to the name of the csv file for possible datasets. The list of datasets is saved as list_2 and converted to a csv file. The list of all possible databases can be used to determine what species will be entered for the species argument of the dataset_retrieve function.

### Find Filters and Attributes For a Specific BioMart Ensembl Database

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

# Query Data From a Specific BioMart Ensembl Database

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


```Python
def dataset_1_final_filter(species, gene_list, file_name):
    dataset = pd.read_csv(species)
    genes = pd.read_csv(gene_list)
    list_1 = list(genes['Gene_name'].unique())
    dataset_query = dataset.query("Gene_name in @list_1")
    dataset_query.to_csv(file_name, index=False)
```


```Python
def gene_ontology_filter(file, go_term, go_name_filter):
    filtered_species = pd.read_csv(file)
    query = filtered_species[filtered_species['GO_term_name'].isin([go_term])]
    query.to_csv(go_name_filter, index=False)
```


```Python
def ref_seq_list(file_name, gene, column_name, name):
    file = pd.read_csv(file_name)
    filtered = file.loc[:, [column_name, "Gene_name"]]
    filtered.drop_duplicates(keep='first', inplace=True)
    list = [gene]
    desired_gene = filtered.query("Gene_name in @list")
    desired_gene.to_csv(name, index=False)
```


```Python
def ref_seq_sequence(email, db_type, id, file_name):
    Entrez.email = email
    net_handle = Entrez.efetch(db=db_type, id=id, rettype='fasta', retmode='text')
    out_handle = open(file_name, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()
```


```Python
def matrix():
    matrix_list = substitution_matrices.load()
    print('The following pre-defined matrices of', ', '.join(matrix_list), 'are available.')
```


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


### Function Arguments

```Python
# Functions are called to access pipeline 
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
