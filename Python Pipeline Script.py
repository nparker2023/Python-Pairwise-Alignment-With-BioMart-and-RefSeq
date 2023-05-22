# Install packages
# This step only applies to installing in Google Collab
# If using another IDE, packages may need to be installed another way
# !pip install pybiomart --quiet
# !pip install biopython --quiet

# Load packages 
from pybiomart import Dataset, Server
import pandas as pd
from Bio import Entrez, Align, SeqIO
from Bio.Align import substitution_matrices


# Finds all possible marts
def mart_finder(file_name_1):
    server = Server(host='http://www.ensembl.org')
    # Mart list is made and saved to file
    list_1 = server.list_marts()
    list_1.to_csv(file_name_1, index=False)


# Finds all BioMart Ensembl databases based on specified mart
# A database represents a species
def database_finder(mart_name, file_name_2):
    server = Server(host='http://www.ensembl.org')
    mart = server[mart_name]
    # Database list is made and saved to file
    list_2 = mart.list_datasets()
    list_2.to_csv(file_name_2, index=False)


# Filters and attributes are found for specified species
# This function is called twice (1 per different species)
def filters_attributes(species, file_1, file_2):
    species_dataset = Dataset(name=species, host='http://www.ensembl.org')
    # List is made for filters and saved to file
    list_1 = species_dataset.filters
    with open(file_1, 'w') as f:
        for item in list_1.keys():
            f.write("%s,%s\n" % (item, list_1[item]))
    # List is made for attributes and saved to file
    list_2 = species_dataset.attributes
    with open(file_2, 'w') as f:
        for item in list_2.keys():
            f.write("%s,%s\n" % (item, list_2[item]))


# Data from BioMart is queried for a specified species
# This function is called twice (1 per different species)
def dataset_retrieve(species, chrom, file_name):
    species_dataset = Dataset(name=species, host='http://www.ensembl.org')
    # Attributes are specified
    species_query = species_dataset.query(
        attributes=['refseq_mrna', 'refseq_peptide', 'ensembl_gene_id', 'external_gene_name', 'description',
                    'start_position', 'end_position', 'strand', 'chromosome_name', 'name_1006'],
        # Filters are specified
        filters={'chromosome_name': [chrom]})
    # Blanks are removed from the dataset
    filtered_set = species_query.dropna()
    # Dataset columns renamed and saved to file
    filtered_set.columns = filtered_set.columns.str.replace(' ', '_')
    filtered_set.to_csv(file_name, index=False)


# Genes are found for a particular species and its corresponding homologs to another specified species
def gene_list(species_1, chrom, species_2_id, species_2_gene_name, file_name):
    # Species 1 is specified
    species_dataset = Dataset(name=species_1, host='http://www.ensembl.org')
    # Attributes include gene name and ID for species_1, and ID and gene name for corresponding homolog species
    gene_list_query = species_dataset.query(
        attributes=['ensembl_gene_id', 'external_gene_name', species_2_id, species_2_gene_name],
        # Genes filtered by chromosomal location
        filters={'chromosome_name': [chrom]})
    # Blanks are removed from the dataset
    filtered_set = gene_list_query.dropna()
    # Dataset columns renamed for easier manipulation and saved to file
    filtered_set.columns = filtered_set.columns.str.replace(' ', '_')
    filtered_set.to_csv(file_name, index=False)


# Queries are filtered out if they don't appear on the gene list for the first species dataset
def gene_list_dataset_1_filter(species, gene_list, species_filter, filter_gene):
    # Species_1 dataset selected
    dataset = pd.read_csv(species)
    # Uses gene list made from gene_list function
    genes = pd.read_csv(gene_list)
    # Removes duplicate gene names from gene list
    list_1 = genes['Gene_name'].unique()
    # Removes genes from species 1 dataset that aren't in list_1 and saves to file
    dataset_query = dataset.query("Gene_name in @list_1")
    dataset_query.to_csv(species_filter, index=False)
    # Removes duplicates from filtered species_1 dataset
    list_2 = dataset_query['Gene_name'].unique()
    # Removes genes from filtered species_1 dataset that aren't in list_2 and saves to file
    genes_filter = genes.query("Gene_name in @list_2")
    genes_filter.to_csv(filter_gene, index=False)


# Queries are filtered out if they don't appear on the gene list for the second species dataset
def gene_list_dataset_2_filter(species, gene_list, column_name, file_name_1, file_name_2):
    # Species_2 dataset is selected
    dataset = pd.read_csv(species)
    # Uses gene list made from gene_list_dataset_1_filter function
    genes = pd.read_csv(gene_list)
    # Removes duplicate gene names from filtered gene list
    list_1 = genes[column_name].unique()
    # Removes genes from species 2 dataset that aren't in list_1 and saves to file
    dataset_query = dataset.query("Gene_stable_ID in @list_1")
    dataset_query.to_csv(file_name_1, index=False)
    # Removes duplicates from filtered species_2 dataset
    list_2 = dataset_query['Gene_stable_ID'].unique()
    # Removes genes from filtered species_2 dataset that aren't in list_2 and saves to file
    dataset_query_2 = genes[genes[column_name].isin(list_2)]
    dataset_query_2.to_csv(file_name_2, index=False)


# First species dataset is updated to reflect filtered second species dataset
def dataset_1_final_filter(species, gene_list, file_name):
    # Use species_1_filter from gene_list_dataset_1_filter function
    dataset = pd.read_csv(species)
    # Uses filter_gene_final from gene_list_dataset_2_filter function
    genes = pd.read_csv(gene_list)
    # Removes duplicate gene names from final filtered gene list
    list_1 = list(genes['Gene_name'].unique())
    # Removes genes from filter species_1 dataset that aren't list_1 and saves to file
    dataset_query = dataset.query("Gene_name in @list_1")
    dataset_query.to_csv(file_name, index=False)


# Finds queries that have a specific gene ontology name
# This function is called twice (1 per different species)
def gene_ontology_filter(file, go_term, go_name_filter):
    # Uses final filter datasets from gene_list_dataset_2_filter and dataset_1_final_filter functions
    filtered_species = pd.read_csv(file)
    # Return dataset entries that have specified GO term and saved to file
    query = filtered_species[filtered_species['GO_term_name'].isin([go_term])]
    query.to_csv(go_name_filter, index=False)


# Returns RefSeqs for a particular gene
# This function is called twice (1 per different species)
def ref_seq_list(file_name, gene, column_name, name):
    # Uses files created from gene_ontology_filter
    file = pd.read_csv(file_name)
    filtered = file.loc[:, [column_name, "Gene_name"]]
    # Duplicate RefSeqs are removed and info saved to file
    filtered.drop_duplicates(keep='first', inplace=True)
    list = [gene]
    desired_gene = filtered.query("Gene_name in @list")
    desired_gene.to_csv(name, index=False)


# Outputs RefSeq based on database type and ID
# This function is called twice (1 per different species)
def ref_seq_sequence(email, db_type, id, file_name):
    # Email address entered in order to use package and access RefSeq
    Entrez.email = email
    # Specified RefSeq is retrieved
    # db can be nucleotide or protein and rettype can be fasta or gb
    net_handle = Entrez.efetch(db=db_type, id=id, rettype='fasta', retmode='text')
    # RefSeq is saved and written to file
    out_handle = open(file_name, "w")
    out_handle.write(net_handle.read())
    out_handle.close()
    net_handle.close()


# Shows all the available substitution matrices
def matrix():
    # List all possible substitution matrices is made
    # Can be used in possible_pairwise_alignment and pairwise_alignment functions
    matrix_list = substitution_matrices.load()
    print('The following pre-defined matrices of', ', '.join(matrix_list), 'are available.')


# Calculates all the possible pairwise alignment variations
def possible_pairwise_alignment(open_gap, extend_gap, matrix, file_1, file_2):
    aligner = Align.PairwiseAligner()
    # Values for open and extended gaps are set
    # Values should be entered as negative integers
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    # Substitution matrix is selected
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    sequence_1 = SeqIO.read(file_1, "fasta")
    sequence_2 = SeqIO.read(file_2, "fasta")
    # Alignment is performed between the two specified sequences
    alignments = aligner.align(sequence_1, sequence_2)
    # Number of possible alignments is reported
    print("There are", len(alignments), "possible alignments.")


# Creates pairwise alignment from two species sequences
def pairwise_alignment(open_gap, extend_gap, matrix, file_1, file_2, file_name, alignment):
    aligner = Align.PairwiseAligner()
    # Values for open and extended gaps are set
    # Values should be entered as negative integers
    aligner.open_gap_score = open_gap
    aligner.extend_gap_score = extend_gap
    # Substitution matrix is selected
    aligner.substitution_matrix = substitution_matrices.load(matrix)
    # Sequences are selected
    sequence_1 = SeqIO.read(file_1, "fasta")
    sequence_2 = SeqIO.read(file_2, "fasta")
    score = aligner.score(sequence_1.seq, sequence_2.seq)
    # Alignment is performed between the two specified sequences
    alignments = aligner.align(sequence_1, sequence_2)
    # File is created to show alignment along with other relevant information
    with open(file_name, 'w') as file:
        file.writelines(["Matrix:", ' ', str(matrix), '\n'])
        file.writelines(["Gap penalty:", ' ', str(abs(open_gap)), '\n'])
        file.writelines(["Extended penalty:", ' ', str(abs(extend_gap)), '\n'])
        file.writelines(["Score:", ' ', str(score), '\n'])
        file.writelines(['\n', str(alignments[alignment])])


if __name__ == '__main__':
    # Functions are called to access pipeline
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
