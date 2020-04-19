## Requirement  
python  
pandas
## Data preprocessing
### 1.gnbr
original gnbr datasets V7:https://zenodo.org/record/1035500#.XlcypZgzZPY  
Treat chemical_disease, gene_disease, gene_gene, gene_disease triples and entities separately  
concat.py:Merge all entities and triples
### 2.drugbank
drugbank.py:Merge the drugbank data set with the gnbr data set (unzip drugbank/drugbank.zip first)
## Train
We used the knowledge graph embedding model Rotate in the DGL framework  
https://github.com/dmlc/dgl/tree/master/apps/kg
