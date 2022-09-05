import KEGGutils as kg
import networkx as nx


def get_disease_gene_graph():
    return kg.KEGGlinkgraph(source_db="disease", target_db="hsa")


def get_drug_disease_graph():
    return kg.KEGGlinkgraph(source_db="disease", target_db="drug")


def get_drug_target_graph():
    return kg.KEGGlinkgraph(source_db="drug", target_db="hsa")


print(kg.__version__)
pathway = kg.KEGGpathway(pathway_id="hsa05215")

print(pathway)
