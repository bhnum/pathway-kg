import networkx as nx

import downloader
from preprocess import process_link_graph, process_pathway


def create_knowledge_graph() -> nx.MultiDiGraph:
    disease_gene = process_link_graph(downloader.get_disease_gene_graph())
    drug_disease = process_link_graph(downloader.get_drug_disease_graph())
    drug_target = process_link_graph(downloader.get_drug_target_graph())
    drug_compound = process_link_graph(downloader.get_drug_compound_graph())

    pathways = map(process_pathway, downloader.get_all_pathways())

    kg = nx.compose_all(
        [*pathways, disease_gene, drug_disease, drug_target, drug_compound]
    )

    return kg


if __name__ == "__main__":
    kg = create_knowledge_graph()
    print(kg)
