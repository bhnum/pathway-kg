from collections.abc import Iterator
from typing import List

import KEGGutils as kg


def get_disease_gene_graph() -> kg.KEGGlinkgraph:
    return kg.KEGGlinkgraph(source_db="disease", target_db="hsa")


def get_drug_disease_graph() -> kg.KEGGlinkgraph:
    return kg.KEGGlinkgraph(source_db="disease", target_db="drug")


def get_drug_target_graph() -> kg.KEGGlinkgraph:
    return kg.KEGGlinkgraph(source_db="drug", target_db="hsa")


def get_drug_compound_graph() -> kg.KEGGlinkgraph:
    return kg.KEGGlinkgraph(source_db="drug", target_db="compound")


def get_all_pathway_names() -> List[str]:
    return kg.keggapi_list("pathway", "hsa")


def get_all_pathways() -> Iterator[kg.KEGGpathway]:
    names = get_all_pathway_names()
    for name in names:
        name_without_path = name[5:]
        yield kg.KEGGpathway(pathway_id=name_without_path)


if __name__ == "__main__":
    print(kg.__version__)

    print(get_disease_gene_graph())
    print(get_drug_disease_graph())
    print(get_drug_target_graph())
    print(get_drug_compound_graph())
    print(get_all_pathway_names())
    print(next(get_all_pathways()))
