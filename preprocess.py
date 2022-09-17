from dataclasses import dataclass
from enum import Enum
from uuid import uuid4

import KEGGutils as kg
import networkx as nx


class NodeType(Enum):
    GeneProduct = "gene product"
    Compound = "compound"

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, eq=True)
class Node:
    name: str
    type: NodeType


class RelationType(Enum):
    ECrel = "enzyme-enzyme relation"
    PPrel = "protein-protein interaction"
    GErel = "gene expression interaction"
    PCrel = "protein-compound interaction"

    @classmethod
    def parse(cls, value: str):
        return cls[value]

    def __str__(self) -> str:
        return self.value


class RelationSubType(Enum):
    Compound = "compound"
    HiddenCompound = "hidden compound"
    Activation = "activation"
    Inhibition = "inhibition"
    Expression = "expression"
    Repression = "repression"
    IndirectEffect = "indirect effect"
    StateChange = "state change"
    BindingAssociation = "binding/association"
    Dissociation = "dissociation"
    MissingInteraction = "missing interaction"
    Phosphorylation = "phosphorylation"
    Dephosphorylation = "dephosphorylation"
    Glycosylation = "glycosylation"
    Ubiquitination = "ubiquitination"
    Methylation = "methylation"

    @classmethod
    def parse(cls, value: str):
        if value == "indirect":
            return cls.IndirectEffect
        return cls(value)

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, eq=True)
class Relation:
    type: RelationType
    subtype: RelationSubType
    pathway: str


@dataclass(frozen=True, eq=True)
class Reaction:
    reversible: bool
    pathway: str


def parse_pathway_node_type(node_type: str):
    if node_type == "gene":
        return NodeType.GeneProduct
    if node_type == "compound":
        return NodeType.Compound
    return None


def process_pathway(pathway: kg.KEGGpathway) -> nx.MultiDiGraph:
    processed_graph = expand_pathway(pathway)
    remove_pathway_groups(processed_graph)
    remove_isolated_nodes(processed_graph)
    return processed_graph


def expand_pathway(pathway: kg.KEGGpathway) -> nx.MultiDiGraph:
    graph = nx.MultiDiGraph()
    for node_id, node_data in pathway.nodes.data():
        if node_data["nodetype"] == "group":
            graph.add_node(node_id, data=Node(node_id, None))
            continue

        if node_data["nodetype"] not in ["gene", "compound"]:
            continue

        if graph.has_node(node_data["name"]):
            continue

        graph.add_node(
            node_data["name"],
            data=Node(
                node_data["name"],
                parse_pathway_node_type(node_data["nodetype"]),
            ),
        )

    for edge_id, edge_data in pathway.relations.items():
        relation_type = edge_data["relation_type"]
        if relation_type == "maplink":
            continue

        source = edge_data["nodes"][0]
        if source == "undefined":
            source = edge_data["node_ids"][0]

        dest = edge_data["nodes"][1]
        if dest == "undefined":
            dest = edge_data["node_ids"][1]

        relation_subtypes = edge_data["subtypes"]

        for relation_subtype, _ in relation_subtypes:
            graph.add_edge(
                source,
                dest,
                uuid4(),
                data=Relation(
                    RelationType.parse(relation_type),
                    RelationSubType.parse(relation_subtype),
                    pathway.name,
                ),
            )

    for edge_id, edge_data in pathway.reactions.items():
        reversible = edge_data["type"] == "reversible"
        for substrate_id, substrate_names in edge_data["substrates"]:
            substrates = (
                substrate_names.split()
                if substrate_names != "undefined"
                else [
                    node for node in graph.nodes if node.startswith(substrate_id + "_")
                ]
            )

            if substrate_names == "undefined":
                continue

            for product_id, product_names in edge_data["products"]:
                products = (
                    product_names.split()
                    if product_names != "undefined"
                    else [
                        node
                        for node in graph.nodes
                        if node.startswith(product_id + "_")
                    ]
                )

                for substrate in substrates:
                    if not graph.has_node(substrate):
                        continue

                    for product in products:
                        if not graph.has_node(product):
                            continue

                        graph.add_edge(
                            substrate,
                            product,
                            uuid4(),
                            data=Reaction(reversible, pathway.name),
                        )

    return graph


def remove_pathway_groups(graph: nx.MultiDiGraph):
    nodes_to_remove = set()

    for node_id, node_data in graph.nodes.data():
        node_data: Node = node_data["data"]
        if node_data.type != None:
            continue

        for source, _, edge_data in graph.in_edges(node_id, data=True):
            for dest in graph.successors(node_id):
                graph.add_edge(source, dest, uuid4(), data=edge_data["data"])

        for _, dest, edge_data in graph.out_edges(node_id, data=True):
            for source in graph.predecessors(node_id):
                graph.add_edge(source, dest, uuid4(), data=edge_data["data"])

        nodes_to_remove.add(node_id)

    graph.remove_nodes_from(nodes_to_remove)


def remove_isolated_nodes(graph: nx.MultiDiGraph):
    graph.remove_nodes_from(list(nx.isolates(graph)))


def process_link_graph(graph: kg.KEGGlinkgraph) -> nx.MultiDiGraph:
    directed_graph = nx.MultiDiGraph()
    directed_graph.add_nodes_from(graph.nodes)

    for source in graph.source_nodes.keys():
        for dest in graph.neighbors(source):
            directed_graph.add_edge(source, dest, uuid4())

    return directed_graph


if __name__ == "__main__":
    pathway = kg.KEGGpathway(pathway_id="hsa00010")
    print(pathway)
    processed = process_pathway(pathway)
    print(processed)

    pathway = kg.KEGGpathway(pathway_id="hsa05215")
    print(pathway)
    processed = process_pathway(pathway)
    print(processed)

    link_graph = kg.KEGGlinkgraph(source_db="disease", target_db="hsa")
    directed = process_link_graph(link_graph)
    print(link_graph)
