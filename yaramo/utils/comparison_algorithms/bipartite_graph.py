from ...model import Topology, Node
import networkx as nx
import math
import itertools
from ..comparison_algorithms.comparison_utils import visualize_topologies
# from ...operations.compare import CompareResult

def calculate_distance_matrix(topology_a: Topology, topology_b: Topology) -> dict[Node, dict[Node, float]]:
    max_dist: float = 0.0

    matrix: dict = {}
    for node_a in topology_a.nodes.values():
        for node_b in topology_b.nodes.values():
            dist = node_a.geo_node.get_distance_to_other_geo_node(node_b.geo_node)
            if dist > max_dist:
                max_dist = dist
            matrix.setdefault(node_a, {})[node_b] = dist
    return matrix, max_dist


def get_distance_cost(node_a: Node, node_b: Node, distance_matrix: dict[Node, dict[Node, float]], max_distance: float) -> float:
    return distance_matrix[node_a][node_b] / max_distance # normalize with max_distance


def get_degree_cost(node_a: Node, node_b: Node) -> float:
    if node_a.is_point() == node_b.is_point():
        return 0.0
    return 1.0


def get_edge_angles(node: Node):
    angles = []
    for edge in node.connected_edges:
        next_geo_node = edge.get_next_geo_node(node)
        x = next_geo_node.x - node.geo_node.x
        y = next_geo_node.y - node.geo_node.y
        angle = math.degrees(math.atan2(y, x))# / 360
        angles.append(angle)
    return sorted(angles)

def angle_difference(angle_1: float, angle_2: float):
    return abs((angle_1 - angle_2 + 180) % 360 - 180)


def get_orientation_cost(node_a: Node, node_b: Node) -> float:
    node_a_angles = get_edge_angles(node_a)
    node_b_angles = get_edge_angles(node_b)

    best_cost = float("inf")
    
    if len(node_a_angles) != len(node_b_angles):
        for angle_a in node_a_angles:
            for angle_b in node_b_angles:
                cost = angle_difference(angle_a, angle_b) / 180
                best_cost = min(cost, best_cost)
        return best_cost / 3.0 # normalize because only one edge can be matched


    for permutation in itertools.permutations(node_b_angles):
        cost = sum(
            angle_difference(a, b)
            for a, b in zip(node_a_angles, permutation)
        )

        cost /= len(node_a_angles)

        cost /= 180.0 # normalize

        best_cost = min(best_cost, cost)

    return best_cost



def get_total_cost(node_a: Node, node_b: Node, distance_matrix: dict[Node, dict[Node, float]], max_distance: float) -> float:
    distance_cost = get_distance_cost(node_a, node_b, distance_matrix, max_distance)
    degree_cost = get_degree_cost(node_a, node_b)
    orientation_cost = get_orientation_cost(node_a, node_b)

    return (
        0.5 * distance_cost + 
        0.2 * degree_cost + 
        0.3 * orientation_cost
    )


def create_bipartite_graph_from_topologies(topology_a: Topology, topology_b: Topology):
    distance_matrix, max_distance = calculate_distance_matrix(topology_a, topology_b)
    graph = nx.Graph()
    graph.add_nodes_from(topology_a.nodes.values(), bipartite=0)
    graph.add_nodes_from(topology_b.nodes.values(), bipartite=1)

    for node_a in topology_a.nodes.values():
        for node_b in topology_b.nodes.values():
            cost = get_total_cost(node_a, node_b, distance_matrix, max_distance)
            graph.add_edge(node_a, node_b, weight=cost)

    return graph


def _calc_bipartite_matching(result: "CompareResult", topology_a: Topology, topology_b: Topology):


    graph = create_bipartite_graph_from_topologies(topology_a, topology_b)

    matching: set[tuple[Node, Node]] = nx.min_weight_matching(graph, weight="weight")

    print(matching)

    for pair in matching:
        node_1 = pair[0]
        node_2 = pair[1]
        node_a = node_1 if node_1 in topology_a.nodes.values() else node_2
        node_b = node_2 if node_2 in topology_b.nodes.values() else node_1
        result.node_matching.element_matching[node_a] = node_b
        # visualize_topologies(topology_a, topology_b, highlight_nodes_a={node_a}, highlight_nodes_b={node_b})

    print(f"element matching: {result.node_matching.element_matching}")

    result.node_matching.not_found_in_b = list(set(topology_a.nodes.values()) - set(result.node_matching.element_matching.keys()))
    result.node_matching.not_found_in_a = list(set(topology_b.nodes.values()) - set(result.node_matching.element_matching.values()))

    print(f"not found in b: {result.node_matching.not_found_in_b}")

    visualize_topologies(topology_a, topology_b, title_a="not found in b", title_b="not found in a", highlight_nodes_a=set(result.node_matching.not_found_in_b), highlight_nodes_b=set(result.node_matching.not_found_in_a))
    # visualize_topologies(topology_a, topology_b, title_b="not found in a", highlight_nodes_b=set(result.node_matching.not_found_in_a))

    return matching
