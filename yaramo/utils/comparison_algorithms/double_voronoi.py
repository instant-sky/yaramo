from ...model import Topology, Node
# from ...operations.compare import CompareResult

import math

def _calc_voronoi_map(
        topology_a: Topology,
        topology_b: Topology,
) -> dict[Node, list[Node]]:
    """
    Maps each node in topology_b to its geographically closes node in 
    topology_a. 

    Returns the inversed map. Keys are the nodes_a of topology_a, mapping
    to a list of nodes_b of topology_b, of which node_a is the closest
    by distance. 
    """

    # TODO:
    # Maybe, not only distance but also orientation should already 
    # matter here?
    # So not just take the closes node, but maybe rather take the closest
    # node, that fullfills other criteria, like "orientation", and 
    
    voronoi_map: dict[Node, list[tuple[Node, float]]] = {}

    for node in topology_a.nodes.values():
        node_coordinates = [node.geo_node.x, node.geo_node.y]
        smallest_distance = math.inf
        matching_reference_node = None

        for reference_node in topology_b.nodes.values():
            reference_node_coordinates = [reference_node.geo_node.x, reference_node.geo_node.y]
            distance = math.dist(node_coordinates, reference_node_coordinates)
            if distance < smallest_distance:
                matching_reference_node = reference_node
                smallest_distance = distance

        voronoi_map.setdefault(matching_reference_node, []).append((node, smallest_distance))

    return voronoi_map



def get_best_candidates(voronoi_map: dict[Node, list[tuple[Node, float]]]) -> dict[Node, Node]:
    for reference_node, candidate_node_tuples in voronoi_map.items():
        candidate_node_tuples.sort(key=lambda t: t[1]) # sort candidate nodes by distance
        voronoi_map[reference_node] = candidate_node_tuples[0][0] # for now, just select the closest node
    return voronoi_map

def double_voronoi_matching(
        topology_a: Topology,
        topology_b: Topology,
):
    """

    """
    voronoi_map_a = _calc_voronoi_map(topology_a, topology_b)
    voronoi_map_b = _calc_voronoi_map(topology_b, topology_a)

    voronoi_map_a = get_best_candidates(voronoi_map_a)
    voronoi_map_b = get_best_candidates(voronoi_map_b)

    if voronoi_map_a == {c: r for r, c in voronoi_map_b.items()}: # Currently, this can only be the case in cases, where both topologies have the same amount of nodes.
        print("same voronoi mapping in both directions")

    print(voronoi_map_a)
    print(voronoi_map_b)