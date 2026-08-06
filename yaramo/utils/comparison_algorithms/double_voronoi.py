from ...model import Topology, Node
from ...operations.compare import CompareResult

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
    
    voronoi_map: dict[Node, list[Node]] = {}

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

        voronoi_map.setdefault(matching_reference_node, [])

    return voronoi_map




def double_voronoi_matching(
        topology_a: Topology,
        topology_b: Topology,
) -> CompareResult:
    """

    """
    voronoi_map_a = _calc_voronoi_map(topology_a, topology_b)
    voronoi_map_b = _calc_voronoi_map(topology_b, topology_a)

    print(voronoi_map_a)
    print(voronoi_map_b)