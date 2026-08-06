from enum import Enum, auto
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Set, Tuple

import networkx as nx
import math
import pyproj

from ..model import Edge, GeoNode, Node, Signal, SignalDirection, Topology
from ..utils.comparison_algorithms.edge_shape import edge_shape_comparison


class CompareMatching:
    def __init__(self):
        self.element_matching = {}
        self.not_found_in_a = []
        self.not_found_in_b = []


class CompareResult:
    def __init__(self):
        self.node_distance: float = -1.0
        self.node_matching: CompareMatching = CompareMatching()
        self.edge_length_difference: float = -1.0
        self.edge_matching: CompareMatching = CompareMatching()
        self.signal_distance: float = -1.0
        self.signal_matching: CompareMatching = CompareMatching()


@dataclass
class GeoTopologyMatchingResult:
    only_in_graph_1: Set[Edge] = field(default_factory=set)
    only_in_graph_2: Set[Edge] = field(default_factory=set)
    contained_in_graph_1: Set[Edge] = field(default_factory=set)
    contained_in_graph_2: Set[Edge] = field(default_factory=set)


class CompareMode(Enum):
    EXACT = auto()
    ISOMORPHIC = auto()
    VORONOI = auto()
    DOUBLE_VORONOI = auto()
    GEOGRAPHIC = auto()
    CONTAINMENT = auto()


class Compare:
    @staticmethod
    def compare(
        topology_a: Topology,
        topology_b: Topology,
        compare_mode: CompareMode,
        given_node_matching: Dict[Node, Node] | None = None,
        exclude_element_list=None,
        skip_signals: bool = False,
    ) -> CompareResult:
        if given_node_matching is None:
            given_node_matching = {}
        if exclude_element_list is None:
            exclude_element_list = []

        result: CompareResult = CompareResult()

        if compare_mode == CompareMode.EXACT:
            Compare._calc_exact_matching(result, topology_a, topology_b)
        if compare_mode == CompareMode.ISOMORPHIC:
            if not Compare._are_topologies_isomorphic(topology_a, topology_b):
                raise ValueError("Both topologies needs to be isomorphic")
            if not given_node_matching:
                raise ValueError(
                    "For isomorphic topologies, at least one mathing node needs to be given."
                )
            Compare._calc_isomorphic_matching(
                result, topology_a, topology_b, given_node_matching, skip_signals
            )
        if compare_mode == CompareMode.VORONOI:
            Compare._calc_voronoi_matching(
                result, topology_a, topology_b, given_node_matching, skip_signals
            )
        if compare_mode == CompareMode.DOUBLE_VORONOI:
            Compare.double_voronoi_matching(
                topology_a, topology_b, given_node_matching, skip_signals
            )
        if compare_mode == CompareMode.GEOGRAPHIC:
            given_node_matching = given_node_matching or {}
            Compare._calc_geographic_matching(
                result, topology_a, topology_b, given_node_matching, skip_signals
            )
        if compare_mode == CompareMode.CONTAINMENT:
            topology_a, topology_b = edge_shape_comparison(topology_a, topology_b)
            Compare._calc_isomorphic_matching(
                result=result, topology_a=topology_a, topology_b=topology_b, given_node_matching=given_node_matching, skip_signals=True
            )
            # result = Compare._geo_based_topology_matching(topology_a, topology_b, buffer_meters=5.0)
            # plot_geo_matching_result(topology_a, topology_b, result)


        result.node_distance = Compare._calc_distance_for_matching(
            result.node_matching, exclude_element_list
        )
        result.edge_length_difference = Compare._calc_distance_for_matching(
            result.edge_matching, exclude_element_list, element_type="edge"
        )
        if not skip_signals:
            result.signal_distance = Compare._calc_distance_for_matching(
                result.signal_matching, exclude_element_list, element_type="signal"
            )
        return result

    @staticmethod
    def _calc_exact_matching(result: CompareResult, topology_a: Topology, topology_b: Topology):
        Compare._calc_exact_element_matching(
            result.node_matching, topology_a.nodes, topology_b.nodes
        )
        Compare._calc_exact_element_matching(
            result.edge_matching, topology_a.edges, topology_b.edges
        )
        Compare._calc_exact_element_matching(
            result.signal_matching, topology_a.signals, topology_b.signals
        )

    @staticmethod
    def _calc_exact_element_matching(
        compare_matching: CompareMatching, element_dict_a: Dict, element_dict_b: Dict
    ):
        for element_uuid_a, element_a in element_dict_a.items():
            if element_uuid_a in element_dict_b:
                compare_matching.element_matching[element_a] = element_dict_b[element_uuid_a]
            else:
                compare_matching.not_found_in_b.append(element_a)
        for element_uuid_b in element_dict_b.keys() - element_dict_a.keys():
            compare_matching.not_found_in_a.append(element_dict_b[element_uuid_b])

    @staticmethod
    def _calc_isomorphic_matching(
        result: CompareResult,
        topology_a: Topology,
        topology_b: Topology,
        given_node_matching: Dict[Node, Node],
        skip_signals: bool,
    ):
        open_nodes: List[Tuple[Node, Node]] = []

        def __add_to_open_nodes(__node_a: Node, __node_b: Node):
            if __node_a is None and __node_b is None:
                raise ValueError(
                    "Graph topology is isomorphic, but railway network graph differs (point is no point)"
                )
            if __node_a is None or __node_b is None:
                raise ValueError(
                    "Graph topology is isomorphic, but railway network graph differs (topology broken)"
                )
            open_nodes.append((__node_a, __node_b))

        def __add_edges_to_matching(__edge_a: Edge, __edge_b: Edge):
            if __edge_a in result.edge_matching.element_matching:
                if result.edge_matching.element_matching[__edge_a] != __edge_b:
                    print(f"bad edge {__edge_a.uuid} {__edge_b.uuid}")
                    raise ValueError(
                        "Graph topology is isomorphic, but railway network graph differs (edge graph broken)"
                    )
            else:
                result.edge_matching.element_matching[__edge_a] = __edge_b

        for node_a, node_b in given_node_matching.items():
            __add_to_open_nodes(node_a, node_b)

        # node and edge matching with a bfs
        while open_nodes:
            current_tuple = open_nodes.pop(0)
            node_a = current_tuple[0]
            node_b = current_tuple[1]
            if node_a in result.node_matching.element_matching:
                if result.node_matching.element_matching[node_a] != node_b:
                    raise ValueError(
                        "Graph topology is isomorphic, but railway network graph differs (connections at points)"
                    )
                continue
            else:
                result.node_matching.element_matching[node_a] = node_b

            __add_to_open_nodes(node_a.connected_on_head, node_b.connected_on_head)
            __add_edges_to_matching(node_a.connected_edge_on_head, node_b.connected_edge_on_head)
            if node_a.is_point():
                __add_to_open_nodes(node_a.connected_on_left, node_b.connected_on_left)
                __add_edges_to_matching(
                    node_a.connected_edge_on_left, node_b.connected_edge_on_left
                )
                __add_to_open_nodes(node_a.connected_on_right, node_b.connected_on_right)
                __add_edges_to_matching(
                    node_a.connected_edge_on_right, node_b.connected_edge_on_right
                )

        if not skip_signals:

            def __add_signal_lists_to_matching(
                signal_list_a: List[Signal], signal_list_b: List[Signal]
            ):
                for i in range(0, len(signal_list_a)):
                    signal_a = signal_list_a[i]
                    signal_b = signal_list_b[i]
                    result.signal_matching.element_matching[signal_a] = signal_b

            # signal matching
            for edge_a in topology_a.edges.values():
                edge_b = result.edge_matching.element_matching[edge_a]
                if not edge_a.signals and not edge_b.signals:
                    continue
                if not edge_a.signals or not edge_b.signals:
                    raise ValueError(
                        f"Signals on edges {edge_a.uuid} and {edge_b.uuid} does not match"
                    )
                in_direction_a = [
                    signal for signal in edge_a.signals if signal.direction == SignalDirection.IN
                ]
                in_direction_a.sort(key=lambda signal: signal.distance_edge)
                other_direction_a = [
                    signal for signal in edge_a.signals if signal.direction == SignalDirection.GEGEN
                ]
                other_direction_a.sort(key=lambda signal: signal.distance_edge)
                in_direction_b = [
                    signal for signal in edge_b.signals if signal.direction == SignalDirection.IN
                ]
                in_direction_b.sort(key=lambda signal: signal.distance_edge)
                other_direction_b = [
                    signal for signal in edge_b.signals if signal.direction == SignalDirection.GEGEN
                ]
                other_direction_b.sort(key=lambda signal: signal.distance_edge)

                if result.node_matching.element_matching[edge_a.node_a] != edge_b.node_a:
                    # edge b is reversed, so switch lists
                    in_direction_b, other_direction_b = list(reversed(other_direction_b)), list(
                        reversed(in_direction_b)
                    )

                if len(in_direction_a) != len(in_direction_b) or len(other_direction_a) != len(
                    other_direction_b
                ):
                    raise ValueError(
                        f"Number of signals on edges {edge_a.uuid} and {edge_b.uuid} differs (per direction)"
                    )
                # TODO: Cases with different amounts of signals on an edge might still be interesting for matching. 
                # TODO: Decide how to treat uneven amount of signals on an edge.
 
                __add_signal_lists_to_matching(in_direction_a, in_direction_b)
                __add_signal_lists_to_matching(other_direction_a, other_direction_b)






#########################################################
# geo based containment matching
#########################################################

    @staticmethod
    def _geo_based_topology_matching(
        topology_1: Topology,
        topology_2: Topology,
        buffer_meters: float = 10.0,
    ) -> GeoTopologyMatchingResult:
        polylines_1 = {
            edge: Compare._get_edge_polyline_wgs84(edge) for edge in topology_1.edges.values()
        }
        polylines_2 = {
            edge: Compare._get_edge_polyline_wgs84(edge) for edge in topology_2.edges.values()
        }

        all_coords = []
        for polyline in list(polylines_1.values()) + list(polylines_2.values()):
            all_coords.extend(polyline)
        if not all_coords:
            return GeoTopologyMatchingResult()
        utm_transformer = Compare._get_utm_transformer(all_coords)

        projected_1 = {
            edge: [utm_transformer.transform(x, y) for x, y in polyline]
            for edge, polyline in polylines_1.items()
        }
        projected_2 = {
            edge: [utm_transformer.transform(x, y) for x, y in polyline]
            for edge, polyline in polylines_2.items()
        }

        ref_polylines_1 = list(projected_1.values())
        only_in_2 = set()
        contained_2 = set()
        for edge, polyline in projected_2.items():
            if Compare._edge_fully_contained(polyline, ref_polylines_1, buffer_meters):
                contained_2.add(edge)
            else:
                only_in_2.add(edge)

        ref_polylines_2 = list(projected_2.values())
        only_in_1 = set()
        contained_1 = set()
        for edge, polyline in projected_1.items():
            if Compare._edge_fully_contained(polyline, ref_polylines_2, buffer_meters):
                contained_1.add(edge)
            else:
                only_in_1.add(edge)

        geo_topology_matching_result = GeoTopologyMatchingResult(
            only_in_graph_1=only_in_1,
            only_in_graph_2=only_in_2,
            contained_in_graph_1=contained_1,
            contained_in_graph_2=contained_2,
        )


        # remove nodes based on the geographical edge matching?
        # There can be nodes, that are part of an "only_in_one" edge, that are in both topologies. For example, former end nodes that were extended by two branches. TODO: What happens if branches are extended?
        # Only remove a node if removing all "only_in_one" edges results in a non-legal node state?

        # This approach might remove nodes that should not be removed:
        # C __                         ___C
        # A __\____ B      <->    A __/_____ B
        #
        # Sollte hier die Weiche gemappt werden oder nicht?

        print("before removal topology 2›: ", len(topology_1.nodes))
        for node in list(topology_1.nodes.values()):
            remaining_edges = {edge for edge in node.connected_edges} - set(only_in_1)
            if len(remaining_edges) != 1 and len(remaining_edges) != 3:
                topology_1.nodes.pop(node.uuid)
        print("after removal topology 1: ", len(topology_1.nodes))

        print("before removal topology 2: ", len(topology_2.nodes))
        for node in list(topology_2.nodes.values()):
            remaining_edges = {edge for edge in node.connected_edges} - set(only_in_2)
            if len(remaining_edges) != 1 and len(remaining_edges) != 3:
                topology_2.nodes.pop(node.uuid)
        print("after removal topology 2: ", len(topology_2.nodes))


        plot_geo_matching_result(topology_1=topology_1, topology_2=topology_2, result=geo_topology_matching_result)

        # Alternatively: Possible to create isomorphic graph by removing the edges that are only present in one of the two graphs?
        # The overall workflow would then look like this:
        # 1. Identify edges that are only present in one topology_a or topology_b based on their geographic shapes. 
        # 2. Remove those edges and their corresponding nodes (this still needs further clarification)
        # 3. Remaining edges should now be present in both graphs (maybe not necessarily, think about two tracks are close to each other, so that the buffer areas overlapt. Switches between those two tracks would not be detected, right?)
        # 4. Either apply double voronoi or isomorphic matching to the network.


        Compare.double_voronoi_matching(topology_a=topology_1, topology_b=topology_2, given_node_matching=dict(), skip_signals=True)





    @staticmethod
    def _get_edge_polyline_wgs84(edge: Edge) -> List[Tuple[float, float]]:
        geo_node_a = edge.node_a.geo_node.to_wgs84()
        geo_node_b = edge.node_b.geo_node.to_wgs84()
        polyline = [(geo_node_a.x, geo_node_a.y)]
        for geo_node in edge.intermediate_geo_nodes:
            gn = geo_node.to_wgs84()
            polyline.append((gn.x, gn.y))
        polyline.append((geo_node_b.x, geo_node_b.y))
        return polyline

    @staticmethod
    def _get_utm_transformer(
        wgs84_coords: List[Tuple[float, float]],
    ) -> pyproj.Transformer:
        mean_lon = sum(x for x, y in wgs84_coords) / len(wgs84_coords)
        mean_lat = sum(y for x, y in wgs84_coords) / len(wgs84_coords)
        utm_zone = int((mean_lon + 180) / 6) + 1
        epsg = 32600 + utm_zone if mean_lat >= 0 else 32700 + utm_zone
        print(epsg)
        return pyproj.Transformer.from_crs(
            "EPSG:4326", f"EPSG:{epsg}", always_xy=True
        )

    @staticmethod
    def _point_segment_distance(
        px: float, py: float, ax: float, ay: float, bx: float, by: float
    ) -> float:
        dx, dy = bx - ax, by - ay
        if dx == 0.0 and dy == 0.0:
            return math.sqrt((px - ax) ** 2 + (py - ay) ** 2)
        t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
        t = max(0.0, min(1.0, t))
        return math.sqrt((px - (ax + t * dx)) ** 2 + (py - (ay + t * dy)) ** 2)

    @staticmethod
    def _point_to_polyline_min_distance(
        px: float, py: float, polyline: List[Tuple[float, float]]
    ) -> float:
        best = float("inf")
        for k in range(len(polyline) - 1):
            d = Compare._point_segment_distance(
                px, py, polyline[k][0], polyline[k][1], polyline[k + 1][0], polyline[k + 1][1]
            )
            if d < best:
                best = d
        return best

    @staticmethod
    def _edge_fully_contained(
        test_polyline: List[Tuple[float, float]],
        ref_polylines: List[List[Tuple[float, float]]],
        buffer_meters: float,
    ) -> bool:
        for px, py in test_polyline:
            min_dist = float("inf")
            for ref_poly in ref_polylines:
                d = Compare._point_to_polyline_min_distance(px, py, ref_poly)
                if d < min_dist:
                    min_dist = d
            if min_dist > buffer_meters:
                return False
        return True










#########################################################
# geo based voronoi matching
#########################################################


    @staticmethod
    def _calc_voronoi_matching(
        result: CompareResult,
        topology_a: Topology,
        topology_b: Topology,
        given_node_matching: Dict[Node, Node],
        skip_signals: bool,
    ):
        print("calculating voronoi distance")
        print("num nodes in topology a: ", len(topology_a.nodes))
        print("num nodes in topology b: ", len(topology_b.nodes))
        # smaller_topology = topology_a if len(topology_a.nodes) < len(topology_b.nodes) else topology_b
        # larger_topology = topology_b if smaller_topology == topology_a else topology_a

        voronoi_map: dict[Node, list[Node]] = {}

        #TOOD: the following implementation can get very slow for large datasets!!
        for node in topology_a.nodes.values():
            print("node to map: ", node)
            node_coordinates = [node.geo_node.x, node.geo_node.y]
            smallest_distance = math.inf
            matching_reference_node = None
            for reference_node in topology_b.nodes.values():
                reference_node_coordinates = [reference_node.geo_node.x, reference_node.geo_node.y]
                distance = math.dist(node_coordinates, reference_node_coordinates)
                if distance < smallest_distance:
                    matching_reference_node = reference_node
                    smallest_distance = distance

            print("matching node: ", matching_reference_node) 
            voronoi_map.setdefault(matching_reference_node, []).append(node)

        print(voronoi_map)
        return voronoi_map


    @staticmethod
    def double_voronoi_matching(
        topology_a: Topology,
        topology_b: Topology,
        given_node_matching: Dict[Node, Node],
        skip_signals: bool,
    ):

        def geo_path_to_ground_truth(correct_nodes: list[Node], node):
            """
            Find path to ground truth node. 
            
            
            """



        # First, calculate Voronoi matching in one direction
        result_1 = CompareResult()
        map_1 = Compare._calc_voronoi_matching(result=result_1, topology_a=topology_a, topology_b=topology_b, given_node_matching=given_node_matching, skip_signals=skip_signals)
        map_2 = Compare._calc_voronoi_matching(result=result_1, topology_a=topology_b, topology_b=topology_a, given_node_matching={v: k for k, v in given_node_matching.items()}, skip_signals=skip_signals)
        locked_nodes_a = given_node_matching.keys()
        locked_nodes_b = given_node_matching.values()

        guaranteed = {
            a: bs[0]
            for a, bs in map_1.items()
            if len(bs) == 1 and map_2.get(bs[0]) == [a]
        }

        # Dieser Ansatz kann natürlich immer noch für falsche Mappings sorgen:
        # Wenn node_a nur in topologie_a und node_b nur in topologie_b ist, beide aber zueinander die nächsten Nodes sind. 
        # Zumindest, wenn es keine anderen Nodes als Kandidaten gibt, würden dann node_a <-> node_b als garantiertes mapping
        # ausgegeben, was aber natürlich trotzdem nicht richtig wäre. Hier bräuchte es dann doch wieder eine topologische 
        # Überprüfung. Vielleicht würde auch eine Kombination mit dem Containment-Fall funktionieren. 

        # for node, candidates in map_1:
                
        print("maps")
        print("map_1: ", map_1)
        print("map_2: ", map_2)
        print("guaranteed: ", {a.name: b.name for a, b in guaranteed.items()})

        # Second, calculate Voronoi matching in other direction

        # Third, check if potentially added nodes lie on existing track. 
        # This can indicate a newly added switch on existing infrastructure. This switch then likely connects to some newly added end node via a new edge.
        # The goal would therefore be, to try matching a geo_node?
        # Also, at the lates then check if there is an edge that corresponds to the split edges created by the new switch. 







# experimental edge based mapping

    @staticmethod
    def _get_edge_polyline(edge: Edge) -> List[Tuple[float, float]]:
        polyline = [(edge.node_a.geo_node.x, edge.node_a.geo_node.y)]
        for geo_node in edge.intermediate_geo_nodes:
            polyline.append((geo_node.x, geo_node.y))
        polyline.append((edge.node_b.geo_node.x, edge.node_b.geo_node.y))
        return polyline

    @staticmethod
    def _frechet_distance(
        poly_a: List[Tuple[float, float]], poly_b: List[Tuple[float, float]]
    ) -> float:
        if not poly_a or not poly_b:
            return float("inf")
        n, m = len(poly_a), len(poly_b)
        ca = [[-1.0] * m for _ in range(n)]
        for i in range(n):
            for j in range(m):
                d = math.dist(poly_a[i], poly_b[j])
                if i == 0 and j == 0:
                    ca[i][j] = d
                elif i == 0:
                    ca[i][j] = max(ca[i][j - 1], d)
                elif j == 0:
                    ca[i][j] = max(ca[i - 1][j], d)
                else:
                    ca[i][j] = max(min(ca[i - 1][j], ca[i][j - 1], ca[i - 1][j - 1]), d)
        return ca[n - 1][m - 1]

    @staticmethod
    def _concatenate_polylines(
        polylines: List[List[Tuple[float, float]]],
    ) -> List[Tuple[float, float]]:
        if not polylines:
            return []
        result = list(polylines[0])
        for poly in polylines[1:]:
            result.extend(poly[1:])
        return result

    @staticmethod
    def _buffer_coverage(
        poly_a: List[Tuple[float, float]],
        poly_b: List[Tuple[float, float]],
        distance: float,
        num_samples: int = 100,
    ) -> float:
        if not poly_a or not poly_b:
            return 0.0

        def _point_segment_distance(px, py, ax, ay, bx, by):
            dx, dy = bx - ax, by - ay
            if dx == 0.0 and dy == 0.0:
                return math.dist((px, py), (ax, ay))
            t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
            t = max(0.0, min(1.0, t))
            return math.dist((px, py), (ax + t * dx, ay + t * dy))

        def _point_to_polyline_dist(px, py, poly):
            best = float("inf")
            for k in range(len(poly) - 1):
                d = _point_segment_distance(px, py, poly[k][0], poly[k][1], poly[k + 1][0], poly[k + 1][1])
                if d < best:
                    best = d
            return best

        def _sample_polyline(poly, n):
            if len(poly) <= 1:
                return poly[:n] if n > 0 else [poly[0]]
            seg_lengths = []
            total = 0.0
            for k in range(len(poly) - 1):
                d = math.dist(poly[k], poly[k + 1])
                seg_lengths.append(d)
                total += d
            if total == 0.0:
                return [poly[0]] * min(n, 1)
            samples = [poly[0]]
            for s in range(1, n - 1):
                target = total * s / (n - 1)
                accumulated = 0.0
                for k in range(len(seg_lengths)):
                    accumulated += seg_lengths[k]
                    if accumulated >= target or abs(accumulated - target) < 1e-10:
                        seg_start = accumulated - seg_lengths[k]
                        t = (target - seg_start) / seg_lengths[k] if seg_lengths[k] > 0 else 0.0
                        x = poly[k][0] + t * (poly[k + 1][0] - poly[k][0])
                        y = poly[k][1] + t * (poly[k + 1][1] - poly[k][1])
                        samples.append((x, y))
                        break
            if n > 1:
                samples.append(poly[-1])
            return samples

        samples_a = _sample_polyline(poly_a, min(num_samples, max(len(poly_a), 1)))
        covered = 0
        for sx, sy in samples_a:
            d = _point_to_polyline_dist(sx, sy, poly_b)
            if d <= distance:
                covered += 1
        return covered / len(samples_a) if samples_a else 0.0

    @staticmethod
    def _build_edge_adjacency(topology: Topology) -> Dict[Node, List[Edge]]:
        adj = defaultdict(list)
        for edge in topology.edges.values():
            adj[edge.node_a].append(edge)
            adj[edge.node_b].append(edge)
        return adj

    @staticmethod
    def _generate_edge_chains_from_node(
        seed_node: Node,
        adj: Dict[Node, List[Edge]],
        used_edges: set,
        max_length: int,
    ) -> List[List[Edge]]:
        chains = []

        def _extend(chain, at_node):
            if len(chain) <= max_length:
                chains.append(list(chain))
            if len(chain) >= max_length:
                return
            for next_edge in adj.get(at_node, []):
                if next_edge in used_edges or next_edge in chain:
                    continue
                next_node = next_edge.get_opposite_node(at_node)
                chain.append(next_edge)
                _extend(chain, next_node)
                chain.pop()

        for start_edge in adj.get(seed_node, []):
            if start_edge in used_edges:
                continue
            next_node = start_edge.get_opposite_node(seed_node)
            _extend([start_edge], next_node)

        return chains

    @staticmethod
    def _polyline_length(poly: List[Tuple[float, float]]) -> float:
        total = 0.0
        for k in range(len(poly) - 1):
            total += math.dist(poly[k], poly[k + 1])
        return total

    @staticmethod
    def _chain_endpoints(chain: List[Edge]) -> Tuple[Node, Node]:
        if len(chain) == 1:
            return chain[0].node_a, chain[0].node_b

        endpoint_nodes = {chain[0].node_a, chain[0].node_b}
        for edge in chain[1:]:
            shared = None
            new = None
            for node in (edge.node_a, edge.node_b):
                if node in endpoint_nodes:
                    shared = node
                else:
                    new = node
            if shared is not None and new is not None:
                endpoint_nodes.remove(shared)
                endpoint_nodes.add(new)

        result = list(endpoint_nodes)
        return result[0], result[1] if len(result) == 2 else (result[0], result[0])

    @staticmethod
    def _calc_geographic_matching(
        result: CompareResult,
        topology_a: Topology,
        topology_b: Topology,
        given_node_matching: Dict[Node, Node],
        skip_signals: bool,
    ):
        # Configurable parameters
        frechet_threshold = 50.0
        max_chain_length = 5
        chain_frechet_threshold = 30.0

        # 1. Extract polylines
        polylines_a = {e: Compare._get_edge_polyline(e) for e in topology_a.edges.values()}
        polylines_b = {e: Compare._get_edge_polyline(e) for e in topology_b.edges.values()}

        # 2. Pairwise Frechet distance
        frechet_scores = {}
        for ea in topology_a.edges.values():
            for eb in topology_b.edges.values():
                d = Compare._frechet_distance(polylines_a[ea], polylines_b[eb])
                frechet_scores[(ea, eb)] = d

        # 3. 1:1 edge matching with length ratio guard
        edge_match_1to1 = {}
        used_edges_b = set()

        def _length_ratio(a, b):
            if b == 0.0:
                return 0.0
            return a / b

        sorted_pairs = sorted(frechet_scores.keys(), key=lambda p: frechet_scores[p])
        for ea, eb in sorted_pairs:
            if ea in edge_match_1to1 or eb in used_edges_b:
                continue
            d_frechet = frechet_scores[(ea, eb)]
            if d_frechet > frechet_threshold:
                continue
            len_a = Compare._polyline_length(polylines_a[ea])
            len_b = Compare._polyline_length(polylines_b[eb])
            ratio = _length_ratio(max(len_a, len_b), min(len_a, len_b))
            if ratio > 1.2:
                continue
            edge_match_1to1[ea] = eb
            used_edges_b.add(eb)

        # 4. Derive node matching from 1:1 edges + given seeds
        node_matching = dict(given_node_matching)

        def _node_distance_xy(na, nb):
            return math.dist(
                (na.geo_node.x, na.geo_node.y),
                (nb.geo_node.x, nb.geo_node.y),
            )

        for ea, eb in edge_match_1to1.items():
            na1, na2 = ea.node_a, ea.node_b
            nb1, nb2 = eb.node_a, eb.node_b
            d_direct = _node_distance_xy(na1, nb1) + _node_distance_xy(na2, nb2)
            d_reverse = _node_distance_xy(na1, nb2) + _node_distance_xy(na2, nb1)
            if d_direct <= d_reverse:
                if na1 not in node_matching and nb1 not in node_matching.values():
                    node_matching[na1] = nb1
                if na2 not in node_matching and nb2 not in node_matching.values():
                    node_matching[na2] = nb2
            else:
                if na1 not in node_matching and nb2 not in node_matching.values():
                    node_matching[na1] = nb2
                if na2 not in node_matching and nb1 not in node_matching.values():
                    node_matching[na2] = nb1

        # 5. Concatenation detection
        adj_b = Compare._build_edge_adjacency(topology_b)
        adj_a = Compare._build_edge_adjacency(topology_a)
        concat_matches = {}  # edge_a -> List[edge_b]

        # 5a. Match unmatched A-edges against chains of B-edges (1:N)
        for ea in topology_a.edges.values():
            if ea in edge_match_1to1:
                continue
            na1, na2 = ea.node_a, ea.node_b
            best_chain = None
            best_frechet = float("inf")

            for na in [na1, na2]:
                matched_nb = node_matching.get(na)
                if matched_nb is None:
                    continue

                chains = Compare._generate_edge_chains_from_node(
                    matched_nb, adj_b, used_edges_b, max_chain_length
                )
                for chain in chains:
                    if not chain:
                        continue
                    if len(chain) < 2:
                        continue
                    chain_poly = Compare._concatenate_polylines([polylines_b[e] for e in chain])
                    d = Compare._frechet_distance(polylines_a[ea], chain_poly)
                    if d < best_frechet:
                        best_frechet = d
                        best_chain = chain

            if best_chain and best_frechet <= chain_frechet_threshold:
                concat_matches[ea] = best_chain

        # 5b. Match unmatched B-edges against chains of A-edges (M:1)
        nb_to_na = {v: k for k, v in node_matching.items()}
        for eb in topology_b.edges.values():
            if eb in used_edges_b:
                continue
            nb1, nb2 = eb.node_a, eb.node_b
            for nb in [nb1, nb2]:
                matched_na = nb_to_na.get(nb)
                if matched_na is None:
                    continue

                chains = Compare._generate_edge_chains_from_node(
                    matched_na, adj_a, set(edge_match_1to1.keys()), max_chain_length
                )
                for chain in chains:
                    if not chain:
                        continue
                    if len(chain) < 2:
                        continue
                    chain_poly = Compare._concatenate_polylines([polylines_a[e] for e in chain])
                    d = Compare._frechet_distance(chain_poly, polylines_b[eb])
                    if d <= chain_frechet_threshold:
                        for ea in chain:
                            if ea not in edge_match_1to1:
                                edge_match_1to1[ea] = eb
                        break

        # 6. Update node matching from concatenation results
        for ea, chain in concat_matches.items():
            na1, na2 = ea.node_a, ea.node_b
            nb_first, nb_last = Compare._chain_endpoints(chain)

            d1 = _node_distance_xy(na1, nb_first) + _node_distance_xy(na2, nb_last)
            d2 = _node_distance_xy(na1, nb_last) + _node_distance_xy(na2, nb_first)
            if d1 <= d2:
                if na1 not in node_matching and nb_first not in node_matching.values():
                    node_matching[na1] = nb_first
                if na2 not in node_matching and nb_last not in node_matching.values():
                    node_matching[na2] = nb_last
            else:
                if na1 not in node_matching and nb_last not in node_matching.values():
                    node_matching[na1] = nb_last
                if na2 not in node_matching and nb_first not in node_matching.values():
                    node_matching[na2] = nb_first

        # 7. Populate result
        for na, nb in node_matching.items():
            result.node_matching.element_matching[na] = nb
        for na in topology_a.nodes.values():
            if na not in node_matching:
                result.node_matching.not_found_in_b.append(na)
        for nb in topology_b.nodes.values():
            if nb not in node_matching.values():
                result.node_matching.not_found_in_a.append(nb)

        used_edges_a = set(edge_match_1to1.keys())
        for ea, eb in edge_match_1to1.items():
            result.edge_matching.element_matching[ea] = eb
        for ea in topology_a.edges.values():
            if ea not in used_edges_a:
                result.edge_matching.not_found_in_b.append(ea)
        for eb in topology_b.edges.values():
            if eb not in used_edges_b:
                result.edge_matching.not_found_in_a.append(eb)

        if not skip_signals:
            signal_dist = 5.0
            for ea, eb in edge_match_1to1.items():
                for sa in ea.signals:
                    for sb in eb.signals:
                        dist_s = abs(sa.distance_edge - sb.distance_edge)
                        if dist_s <= signal_dist and sa.direction == sb.direction:
                            result.signal_matching.element_matching[sa] = sb
                            break
            for sa in topology_a.signals.values():
                if sa not in result.signal_matching.element_matching:
                    result.signal_matching.not_found_in_b.append(sa)
            for sb in topology_b.signals.values():
                if sb not in result.signal_matching.element_matching.values():
                    result.signal_matching.not_found_in_a.append(sb)









    @staticmethod
    def _are_topologies_isomorphic(topology_a: Topology, topology_b: Topology):
        if len(topology_a.nodes) != len(topology_b.nodes) or len(topology_a.edges) != len(
            topology_b.edges
        ):
            # Catch easy case before running expensive network x lib
            return False
        graph_a = topology_a.to_networkx_graph()
        graph_b = topology_b.to_networkx_graph()
        return nx.is_isomorphic(graph_a, graph_b)

    @staticmethod
    def _calc_distance_for_matching(
        matching: CompareMatching,
        exclude_element_list,
        element_type: str = "node",
    ):
        if not matching.element_matching:
            return -1.0

        distance_sum: float = 0.0
        for element_a in matching.element_matching:
            element_b = matching.element_matching[element_a]
            if element_a in exclude_element_list or element_b in exclude_element_list:
                continue

            if element_type == "node":
                geo_node_a: GeoNode = element_a.geo_node
                geo_node_b: GeoNode = element_b.geo_node
                distance = abs(geo_node_a.get_distance_to_other_geo_node(geo_node_b))
                print(f"From {element_a.uuid} to {element_b.uuid}: {distance}")
                distance_sum += distance
            elif element_type == "edge":
                print(
                    f"Edge {element_a.uuid} compared to {element_b.uuid}: {abs(element_a.length - element_b.length)}"
                )
                distance_sum += abs(element_a.length - element_b.length)
            elif element_type == "signal":
                x_a, y_a = element_a.get_calculated_coordinates()
                geo_node_a = GeoNode.get_new_geo_node_same_type(
                    element_a.edge.node_a.geo_node, x_a, y_a
                )
                x_b, y_b = element_b.get_calculated_coordinates()
                geo_node_b = GeoNode.get_new_geo_node_same_type(
                    element_b.edge.node_a.geo_node, x_b, y_b
                )
                distance_sum += geo_node_a.get_distance_to_other_geo_node(geo_node_b)

        return distance_sum


def plot_geo_matching_result(
    topology_1: Topology,
    topology_2: Topology,
    result: GeoTopologyMatchingResult,
    save_path: str = None,
):
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    def _plot_polyline(poly, color, linewidth=1.0, alpha=1.0):
        xs = [p[0] for p in poly]
        ys = [p[1] for p in poly]
        ax.plot(xs, ys, color=color, linewidth=linewidth, alpha=alpha)

    print(result)
    for edge in result.contained_in_graph_1:
        _plot_polyline(Compare._get_edge_polyline(edge), "gray", 1.0, 0.5)
    for edge in result.contained_in_graph_2:
        _plot_polyline(Compare._get_edge_polyline(edge), "gray", 1.0, 0.5)

    for edge in result.only_in_graph_1:
        _plot_polyline(Compare._get_edge_polyline(edge), "red", 2.0)
    for edge in result.only_in_graph_2:
        _plot_polyline(Compare._get_edge_polyline(edge), "blue", 2.0)

    legend_elements = [
        Line2D([0], [0], color="gray", linewidth=1, alpha=0.5, label="In both graphs"),
        Line2D([0], [0], color="red", linewidth=2, label="Only in Graph 1"),
        Line2D([0], [0], color="blue", linewidth=2, label="Only in Graph 2"),
    ]
    ax.legend(handles=legend_elements)
    ax.set_aspect("equal")
    ax.set_title("Geo-based Topology Matching")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()
