from ...model import Topology, Edge, Node, GeoNode, Wgs84GeoNode, EdgeConnectionDirection
from .comparison_utils import visualize_topologies, visualize_topology, get_edge_polyline

import math
import networkx as nx
import matplotlib.pyplot as plt
import pyproj



def point_segment_distance(
        px: float, py: float, ax: float, ay: float, bx: float, by: float
) -> float:
    """
    Return the distance of a point p(px, py) to the segment (line) 
    between point a(ax, ay) and b(bx, by).
    """
    dx, dy = bx - ax, by - ay
    if dx == 0.0 and dy == 0.0: # points a and b have the same coordinates
        return math.sqrt((px - ax)**2 + (py - ay)**2)
    t = ((px - ax) * dx + (py - ay) * dy) / (dx * dx + dy * dy)
    t = max(0.0, min(1.0 , t)) # make sure, that distance to line segment and not line is calculated
    return math.sqrt((px - (ax + t * dx)) ** 2 + (py - (ay + t * dy)) ** 2) # distance between p and projected point on segment (or a/b if those are closer than other points on line segment)

def point_to_polyline_min_distance(px: float, py: float, reference_polyline: list[tuple[float, float]]) -> float:
    best = float("inf")
    for k in range(len(reference_polyline) - 1):
        d = point_segment_distance(
            float(px), float(py), float(reference_polyline[k][0]), float(reference_polyline[k][1]), float(reference_polyline[k+1][0]), float(reference_polyline[k+1][1])
        )
        if d < best:
            best = d
    return best

def node_to_closest_other_node_min_distance(node: Node, topology: Topology) -> float:
    min_distance = float("inf")
    for reference_node in topology.nodes.values():
        if reference_node is node:
            continue
        distance = abs(node.geo_node.get_distance_to_other_geo_node(reference_node.geo_node))#math.dist((node.geo_node.x, node.geo_node.y), (reference_node.geo_node.x, reference_node.geo_node.y))
        if distance < min_distance:
            min_distance = distance
    return min_distance

def _polyline_length(polyline: list[tuple[float, float]]) -> float:
    return sum(
        math.hypot(bx - ax, by - ay)
        for (ax, ay), (bx, by) in zip(polyline, polyline[1:])
    )


def _covered_length(
    polyline: list[tuple[float, float]],
    reference_polylines: list[list[tuple[float, float]]],
    buffer: float = 10,
    step: float | None = None,
) -> float:
    """Total length of ``polyline`` that lies within ``buffer`` of any reference polyline.

    Points are sampled along every segment (default every ``buffer / 2`` meters) so
    that long segments without intermediate vertices are measured correctly.
    """
    if step is None:
        step = buffer / 2
    covered = 0.0
    for (ax, ay), (bx, by) in zip(polyline, polyline[1:]):
        seg_length = math.dist((ax, ay), (bx, by))#math.hypot(bx - ax, by - ay)
        if seg_length == 0.0:
            continue
        n = max(1, math.ceil(seg_length / step))
        for i in range(n):
            t = (i + 0.5) / n
            px = ax + (bx - ax) * t
            py = ay + (by - ay) * t
            if any(
                point_to_polyline_min_distance(px, py, reference_polyline) <= buffer
                for reference_polyline in reference_polylines
            ):
                covered += seg_length / n
    return covered


def edge_overlap_ratio(
    polyline_a: list[tuple[float, float]],
    polyline_b: list[tuple[float, float]],
    buffer: float = None,
    step: float | None = None,
) -> float:
    """Fraction of the shorter of two polylines that lies within ``buffer`` of the other.

    Returns ~1.0 if one polyline is (almost) completely contained in the other, no
    matter how much longer the containing one is. Returns ~0.0 for polylines that
    take completely different routes.
    """
    len_a = _polyline_length(polyline_a)
    len_b = _polyline_length(polyline_b)
    min_len = min(len_a, len_b)
    if min_len == 0.0:
        return 0.0
    covered_a_in_b = _covered_length(polyline_a, [polyline_b], buffer, step)
    covered_b_in_a = _covered_length(polyline_b, [polyline_a], buffer, step)
    return min(1.0, max(covered_a_in_b, covered_b_in_a) / min_len)


def _get_utm_transformer(
    polylines: list[list[tuple[float, float]]],
) -> pyproj.Transformer:
    all_coords = [coord for polyline in polylines for coord in polyline]
    mean_lat = sum(x for x, y in all_coords) / len(all_coords)
    mean_lon = sum(y for x, y in all_coords) / len(all_coords)
    utm_zone = int((mean_lon + 180) / 6) + 1
    epsg = 32600 + utm_zone if mean_lat >= 0 else 32700 + utm_zone
    return pyproj.Transformer.from_crs(
        "EPSG:4326", f"EPSG:{epsg}", always_xy=True
    )


def _project_polyline(
    polyline: list[tuple[float, float]],
    transformer: pyproj.Transformer,
) -> list[tuple[float, float]]:
    return [transformer.transform(y, x) for x, y in polyline]




def _intermediate_geo_nodes_from_node(edge: Edge, from_node: Node) -> list[GeoNode]:
    """Return the intermediate geo nodes of ``edge`` ordered from ``from_node`` towards the other end."""
    if edge.node_a is from_node:
        return list(edge.intermediate_geo_nodes)
    return list(reversed(edge.intermediate_geo_nodes))

def expand_edges_to_remove(edges_to_remove: list[Edge]):
    edges_to_expand = edges_to_remove.copy()
    print(edges_to_expand)
    print(len(edges_to_expand))
    i = 0
    while i < len(edges_to_expand):
        edge = edges_to_expand[i]
        node_a = edge.node_a
        node_b = edge.node_b
        node_a.calc_anschluss_of_all_edges()
        node_b.calc_anschluss_of_all_edges()
        for node in [node_a, node_b]:
            if node.is_point() and edge == node.connected_edge_on_head:
                # print(node.connected_edges)
                left_edge = node.connected_edge_on_left
                right_edge = node.connected_edge_on_right
                print(left_edge)
                print(right_edge)
                for lr_edge in [left_edge, right_edge]:
                    if lr_edge not in edges_to_remove:
                        print("added edge to expandables")
                        edges_to_remove.append(lr_edge)
                        edges_to_expand.append(lr_edge)
        print(len(edges_to_expand))
        print(edges_to_expand)
        i += 1

    return edges_to_remove
        


def clean_topology(topology: Topology, edges_to_remove: list[Edge]):
    print(f"len edges to remove before expanding: {len(edges_to_remove)}")

    # edges_to_remove = expand_edges_to_remove(edges_to_remove)

    print(f"len edges to remove after expanding: {len(edges_to_remove)}")


    for edge in edges_to_remove:
        topology.edges.pop(edge.uuid)
    
    edges_to_remove_per_node: dict[Node, set[Edge]] = {}
    for edge in edges_to_remove:
        edges_to_remove_per_node.setdefault(edge.node_a, set()).add(edge)
        edge.node_a.calc_anschluss_of_all_edges()
        if edge == edge.node_a.connected_edge_on_head:
            edges_to_remove_per_node[edge.node_a].update([edge.node_a.connected_edge_on_left, edge.node_a.connected_edge_on_right])

        edges_to_remove_per_node.setdefault(edge.node_b, set()).add(edge)
        edge.node_b.calc_anschluss_of_all_edges()
        if edge == edge.node_b.connected_edge_on_head:
            edges_to_remove_per_node[edge.node_b].update([edge.node_b.connected_edge_on_left, edge.node_b.connected_edge_on_right])

        
    for node, edges in list(edges_to_remove_per_node.items()):
        # node.calc_anschluss_of_all_edges()
        for edge in edges:
            if edge in node.connected_edges:
                node.remove_edge(edge)
        if len(node.connected_edges) == 3:
            raise ValueError("no edges removed but node was connected to removed edge.")
        if len(node.connected_edges) == 2: # illegal node state -> node needs to be removed, adjacent edges form union
            print(f"removing node and checking edge")
            topology.nodes.pop(node.uuid)
            
            edge_a: Edge = node.connected_edges[0]
            edge_b: Edge = node.connected_edges[1]
            print("remaining edges: ", set([edge_a, edge_b]))
            print("edges connected left, right: ", set([node.connected_edge_on_left, node.connected_edge_on_right]))
            # if set([edge_a, edge_b]) == set([node.connected_edge_on_left, node.connected_edge_on_right]):
            #     raise ValueError("Edge connected on head removed!")
            if edge_a.uuid in topology.edges:
                topology.edges.pop(edge_a.uuid)
            if edge_b.uuid in topology.edges:
                topology.edges.pop(edge_b.uuid)

            node_a: Node = edge_a.get_opposite_node(node)
            node_b: Node = edge_b.get_opposite_node(node)

            new_edge: Edge = Edge(node_a, node_b, intermediate_geo_nodes=_intermediate_geo_nodes_from_node(edge_a, node_a)
                                   + [node.geo_node] 
                                   + _intermediate_geo_nodes_from_node(edge_b, node))
            topology.add_edge(new_edge)

            node_a.replace_edge(edge_a, new_edge)
            node_b.replace_edge(edge_b, new_edge)
            # visualize_topology(
            #     topology,
            #     title=f"merged edges at node {_short_id(node.uuid, False)}",
            #     highlight_nodes={node},
            #     highlight_edges={edge_a, edge_b},
            # )
        if len(node.connected_edges) == 1: # switch -> end, legal
            print("legal node remaining")
            continue
        if len(node.connected_edges) == 0: # all edges removed, just remove node
            print(f"removing node {node}")
            topology.nodes.pop(node.uuid)
            # visualize_topology(topology, title=f"removed node {_short_id(node.uuid, False)}", highlight_nodes={node})

    topology.update_edge_lengths()
    return topology


def edge_shape_comparison(
    topology_a: Topology,
    topology_b: Topology,
    buffer: float = 2, # In Deutschland Mindestgleisabstand 4 Meter. Wenn also buffer ungefähr < 2 m sollten Überleitstellen ausgeschlossen werden.
    overlap_threshold: float = 0.6,
):
    # get polylines of topologies a and b
    raw_polylines_a: dict[Edge, list[tuple[float, float]]] = {
        edge: get_edge_polyline(edge) for edge in topology_a.edges.values()
    }

    raw_polylines_b: dict[Edge, list[tuple[float, float]]] = {
        edge: get_edge_polyline(edge) for edge in topology_b.edges.values()
    }

    polylines_a = raw_polylines_a
    polylines_b = raw_polylines_b

    is_wgs84 = any(
        isinstance(node.geo_node, Wgs84GeoNode)
        for topology in (topology_a, topology_b)
        for node in topology.nodes.values()
    )
    if is_wgs84:
        transformer = _get_utm_transformer(
            list(raw_polylines_a.values()) + list(raw_polylines_b.values())
        )
        polylines_a = {
            edge: _project_polyline(polyline, transformer)
            for edge, polyline in raw_polylines_a.items()
        }
        polylines_b = {
            edge: _project_polyline(polyline, transformer)
            for edge, polyline in raw_polylines_b.items()
        }


    # Get edges that are only present in one topology based on their geographical overlap with the edges of the other topology 
    only_in_a: list[Edge] = []
    only_in_b: list[Edge] = []

    for edge, polyline in polylines_a.items():
        if not any(
            edge_overlap_ratio(polyline, reference_polyline, buffer=buffer) >= overlap_threshold 
            for reference_polyline in polylines_b.values()
        ):
            only_in_a.append(edge)

    for edge, polyline in polylines_b.items():
        if not any(
            edge_overlap_ratio(polyline, reference_polyline, buffer=buffer) >= overlap_threshold
            for reference_polyline in polylines_a.values()
        ):
            only_in_b.append(edge)


    only_in_a = expand_edges_to_remove(only_in_a)
    only_in_b = expand_edges_to_remove(only_in_b)

    print(f"only in a: {only_in_a}")
    print(f"only in b: {only_in_b}")

    # show both topologies, highlight the edges to be removed. 
    visualize_topologies(topology_a, topology_b, highlight_edges_a=only_in_a, highlight_edges_b=only_in_b)

    print("continuing3")

    # remove edges that are only present in one of the two topologies
    # propagate changes in case of illegal states (illegal nodes, geographically impossible edges)
    topology_a = clean_topology(topology_a, only_in_a)
    topology_b = clean_topology(topology_b, only_in_b)

    # TODO: What if I would also give the start/end nodes some greateer weight?

    visualize_topologies(topology_a, topology_b)

    print(f"topology_1 nodes {len(topology_a.nodes)}")
    print(f"topology_2 nodes {len(topology_b.nodes)}")
    print(f"topology_1 edges {len(topology_a.edges)}")
    print(f"topology_2 edges {len(topology_b.edges)}")

    return topology_a, topology_b