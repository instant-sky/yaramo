from ...model import Topology, Edge, Node, GeoNode, Wgs84GeoNode, EdgeConnectionDirection

import math
import networkx as nx
import matplotlib.pyplot as plt
import pyproj


def get_edge_polyline(edge: Edge) -> list[tuple[float, float]]:
    coordinates = [(geo_node.x, geo_node.y) for geo_node in [edge.node_a.geo_node] + edge.intermediate_geo_nodes + [edge.node_b.geo_node]]
    return coordinates

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

def edge_fully_contained(
    polyline: list[tuple[float, float]],
    reference_polylines: list[list[tuple[float, float]]],
    buffer: float = 10,
) -> bool:
    for px, py in polyline:
        min_dist = float("inf")
        for reference_polyline in reference_polylines:
            d = point_to_polyline_min_distance(px, py, reference_polyline)
            if d < min_dist:
                min_dist = d
        if min_dist > buffer:
            return False
    return True


def _get_utm_transformer(
    polylines: list[list[tuple[float, float]]],
) -> pyproj.Transformer:
    all_coords = [coord for polyline in polylines for coord in polyline]
    mean_lon = sum(x for x, y in all_coords) / len(all_coords)
    mean_lat = sum(y for x, y in all_coords) / len(all_coords)
    utm_zone = int((mean_lon + 180) / 6) + 1
    epsg = 32600 + utm_zone if mean_lat >= 0 else 32700 + utm_zone
    return pyproj.Transformer.from_crs(
        "EPSG:4326", f"EPSG:{epsg}", always_xy=True
    )


def _project_polyline(
    polyline: list[tuple[float, float]],
    transformer: pyproj.Transformer,
) -> list[tuple[float, float]]:
    return [transformer.transform(x, y) for x, y in polyline]


def _short_id(uuid: str, full: bool) -> str:
    return uuid if full else uuid[-4:]


def _to_plot_coordinates(geo_node: GeoNode) -> tuple[float, float]:
    """Return the (horizontal, vertical) plot coordinates of a geo node.

    ``Wgs84GeoNode`` stores latitude in ``x`` and longitude in ``y`` (the order
    expected by the haversine library), so the values are swapped here to get the
    usual map orientation with longitude east-west on the x-axis.
    """
    if isinstance(geo_node, Wgs84GeoNode):
        return (float(geo_node.y), float(geo_node.x))
    return (float(geo_node.x), float(geo_node.y))


def visualize_topology(
    topology: Topology,
    title: str | None = None,
    highlight_nodes: set[Node] | None = None,
    highlight_edges: set[Edge] | None = None,
    label_nodes: bool = True,
    label_edges: bool = False,
    full_labels: bool = False,
    show: bool = True,
    block: bool = True,
    save_path: str | None = None,
    ax=None,
):
    """Draw a Topology using the geo positions of its Nodes and the polylines of
    its Edges (node_a -> intermediate_geo_nodes -> node_b).

    Nodes/edges in ``highlight_nodes``/``highlight_edges`` are drawn in red. They
    are also drawn if they were already removed from the topology (but still carry
    their geo information), which is useful right after a node got merged away.

    ``show=False`` can be used to draw into a provided ``ax`` without blocking
    (e.g. to plot two topologies next to each other).
    """
    highlight_nodes = highlight_nodes or set()
    highlight_edges = highlight_edges or set()

    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 8))

    # positions for all current nodes, plus highlighted nodes that are no longer
    # part of the topology (e.g. a node that was just removed)
    positions: dict[Node, tuple[float, float]] = {}
    for node in topology.nodes.values():
        if node.geo_node is not None:
            positions[node] = _to_plot_coordinates(node.geo_node)
    for node in highlight_nodes:
        if node not in positions and node.geo_node is not None:
            positions[node] = _to_plot_coordinates(node.geo_node)

    # edges as polylines, plus highlighted edges no longer in the topology
    edges = list(topology.edges.values())
    for edge in highlight_edges:
        if edge not in edges:
            edges.append(edge)
    for edge in edges:
        coords = get_edge_polyline(edge)
        if isinstance(edge.node_a.geo_node, Wgs84GeoNode):
            coords = [(c[1], c[0]) for c in coords]
        xs, ys = [c[0] for c in coords], [c[1] for c in coords]
        highlighted = edge in highlight_edges
        ax.plot(
            xs,
            ys,
            color="tab:red" if highlighted else "0.6",
            lw=3.0 if highlighted else 1.2,
            solid_capstyle="round",
            zorder=2,
        )
        for geo_node in edge.intermediate_geo_nodes:
            px, py = _to_plot_coordinates(geo_node)
            ax.scatter([px], [py], color="0.3", s=10, zorder=2.5)
        if label_edges or highlighted:
            mid = coords[len(coords) // 2]
            ax.annotate(
                _short_id(edge.uuid, full_labels),
                (mid[0], mid[1]),
                fontsize=7,
                color="tab:red" if highlighted else "0.4",
                zorder=4,
            )

    # nodes at their geo positions
    for node, (x, y) in positions.items():
        highlighted = node in highlight_nodes
        ax.scatter(
            [x],
            [y],
            color="tab:red" if highlighted else "tab:blue",
            s=90 if highlighted else 40,
            marker="^" if highlighted else "o",
            zorder=3,
        )
        if label_nodes:
            ax.annotate(
                _short_id(node.uuid, full_labels),
                (x, y),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
                zorder=4,
            )

    ax.set_title(title or f"{len(topology.nodes)} nodes, {len(topology.edges)} edges")

    # WGS84 uses (longitude, latitude): stretch x by 1/cos(lat) so the plot is not distorted
    wgs84_lats = [pos[1] for node, pos in positions.items() if isinstance(node.geo_node, Wgs84GeoNode)]
    if wgs84_lats:
        mean_lat = sum(wgs84_lats) / len(wgs84_lats)
        ax.set_aspect(1.0 / math.cos(math.radians(mean_lat)))
    else:
        ax.set_aspect("equal")

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show(block=block)
    return ax


def visualize_topologies(
    topology_a: Topology,
    topology_b: Topology,
    title_a: str = "Topology A",
    title_b: str = "Topology B",
    highlight_nodes_a: set[Node] | None = None,
    highlight_nodes_b: set[Node] | None = None,
    highlight_edges_a: set[Edge] | None = None,
    highlight_edges_b: set[Edge] | None = None,
    show: bool = True,
    block: bool = True,
    save_path: str | None = None,
):
    """Draw two Topologies side by side on shared axes for debugging comparisons."""
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(18, 8), sharex=True, sharey=True)
    visualize_topology(
        topology_a,
        title=title_a,
        highlight_nodes=highlight_nodes_a,
        highlight_edges=highlight_edges_a,
        show=False,
        ax=ax_a,
    )
    visualize_topology(
        topology_b,
        title=title_b,
        highlight_nodes=highlight_nodes_b,
        highlight_edges=highlight_edges_b,
        show=False,
        ax=ax_b,
    )
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    if show:
        plt.show(block=block)
    return fig, (ax_a, ax_b)


def clean_topology(topology: Topology, edges_to_remove: list[Edge]):
    for edge in edges_to_remove:
        topology.edges.pop(edge.uuid)
    
    edges_to_remove_per_node_a: dict[Node, list[Edge]] = {}
    for edge in edges_to_remove:
        edges_to_remove_per_node_a.setdefault(edge.node_a, []).append(edge)
        edges_to_remove_per_node_a.setdefault(edge.node_b, []).append(edge)
        
    for node, edges in edges_to_remove_per_node_a.items():
        for edge in edges:
            node.remove_edge(edge)
        if len(node.connected_edges) == 3:
            raise ValueError("no edges removed but node was connected to removed edge.")
        if len(node.connected_edges) == 2: # illegal node state -> node needs to be removed, adjacent edges form union
            print(f"removing node and checking edge")
            topology.nodes.pop(node.uuid)
            
            edge_a: Edge = node.connected_edges[0]
            edge_b: Edge = node.connected_edges[1]
            if edge_a.uuid in topology.edges:
                topology.edges.pop(edge_a.uuid)
            if edge_b.uuid in topology.edges:
                topology.edges.pop(edge_b.uuid)

            node_a: Node = edge_a.get_opposite_node(node)
            node_b: Node = edge_b.get_opposite_node(node)

            new_edge: Edge = Edge(node_a, node_b, intermediate_geo_nodes=node.connected_edges[0].intermediate_geo_nodes + [node.geo_node] + node.connected_edges[1].intermediate_geo_nodes)
            topology.add_edge(new_edge)

            node_a.replace_edge(edge_a, new_edge)
            node_b.replace_edge(edge_b, new_edge)
            visualize_topology(
                topology,
                title=f"merged edges at node {_short_id(node.uuid, False)}",
                highlight_nodes={node},
                highlight_edges={edge_a, edge_b},
            )
        if len(node.connected_edges) == 1: # switch -> end, legal
            print("legal node remaining")
            continue
        if len(node.connected_edges) == 0: # all edges removed, just remove node
            print(f"removing node {node}")
            topology.nodes.pop(node.uuid)
            visualize_topology(topology, title=f"removed node {_short_id(node.uuid, False)}", highlight_nodes={node})
    return topology


def edge_shape_comparison(
    topology_a: Topology,
    topology_b: Topology,
):
    raw_polylines_a: dict[Edge, list[tuple[float, float]]] = {
        edge: get_edge_polyline(edge) for edge in topology_a.edges.values()
    }

    raw_polylines_b = {
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

    only_in_a: dict[Edge, tuple[float, float]] = {}
    only_in_b: dict[Edge, tuple[float, float]] = {}


    for edge, polyline in polylines_a.items():
        if not edge_fully_contained(polyline, polylines_b.values()):
            only_in_a[edge] = raw_polylines_a[edge]

    for edge, polyline in polylines_b.items():
        if not edge_fully_contained(polyline, polylines_a.values()):
            only_in_b[edge] = raw_polylines_b[edge]

    print(f"only in a: {only_in_a}")
    print(f"only in b: {only_in_b}")

    visualize_topologies(topology_a, topology_b)

    print("continuing3")

    topology_a = clean_topology(topology_a, only_in_a)
    topology_b = clean_topology(topology_b, only_in_b)

    visualize_topologies(topology_a, topology_b)

    topology_a.update_edge_lengths()
    topology_b.update_edge_lengths()

    return topology_a, topology_b

    # remove edges from topology
    # for edge in only_in_a:
    #     topology_a.edges.pop(edge.uuid)
    
    # edges_to_remove_per_node_a: dict[Node, list[Edge]] = {}
    # for edge in only_in_a:
    #     edges_to_remove_per_node_a.setdefault(edge.node_a, []).append(edge)
    #     edges_to_remove_per_node_a.setdefault(edge.node_b, []).append(edge)
        
    # for node, edges in edges_to_remove_per_node_a.items():
    #     for edge in edges:
    #         node.remove_edge(edge)
    #     if len(node.connected_edges) == 3:
    #         raise ValueError("no edges removed but node was connected to removed edge.")
    #     if len(node.connected_edges) == 2: # illegal node state -> node needs to be removed, adjacent edges form union
    #         topology_a.nodes.pop(node.uuid)
            
    #         edge_a: Edge = node.connected_edges[0]
    #         edge_b: Edge = node.connected_edges[1]
    #         if edge_a.uuid in topology_a.edges:
    #             topology_a.edges.pop(edge_a.uuid)
    #         if edge_b.uuid in topology_a.edges:
    #             topology_a.edges.pop(edge_b.uuid)

    #         node_a: Node = edge_a.get_opposite_node(node)
    #         node_b: Node = edge_b.get_opposite_node(node)

    #         new_edge: Edge = Edge(node_a, node_b, intermediate_geo_nodes=node.connected_edges[0].intermediate_geo_nodes + [node.geo_node] + node.connected_edges[1].intermediate_geo_nodes)
    #         topology_a.add_edge(new_edge)

    #         node_a.replace_edge(edge_a, new_edge)
    #         node_b.replace_edge(edge_b, new_edge)
    #         visualize_topology(
    #             topology_a,
    #             title=f"merged edges at node {_short_id(node.uuid, False)}",
    #             highlight_nodes={node},
    #             highlight_edges={edge_a, edge_b},
    #         )
    #     if len(node.connected_edges) == 1: # switch -> end, legal
    #         continue
    #     if len(node.connected_edges) == 0: # all edges removed, just remove node
    #         topology_a.nodes.pop(node.uuid)

    visualize_topology(
        topology=topology_a,
        title="topology_a after removal",
    )
    print(f"nodes: {len(topology_a.nodes)}")
    print(f"edges: {len(topology_a.edges)}")

    # remove edges with corresponding nodes (using combine/split methods by Arne?)
    # to create isomorph graphs

    # use isomorphic matching