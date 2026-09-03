from ...model import Topology, Node, Edge, Wgs84GeoNode, GeoNode

import matplotlib.pyplot as plt
import math


def get_edge_polyline(edge: Edge) -> list[tuple[float, float]]:
    coordinates = [(geo_node.x, geo_node.y) for geo_node in [edge.node_a.geo_node] + edge.intermediate_geo_nodes + [edge.node_b.geo_node]]
    return coordinates

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

