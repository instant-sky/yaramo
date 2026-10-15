from yaramo.model import Topology

def edge_ends_comparison(
        topology_a: Topology,
        topology_b: Topology,
        buffer: int,
):
   """
   This approach is based, but only considers the nodes connected to edges.
   So instead of comparing overlap of edges as is the case for the edge_shape comparison, 
   this approach compares the positions of the start/end nodes of edges.
   """ 

   pass



    # Things to consider:
    # Funktioniert der Ansatz überhaupt? Wenn eine Kante unterbrochen wird, zum Beispiel durch mehrere Überleitstellen, 
    # dann gibt es keine Anhaltspunkte mehr für diese Kante