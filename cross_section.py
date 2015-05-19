##
import pylab as pl
import mayavi.mlab as mlab
import numpy as np
import numpy.linalg as la
import collections
import itertools
import utils
import ply
##
def make_edge(v1, v2):
    """
    We store edges as tuple where the vertex indices are sorted (so
    the edge going from v1 to v2 and v2 to v1 is the same)
    """
    return tuple(sorted((v1, v2)))


class TriangleMesh(object):
    def __init__(self, verts, tris):
        """
        Args:
            verts: The 3D vertex positions
            tris: A list of triplet containing vertex indices for each triangle
        """
        self.verts = verts
        # For each edge, contains the list of triangles it belongs to
        # If the mesh is closed, each edge belongs to 2 triangles
        self.edges_to_tris = collections.defaultdict(lambda : [])
        # For each triangle, contains the edges it contains
        self.tris_to_edges = {}
        # For each vertex, the list of triangles it belongs to
        self.verts_to_tris = collections.defaultdict(lambda : [])

        self.tris = faces

        # Fill data structures
        for tid, f in enumerate(faces):
            tri_edges = []
            for i in xrange(3):
                v1 = f[i]
                v2 = f[(i + 1) % 3]
                e = make_edge(v1, v2)
                self.edges_to_tris[e].append(tid)
                tri_edges.append(e)
                self.verts_to_tris[f[i]].append(tid)
            self.tris_to_edges[tid] = tri_edges

    def edges_for_triangle(self, tidx):
        """Returns the edges forming triangle with given index"""
        return self.tris_to_edges[tidx]

    def triangles_for_edge(self, edge):
        return self.edges_to_tris[edge]

    def triangles_for_vert(self, vidx):
        """Returns the triangles `vidx` belongs to"""
        return self.verts_to_tris[vidx]


class Plane(object):
    def __init__(self, orig, normal):
        self.orig = orig
        self.n = normal / la.norm(normal)


def point_to_plane_dist(p, plane):
    return np.dot((p - plane.orig), plane.n)


INTERSECT_EDGE = 0
INTERSECT_VERTEX = 1


def compute_edge_plane_intersection(mesh, e, plane):
    """
    Compute the intersection between and edge and a plane
    """
    d1 = point_to_plane_dist(mesh.verts[e[0]], plane)
    d2 = point_to_plane_dist(mesh.verts[e[1]], plane)

    if d1 * d2 < 0:
        # intersection factor (between 0 and 1)
        # here is a nice drawing :
        # https://ravehgonen.files.wordpress.com/2013/02/slide8.png
        s = d1 / (d1 - d2)
        vdir = mesh.verts[e[1]] - mesh.verts[e[0]]
        ip = mesh.verts[e[0]] + vdir * s

        return (INTERSECT_EDGE, ip, e)
    elif np.fabs(d1) < 1e-5:
        # point on plane
        return (INTERSECT_VERTEX, mesh.verts[e[0]], e[0])
    elif np.fabs(d2) < 1e-5:
        # point on plane
        return (INTERSECT_VERTEX, mesh.verts[e[1]], e[1])
    else:
        return None


def cross_section_from(mesh, tid, plane):
    """
    Create a cross-section of the mesh starting at the given triangle

    This returns p = [] if the plane does not intersect the given triangle

    Args:
        mesh: A TriangleMesh instance
        tid: the triangle index at which to start the cross-section polyline
        plane: the cutting plane
    Returns:
        p: a list of consecutive points forming a polyline
        visited_tris: the set of visited triangles
    """
    # We keep track of visited triangles. This is not required by this
    # function but can be useful for callers
    visited_tris = set([tid])

    # Take a pen an draw a plane cutting some triangles
    #
    # To create the cutting polyline, we explore the edges that are cut by
    # the plane. At any given point in time, the next point in our line can
    # to be found in a "to visit" edges set. Amongst the edges to visit,
    # we will exactly one intersection and from this intersection, we'll
    # create a new "to visit" set and continue.
    #
    # The set of edges "to visit" is determined by the type of intersection :
    # - If we found that our plane intersects an edge (cuts it), the "to visit"
    #   set is simply the list of edges of the neighboring triangle to the
    #   cut edge.
    # - If our plane intersects a vertex, the "to visit" set is all the
    #   edges of the triangles that contain the intersected vertex (again,
    #   draw it)
    #
    # To avoid going back on our steps, we maintain a set of explored edges
    #
    # We start by visiting the edges in the 'initial' triangle (tid)

    p = []
    edges_to_visit = set(mesh.edges_for_triangle(tid))
    visited_edges = set()

    while len(edges_to_visit) > 0:
        e = edges_to_visit.pop()
        visited_edges.add(e)

        intersection = compute_edge_plane_intersection(mesh, e, plane)
        if intersection is None:
            continue
        elif intersection[0] == INTERSECT_EDGE:
            p.append(intersection[1])

            next_tids = mesh.triangles_for_edge(e)

        elif intersection[0] == INTERSECT_VERTEX:
            p.append(intersection[1])

            next_tids = mesh.triangles_for_vert(intersection[2])

        # We will not find any more intersections in current edges
        visited_edges.update(edges_to_visit)

        # Update the list of edges to visit from the triangles in next_tids,
        # discarding already visited edges
        edges_to_visit = []
        for tid in mesh.triangles_for_edge(e):
            visited_tris.add(tid)
            for new_edge in mesh.edges_for_triangle(tid):
                if new_edge not in visited_edges:
                    edges_to_visit.append(new_edge)
    return np.array(p), visited_tris


def cross_section(mesh, plane):
    """
    Compute a cross-section of the mesh using the provided cut-plane.

    Returns:
        P: a list of polylines (if the mesh is non-convex, the cross-section
           can consist of multiple disjoints polylines)
    """

    # the set of all triangles
    T = set(range(len(mesh.tris)))

    # the list of all cross-section polylines
    P = []

    while len(T) > 0:
        tid = T.pop()
        p, visited_tris = cross_section_from(mesh, tid, plane)
        if len(p) > 0:
            P.append(p)
        T.difference_update(visited_tris)

    return P


##
if __name__ == '__main__':
    ##
    with open('mesh.ply') as f:
        verts, faces, _ = ply.load_ply(f)

    mesh = TriangleMesh(verts, faces)
    ##
    def show(plane):
        P = cross_section(mesh, plane)
        colors = [
            (0, 1, 1),
            (1, 0, 1),
            (0, 0, 1)
        ]
        print len(P)

        if True:
            utils.trimesh3d(mesh.verts, mesh.tris, color=(1,1,1))
            utils.show_plane(plane.orig, plane.n, scale=1, color=(1,0,0),
                             opacity=0.5)

            for p, color in zip(P, itertools.cycle(colors)):
                p = np.array(p)
                #points3d(np.array(p), point_size=2, color=color)
                mlab.plot3d(p[:,0], p[:,1], p[:,2], tube_radius=None,
                            line_width=3.0, color=color)

    ##
    plane_orig = (1, 0, 0)
    plane_norm = (1, 0, 0)

    plane = Plane(plane_orig, plane_norm)
    show(plane)
    ##
    plane_orig = (0.7, 0, 0)
    plane_norm = (0.2, 0.5, 0.3)

    plane = Plane(plane_orig, plane_norm)
    show(plane)
    ##
    mlab.show()
    ##
