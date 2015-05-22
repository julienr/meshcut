##
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
        self.edges_to_tris = collections.defaultdict(lambda: [])
        # For each triangle, contains the edges it contains
        self.tris_to_edges = {}
        # For each vertex, the list of triangles it belongs to
        self.verts_to_tris = collections.defaultdict(lambda: [])

        self.tris = tris

        # Fill data structures
        for tid, f in enumerate(tris):
            tri_edges = []
            for i in xrange(3):
                v1 = f[i]
                v2 = f[(i + 1) % 3]
                e = make_edge(v1, v2)
                self.edges_to_tris[e].append(tid)
                tri_edges.append(e)
                self.verts_to_tris[f[i]].append(tid)
            self.tris_to_edges[tid] = tri_edges

        # Sanity check : max 2 faces per edge
        for e, tris in self.edges_to_tris.items():
            assert len(tris) <= 2

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


def triangle_intersect_plane(mesh, tid, plane):
    dists = [point_to_plane_dist(mesh.verts[vid], plane)
             for vid in mesh.tris[tid]]
    side = np.sign(dists)
    return not (side[0] == side[1] == side[2])


INTERSECT_EDGE = 0
INTERSECT_VERTEX = 1


def compute_triangle_plane_intersections(mesh, tid, plane):
    """
    Compute the intersection between a triangle and a plane

    Returns a list of intersections in the form
        (INTERSECT_EDGE, <intersection point>, <edge>) for edges intersection
        (INTERSECT_VERTEX, <intersection point>, <vertex index) for vertices

    This return between 0 and 2 intersections :
    - 0 : the plane does not intersect the plane
    - 1 : one of the triangle's vertices lies on the plane (so it just
          "touches" the plane without really intersecting)
    - 2 : the plane slice the triangle in two parts (either vertex-edge,
          vertex-vertex or edge-edge)
    """
    # TODO: Use a distance cache
    dists = {vid:point_to_plane_dist(mesh.verts[vid], plane)
             for vid in mesh.tris[tid]}

    # Iterate through the edges, cutting the ones that intersect
    intersections = []
    for e in mesh.edges_for_triangle(tid):
        v1 = mesh.verts[e[0]]
        d1 = dists[e[0]]
        v2 = mesh.verts[e[1]]
        d2 = dists[e[1]]

        if d1 * d2 < 0:
            # intersection factor (between 0 and 1)
            # here is a nice drawing :
            # https://ravehgonen.files.wordpress.com/2013/02/slide8.png
            s = d1 / (d1 - d2)
            vdir = v2 - v1
            ipos = v1 + vdir * s
            intersections.append((INTERSECT_EDGE, ipos, e))
        elif np.fabs(d1) < 1e-10:
            # point on plane
            intersections.append((INTERSECT_VERTEX, v1, e[0]))

    assert len(intersections) <= 2
    return intersections


def get_next_triangle(mesh, from_tid, plane, intersection):
    """
    Returns the next triangle to visit given the intersection and
    the triangle we're coming from
    """
    if intersection[0] == INTERSECT_EDGE:
        tris = mesh.triangles_for_edge(intersection[2])
        for tid in tris:
            if tid != from_tid:
                return tid
        return None
    elif intersections[0] == INTERSECT_VERTEX:
        tris = mesh.triangles_for_vert(intersection[2])
        for tid in tris:
            if tid != from_tid and triangle_intersect_plane(mesh, tid, plane):
                return tid
        return None


def cross_section(mesh, plane, dist_tol=1e-8):
    """
    Args:
        dist_tol:
    """
    # Set of all triangles
    T = set(range(len(mesh.tris)))
    # List of all cross-section polylines
    P = []

    while len(T) > 0:
        tid = T.pop()
        if triangle_intersect_plane(mesh, tid, plane):
            # We found a starting triangle for a new polyline
            p = []
            intersections = compute_triangle_plane_intersections(
                mesh, tid, plane)
            assert len(intersections) == 2

            # We can start in either direction (intersections[0] or [1]), this
            # is arbitrary for the first triangle
            p.append(intersections[0][1])
            tid = get_next_triangle(mesh, tid, plane, intersections[0])

            # Loop until we have explored all the triangles for the current
            # polyline
            while tid in T:
                T.remove(tid)
                intersections = compute_triangle_plane_intersections(
                        mesh, tid, plane)
                assert len(intersections) == 2
                # Of the two returned intersections, one should have the
                # intersection point equal to p[-1]
                if la.norm(intersections[0][1] - p[-1]) < 1e-8:
                    intersect = intersections[1]
                else:
                    assert la.norm(intersections[1][1] - p[-1]) < 1e-8
                    intersect = intersections[0]

                p.append(intersect[1])
                tid = get_next_triangle(mesh, tid, plane, intersect)

                if tid is None:
                    print 'Degenerate case (probably non-closed mesh)'
                    break

            P.append(p)

    return P

if False:
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (1.28380000591278076172, -0.12510000169277191162, 0)
    plane_norm = (1, 0, 0)

    plane = Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=3)

##
if __name__ == '__main__':
    ##
    with open('mesh.ply') as f:
        verts, faces, _ = ply.load_ply(f)

    mesh = TriangleMesh(verts, faces)
    ##

    def show(plane, expected_n_contours):
        P = cross_section(mesh, plane)
        colors = [
            (0, 1, 1),
            (1, 0, 1),
            (0, 0, 1)
        ]
        print "num contours : ", len(P), ' expected : ', expected_n_contours

        if True:
            utils.trimesh3d(mesh.verts, mesh.tris, color=(1, 1, 1),
                            opacity=0.5)
            utils.show_plane(plane.orig, plane.n, scale=1, color=(1, 0, 0),
                             opacity=0.5)

            for p, color in zip(P, itertools.cycle(colors)):
                p = np.array(p)
                # points3d(np.array(p), point_size=2, color=color)
                mlab.plot3d(p[:, 0], p[:, 1], p[:, 2], tube_radius=None,
                            line_width=3.0, color=color)


    ##
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (1.28380000591278076172, -0.12510000169277191162, 0)
    plane_norm = (1, 0, 0)

    plane = Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=3)
    ##
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (0.93760002, -0.12909999, 0)
    plane_norm = (1, 0, 0)

    plane = Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=1)

    ##
    plane_orig = (1, 0, 0)
    plane_norm = (1, 0, 0)

    plane = Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=3)
    ##
    plane_orig = (0.7, 0, 0)
    plane_norm = (0.2, 0.5, 0.3)

    plane = Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=2)
    ##
    mlab.show()
    ##
