##
import pylab as pl
import mayavi.mlab as mlab
import numpy as np
import numpy.linalg as la
import collections
import itertools
import ply
%matplotlib qt
##
with open('mesh.ply') as f:
    verts, faces, _ = ply.load_ply(f)

##
# plane defined by origin and normal
plane_orig = (1, 0, 0)
plane_norm = (1, 0, 0)

def points3d(verts, point_size=3, **kwargs):
    if 'mode' not in kwargs:
        kwargs['mode'] = 'point'
    p = mlab.points3d(verts[:,0], verts[:,1], verts[:,2], **kwargs)
    p.actor.property.point_size = point_size

def trimesh3d(verts, faces, **kwargs):
    mlab.triangular_mesh(verts[:,0], verts[:,1], verts[:,2], faces,
            **kwargs)

def orthogonal_vector(v):
    """Return an arbitrary vector that is orthogonal to v"""
    if v[1] != 0 or v[2] != 0:
        c = (1, 0, 0)
    else:
        c = (0, 1, 0)
    return np.cross(v, c)

def show_plane(orig, n, scale=1.0, **kwargs):
    """
    Show the plane with the given origin and normal. scale give its size
    """
    b1 = orthogonal_vector(n)
    b1 /= la.norm(b1)
    b2 = np.cross(b1, n)
    b2 /= la.norm(b2)
    verts = [orig + scale*(-b1 - b2),
             orig + scale*(b1 - b2),
             orig + scale*(b1 + b2),
             orig + scale*(-b1 + b2)]
    faces = [(0, 1, 2), (0, 2, 3)]
    trimesh3d(np.array(verts), faces, **kwargs)

trimesh3d(verts, faces, color=(1,1,1))
show_plane(plane_orig, plane_norm, scale=0.5, color=(1,0,0), opacity=0.5)
##
def point_to_plane_dist(p, plane_orig, plane_norm):
    return np.dot((p - plane_orig), plane_norm)

def classify_faces(verts, faces, plane_orig, plane_norm):
    faces_pos = []
    faces_mid = []
    faces_neg = []

    for f in faces:
        sides = [point_to_plane_dist(p, plane_orig, plane_norm)
                 for p in verts[f]]
        sides = map(np.sign, sides)
        if len(np.unique(sides)) > 1:
            faces_mid.append(f)
        else:
            if sides[0] < 0:
                faces_neg.append(f)
            else:
                faces_pos.append(f)
    return faces_pos, faces_mid, faces_neg

faces_pos, faces_mid, faces_neg = classify_faces(verts, faces, plane_orig,
                                                 plane_norm)

if True:
    trimesh3d(verts, faces_pos, color=(1,0,0))
    trimesh3d(verts, faces_mid, color=(0,1,0))
    trimesh3d(verts, faces_neg, color=(0,0,1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1,0,0), opacity=0.5)

##
def slice_triangle_plane(verts, tri, plane_orig, plane_norm):
    """
    Args:
        verts : the vertices of the mesh
        tri: the face to cut
        plane_orig: origin of the plane
        plane_norm: normal to the plane
    """
    dists = [point_to_plane_dist(p, plane_orig, plane_norm)
             for p in verts[tri]]

    if np.sign(dists[0]) == np.sign(dists[1]) \
            and np.sign(dists[1]) == np.sign(dists[2]):
        # Triangle is on one side of the plane
        return []

    # Iterate through the edges, cutting the ones that intersect
    intersect_points = []
    for fi in xrange(3):
        v1 = verts[tri[fi]]
        d1 = dists[fi]
        v2 = verts[tri[(fi + 1) % 3]]
        d2 = dists[(fi + 1) % 3]

        if d1 * d2 < 0:
            # intersection factor (between 0 and 1)
            # here is a nice drawing :
            # https://ravehgonen.files.wordpress.com/2013/02/slide8.png
            s = d1 / (d1 - d2)
            vdir = v2 - v1
            intersect_points.append(v1 + vdir * s)
        elif np.fabs(d1) < 1e-5:
            # point on plane
            intersect_points.append(v1)

    return intersect_points

inter_points = []
for tri in faces:
    inter_points += slice_triangle_plane(verts, tri, plane_orig, plane_norm)
inter_points = np.array(inter_points)

if True:
    trimesh3d(verts, faces_pos, color=(1,0,0))
    trimesh3d(verts, faces_mid, color=(0,1,0))
    trimesh3d(verts, faces_neg, color=(0,0,1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1,0,0), opacity=0.5)

    points3d(inter_points, point_size=2, color=(1,1,1))
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


#p = cross_section_from(mesh, 2, T, plane_orig, plane_norm)
plane = Plane(plane_orig, plane_norm)
P = cross_section(mesh, plane)

##
def project_on_plane(p, plane):
    v = (p - plane.orig)
    dist = np.dot(v, plane.n)
    return p - dist * plane.n
    #d = point_to_plane_dist(p, plane)
    #return (p - plane.orig) - plane.n * d) - plane.orig

p2d = np.array([project_point_on_plane(p, plane) for p in P[0]])
##
p = np.array(p)
pl.plot(p[:,1], p[:,2], '*-')
##
colors = [
    (0, 1, 1),
    (1, 0, 1),
    (0, 0, 1)
]
print len(P)

if True:
    #trimesh3d(verts, faces_pos, color=(1,0,0))
    #trimesh3d(verts, faces_mid, color=(0,1,0))
    #trimesh3d(verts, faces_neg, color=(0,0,1))
    trimesh3d(verts, faces, color=(1,1,1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1,0,0), opacity=0.5)

    for p, color in zip(P, itertools.cycle(colors)):
        p = np.array(p)
        #points3d(np.array(p), point_size=2, color=color)
        mlab.plot3d(p[:,0], p[:,1], p[:,2], tube_radius=None,
                    line_width=3.0, color=color)
##
if True:
    #trimesh3d(verts, faces_pos, color=(1,0,0))
    #trimesh3d(verts, faces_mid, color=(0,1,0))
    #trimesh3d(verts, faces_neg, color=(0,0,1))
    trimesh3d(verts, faces, color=(1,1,1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1,0,0), opacity=0.5)

    mlab.plot3d(p[:,0], p[:,1], p[:,2], tube_radius=None,
                line_width=3.0, color=(1,0,0))
    #points3d(np.array(p), point_size=2, color=color)

##


def slice_triangle(mesh, tidx, plane_orig, plane_norm, prev_edge=None):
    """
    Slice the triangle with index `tidx` using the given plane

    Args:
        prev_edge: The edge we came from

    Returns
        intersect_points: 2 points defining the line segment that is
            the intersection between the plane and this triangle
        next_tids: A list of the next triangles to visit (This will only
            contain more than one element if we intersected a vertex).
            Note that this can contain the current tidx
    """
    # distances keyed by vertex id
    dists = {}
    for vid in mesh.tris[tidx]:
        dists[vid] = dist_to_point(mesh.verts[vid], plane_orig, plane_norm)

    if on_same_side(*dists.values()):
        # Plane does not cut triangle
        return [], []

    edges = mesh.edges_for_triangle(tidx)
    intersect_points = []

    next_tids = []
    for e in edges:
        d1 = dists[e[0]]
        d2 = dists[e[1]]
        if d1 * d2 < 0:
            # intersection factor (between 0 and 1)
            # here is a nice drawing :
            # https://ravehgonen.files.wordpress.com/2013/02/slide8.png
            s = d1 / (d1 - d2)
            vdir = mesh.verts[e[1]] - mesh.verts[e[0]]
            intersect_points.append(mesh.verts[e[0]] + vdir * s)

            # Determine if we should advance on the next neighboring triangle
            # on this edge
            next_tids += mesh.triangles_for_edge(e)

        elif np.fabs(d1) < 1e-5:
            # point on plane
            intersect_points.append(mesh.verts[e[0]])

            next_tids += mesh.triangles_for_vert(e[0])

    return intersect_points, set(next_tids)

#
# http://stackoverflow.com/questions/2797431/generate-2d-cross-section-polygon-from-3d-mesh

mesh = TriangleMesh(verts, faces)

iters = 0

# the set of all triangles
T = set(range(len(mesh.tris)))
# output points
P = []
while len(T) > 0:
    tid = T.pop()
    l, next_tids = slice_triangle(mesh, tid, plane_orig, plane_norm)
    # If this triangle intersects the plane, we have a starting triangle
    if len(l) > 0:
        print 'starting new polyline at ', tid
        s = l[0]
        p = list(l)

        while not np.allclose(p[-1], s, 1e-5):
            assert len(next_tids) > 0, "Non-closed mesh ? "

            for tid in next_tids:
                iters += 1
                if tid not in T:
                    # Do not re-explore already explored points
                    continue

                #print tid
                T.remove(tid)
                l, next_tids = slice_triangle(mesh, tid, plane_orig,
                                              plane_norm)
                if len(l) > 0:
                    break
            else:
                # This is a non-closed mesh
                print 'Couldn\'t find how to continue'
                break

            if la.norm(p[-1] - l[0]) < la.norm(p[-1] - l[1]):
                p.append(l[0])
            else:
                p.append(l[1])
        P.append(p)
##
p = np.array(P[1])
pl.plot(p[:,1], p[:,2], '*-')
##
colors = [
    (0, 1, 1),
    (1, 0, 1),
    (0, 0, 1)
]
print len(P)

if True:
    #trimesh3d(verts, faces_pos, color=(1,0,0))
    #trimesh3d(verts, faces_mid, color=(0,1,0))
    #trimesh3d(verts, faces_neg, color=(0,0,1))
    trimesh3d(verts, faces, color=(1,1,1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1,0,0), opacity=0.5)

    for p, color in zip(P, itertools.cycle(colors)):
        p = np.array(p)
        #mlab.plot3d(p[:,0], p[:,1], p[:,2], representation='wireframe')
        points3d(np.array(p), point_size=2, color=color)

##
