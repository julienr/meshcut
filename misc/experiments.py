"""Experimental stuff. Run this in ipython from the 'examples' directory"""
##
import numpy as np
import ply
import mayavi.mlab as mlab  # noqa
from utils import points3d, trimesh3d, show_plane
# %matplotlib qt
##
with open('data/mesh.ply') as f:
    verts, faces, _ = ply.load_ply(f)

##
# plane defined by origin and normal
plane_orig = (1.0, 0.0, 0.0)
plane_norm = (1.0, 0.0, 0.0)

trimesh3d(verts, faces, color=(1, 1, 1))
show_plane(plane_orig, plane_norm, scale=0.5, color=(1, 0, 0), opacity=0.5)
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
    trimesh3d(verts, faces_pos, color=(1, 0, 0))
    trimesh3d(verts, faces_mid, color=(0, 1, 0))
    trimesh3d(verts, faces_neg, color=(0, 0, 1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1, 0, 0), opacity=0.5)

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
    for fi in range(3):
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
    trimesh3d(verts, faces_pos, color=(1, 0, 0))
    trimesh3d(verts, faces_mid, color=(0, 1, 0))
    trimesh3d(verts, faces_neg, color=(0, 0, 1))
    show_plane(plane_orig, plane_norm, scale=0.5, color=(1, 0, 0), opacity=0.5)

    points3d(inter_points, point_size=2, color=(1, 1, 1))
##
