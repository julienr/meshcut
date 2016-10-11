"""
This demonstrates cutting on a .stl mesh, which requires some pre-processing
because STL stores 3 vertices for each face (so vertices are duplicated).

The merge_close_vertices is a naive implementation of this pre-processing
that will probably take a long time on large meshes. libigl provides
a remove_duplicates function that should be better :
  https://github.com/libigl/libigl/blob/master/include/igl/remove_duplicates.h

Thanks to @superzanti for providing the model
"""
##
import os
import numpy as np
import numpy.linalg as la
import mayavi.mlab as mlab
import itertools
import utils
import ply
import scipy.spatial.distance as spdist

import meshcut
##
def merge_close_vertices(verts, faces, close_epsilon=1e-5):
    """
    Will merge vertices that are closer than close_epsilon.

    Warning, this has a O(n^2) memory usage because we compute the full
    vert-to-vert distance matrix. If you have a large mesh, might want
    to use some kind of spatial search structure like an octree or some fancy
    hashing scheme

    Returns: new_verts, new_faces
    """
    # Pairwise distance between verts
    D = spdist.cdist(verts, verts)

    # Compute a mapping from old to new : for each input vert, store the index
    # of the new vert it will be merged into
    close_epsilon = 1e-5
    old2new = np.zeros(D.shape[0], dtype=np.int)
    # A mask indicating if a vertex has already been merged into another
    merged_verts = np.zeros(D.shape[0], dtype=np.bool)
    new_verts = []
    for i in range(D.shape[0]):
        if merged_verts[i]:
            continue
        else:
            # The vertices that will be merged into this one
            merged = np.flatnonzero(D[i, :] < close_epsilon)
            old2new[merged] = len(new_verts)
            new_verts.append(verts[i])
            merged_verts[merged] = True

    new_verts = np.array(new_verts)

    # Recompute face indices to index in new_verts
    new_faces = np.zeros((len(faces), 3), dtype=np.int)
    for i, f in enumerate(faces):
        new_faces[i] = (old2new[f[0]], old2new[f[1]], old2new[f[2]])

    # again, plot with utils.trimesh3d(new_verts, new_faces)
    return new_verts, new_faces


def load_stl(stl_fname):
    import stl
    m = stl.mesh.Mesh.from_file(stl_fname)

    # Flatten our vert array to Nx3 and generate corresponding faces array
    verts = m.vectors.reshape(-1, 3)
    faces = np.arange(len(verts)).reshape(-1, 3)

    verts, faces = merge_close_vertices(verts, faces)
    return verts, faces

##
if __name__ == '__main__':
    ##
    example_dir = os.path.join(os.path.dirname(meshcut.__file__), 'examples')
    example_fname = os.path.join(example_dir, 'data', 'sphere.stl')
    verts, faces = load_stl(example_fname)

    mesh = meshcut.TriangleMesh(verts, faces)
    ##
    def show(plane):
        P = meshcut.cross_section_mesh(mesh, plane)
        colors = [
            (0, 1, 1),
            (1, 0, 1),
            (0, 0, 1)
        ]
        print "num contours : ", len(P)

        if True:
            utils.trimesh3d(mesh.verts, mesh.tris, color=(1, 1, 1),
                            opacity=0.5)
            utils.show_plane(plane.orig, plane.n, scale=1, color=(1, 0, 0),
                             opacity=0.5)

            for p, color in zip(P, itertools.cycle(colors)):
                p = np.array(p)
                #utils.points3d(np.array(p), point_size=3, color=(1,1,1))
                mlab.plot3d(p[:, 0], p[:, 1], p[:, 2], tube_radius=None,
                            line_width=3.0, color=color)
        return P

    ##
    plane_orig = (0, 0.5, 0)
    plane_norm = (0, 1, 0)
    plane = meshcut.Plane(plane_orig, plane_norm)
    P = show(plane)
    ##
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (0.6, -0.12510000169277191162, 0)
    plane_norm = (0.3, 0.7, -0.2)
    plane_norm /= la.norm(plane_norm)

    plane = meshcut.Plane(plane_orig, plane_norm)
    P = show(plane)
    ##
    mlab.show()
    ##
