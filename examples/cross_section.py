##
import os
import numpy as np
import mayavi.mlab as mlab
import itertools
import utils
import ply

import meshcut

##

if False:
    reload(meshcut)
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (1.28380000591278076172, -0.12510000169277191162, 0)
    plane_norm = (1, 0, 0)

    plane = meshcut.Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=3)

##
if __name__ == '__main__':
    ##
    example_dir = os.path.join(os.path.dirname(meshcut.__file__), 'examples')
    example_fname = os.path.join(example_dir, 'data', 'mesh.ply')
    with open(example_fname) as f:
        verts, faces, _ = ply.load_ply(f)

    mesh = meshcut.TriangleMesh(verts, faces)
    ##

    def show(plane, expected_n_contours):
        P = meshcut.cross_section_mesh(mesh, plane)
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
                utils.points3d(np.array(p), point_size=3, color=(1,1,1))
                mlab.plot3d(p[:, 0], p[:, 1], p[:, 2], tube_radius=None,
                            line_width=3.0, color=color)
        return P

    ##
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (1.28380000591278076172, -0.12510000169277191162, 0)
    plane_norm = (1, 0, 0)

    plane = meshcut.Plane(plane_orig, plane_norm)
    P = show(plane, expected_n_contours=3)
    ##
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (0.93760002, -0.12909999, 0)
    plane_norm = (1, 0, 0)

    plane = meshcut.Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=1)

    ##
    plane_orig = (1, 0, 0)
    plane_norm = (1, 0, 0)

    plane = meshcut.Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=3)
    ##
    plane_orig = (0.7, 0, 0)
    plane_norm = (0.2, 0.5, 0.3)

    plane = meshcut.Plane(plane_orig, plane_norm)
    show(plane, expected_n_contours=2)
    ##
    mlab.show()
    ##
