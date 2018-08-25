from __future__ import print_function
import os
import sys
import meshcut


EXAMPLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'examples')
DATA_DIR = os.path.join(EXAMPLES_DIR, 'data')


def load_human():
    # TODO: This is ugly, make ply a submodule of meshcut (requires
    # some refactoring)
    print('examples :', EXAMPLES_DIR)
    sys.path.append(EXAMPLES_DIR)
    import ply

    fname = os.path.join(DATA_DIR, 'mesh.ply')
    with open(fname) as f:
        verts, faces, _ = ply.load_ply(f)

    return meshcut.TriangleMesh(verts, faces)


def test_plane_aligned_with_edges():
    mesh = load_human()
    # This will align the plane with some edges, so this is a good test
    # for vertices intersection handling
    plane_orig = (1.28380000591278076172, -0.12510000169277191162, 0)
    plane_norm = (1, 0, 0)
    plane = meshcut.Plane(plane_orig, plane_norm)
    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 3


def test_plane_not_axis_aligned():
    mesh = load_human()
    plane_orig = (0.7, 0, 0)
    plane_norm = (0.2, 0.5, 0.3)
    plane = meshcut.Plane(plane_orig, plane_norm)

    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 2
