from __future__ import print_function
import os
import sys
import stl
import meshcut
import numpy as np


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


def load_sphere():
    fname = os.path.join(DATA_DIR, 'sphere.stl')
    m = stl.mesh.Mesh.from_file(fname)

    # Flatten our vert array to Nx3 and generate corresponding faces array
    verts = m.vectors.reshape(-1, 3)
    faces = np.arange(len(verts)).reshape(-1, 3)

    verts, faces = meshcut.merge_close_vertices(verts, faces)
    return meshcut.TriangleMesh(verts, faces)


def load_non_manifold_1():
    fname = os.path.join(DATA_DIR, 'non_manifold_1.stl')
    m = stl.mesh.Mesh.from_file(fname)

    # Flatten our vert array to Nx3 and generate corresponding faces array
    verts = m.vectors.reshape(-1, 3)
    faces = np.arange(len(verts)).reshape(-1, 3)

    verts, faces = meshcut.merge_close_vertices(verts, faces)
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


def test_plane_contain_edge():
    """
    Test with a plane that contains the bottom edge of the mesh
    """
    # Flattened view
    # 1 is connected to 7 and 6 to 0
    #
    #  -- 1 ---- 3 ---- 5 ---- 7 --
    #   / |    / |   /  |   /  |
    #     |  /   | /    | /    |  /
    #  -- 0 ---- 2 ---- 4 ---- 6 --
    #
    # Top view
    #     6 - 4
    #     |   |
    #     0 - 2
    #

    verts = [
        (0, 0, 0),
        (0, 1, 0),
        (1, 0, 0),
        (1, 1, 0),
        (1, 0, 1),
        (1, 1, 1),
        (0, 0, 1),
        (0, 1, 1)
    ]
    faces = [
        (0, 1, 3),
        (0, 3, 2),
        (2, 3, 5),
        (2, 5, 4),
        (4, 5, 7),
        (4, 7, 6),
        (6, 7, 1),
        (6, 1, 0),
    ]

    mesh = meshcut.TriangleMesh(verts, faces)
    plane_orig = verts[2]
    plane_orig = (plane_orig[0], plane_orig[1], plane_orig[2])
    print('plane_orig', plane_orig)
    # Align exactly with the 0 - 2 - 4 line
    plane_norm = (0, 1, 0)
    plane = meshcut.Plane(plane_orig, plane_norm)

    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 1
    # We should have the four bottom points in our slice
    assert len(p[0]) == 4
    assert np.all(p[0][:, 1] == 0)


def test_plane_triangles_common_edge_on_plane():
    """
    Test with a plane that contains the middle edge of triangles
    """
    # Flattened view
    # 1 is connected to 8 and 10
    #
    #     2     5    8
    #   /   \ /  \ /   \
    # -1-----3----6----- (back to 1)
    #   \   / \  / \   /
    #     4     7    10
    #
    # Top view
    #      1 - 3
    #      | /
    #      6
    verts = [
        (0, 0, 0),
        (0, 1, 0),
        (1, 0, 0),
        (1, 1, 0),
        (1, 0, 1),
        (1, 1, 1),
        (0, 0, 1),
        (0, 1, 1)
    ]
    faces = [
        (0, 1, 3),
        (0, 3, 2),
        (2, 3, 5),
        (2, 5, 4),
        (4, 5, 7),
        (4, 7, 6),
        (6, 7, 1),
        (6, 1, 0),
    ]

    mesh = meshcut.TriangleMesh(verts, faces)
    plane_orig = verts[2]
    plane_orig = (plane_orig[0], plane_orig[1], plane_orig[2])
    print('plane_orig', plane_orig)
    # Align exactly with the 0 - 2 - 4 line
    plane_norm = (0, 1, 0)
    plane = meshcut.Plane(plane_orig, plane_norm)

    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 1
    # We should have the four bottom points in our slice
    assert len(p[0]) == 4
    assert np.all(p[0][:, 1] == 0)


def test_sphere():
    mesh = load_sphere()

    plane_orig = (0, 0.7, 0)
    plane_norm = (0, 0.5, 0.5)
    plane = meshcut.Plane(plane_orig, plane_norm)
    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 1
    assert len(p[0]) > 10

    plane_orig = (0, 0.75, 0)
    plane_norm = (0, 1, 0)
    plane = meshcut.Plane(plane_orig, plane_norm)
    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 1
    assert len(p[0]) > 10


def test_non_manifold_join_disjoint_lines():
    """
    Regression test for #20

    On a non-manifold mesh, if the starting trial is not on the boundary,
    we'll end up walking polylines in two directions. But instead of returning
    those as two cuts, we should join them because they originate from the
    same triangle
    """
    mesh = load_non_manifold_1()

    plane_orig = (387.224, 232.544, 6.287)
    plane_norm = (0, 1, 0)
    plane = meshcut.Plane(plane_orig, plane_norm)

    p = meshcut.cross_section_mesh(mesh, plane)
    assert len(p) == 1
    assert len(p[0]) > 10
