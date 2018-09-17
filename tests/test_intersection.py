from __future__ import print_function
import meshcut
import numpy as np


def test_triangle_plane_intersection_edge_edge():
    verts = [
        (0, 0, 0),
        (0, 1, 0),
        (1, 0, 0),
    ]
    faces = [
        (0, 1, 2),
    ]

    mesh = meshcut.TriangleMesh(verts, faces)
    plane = meshcut.Plane(
        orig=(0, 0.5, 0),
        normal=(0, 1, 0)
    )

    intersections = meshcut.compute_triangle_plane_intersections(
        mesh, 0, plane)

    assert len(intersections) == 2
    assert intersections[0][0] == meshcut.INTERSECT_EDGE
    assert np.allclose(intersections[0][1][1], 0.5)
    assert intersections[1][0] == meshcut.INTERSECT_EDGE
    assert np.allclose(intersections[1][1][1], 0.5)


def test_triangle_plane_intersection_edge_on_plane():
    verts = [
        (0, 0, 0),
        (0, 1, 0),
        (1, 0, 0),
    ]
    faces = [
        (0, 1, 2),
    ]

    mesh = meshcut.TriangleMesh(verts, faces)
    plane = meshcut.Plane(
        orig=(0, 0, 0),
        normal=(0, 1, 0)
    )

    intersections = meshcut.compute_triangle_plane_intersections(
        mesh, 0, plane)

    assert len(intersections) == 2
