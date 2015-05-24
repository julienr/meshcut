========================================================
meshcut - Python utilities to slice 3D triangular meshes
========================================================

For now, this computes the planar cross-section of a 3D mesh.

Given a 3D mesh and a cut plane, this computes a (set of) polyline(s) that
results from cutting the mesh with the plane.

Requires python and numpy. Mayavi is used for visualisation in the examples.

Run with ::

    python examples/cross_section.py


.. image:: screenshot.png
   :width: 50%
