===============================================
Python code to compute cross-section of 3D mesh
===============================================

Simple (and not thoroughly tested) implementation of a 3D mesh cross-section
algorithm.

Note that this works on some tests, but I'm not 100% confident the algorithm
handles special cases (if an edge lies exactly on the plane) correctly.

Given a 3D mesh and a cut plane, this computes a (set of) polyline(s) that
results from cutting the mesh with the plane.

Requires python and numpy. Mayavi is used for visualisation.

Run with ::

    python cross_section.py


.. image:: screenshot.png
   :width: 50%
