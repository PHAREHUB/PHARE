# FieldData tests

These tests assess the FieldData copy and streaming features work as equivalent
SAMRAI components. Assuming a Yee grid layout in 1D, Ex and Ey are respectively
dual and primal, thus copy and streaming operation should be the same as those
of SAMRAI CellData and NodeData specific PatchData.

There is no SAMRAI PatchData able to describe the 2D and 3D Yee layout so these
tests are not easily generalized to 2D and 3D.

However FieldData are used in higher level tests.

