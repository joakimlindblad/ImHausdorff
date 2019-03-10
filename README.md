# ImHausdorff
Matlab code for computing Hausdorff distance between binary images or segmentation maps.

---

Noticing that all Hausdorff distance functions I could find on MATLAB File Exchange
seem to compute all pairwise point-to-point distances making them very slow for large
images, I felt compelled to write a version `imhausdorff` which uses `bwdist` and
thus scale linearly for large image data. The code should work for images of any
dimension. Do try the examples in the help text.

Including point cloud version `hausdorff` as well. Functional for small data sets.

/ Joakim
