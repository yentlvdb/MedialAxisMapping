# MedialAxisMapping

This code implements the method that is described in the work entitled 'Corelative light and electron microscopy (CLEM)
image registration by blood vessel axes alignment'

This work proposes a method to align a subset of challenging CLEM alignment problems where an initial alignment is potentially difficult to find. Our method bases its alignment on a structure nicely spread throughout tissues and that is stainable in LM: the blood vessels. It focuses on aligning the blood vessel medial axes rather than the contour. It combines a feature-based and intensity-based approach: first a set of potential global overlays is selected based on the blood vessel bifurcations present in the volume. From this set a final overlay is selected based on an intensity similarity metric.

All code can be found back in algorithm.py (supplementary info folder).

Data to test the algorithm is provided in the data folder. In this folder the EM volume and LM volume (blood vessel and nucleus staining) of a sample containing choroid plexus tissue can be found back.

A notebook describing the algorithm using a manual input step by step can be found back in the notebooks folder.
