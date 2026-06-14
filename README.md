# MarchenCode


This repo contains GPU-optimized Python code to compute the Green's functions from single-sided surface reflectivity data. Learners will also find code to perform Marchenko-based imaging using decomposed Green’s functions (upgoing). I embarked on this project in an effort to understand the Marchenko methods. For learners proficient in MATLAB, you can download the original version of the code in MATLAB [here](https://wiki.seg.org/wiki/Software:Marchenko_for_imaging) by [`Angus Lomas and Andrew Curtis`](https://www.geos.ed.ac.uk/~acurtis/assets/Lomas_Curtis_Geop_March_2019.pdf)

marchenko_gpu_optimized.ipynb computes the **Marchenko Green's functions** and **focusing functions** from the **input reflectivity**, **direct arrival**, **filter (theta)**, and **ground-truth Green’s functions** (for comparison).

