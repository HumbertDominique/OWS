This document does not cover the conda environment setup. The requierments for use with a GPU are located in "requirements-TF-GPU.txt". For use without GPU, install "requirements.txt".

# Data synthetization

1. Input the variable values in the "Inputs" section,
1. Setup the data to save in the "Phase screen & PSF" section,
1. Run the required sections.

All section preceding the desired sections must be run at least once before in every session.


# Machine learning

1. Input the variable values in the "Inputs" section,
1. Setup the dataset assembly in the "Assemble dataset" section,
1. Setup the required dataset in the "Prepare the dataset" section,
1. Setup the required dataset in the "Inference" section. Note that manipulation of the dataset at the revious point must be repeated here.
1. Run the required sections.

All section preceding the desired sections must be run at least once before in every session. Only "Assemble dataset" can be run only once as long as the .npz file are not erased.
