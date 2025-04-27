This is a repository containing the code to generate the results in the paper **"Depth-Dependent Stress Sensitivity and Its Influence on Temporal Variations in Low-Frequency Rayleigh Wave Velocities"**, by *Linxuan Li1, Jiangtao Li1, and Xiaodong Song*.

**Data files:**

"Bao_TP_input.txt" and "Bao_BB_input.txt": S velocity models

"Bao_TP_sensitivity.txt" and "Bao_BB_sensitivity.txt": the computed sensitivity kernels (K_V)


**Main scripts:**

"fit_in_situ.m" is used to predict the in-situ measurements in Table 1.

"z_h.m" is used to plot Figures 3 and S3–S5 (one has to manually change the input stress perurbations).

"nominal_eta" is used to plot Figure 4.

"velocity_model.m" is used to plot Figures S1b and S1c.

"plt_EQ_crosssection.m" is used to plot Figure S2c.

**Folder "Lab"**

Contains the the scripts and original data (.xlsx) for reproducing Figure 1. One can use ./(Un)Conolidated_Lab/Out_sensitivity/(un)consolidated_plot.m to generate Figure 1.

**Reference:**

Li, L., Li, J., & Song, X. (2025). Depth-Dependent Stress Sensitivity and Its Influence on Temporal Variations in Low-Frequency Rayleigh Wave Velocities.

**Contact:**

Linxuan Li: lxli@caltech.edu or lucas.linxuan.li@gmail.com

