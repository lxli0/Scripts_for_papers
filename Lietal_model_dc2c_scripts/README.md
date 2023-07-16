This is a repository containing the code to generate the results in the paper **"Modeling seismic velocity changes caused by seasonal hydrosphere variations"**, by *Linxuan Li, Jiangtao Li, and Xiaodong Song*.

**Data files:**
"Bao_TP_velocity.txt" and "Bao_BB_velocity.txt": S velocity models
"Bao_TP_sensitivity.txt" and "Bao_BB_sensitivity.txt": the computed sensitivity kernels

**Main scripts:**
"velocity_KV.m" is used to plot Figure 3
"main.m" is used to calculate the dc/c under different GWS assumptions and hydraulic coefficients
"out_model_result_maintextt.m" is used to generate Figures 4, 5, and S2.
"test_stress_sensitivity.m" is used to generate the results in Section 4.2.
"out_test_stress_sensitivity" is used to generate Figures 6 and S5.
"out_model_result_SI.m" is used to generate Figures S3 and S4.
"Map" is used to generate Figure S1 using GMT.

**Reference:**
Li, L., Li, J., & Song, X.. (in preparation). Modeling seismic velocity changes caused by seasonal hydrosphere variations.
