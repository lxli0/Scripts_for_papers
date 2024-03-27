This is a repository containing the code to generate the results in the paper **"Can we obtain reliable seismic b-values for real-time catalogs?"**, by *Linxuan Li and Gang Luo*.

You will find:

"benchmark_three_method": 
    Use "check_random_number_for_KMS.m" to evaluate the rationality of the number of random catalogs (10) to obtain a credible b-value for the KMS method;
    Use "Benchmark.m" to generate Fig. 1.
    
"Incompleteness":
    Use "sensitivity_discrete.m" and "sensitivity_continous.m" to generate Figs 2 and S1;
    Use "Test_MBS_discrete.m" and "Test_MBS_continuous.m" to generate Figs 3 and S3;
    Use "Test_MBS_KMS.m" to generate Fig. S2.
    
"Typpical_case_configuration": for Fig. 4.

"Catalog_Size":
    Use "catalog_size_discrete.m" and "catalog_size_continuous.m" to generate Figs 5 and S4.
    
"STAI":
    Use "gen_catalog_main.m" to generate random catalogs used for test (stored in the two folders);
    Use "STAI_discrete.m" and "STAI_continuous.m" to generate Figs 6 and S5.
    
"Magnitude_uncertainty":
    Use "magnitude_uncertainty_discrete.m" and "magnitude_uncertainty_continuous.m" to generate Figs 7 and S6.
    
"Italy" and "Turkey":
    Use "Plot_magnitude_distribution.m" and "Plot_magnitude_difference_distribution.m" to generate Figs 8a, 8b, 9a, and 9b.
    Use "main.m" (Italy) and "different_Mc_Mc_pirme" (Turkey) to generate Figs 8c, 9c.
    Use "main.m" (Turkey) to generate Fig S7.
    
"Feasibility":
    Use "gen_catalog_main.m" to generate random catalogs used for test (stored in the two folders);
    Use "Feasibility_discrete.m" and "Feasibility_continuous.m" to generate the results in Tables 1 and S1.



We are pleased to answer any questions (email: lxli@caltech.edu or lucas.linxuan.li@gmail.com).


Reference: Linxuan Li and Gnag Luo (2024). Can we obtain reliable seismic b-values for real-time catalogs?.
