%% generate catalogs
clc,clear
Sim_num=1e3;
p=parpool(6);    
parfor i=1:Sim_num
    mp=7;
    [T0,M0]=gen_catalog0(mp,i);
end
delete(gcp('nocreate'))