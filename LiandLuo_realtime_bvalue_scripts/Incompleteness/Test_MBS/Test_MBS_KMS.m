%% Test whether MBS can obtain a credible Mc for the KMS method
clc,clear
Sim = 1e4; % number of simulated catalogs
L = 1e4; % number of events in each simulated catalog (initial complete catalog)
b=1;
delta=0.05;
Magn_Start = -delta;
mu=0.8;
sigma0=0.2;
p=parpool(6);    
parfor i = 1 : Sim 
    Pcm=rand(1,L);
    Magn0=-1/b.*log10(Pcm)+ Magn_Start;
    cri=normcdf(Magn0,mu,sigma0);
    ran_cri=rand(1,L);
    jkf=find(ran_cri<=cri);
    Magn1=Magn0(jkf);
    Magn2=roundn(Magn1,log10(2*delta));
    [Mc_kms(i),b_kms_con(i),number_kms(i)]=MBS_KMS(Magn1);
    [Mc_kms(i),b_kms_dis(i),number_kms(i)]=MBS_KMS(Magn2);
    Sim-i
end
delete(gcp('nocreate'))
%% plot
figure('units','normalized','position',[0.1,0.1,0.7,0.3])
edges = [0.7:0.02:1.3];
subplot(1,2,1)
h3 = histogram( b_kms_con,edges,'Normalization','probability') ;
h3.EdgeColor=[1 1 1];
h3.FaceColor = 1/255*[250 192 61];
hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
xlabel( 'b-value' )
ylabel( 'Frequency' )
set(gca,'fontsize',16);
grid on;box on;
xlim([0.7 1.3])
ylim([0 0.2])
title('KMS (continuous)');

subplot(1,2,2)
h3 = histogram( b_kms_dis,edges,'Normalization','probability') ;
h3.EdgeColor=[1 1 1];
h3.FaceColor = 1/255*[250 192 61];
hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
xlabel( 'b-value' )
ylabel( 'Frequency' )
set(gca,'fontsize',16);
grid on;box on;
xlim([0.7 1.3])
ylim([0 0.2])
title('KMS (discrete)');

mean(b_kms_con)
std(b_kms_con)
mean(b_kms_dis)
std(b_kms_dis)
