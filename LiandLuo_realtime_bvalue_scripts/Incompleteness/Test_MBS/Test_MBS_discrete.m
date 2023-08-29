%% Test whether the Mc/Mc' obtained can help estimate credible b-value (discrete magnitudes)
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
    Magn=Magn0(jkf);
    Magn=roundn(Magn,log10(2*delta));
    num0(i)=length(Magn);
    [Mc_mle(i),b_mle(i),number_mle(i)]=MBS_MLE_discrete(Magn);
    [Mc_pos(i),b_pos(i),number_pos(i)]=MBS_pos_discrete(Magn);
    Mc_kms(i)=Mc_mle(i);
    Magn_kms=Magn((Magn-Mc_kms(i))>=-1e-6);
    number_kms(i)=length(Magn_kms);
    b_kms(i)=KMS(Magn_kms');
    Sim-i
end
delete(gcp('nocreate'))
%% plot
figure('units','normalized','position',[0.1,0.1,1,0.3])

subplot(1,3,1)
edges = [0.7:0.02:1.3];
h1 = histogram( b_mle,edges,'Normalization','probability') ;
h1.EdgeColor=[1 1 1];
h1.FaceColor = 1/255*[46 89 167];
hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
ylabel( 'Frequency' )
xlabel( 'b-value' )
set(gca,'fontsize',16);
grid on;box on;
xlim([0.7 1.3])
ylim([0 0.2])
title('MLE');

subplot(1,3,2)
h2 = histogram( b_pos,edges,'Normalization','probability') ;
h2.EdgeColor=[1 1 1];
h2.FaceColor = 1/255*[209 41 32];
hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
ylabel( 'Frequency' )
xlabel( 'b-value' )
set(gca,'fontsize',16);
grid on;box on;
xlim([0.7 1.3])
ylim([0 0.2])
title('b-positive');

subplot(1,3,3)
h3 = histogram( b_kms,edges,'Normalization','probability') ;
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
title('KMS');
%%
mean(b_mle)
std(b_mle)
mean(b_pos)
std(b_pos)
mean(b_kms)
std(b_kms)