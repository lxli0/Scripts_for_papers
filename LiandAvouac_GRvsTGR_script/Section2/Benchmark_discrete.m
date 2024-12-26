clc,clear
b=1;
L=1e3;
mcorner_TGR=2;
M_corner=10^(1.5*mcorner_TGR+9.1);   
delta=0.05;
mmin_set=0;
mmin_true=mmin_set-delta;
Mmin_set=10^(1.5*mmin_true+9.1);
parpool(8);
parfor i=1:1e3
    Dm_GR=-1/b.*log10(rand(1,L))+mmin_true;  
    Dm_GR=roundn(Dm_GR,log10(2*delta)); 
    [b_GR_est(i),loglikelihood_GR_rand]=Estimation_GR_discrete(Dm_GR,mmin_set,delta);

    Dm_TGR=gentgr(2/3*b,L,Mmin_set,M_corner);
    Dm_TGR=roundn(Dm_TGR,log10(2*delta)); 
    [b_TGR_est(i),mcorner_TGR_est(i),loglikelihood_TGR_rand]=Estimation_TGR_discrete(Dm_TGR,mmin_set,delta);
end
delete(gcp('nocreate'))

figure('units','normalized','position',[0.1,0.1,0.6,0.3])
subplot(1,3,1)
h=histogram(b_GR_est,'Normalization','probability');
h.FaceColor = 1/255*[46 89 167];
h.EdgeColor=1/255*[46 89 167];
hold on;
plot([1 1],[0 0.15],'--k')
mean_b_GR=mean(b_GR_est)
std_b_GR=std(b_GR_est)
set(gca,'fontsize',16)
xlabel('Estimated b-value (GR)')
ylabel('Frequency')
grid on;box on;

subplot(1,3,2)
h=histogram(b_TGR_est,'Normalization','probability');
h.FaceColor = 1/255*[46 89 167];
h.EdgeColor=1/255*[46 89 167];
hold on;
plot([1 1],[0 0.15],'--k')
mean_b_TGR=mean(b_TGR_est)
std_b_TGR=std(b_TGR_est)
xlabel('Estimated b-value (TGR)')
ylabel('Frequency')
set(gca,'fontsize',16)
grid on;box on;

subplot(1,3,3)
h=histogram(mcorner_TGR_est,'Normalization','probability');
h.FaceColor = 1/255*[46 89 167];
h.EdgeColor=1/255*[46 89 167];
hold on;
plot([2 2],[0 0.2],'--k')
mean_mcorner_TGR=mean(mcorner_TGR_est)
std_mcorner_TGR=std(mcorner_TGR_est)
xlabel('Estimated m_{corner} (TGR)')
ylabel('Frequency')
set(gca,'fontsize',16)
grid on;box on;


function [Dm]=gentgr(beta,L,Mmin,M_corner)
        U10=rand(1,L);
        U20=rand(1,L);
        U1=-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
end