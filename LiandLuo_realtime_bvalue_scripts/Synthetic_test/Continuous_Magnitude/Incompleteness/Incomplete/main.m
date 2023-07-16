%% Test whether the Mc/Mc' obtained can help estimate credible b-value (continuous magnitudes)
clc,clear
Sim = 1e3;% number of simulated catalogs
L = 6e3;% number of events in each simulated catalog
b=1;
Magn_Start = 0;
mu=0.8;
sigma0=0.2;
p=parpool(8);    
parfor i = 1 : Sim 
    Pcm=rand(1,L);
    Magn0=-1/b.*log10(Pcm)+ Magn_Start;
    cri=normcdf(Magn0,mu,sigma0);
    ran_cri=rand(1,L);
    jkf=find(ran_cri<=cri);
    Magn=Magn0(jkf);
    num0(i)=length(Magn);
    [Mc_mle(i),b_mle(i),number_mle(i)]=MBS_MLE(Magn);
    [Mc_pos(i),b_pos(i),number_pos(i)]=MBS_pos(Magn);  
    Mc_kms(i)=Mc_mle(i);
    number_kms(i)=length(Magn(Magn>=Mc_kms(i)));
    b_kms(i)=KMS(Magn,Mc_kms(i));
    Sim-i
end
delete(gcp('nocreate'))
figure('units','normalized','position',[0.1,0.1,0.9,0.3])
subplot(1,3,1)
edges = [0.7:0.02:1.3];
h1 = histogram( b_mle,edges,'Normalization','probability') ;
h1.EdgeColor=[1 1 1];
h1.FaceColor = 1/255*[46 89 167];
hold on;
h2 = histogram( b_pos,edges,'Normalization','probability') ;
h2.EdgeColor=[1 1 1];
h2.FaceColor = 1/255*[209 41 32];
h3 = histogram( b_kms,edges,'Normalization','probability') ;
h3.EdgeColor=[1 1 1];
h3.FaceColor = 1/255*[250 192 61];
xlabel( 'b-value' )
ylabel( 'Frequency' )
set(gca,'fontsize',16);
grid on;box on;
xlim([0.7 1.3])
ylim([0 0.15])
subplot(1,3,2)
edges = [0:0.1:1.6];
h1 = histogram( Mc_mle,edges,'Normalization','probability');
h1.EdgeColor=[1 1 1];
h1.FaceColor = 1/255*[46 89 167];
hold on;
h2 = histogram( Mc_pos,edges,'Normalization','probability');
h2.EdgeColor=[1 1 1];
h2.FaceColor = 1/255*[209 41 32];
hold on;
h3 = histogram( Mc_kms,edges,'Normalization','probability');
h3.EdgeColor=[1 1 1];
h3.FaceColor = 1/255*[250 192 61];
xlabel('Mc or Mc''');
ylabel('Frequency');
xlim([-0.1 1.5])
ylim([0 0.5])
set(gca,'fontsize',16);
grid on;box on;
subplot(1,3,3)
edges = [100:25:800];
h1 = histogram( number_mle,edges,'Normalization','probability');
h1.EdgeColor=[1 1 1];
h1.FaceColor = 1/255*[46 89 167];
hold on;
h2 = histogram( number_pos,edges,'Normalization','probability');
h2.EdgeColor=[1 1 1];
h2.FaceColor = 1/255*[209 41 32];
hold on;
h3 = histogram( number_kms,edges,'Normalization','probability');
h3.EdgeColor=[1 1 1];
h3.FaceColor = 1/255*[250 192 61];
legend([h1,h2,h3],'MLE','b-positive','KMS');
xlim([100,800])
ylim([0 0.2])
xlabel('Sample Number');
ylabel('Frequency');
set(gca,'fontsize',16);
grid on;box on;


mean(b_mle)
mean(b_pos)
mean(b_kms)
