% The script describes how nonuniform magnitude errors bias the b-value estimation for the three methods (continuous magnitudes)
clc,clear
Sim=1e4; % number of random catalogs
L=1e4; % number of events in each simulated catalog
Magn_Start = 0 ;
Mc = 1.0 ;
Magn_Thr=1.05;
Mc_prime=0;
Error=[0.25 0.05];
b=1; % set b-value, we should expect the methods provide estimation close to this value
parpool(6);
parfor i=1:Sim
    Pcm=rand(1,L);
    Magn0=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitude
    error1=normrnd(0,Error(1),[1 L]);
    error2=normrnd(0,Error(2),[1 L]);
    jkf1=Magn0<Magn_Thr;
    jkf2=Magn0>=Magn_Thr;
    Magn=[Magn0(jkf1)+error1(jkf1),Magn0(jkf2)+error2(jkf2)];
    randIndex_A = randperm(L);
    Magn = Magn(randIndex_A);
    jkf3=find((Magn-Mc)>=-1e-6);
    Magn=Magn(jkf3); % filter earthquakes by Mc
  
    %MLE   
    b_mle(i)=1/(log(10)*(mean(Magn)-Mc));

    %b-positive
    Dm0=diff(Magn);
    jkf4=find((Dm0-Mc_prime)>=-1e-6);
    Dm=Dm0(jkf4);% filter earthquakes by Mc'
    b_pos(i)=1/(log(10)*(mean(Dm)-Mc_prime));

    %KMS
    b_kms(i)=KMS(Magn);
    
    Sim-i
end
delete(gcp('nocreate'))

%% plot
figure('units','normalized','position',[0.1,0.1,1,0.3])

subplot(1,3,1)
edges = 0.9:0.02:1.4;
h1 = histogram( b_mle,edges,'Normalization','probability') ;
h1.EdgeColor=[1 1 1];
h1.FaceColor = 1/255*[46 89 167];
hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
ylabel( 'Frequency' )
xlabel( 'b-value' )
set(gca,'fontsize',16);
grid on;box on;
ylim([0 0.3])
xlim([0.9 1.4])
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
ylim([0 0.3])
xlim([0.9 1.4])
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
ylim([0 0.3])
xlim([0.9 1.4])
title('KMS');

%%
mean(b_mle)
std(b_mle)
mean(b_pos)
std(b_pos)
mean(b_kms)
std(b_kms)
