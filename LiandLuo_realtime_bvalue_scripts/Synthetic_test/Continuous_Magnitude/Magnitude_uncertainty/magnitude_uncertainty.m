% The script describes how nonuniform magnitude errors bias the b-value estimation for the three methods (continuous magnitudes)
clc,clear
Sim=1e3; % number of random catalogs
L=1e4; % number of events in each simulated catalog
Magn_Start = 0 ;
Mc = 1.0 ;
Magn_Thr=1.15;
Mc_prime=0;
Error=[0.25 0.05];
b=1; % set b-value, we should expect the methods provide estimation close to this value
parpool(8);
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
    size=length(Magn);
    kms0=-15.15*(log10(size))^(-2.14)+11.85;
    kms=zeros(1,10);
    for l=1:10
        dt_Poisson=exprnd(1,1,size-1);
        t_Poisson=cumsum(dt_Poisson);
        t_Poisson=[0,t_Poisson];
        randIndex_A = randperm(size);
        Magn_cal= Magn(randIndex_A);
        kms(l)=VGA(t_Poisson,Magn,Mc);
    end
    im_kms=mean(kms);
    b_kms(i)=im_kms/kms0;
    
    Sim-i
end
delete(gcp('nocreate'))

% plot
figure('units','normalized','position',[0.1,0.1,0.3,0.7])
edges = 0.9:0.02:1.4;

subplot(3,1,1)
h1=histogram( b_mle,edges,'Normalization','probability');
h1.EdgeColor=[1 1 1];
h1.FaceColor=1/255*[48 80 149];
set(gca,'fontsize',16)
ylabel( 'Frequency' )
set(gca,'xticklabel',[]);
ylim([0 0.3])
xlim([0.9 1.4])
grid on;box on;hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
set(gca,'position',[0.15 0.69 0.75 0.25])

subplot(3,1,2)
h2=histogram( b_pos,edges,'Normalization','probability');
h2.EdgeColor=[1 1 1];
h2.FaceColor=1/255*[168 88 43];
set(gca,'fontsize',16)
ylabel( 'Frequency' )
set(gca,'xticklabel',[]);
ylim([0 0.3])
xlim([0.9 1.4])
grid on;box on;hold on;
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
set(gca,'position',[0.15 0.40 0.75 0.25])

subplot(3,1,3)
h3=histogram( b_kms,edges,'Normalization','probability');
h3.EdgeColor=[1 1 1];
h3.FaceColor=1/255*[179 142 72];
hold on;
xlabel( 'b-value' )
ylabel( 'Frequency' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
grid on;
set(gca,'fontsize',16)
ylim([0 0.3])
xlim([0.9 1.4])
set(gca,'position',[0.15 0.11 0.75 0.25])

% calculate KMS for a particular catalog
function [Slope] =VGA(Dt0,Dm0,Mc)
    ij=find((Dm0-Mc)>=-1e-6);
    Dt=Dt0(ij);
    Dm=Dm0(ij);
    L=length(Dm);
    VG=2*ones(1,L);
    VG(1)=1;
    VG(end)=1;
    for i=1:L-2
        for j=i+2:L
            tem_Dm=Dm(i+1:j-1);
            tem_Dt=Dt(i+1:j-1);
            cri_Dm=Dm(i)+(Dm(j)-Dm(i))*(tem_Dt-Dt(i))/(Dt(j)-Dt(i));
            cri=tem_Dm-cri_Dm;
            if max(cri)<0
                VG(i)=VG(i)+1;
                VG(j)=VG(j)+1;
            end
        end
    end
    p=polyfit(Dm,VG,1);
    Slope=p(1);
end
