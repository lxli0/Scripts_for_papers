% The script analyze the standard error for different catalog sizes (continuous magnitudes)
clc,clear
figure('units','normalized','position',[0.1,0.1,0.3,0.7])
load('lajolla')
color=flip(lajolla);

L=[50,100,200,350,500,750,1000,1250,1500,1750,2000]; % number of events in each simulated catalog
Sim=floor(1e5./L); % number of random catalogs
b0=[0.6,0.8,1.0,1.2,1.4]; % set b-value
sigma_MLE=1./sqrt(L); % equation from Aki 1965
   
for i=1:length(b0)
    b=b0(i);
    parpool(6); 
    parfor j=1:length(Sim)
        b_mle=zeros(1,Sim(j));
        b_pos=zeros(1,Sim(j));
        b_kms=zeros(1,Sim(j));
        for k=1:Sim(j)
            Pcm=rand(1,3*L(j));
            Magn0=-1/b.*log10(Pcm); % generate random magnitude
            Mmin=0;
            %MLE
            Magn=Magn0(1:L(j));
            b_mle(k)=1/(log(10)*(mean(Magn)-Mmin));
            
            %b-positive
            Dm00=diff(Magn0);
            jkf=Dm00>=-1e-6;
            Dm0=Dm00(jkf);
            Dm=Dm0(1:L(j));
            b_pos(k)=1/(log(10)*(mean(Dm)-Mmin));
           
            %KMS
            kms0=-15.15*(log10(L(j)))^(-2.14)+11.85; % KMS/b ratio
            kms=zeros(1,10);
            for l=1:10
                dt_Poisson=exprnd(1,1,L(j)-1);
                t_Poisson=cumsum(dt_Poisson);
                t_Poisson=[0,t_Poisson]; % generate random occurrence time
                randIndex_A = randperm(L(j));
                Magn = Magn(randIndex_A); % randomly arrange the occurrence order
                kms(l)=VGA(t_Poisson,Magn,Mmin); % calculate KMS for the generated catalog
            end
            im_kms=mean(kms);
            b_kms(k)=im_kms/kms0;
        end
        std_mle(i,j)=std(b_mle)/b;
        std_pos(i,j)=std(b_pos)/b;
        std_kms(i,j)=std(b_kms)/b;        
    end
    delete(gcp('nocreate'))
end

% plot    
subplot(3,1,1)
for i=1:length(b0)
    h(i)=scatter(L,std_mle(i,:),40,color(214-140/length(b0)*(i-1),:),'filled');
    hold on;
end
h(length(b0)+1)=plot(L,sigma_MLE,'linewidth',1.5);
grid on;
box on;
xlabel('Number of Magnitudes')
ylabel('\sigma/\itb','Fontname', 'Times New Roman')
ylim([0 0.22])
set(gca,'fontsize',16)
legend(h,'b=0.6','b=0.8','b=1.0','b=1.2','b=1.4','Aki, 1965','NumColumns',2,'Location','Northeast')

subplot(3,1,2)
for i=1:length(b0)
    scatter(L,std_pos(i,:),40,color(214-140/length(b0)*(i-1),:),'filled');
    hold on;
end
plot(L,sigma_MLE,'linewidth',1.5)
grid on;
box on;
xlabel('Number of Magnitude Differences')
ylabel('\sigma/\itb','Fontname', 'Times New Roman')
ylim([0 0.22])
set(gca,'fontsize',16)

subplot(3,1,3)
for i=1:length(b0)
    scatter(L,std_kms(i,:),40,color(214-140/length(b0)*(i-1),:),'filled');
    hold on;
end
mean_std_kms=mean(std_kms,1);
nihe_y=log(mean_std_kms);
nihe_x=log(L);
p=polyfit(nihe_x,nihe_y,1); 
w1=exp(p(2));w2=p(1);
hh=plot(L,w1.*L.^(w2),'linewidth',1.5,'color','k');
legend(hh,'fitting curve')
grid on;
box on;
xlabel('Number of Magnitudes')
ylabel('\sigma/\itb','Fontname', 'Times New Roman')
ylim([0 0.22])
set(gca,'fontsize',16)

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

