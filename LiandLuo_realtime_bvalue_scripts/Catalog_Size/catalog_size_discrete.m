% The script analyze the standard error for different catalog sizes (discrete magnitudes)
clc,clear
L=[50,100,200,350,500,750,1000,1250,1500,1750,2000]; % number of events in each simulated catalog
Sim=floor(1e7./L); % number of random catalogs
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
            Pcm=rand(1,10*L(j));
            delta=0.05; % magnitudes are discretized by 2*delta
            Magn_Start = -delta;% starting magnitude for simulations
            Magn00=-1/b.*log10(Pcm)+Magn_Start; % generate random magnitude
            Magn0=roundn(Magn00,log10(2*delta)); % magnitude discretization
            Mmin=0;
            %MLE
            Magn=Magn0(1:L(j));
            b_mle(k)=1/(log(10)*(mean(Magn)-Mmin+delta));
            
            %b-positive
            Dm00=diff(Magn0);
            Mc_prime=2*delta;
            jkf=find((Dm00-Mc_prime)>=-1e-6);
            Dm0=Dm00(jkf);
            Dm=Dm0(1:L(j));
            b_pos(k)=1/delta/log(10)*acoth((mean(Dm)-Mc_prime+delta )/delta); 
           
            %KMS
            b_kms(k)=KMS(Magn);
        end
        std_mle(i,j)=std(b_mle)/b;
        std_pos(i,j)=std(b_pos)/b;
        std_kms(i,j)=std(b_kms)/b;        
    end
    delete(gcp('nocreate'))
end

% plot
figure('units','normalized','position',[0.1,0.1,1,0.3])
str=['o','x','+','s','d'];
color=1/255*[46 89 167;209 41 32;250 192 61];

subplot(1,3,1)
for i=1:length(b0)
    h(i)=plot(L,std_pos(i,:),str(i),'color',color(1,:));
    hold on;
end
h(length(b0)+1)=plot(L,sigma_MLE,'linewidth',1.5,'color',color(2,:));
grid on;
box on;
xlabel('\itn')
ylabel('\sigma/\itb','Fontname', 'Times New Roman')
ylim([0 0.22])
set(gca,'fontsize',16)
legend(h,'b=0.6','b=0.8','b=1.0','b=1.2','b=1.4','Aki, 1965','NumColumns',2,'Location','Northeast')
title('MLE');
subplot(1,3,2)
for i=1:length(b0)
    plot(L,std_pos(i,:),str(i),'color',color(1,:));
    hold on;
end
plot(L,sigma_MLE,'linewidth',1.5,'color',color(2,:))
grid on;
box on;
xlabel('\itn')
ylabel('\sigma/\itb','Fontname', 'Times New Roman')
ylim([0 0.22])
set(gca,'fontsize',16)
title('b-positive');

subplot(1,3,3)
for i=1:length(b0)
    plot(L,std_kms(i,:),str(i),'color',color(1,:));
    hold on;
end
mean_std_kms=mean(std_kms,1);
nihe_y=log(mean_std_kms);
nihe_x=log(L);
p=polyfit(nihe_x,nihe_y,1); 
w1=exp(p(2));w2=p(1);
hh=plot(L,w1.*L.^(w2),'linewidth',1.5,'color',color(3,:));
legend(hh,'fitting curve')
grid on;
box on;
xlabel('\itn')
ylabel('\sigma/\itb','Fontname', 'Times New Roman')
ylim([0 0.22])
set(gca,'fontsize',16)
title('KMS');

