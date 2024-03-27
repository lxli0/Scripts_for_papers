% Ths script provides an illustration of inconsistent magnitude errors
clc,clear
figure('units','normalized','position',[0.1,0.1,0.25,0.25])
L=1e4*1e4; 
Magn_Start =0; % starting magnitude for simulations
Mc = 1.0 ;
Magn_Thr=1.15;
Mc_prime=0;
Error=[0.25 0.05];
b=1; % set b-value
Pcm=rand(1,L);
Magn0=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitude
error1=normrnd(0,Error(1),[1 L]);
error2=normrnd(0,Error(2),[1 L]);
jkf1=Magn0<Magn_Thr;
jkf2=Magn0>=Magn_Thr;
Magn=[Magn0(jkf1)+error1(jkf1),Magn0(jkf2)+error2(jkf2)]; % add magnitude uncertainty
randIndex_A = randperm(L);
Dm = Magn(randIndex_A); %the generated catalog

% plot
dm=0.1;
m=0:dm:max(Dm);
n0=hist(Dm,m);   
cn0(1)=L;
for i=2:length(m)
    cn0(i)=L-sum(n0(1:i-1)); 
end
n=log10(n0);
cn=log10(cn0);
color=1/255*[132 94 194;178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 10, 'markerfacecolor', color(3,:), ...
            'MarkerEdgeColor', color(3,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 10, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
grid on;
box on;
xlabel('Magnitude');
ylabel('Number');
set(gca,'fontsize',16);
xlim([1 2.5])
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');