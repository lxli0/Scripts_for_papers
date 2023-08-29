clc,clear
D=load('Ding_Turekey.txt');
Dt0=D(:,1);
jkf=find(Dt0>=2023.1013004791);
Dt=Dt0(jkf);
Dm=D(jkf,6);

L=length(Dm);
dm=0.1;
m=0:dm:floor(max(Dm)/dm)*dm;
n0=hist(Dm,m);   
cn0(1)=L;
for i=2:floor(max(Dm)/dm)+1
    cn0(i)=L-sum(n0(1:i-1));
end
n=log10(n0);
cn=log10(cn0);
%% Calculation by MBS
[Mc_MBS,b,number_mle]=MBS_MLE_discrete(Dm);
ind=find(abs(m-Mc_MBS)<=1e-6);
jkf=Dm>=(Mc_MBS-1e-6);
Dm=Dm(jkf);
a=log10(cn0(ind))+b*m(ind);

%% plot
figure('units','normalized','position',[0.1,0.1,0.3,0.3])
color=1/255*[178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 6, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 6, 'markerfacecolor', color(1,:), ...
            'MarkerEdgeColor', color(1,:));
hold on;
scatter(Mc_MBS,cn0(ind-2),180,[0 0 0],'v','filled');hold on;
hold on;
hhf=find(cn0>=10);
plot_range=Mc_MBS:0.1:m(hhf(end));
y=a-b*plot_range;
plot(plot_range,10.^y,'-','color',[0 0 0],'linewidth',2.5);
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');

grid on;
box on;
xlabel('Magnitude');
ylabel('Number');
set(gca,'fontsize',16);



