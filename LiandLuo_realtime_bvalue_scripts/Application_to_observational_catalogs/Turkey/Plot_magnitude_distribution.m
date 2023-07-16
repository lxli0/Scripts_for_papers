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
%% Mc
%maximum curvatrue
[~,ind]=max(n0);
Mc_cur=m(ind)+0.2; 
ind1=ind+2;
%MBS
[Mc_MBS,b_mle00,number_mle]=MBS_MLE(Dm);
ind2=find(abs(m-Mc_MBS)<=1e-6);

% b
jkf1=Dm>=(Mc_cur-1e-6);
Dm1=Dm(jkf1);
b1=1/(log(10)*(mean(Dm1)-Mc_cur+0.05));
a1=log10(cn0(ind1))+b1*m(ind1);

jkf2=Dm>=(Mc_MBS-1e-6);
Dm2=Dm(jkf2);
b2=1/(log(10)*(mean(Dm2)-Mc_MBS+0.05));
a2=log10(cn0(ind2))+b2*m(ind2);
%% plot
figure('units','normalized','position',[0.1,0.1,0.35,0.3])
color=1/255*[178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 6, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 6, 'markerfacecolor', color(1,:), ...
            'MarkerEdgeColor', color(1,:));
hold on;
scatter(Mc_cur,cn0(ind),180,[0.5 0.5 0.5],'v','filled');hold on;
hold on;
scatter(Mc_MBS,cn0(ind2-2),180,[0 0 0],'v','filled');hold on;
hold on;
hhf=find(cn0>=10);
plot_range1=Mc_cur:0.1:m(hhf(end))-1;
y1=a1-b1*plot_range1;
plot(plot_range1,10.^y1,'-','color',[0.5 0.5 0.5],'linewidth',2.5);
hold on;
plot_range2=Mc_MBS:0.1:m(hhf(end));
y2=a2-b2*plot_range2;
plot(plot_range2,10.^y2,'-','color',[0 0 0],'linewidth',2.5);
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');

grid on;
box on;
xlabel('Magnitude');
ylabel('Number');
set(gca,'fontsize',16);



