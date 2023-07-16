% Ths script provides an illustration of small catalog size
clc,clear
figure('units','normalized','position',[0.1,0.1,0.25,0.25])
L=50;% number of events in each simulated catalog
b=1; % set b-value
Pcm=rand(1,L);
Dm=-1/b.*log10(Pcm)-0.05; % generate random magnitude
% plot
dm=0.1;
m=0:dm:floor(max(Dm)/dm)*dm;
n0=hist(Dm,m);   
cn0(1)=L;
for i=2:floor(max(Dm)/dm)+1
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
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
hold on;


