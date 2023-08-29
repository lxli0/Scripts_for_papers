clc,clear
%% Load Data
Cat_Raw=load( 'Cat_16Nov_ZMAP.txt' ) ;
% select the events with the previous threshold criteria
Center   = [ 13.324 , 44.013 ] ;  % Long and Lat
DistMax  = 30;                   % Km
DepthMax = 30;                   % Km
Cat = Cat_Raw( Cat_Raw( : , 7 ) <= DepthMax & ...
      distance( Cat_Raw( : , 2 ) , Cat_Raw( : , 1 ) , Center(2) , Center(1) ) ...
      .* pi/180*6371 <= DistMax , : ) ;
% select events based on STAI
DeltaT = 1/6 ; %in days
Time = datenum( Cat( : , [ 3 : 5 , 8 : 10 ]) ) - ...
        datenum( Cat( 1 , [ 3 : 5 , 8 : 10 ]) ) ;
Cat_NoSTAI = Cat( Time >= DeltaT , : ) ;
Magn=Cat_NoSTAI(:,6);
Dm=diff(Magn);
L=length(Dm);
dm=0.1;
m=floor(min(Dm)/dm)*dm:dm:floor(max(Dm)/dm)*dm;
n0=hist(Dm,m);
jkf=find(m>=-1e-6);
m=m(jkf);
n0=n0(jkf);
cn0(1)=sum(n0);
for i=2:length(m)
    cn0(i)=sum(n0)-sum(n0(1:i-1));
end
n=log10(n0);
cn=log10(cn0);

%% plot
figure('units','normalized','position',[0.1,0.1,0.3,0.3])
color=1/255*[178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 8, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 8, 'markerfacecolor', color(1,:), ...
            'MarkerEdgeColor', color(1,:));
hold on;
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
grid on;
box on;
xlabel('Magnitude Differnece');
ylabel('Number');
set(gca,'fontsize',16);
xlim([0 inf])
ylim([0 1e3])
