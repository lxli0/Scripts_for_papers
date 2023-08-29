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
Dm=Cat_NoSTAI(:,6);

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

%% plot
figure('units','normalized','position',[0.1,0.1,0.35,0.3])
color=1/255*[178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 7, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 7, 'markerfacecolor', color(1,:), ...
            'MarkerEdgeColor', color(1,:));
hold on;
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
grid on;
box on;
xlabel('Magnitude');
ylabel('Number');
set(gca,'fontsize',16);
xlim([0.9 inf])

