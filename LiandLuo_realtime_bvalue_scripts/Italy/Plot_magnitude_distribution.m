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
Dm=Cat_NoSTAI(:,6)+1e-12;

L=length(Dm);
dm=0.1;
m=0:dm:max(Dm);
n0=hist(Dm,m);   
cn0(1)=L;
for i=2:length(m)
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
figure('units','normalized','position',[0.1,0.1,0.25,0.25])
color=1/255*[178 91 0;0 139 200];
plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 6, 'markerfacecolor', color(2,:), ...
            'MarkerEdgeColor', color(2,:));
hold on;
plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
            'markersize', 6, 'markerfacecolor', color(1,:), ...
            'MarkerEdgeColor', color(1,:));
hold on;
scatter(Mc_MBS,cn0(ind)*1.4,180,[0 0 0],'v','filled'); 
%"time 1.4" is to make sure that the symbol does not overlap with the MFD 
hhf=find(cn0>=1);
plot_range=Mc_MBS:0.1:m(hhf(end));
y=a-b*plot_range;
plot(plot_range,10.^y,'-','color',[0 0 0],'linewidth',2.5);
legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
grid on;
box on;
xlabel('Magnitude');
ylabel('Number');
set(gca,'fontsize',16);
xlim([0.9 inf])
ylim([1 1e3])
