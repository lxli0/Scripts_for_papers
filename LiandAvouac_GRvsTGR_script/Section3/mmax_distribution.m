clc,clear
b=1;
N=1e3;
m0_exp=log10(N)/b;
for i=1:1e4
    Dm=-1/b.*log10(rand(1,N)); % generate random magnitude
    [Mmax(i),ind(i)]=max(Dm);
end
%
m=1:0.1:7;
f=b*log(10)*10.^(-b*m);
F=1-10.^(-b.*m);
ff=N.*f.*F.^(N-1);
% non-bias
dm=-1:0.1:4;
for i=1:length(dm)
    Nmin(i)=ceil(10^(-b*dm(i)));
    N=Nmin:1:1e3;
    somefun=1./N.*(1-1./N*10^(-b*dm(i))).^(max(N)-1);
    fff(i)=b*log(10)*10^(-b*dm(i))*sum(somefun);
end
%%
color=1/255*[132 94 194;178 91 0;0 139 200];

figure('units','normalized','position',[0.1,0.1,0.6,0.3])
subplot(1,2,1)
% total
h10=plot(m-m0_exp,ff,'linewidth',2.5,'color',color(1,:));
hold on;
delta_M_total=Mmax-m0_exp;
h11=histogram(delta_M_total,'Normalization','pdf', 'FaceColor',color(1,:), 'EdgeColor', 'none');
set(gca,'fontsize',16)
hold on;

%prior
h20=plot(dm,fff,'linewidth',2.5,'color',color(2,:));
Mmax_exp=log10(ind-1)/b;
delta_M_prior=Mmax-Mmax_exp;
h21=histogram(delta_M_prior,'Normalization','pdf', 'FaceColor',color(2,:), 'EdgeColor', 'none');
xlabel('$m_{max}-\hat{m}_{max}$', 'Interpreter', 'latex');
ylabel('Frequency');
set(gca,'fontsize',16)
grid on;box on;grid minor;
legend([h10,h11,h20,h21],'All (Theoretic)','All (Numerical)','Prior (Theoretic)','Prior (Numerical)');


subplot(1,2,2)
for i=1:length(Mmax)
    p1(i)=length(find(delta_M_total<=delta_M_total(i)))/length(delta_M_total); 
    p2(i)=length(find(delta_M_prior<=delta_M_prior(i)))/length(delta_M_prior); 
end
density_2D=density2C(p1',p2',0:0.1:1,0:0.1:1);
scatter(p1,p2,40,density_2D,'filled');
colormap('autumn')
grid on;box on;grid minor;
xlabel('$p$ (total)', 'Interpreter', 'latex');
ylabel('$p$ (prior)', 'Interpreter', 'latex');
set(gca,'fontsize',16)

function [h]=density2C(X,Y,XList,YList)
    [XMesh,YMesh]=meshgrid(XList,YList);
    XYi=[XMesh(:) YMesh(:)];
    F=ksdensity([X,Y],XYi);
    ZMesh=zeros(size(XMesh));
    ZMesh(1:length(F))=F;
    h=interp2(XMesh,YMesh,ZMesh,X,Y);
end