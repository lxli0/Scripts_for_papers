clc,clear
b=1;
N=1e3;
m0_exp=log10(N)/b;
for i=1:1e4
    Dm=-1/b.*log10(rand(1,N)); % generate random magnitude
    dDm=diff(Dm);
    dDm_pos=dDm(dDm>=0);
    
    b_mle(i)=1/log(10)/mean(Dm);
    b_pos(i)=1/log(10)/mean(dDm_pos);
    
    mmax(i)=max(Dm);
    dmmax(i)=max(dDm_pos);
   
    F_m=1-10.^(-b_mle(i)*mmax(i));
    F_dm_pos=1-10.^(-b_pos(i)*dmmax(i));
    
    P_m(i)=F_m.^N;
    P_dm_pos(i)=F_dm_pos^(length(dDm_pos));
    
end
figure('units','normalized','position',[0.1,0.1,0.5,0.3])
subplot(1,2,1)
density_2D=density2C(b_mle',b_pos',0.6:0.05:1.4, 0.6:0.05:1.4);
scatter(b_mle,b_pos,40,density_2D,'filled');
colormap('autumn')
grid on;box on;grid minor;
xlabel('$b_{MLE}$', 'Interpreter', 'latex');
ylabel('$b_{pos}$', 'Interpreter', 'latex');
set(gca,'fontsize',16)
xlim([0.9 1.1])
ylim([0.9 1.1])

subplot(1,2,2)
density_2D=density2C(P_m',P_dm_pos',0:0.1:1,0:0.1:1);
scatter(P_m,P_dm_pos,40,density_2D,'filled');
colormap('autumn')
grid on;box on;grid minor;
xlabel('$p(m_{max})$', 'Interpreter', 'latex');
ylabel('$p(\delta m_{max})$', 'Interpreter', 'latex');
set(gca,'fontsize',16)

function [h]=density2C(X,Y,XList,YList)
    [XMesh,YMesh]=meshgrid(XList,YList);
    XYi=[XMesh(:) YMesh(:)];
    F=ksdensity([X,Y],XYi);
    ZMesh=zeros(size(XMesh));
    ZMesh(1:length(F))=F;
    h=interp2(XMesh,YMesh,ZMesh,X,Y);
end
