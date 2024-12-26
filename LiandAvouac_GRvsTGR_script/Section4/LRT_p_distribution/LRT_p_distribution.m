%{
clc,clear
b0=[0.6 0.8 1.0 1.2 1.4];
L0=[100 200 500 1000 2000];
zeta0=-0.5:0.1:1.4;
correction_LRT=zeros(length(b0),length(L0),length(zeta0));
times=500;
R=zeros(length(b0),length(L0),length(zeta0),times);
for k=1:length(b0)
    b=b0(k);
    for j=1:length(L0)
        L=L0(j);
        Target_LRT=zeros(length(zeta0),times);
        parpool(20);
        parfor jk=1:length(zeta0)
            zeta=zeta0(jk);         
            for i=1:times
                %% generate catalogs
                beta=2/3*b;
                mmin=0;
                Mmin=10^(1.5*(mmin+6.07));
                mcorner_TGR=log10(L)/b-zeta;
                M_corner=10^(1.5*(mcorner_TGR+6.07));

                U10=rand(1,L);
                U20=rand(1,L);
                U1=Mmin-M_corner*log(U10);
                U2=Mmin*U20.^(-1/beta);
                DM=min(U1,U2);
                Dm=2/3*log10(DM)-6.07;    

                %% Estimation
                [b_GR,loglikelihood_GR]=Estimation_GR(Dm,mmin);
                [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR(Dm,mmin);
                
                % Likelihood-ratio test
                R(k,j,jk,i)=loglikelihood_TGR-loglikelihood_GR;
                if 2*R(k,j,jk,i)>3.84
                    Target_LRT(jk,i)=2;
                else
                    Target_LRT(jk,i)=1;
                end

            end
        end
        delete(gcp('nocreate'))  
        data = R(k,j,:,:);
        data = squeeze(data);
        filename = sprintf('LRT_p_%d_%d.txt', k,j);
        dlmwrite(filename, data, 'delimiter', '\t');

    end
end
  
save('LRT_condition.mat')
%}

%% plot
clc,clear
b0=[0.6 0.8 1.0 1.2 1.4];
L0=[100 200 500 1000 2000];
zeta0=-0.5:0.1:1.4;
redShades = [
    255, 0, 0;
    255, 51, 51;
    255, 102, 102;
    255, 153, 153;
    255, 204, 204
];
blueShades = [
    0, 0, 255;
    51, 51, 255;
    102, 102, 255;
    153, 153, 255;
    204, 204, 255
];
blackShades = [
    50, 50, 50;
    100, 100, 100;
    150, 150, 150;
    200, 200, 200;
    225, 225, 225
];
purpleShades = [
    128, 0, 128;
    153, 51, 153;
    178, 102, 178;
    204, 153, 204;
    229, 204, 229
];
orangeShades = [
    255, 165, 0;
    255, 178, 51;
    255, 191, 102;
    255, 204, 153;
    255, 217, 204
];

colors = 1/255*[
    blackShades;
    purpleShades;
    blueShades;
    orangeShades;
    redShades;
];
%%
figure('units','normalized','position',[0.1,0.1,0.5,0.4])
for k=1:length(b0)
    b=b0(k);
    for j=1:length(L0)
        L=L0(j);
        filename =['LRT_p_',num2str(k),'_',num2str(j),'.txt'];
        D=load(filename);
        subplot(length(b0),length(L0),(k-1)*length(L0)+j)
        for jk=1:length(zeta0)
            color=summer(length(zeta0));
            P_LR=1-chi2cdf(2*D(jk,:),1);
            [h, stats] = ecdf(P_LR);
            plot(stats, h, 'Color', color(jk,:));hold on;
            legzeta{jk} = num2str(zeta0(jk));
            correction_LRT(k,j,jk)=length(find(P_LR<0.05))/length(P_LR);
        end
    fortitle=['$$b=',num2str(b),', n=',num2str(L), '$$'];
    title(fortitle, 'Interpreter', 'latex');
    xlim([0 1])
    xlabel('$$p$$', 'Interpreter', 'latex')
    ylabel('$$CDF$$', 'Interpreter', 'latex')
    grid on;box on;grid minor;
    set(gca,'fontsize',16)
    end
end

leg = legend(legzeta, 'Interpreter', 'latex', 'NumColumns', 1, 'Location', 'eastoutside');
title(leg, 'placeholder');
set(leg.Title, 'String', '$a/b-m_{corner}$', 'Interpreter', 'latex');

figure('units','normalized','position',[0.1,0.1,0.5,0.4])
for k=1:length(b0)
    b=b0(k);
    D=squeeze(correction_LRT(k,:,:));
    for j=1:length(L0)
        L=L0(j);
        plot(zeta0,D(j,:),'color',colors(5*(k-1)+j, :),'linewidth',2.5);
        hold on;
        legendEntries{length(L0)*(k-1)+j} = sprintf('$b=%.1f$, $n=%d$', b, L);
    end
end
set(gca, 'XDir','reverse'); 
grid on; box on;grid minor;
set(gca,'fontsize',16)
xlabel('$a/b-m_{corner}$', 'Interpreter', 'latex');
ylabel('$R$', 'Interpreter', 'latex');
legend(legendEntries, 'Interpreter', 'latex','Numcolumns',2,'location','eastoutside')


figure('units','normalized','position',[0.1,0.1,0.5,0.4])
for jk=1:length(zeta0)
    subplot(4,5,jk)
    for k=1:length(b0)
        b=b0(k);
        for j=1:length(L0)
            L=L0(j);
            filename =['LRT_p_',num2str(k),'_',num2str(j),'.txt'];
            D=load(filename);
            P_LR=1-chi2cdf(2*D(jk,:),1);
            [h, stats] = ecdf(P_LR);
            plot(stats, h, 'Color', colors(5*(k-1)+j, :));hold on;
        end
    end
    fortitle = ['$$a/b - m_{corner} = ', num2str(zeta0(jk)), '$$'];
    title(fortitle, 'Interpreter', 'latex');
    xlim([0 1])
    xlabel('$$p$$', 'Interpreter', 'latex')
    ylabel('$$CDF$$', 'Interpreter', 'latex')
    grid on;box on;grid minor;
    set(gca,'fontsize',16)
end
legend(legendEntries, 'Interpreter', 'latex','Numcolumns',1,'location','eastoutside')
