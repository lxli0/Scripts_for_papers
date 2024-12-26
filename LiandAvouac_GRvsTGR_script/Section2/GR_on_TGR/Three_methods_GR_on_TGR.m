%{
clc,clear
b0=[0.6 1 1.5];
L0=[50 100 150 200 300 400 600 800 1e3];
zeta0=0.1:0.1:1;
mmin_set=0;
mmin_true=mmin_set;
Mmin_set=10^(1.5*(mmin_true+6.07));
for s=1:length(b0)
    b=b0(s);
    for j=1:length(zeta0)
        zeta=zeta0(j);
        parpool(9);
        parfor k=1:length(L0)
            L=L0(k);
            mcorner_TGR=zeta*log10(L)/b;
            M_corner=10^(1.5*(mcorner_TGR+6.07));   
            b_mle=zeros(1,1e3);
            b_pos=zeros(1,1e3);
            b_kms=zeros(1,1e3);
            ll=zeros(1,1e3);
            for i=1:1e3
                Dm_TGR=gentgr(2/3*b,L,Mmin_set,M_corner);
                Mmin=0;
                Dm_TGR=Dm_TGR(Dm_TGR>=Mmin)
                ll(i)=length(Dm_TGR);
                %MLE
                b_mle(i)=1/(log(10)*(mean(Dm_TGR)-Mmin));

                %b-positive
                Dm00=diff(Dm_TGR);
                jkf=Dm00>=0;
                Dm=Dm00(jkf);
                b_pos(i)=1/(log(10)*(mean(Dm)-Mmin));

                %KMS
                b_kms(i)=KMS(Dm_TGR);
            end
            mean_mle(s,j,k)=mean(b_mle);
            mean_pos(s,j,k)=mean(b_pos);
            mean_kms(s,j,k)=mean(b_kms);   
            std_mle(s,j,k)=std(b_mle);
            std_pos(s,j,k)=std(b_pos);
            std_kms(s,j,k)=std(b_kms);
            mean_ll(s,j,k)=mean(ll);
        end
        delete(gcp('nocreate'))
    end
end
save('GRonTGR.mat')

%}
clc,clear
load('GRonTGR.mat');
fig=figure('units','normalized','position',[0.1,0.1,0.6,0.6]);
load('batlow.mat');
colormap_input=batlow;
original_indices = linspace(1, size(colormap_input, 1), size(colormap_input, 1));
target_indices = linspace(1, size(colormap_input, 1), length(zeta0));
color_plt = interp1(original_indices, colormap_input, target_indices);
%color_plt=jet(length(zeta0));
for s=1:length(b0)
    b=b0(s);
    subplot(3,3,s)
    for j=1:length(zeta0)
        scatter(L0,mean_mle(s,j,:)./b,80,color_plt(j,:),'filled');
        hold on;
    end
    ylabel('$b_{MLE}/b_{set}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
    ylim([0.9 4])
    set(gca,'XTicklabel', []);
    %fortitle = ['\textbf{\textit{b}=', num2str(b), '}'];
    %title(fortitle, 'Interpreter', 'latex');
    axes('Position',[0.24+0.28*(s-1) 0.89 0.1 0.1])
    box on
    for j=1:length(zeta0)
        scatter(L0,std_mle(s,j,:)./mean_mle(s,j,:),80,color_plt(j,:),'filled');
        hold on;
    end
    hh=plot(L0,1./sqrt(L0),'k','linewidth',2);
    if s==1
        legend(hh,'$1/\sqrt{n}$', 'Interpreter', 'latex');
    end
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
    ylim([0 0.3])

    
    subplot(3,3,3+s)
    for j=1:length(zeta0)
        scatter(L0,mean_pos(s,j,:)/b,80,color_plt(j,:),'filled');
        hold on;
    end
    ylabel('$b_{pos}/b_{set}$', 'Interpreter', 'latex');set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
    ylim([0.9 4])
    set(gca,'XTicklabel', []);
    axes('Position',[0.24+0.28*(s-1) 0.59 0.1 0.1])
    box on
    for j=1:length(zeta0)
        scatter(L0,std_pos(s,j,:)./mean_pos(s,j,:),80,color_plt(j,:),'filled');
        hold on;
    end
    hh=plot(L0,1./sqrt(L0),'k','linewidth',2);
    set(gca, 'fontsize', 16);
    grid on;box on;
    ylim([0 0.3])


    subplot(3,3,6+s)
    for j=1:length(zeta0)
        hh_all(j)=scatter(L0,mean_kms(s,j,:)/b,80,color_plt(j,:),'filled');
        hold on;
    end
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$b_{KMS}/b_{set}$', 'Interpreter', 'latex');set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
    ylim([0.9 4])
    axes('Position',[0.24+0.28*(s-1) 0.29 0.1 0.1])
    box on
    for j=1:length(zeta0)
        scatter(L0,std_kms(s,j,:)./mean_kms(s,j,:),80,color_plt(j,:),'filled');
        hold on;
    end
    hh=plot(L0,1./sqrt(L0),'k','linewidth',2);
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
    ylim([0 0.3])
end
figure(fig);
forleg = arrayfun(@(x) ['$\zeta=' num2str(x) '$'], zeta0, 'UniformOutput', false);
hL=legend(hh_all, forleg{:}, 'NumColumns', 1, 'Location', 'westoutside', 'Interpreter', 'latex');
newPosition = [0.03 0.42 0.05 0.2]; % left, bottom, width, height
newUnits = 'normalized';
set(hL, 'Position', newPosition);

function [Dm]=gentgr(beta,L,Mmin,M_corner)
        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
end
