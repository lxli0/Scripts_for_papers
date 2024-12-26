%{
%%
clc,clear
beta=1;
b=1;
gamma=[1 0.5];
L0=[100 200 400 800];
mmin_set=0;
mmin_true=mmin_set;
Mmin_set=10^(1.5*(mmin_true+6.07));

rep_times=1000;
size1=length(gamma);
size2=length(L0);
mcorner_TGR_est=zeros(size1,size2,rep_times);

for k=1:size1
    for j=1:size2
        L=L0(j);
        mcorner_TGR=zeros(L,1);
        for i = 1:L
            mcorner_TGR(i)=log10(L)/b-beta-gamma(k)*(L-i)/L;
        end
        mean_mco(k,j)=mean(mcorner_TGR);

        parpool(15);
        parfor jk=1:rep_times
            x = zeros(L, 1);
            for i = 1:L 
                M_corner=10^(1.5*(mcorner_TGR(i)+6.07));                   
                x(i)=gentgr(2/3*b,1,Mmin_set,M_corner);
            end
            M_corner=10^(1.5*(mean_mco(k,j)+6.07));  
            x_mean=gentgr(2/3*b,L,Mmin_set,M_corner);

            [b_TGR_est(k,j,jk),mcorner_TGR_est(k,j,jk),loglikelihood_TGR_rand]=Estimation_TGR_NR_continuous(x,mmin_set);
            [b_TGR_est_mean(k,j,jk),mcorner_TGR_est_mean(k,j,jk),loglikelihood_TGR_rand]=Estimation_TGR_NR_continuous(x_mean,mmin_set);
        end
        delete(gcp('nocreate'))
    end    
end

save('TGR2.mat');
%}
%%
clc,clear
load('TGR2.mat');
load('batlow.mat');
colormap_input=batlow;
original_indices = linspace(1, size(colormap_input, 1), size(colormap_input, 1));
target_indices = linspace(1, size(colormap_input, 1), size1);
color_plt = interp1(original_indices, colormap_input, target_indices);

figure('units','normalized','position',[0.1,0.1,0.3,0.8])
subplot('Position',[0.1 0.7 0.8 0.25]);
mainColors = colormap(lines(size1));

for k=1:size1
    for j=1:size2
        L=L0(j);
        mcorner_TGR=[];
        blendFactor = 0.9*(1 - j/size2);  % Blending towards white, less so for larger indices
        gradientColor = mainColors(k, :) + (1 - mainColors(k, :)) * blendFactor;

        for i = 1:L            
            mcorner_TGR(i)=log10(L)/b-beta-gamma(k)*(L-i)/L;
        end
        h((k-1)*size2+j)=plot((1:1:L)/L,mcorner_TGR,'color',gradientColor,'linewidth',2.5);hold on;
        [~,ind]=min(abs(mcorner_TGR-mean(mcorner_TGR)));
        scatter(ind/L,mean(mcorner_TGR),30,gradientColor,'filled');hold on;
        forleg((k-1)*size2+j) = arrayfun(@(x,y) ['$n=' num2str(x) ', \gamma=' num2str(y) '$'], L,gamma(k), 'UniformOutput', false);
    end
end
grid on;box on; grid minor;
set(gca, 'fontsize', 16);
ylim([-1 3]);
legend(h,forleg, 'NumColumns', 2,'location','southeast', 'Interpreter', 'latex');
ylabel('$m_{corner}$', 'Interpreter', 'latex');
xlabel('$i/n$', 'Interpreter', 'latex');
syms n b gamma i
equation = sprintf('$m_{corner}(i) = \\log_{10}(\\frac{n}{b}) - 1- \\gamma\\frac{(n-i)}{n}$');
title(equation,'Interpreter','latex');


color_plt=[];
for j=1:size2
    for k=1:size1
        values_mco(:,(j-1)*(2+1)+k) = reshape(mcorner_TGR_est(k,j,:),[],1)-mean_mco(k,j);
        values_b(:,(j-1)*(2+1)+k)= reshape(b_TGR_est(k,j,:),[],1);
    end
    values_mco(:,j*(2+1))=nan;
    values_b(:,j*(2+1))=nan;
    c=colormap(lines(size1));
    color_plt=[color_plt;c;ones(1,3)];
end
values_mco=values_mco(:,1:end-1);
values_b=values_b(:,1:end-1);

subplot('Position',[0.1 0.35 0.8 0.25])
h=boxplot(values_b,'color','k');
set(findobj(gca,'type','line'),'linew',2)
h_box=findobj(gca,'Tag','Box');
for j0=1:length(h_box)
    j=length(h_box)-j0+1; % findobj ordered from the most recently added to the oldest
    patch(get(h_box(j0),'XData'),get(h_box(j0),'YData'),color_plt(j,:),'FaceAlpha',0.8,'EdgeAlpha',0.3);
end
hold on;
for ii = 1:size1
    hh(ii)=plot(NaN,1,'color', color_plt(ii,:), 'LineWidth', 20);
end
legend(hh,{'$\gamma=1$', '$\gamma=0.5$'}, 'Interpreter', 'latex');
xticks([1.5 4.5 7.5 10.5])
xticklabels({'100', '200', '400', '800'})
set(gca, 'fontsize', 16);
ylim([0.5 2]);
ylabel('$b_{est}/b_{0}$', 'Interpreter', 'latex');xlabel('$n$', 'Interpreter', 'latex'); 
grid on;box on; grid minor;



subplot('Position',[0.1 0.05 0.8 0.25])
h=boxplot(values_mco,'color','k');
set(findobj(gca,'type','line'),'linew',2)
h_box=findobj(gca,'Tag','Box');
for j0=1:length(h_box)
    j=length(h_box)-j0+1; % findobj ordered from the most recently added to the oldest
    patch(get(h_box(j0),'XData'),get(h_box(j0),'YData'),color_plt(j,:),'FaceAlpha',0.8,'EdgeAlpha',0.3);
end
hold on;
xticks([1.5 4.5 7.5 10.5])
xticklabels({'100', '200', '400', '800'})
set(gca, 'fontsize', 16);
ylim([-0.5 1.5]);
ylabel('$m_{corner}-m_{corner_{0}}$', 'Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex'); 
grid on;box on; grid minor;
%}
%%
clc,clear
load('TGR2.mat');
load('batlow.mat');
colormap_input=batlow;
original_indices = linspace(1, size(colormap_input, 1), size(colormap_input, 1));
target_indices = linspace(1, size(colormap_input, 1), size1);
color_plt = interp1(original_indices, colormap_input, target_indices);

figure('units','normalized','position',[0.1,0.1,0.2,0.5])
mainColors = colormap(lines(size1));

color_plt=[];
for j=1:size2
    for k=1:size1
        values_mco(:,(j-1)*(2+1)+k) = reshape(mcorner_TGR_est_mean(k,j,:),[],1)-mean_mco(k,j);
        values_b(:,(j-1)*(2+1)+k)= reshape(b_TGR_est_mean(k,j,:),[],1);
    end
    values_mco(:,j*(2+1))=nan;
    values_b(:,j*(2+1))=nan;
    c=colormap(lines(size1));
    color_plt=[color_plt;c;ones(1,3)];
end
values_mco=values_mco(:,1:end-1);
values_b=values_b(:,1:end-1);

subplot(2,1,1)
h=boxplot(values_b,'color','k');
set(findobj(gca,'type','line'),'linew',2)
h_box=findobj(gca,'Tag','Box');
for j0=1:length(h_box)
    j=length(h_box)-j0+1; % findobj ordered from the most recently added to the oldest
    patch(get(h_box(j0),'XData'),get(h_box(j0),'YData'),color_plt(j,:),'FaceAlpha',0.8,'EdgeAlpha',0.3);
end
hold on;
for ii = 1:size1
    hh(ii)=plot(NaN,1,'color', color_plt(ii,:), 'LineWidth', 20);
end
legend(hh,{'$\lambda=1$', '$\lambda=0.5$'}, 'Interpreter', 'latex');
xticks([1.5 4.5 7.5 10.5])
xticklabels({'100', '200', '400', '800'})
set(gca, 'fontsize', 16);
ylim([0.5 2]);
ylabel('$b_{est}/b_{0}$', 'Interpreter', 'latex');xlabel('$n$', 'Interpreter', 'latex'); 
grid on;box on; grid minor;


syms n b lambda i
equation = sprintf('$m_{corner} = \\frac{\\sum{\\log_{10}(\\frac{n}{b}) - 1- \\gamma\\frac{(n-i)}{n}}}{n}$');

sgtitle(equation,'Interpreter','latex');

subplot(2,1,2)
h=boxplot(values_mco,'color','k');
set(findobj(gca,'type','line'),'linew',2)
h_box=findobj(gca,'Tag','Box');
for j0=1:length(h_box)
    j=length(h_box)-j0+1; % findobj ordered from the most recently added to the oldest
    patch(get(h_box(j0),'XData'),get(h_box(j0),'YData'),color_plt(j,:),'FaceAlpha',0.8,'EdgeAlpha',0.3);
end
hold on;
xticks([1.5 4.5 7.5 10.5])
xticklabels({'100', '200', '400', '800'})
set(gca, 'fontsize', 16);
ylim([-0.5 1.5]);
ylabel('$m_{corner}-m_{corner_{0}}$', 'Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex'); 
grid on;box on; grid minor;
