%{
clc,clear
b0=0.6:0.1:1.5;
L0=[50 100 150 200 300 400 600 800 1e3];
mmin_set=0;
mmin_true=mmin_set;

rep_times=1000;
size1=length(b0);
size2=length(L0);
b_GR_est=zeros(rep_times,size2,size1);
b_TGR_est=zeros(rep_times,size2,size1);
mcorner_TGR_est=zeros(rep_times,size2,size1);
mean_b_GR=zeros(size1,size2);
std_b_GR=zeros(size1,size2);
mean_b_TGR=zeros(size1,size2);
std_b_TGR=zeros(size1,size2);
mean_mcorner=zeros(size1,size2);
std_mcorner=zeros(size1,size2);
parpool(10);
parfor k=1:size1
    b=b0(k);
    for j=1:size2
        j
        L=L0(j);  
        for i=1:rep_times
            Dm_GR=-1/b.*log10(rand(1,L))+mmin_true;  
            [b_GR_est(i,j,k),loglikelihood_GR_rand]=Estimation_GR_continuous(Dm_GR,mmin_set);
            [b_TGR_est(i,j,k),mcorner_TGR_est(i,j,k),loglikelihood_TGR_rand]=Estimation_TGR_NR_continuous(Dm_GR,mmin_set);
        end
    end    
end
delete(gcp('nocreate'))

for k=1:size1
    for j=1:size2
        mean_b_GR(k,j)=mean(b_GR_est(:,j,k));
        std_b_GR(k,j)=std(b_GR_est(:,j,k));
        mean_b_TGR(k,j)=mean(b_TGR_est(:,j,k));
        std_b_TGR(k,j)=std(b_TGR_est(:,j,k));
        mean_mcorner(k,j)=mean(mcorner_TGR_est(:,j,k));
        std_mcorner(k,j)=std(mcorner_TGR_est(:,j,k));
    end
end

save('TGR_on_GR.mat');
%}
clc,clear
close all
load('TGR_on_GR.mat');

for k=1:size1
    x=L0;
    y=mean_b_TGR(k,:)/b0(k);
    fittingFunction = fittype('1-a * x^(-1)', 'independent', 'x', 'coefficients', {'a'});
    [fitResult, gof] = fit(x',y', fittingFunction);    
    fita(k)=fitResult.a;
end
fittingFunction = fittype('a * x+b', 'independent', 'x', 'coefficients', {'a','b'});
[fitResult, gof] = fit(b0',fita', fittingFunction);
p1=fitResult.a
p2=fitResult.b

load('batlow.mat');
colormap_input=batlow;
original_indices = linspace(1, size(colormap_input, 1), size(colormap_input, 1));
target_indices = linspace(1, size(colormap_input, 1), size1);
color_plt = interp1(original_indices, colormap_input, target_indices);


figure('units','normalized','position',[0.1,0.1,0.98,0.9])
for k=1:size1
    subplot(2,2,1)
    x=L0;
    y=mean_b_TGR(k,:)/b0(k);
    plot(L0, 1- (3.38*b0(k)-1.20)./L0, 'Color', color_plt(k,:), 'LineWidth', 2);
    hold on;
    h(k)=scatter(L0, mean_b_TGR(k,:)/b0(k),30,color_plt(k,:),'filled');
    forleg{k} = ['\textit{b=}', num2str(b0(k))];
    hold on;
    hold on;
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$b_{TGR}/b_{set}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;

    subplot(2,2,2)
    scatter(L0, std_b_TGR(k,:)./mean_b_TGR(k,:),30,color_plt(k,:),'filled');
    hold on;
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$\sigma_{b_{TGR}}/b_{TGR}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
end
subplot(2,2,1)
legend(h,forleg, 'NumColumns', 2,'location','southeast', 'Interpreter', 'latex');

subplot(2,2,2)
hh=plot(L0,1./sqrt(L0),'-k','linewidth',2);
legend(hh,'$1/\sqrt{n}$', 'Interpreter', 'latex');

% Define the offset and size for the smaller subplots
offsetX = 0.13; % X offset for the bottom left corner
offsetY = 0.1; % Y offset for the bottom left corner
width = 0.775 / (size1/2); % Each subplot's width (half of the figure width divided by 4)
height = 1 / (size1/2); % Each subplot's height (half of the figure height divided by 4)

% Smaller subplots in the bottom left, without titles and space
for j = 1:size1
    if j<=(size1/2)
        i=j+(size1/2);
    else
        i=j-(size1/2);
    end
    col = mod(i-1, (size1/2));
    row = floor((i-1) / (size1/2));
    
    left = offsetX + col * width;
    bottom = offsetY + row * height;
    
    subplot('Position', [left bottom width height]);
    values=zeros(rep_times,size2);
    for k=1:size2
    values(:,k) = reshape(mcorner_TGR_est(:,k,j)-log10(L0(k))/b0(j),[],1);
    end

    boxplot(values, 'Labels',string(L0),'Colors', color_plt(j,:)); 
    set(findobj(gca,'type','line'),'linew',2)
    h = findobj(gca,'Tag','Box');
    for jkf=1:length(h)
        patch(get(h(jkf),'XData'),get(h(jkf),'YData'),color_plt(j,:),'FaceAlpha',.5,'EdgeAlpha',0.3); 
    end
  
    if mod(j,(size1/2))~=1
        set(gca, 'YTickLabel', []);
    else
        ylabel('$m_{corner}-a/b$', 'Interpreter', 'latex');
    end
    if j<=(size1/2)
        set(gca, 'XTickLabel', []);
    else
        xlabel('$n$', 'Interpreter', 'latex');            
    end
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
    yticks(-1:1:2);
    ylim([-2 3]);
end
%}