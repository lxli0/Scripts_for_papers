figure('units','normalized','position',[0.1,0.1,0.45,0.65]);
%% linear change
clc,clear
L0=[50 100 150 200 300 400 600 800 1e3];
do_time=2e4;
alpha0=[0.05 0.1 0.2 0.4];
for j=1:length(L0)
    L=L0(j);
    for k=1:length(alpha0)
        alpha=alpha0(k);
        b=zeros(L,1);
        for i=1:L
            b(i)=1+alpha*(2*i-L-1)/2/(L-1);
        end
        b_true(j,k)=mean(b);
        b_MLE=zeros(1,do_time);
        for jk=1:do_time
            x = zeros(L, 1);
            for i = 1:L
                x(i) = -1/b(i).*log10(rand(1));  
            end
            b_MLE(jk)=1/log(10)/mean(x);
        end
        mean_b(j,k)=mean(b_MLE);
        std_b(j,k)=std(b_MLE);
    end
end
load('batlow.mat');
colormap_input=batlow;
original_indices = linspace(1, size(colormap_input, 1), size(colormap_input, 1));
target_indices = linspace(1, size(colormap_input, 1), length(alpha0));
color_plt = interp1(original_indices, colormap_input, target_indices);
subplot(3,2,1)
for k=1:length(alpha0)
        alpha=alpha0(k);
        L=L0(end);
        b=zeros(L,1);
        for i=1:L
            b(i)=1+alpha*(2*i-L-1)/2/(L-1);
        end
        h(k)=plot((1:1:L)/L,b,'color',color_plt(k,:),'linewidth',2.5);hold on;
end
grid on;box on;grid minor;
xlabel('$i/n$', 'Interpreter', 'latex');
ylabel('$b$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
forleg = arrayfun(@(x) ['$\alpha=' num2str(x) '$'], alpha0, 'UniformOutput', false);
legend(h,forleg, 'NumColumns', 2,'location','southeast', 'Interpreter', 'latex');
ylim([1-max(alpha0)/2, 1+max(alpha0)/2]);
for k=1:length(alpha0)
    subplot(3,2,3)
    scatter(L0, mean_b(:,k)./b_true(:,k),40,color_plt(k,:),'filled');

    hold on;
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$b_{est}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;

    subplot(3,2,5)
    scatter(L0, std_b(:,k)./mean_b(:,k),40,color_plt(k,:),'filled');
    hold on;
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$\sigma/b_{est}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
end
subplot(3,2,3)
ylim([0.95 1.05])
grid on;box on;grid minor;
subplot(3,2,5)
hh=plot(L0,1./sqrt(L0),'-k','linewidth',2);
legend(hh,'$1/\sqrt{n}$', 'Interpreter', 'latex');
grid on;box on;grid minor;


%% step change
clc,clear
L0=[50 100 150 200 300 400 600 800 1e3];
do_time=2e4;
alpha0=[0.05 0.1 0.2 0.4];
for j=1:length(L0)
    L=L0(j);
    for k=1:length(alpha0)
        alpha=alpha0(k);
        b=zeros(L,1);
        b(1:L/2)=1-alpha/2;
        b(L/2+1:L)=1+alpha/2;
        b_true(j,k)=mean(b);
        b_MLE=zeros(1,do_time);
        for jk=1:do_time
            % Generate synthetic data
            x = zeros(L, 1);
            for i = 1:L
                x(i) = -1/b(i).*log10(rand(1));  
            end
            b_MLE(jk)=1/log(10)/mean(x);
        end
        mean_b(j,k)=mean(b_MLE);
        std_b(j,k)=std(b_MLE);
    end
end
load('batlow.mat');
colormap_input=batlow;
original_indices = linspace(1, size(colormap_input, 1), size(colormap_input, 1));
target_indices = linspace(1, size(colormap_input, 1), length(alpha0));
color_plt = interp1(original_indices, colormap_input, target_indices);

subplot(3,2,2)
for k=1:length(alpha0)
        alpha=alpha0(k);
        L=L0(end);
        b=zeros(L,1);
        b(1:L/2)=1-alpha/2;
        b(L/2+1:L)=1+alpha/2;
        h(k)=plot((1:1:L)/L,b,'color',color_plt(k,:),'linewidth',2.5);hold on;
end
grid on;box on;grid minor;
xlabel('$i/n$', 'Interpreter', 'latex');
ylabel('$b$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
ylim([1-max(alpha0)/2, 1+max(alpha0)/2]);
for k=1:length(alpha0)
    subplot(3,2,4)
    scatter(L0, mean_b(:,k)./b_true(:,k),40,color_plt(k,:),'filled');

    hold on;
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$b_{est}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;

    subplot(3,2,6)
    scatter(L0, std_b(:,k)./mean_b(:,k),40,color_plt(k,:),'filled');
    hold on;
    xlabel('$n$', 'Interpreter', 'latex');
    ylabel('$\sigma/b_{est}$', 'Interpreter', 'latex');
    set(gca, 'fontsize', 16);
    grid on; box on; grid minor;
end
subplot(3,2,4)
ylim([0.95 1.05])
grid on;box on;grid minor;
subplot(3,2,6)
hh=plot(L0,1./sqrt(L0),'-k','linewidth',2);
legend(hh,'$1/\sqrt{n}$', 'Interpreter', 'latex');
grid on;box on;grid minor;


%%
fig2=figure('units','normalized','position',[0.1,0.1,0.45,0.2]);
for j=1:length(L0)
    L=L0(j);   
    for jk=1:do_time        
        x_mean=-log10(rand(L,1));  
        b_MLE_mean(jk)=1/log(10)/mean(x_mean);
    end
    mean_b_mean(j)=mean(b_MLE_mean);
    std_b_mean(j)=std(b_MLE_mean);
end
subplot(1,2,1)
scatter(L0, mean_b_mean/1,40,'filled');
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$b_{est}/b_{0}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
ylim([0.95 1.05])
grid on;box on;grid minor;
subplot(1,2,2)
scatter(L0, std_b_mean/1,40,'filled');
hold on;
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$\sigma/b_{est}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
hh=plot(L0,1./sqrt(L0),'-k','linewidth',2);
legend(hh,'$1/\sqrt{n}$', 'Interpreter', 'latex');