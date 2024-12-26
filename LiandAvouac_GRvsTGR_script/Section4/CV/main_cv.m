%{
clc,clear
L0=[100 200 500 1000];
times=1000;
Bingo_matrix=zeros(length(L0),times);
for j=1:length(L0)
    L=L0(j);
    parpool(16);
    parfor i=1:times
        %% Generate Catalog 
        b=1; % set b-value
        Pcm=rand(1,L);
        Mmin_true=0;
        Dm0=-1/b.*log10(Pcm)+Mmin_true; % generate random magnitude
        Mmin=Mmin_true; % I assume the Mc is already known
     i   
        for k=1:100
            tk=randi([1, length(Dm0)], 1, length(Dm0));
            Dm=Dm0(tk);
           
            % Estimation for TGR
            [b_TGR(i,k),Mcorner_TGR(i,k),loglikelihood_TGR(i,k)]=Estimation_TGR_continuous(Dm,Mmin);
        end
    end
    delete(gcp('nocreate'))
    for i=1:times
        eta(i,j)=mean(1./(Mcorner_TGR(i,:)));
        std_eta(i,j)=std((1./(Mcorner_TGR(i,:))));
        CV(i,j)=std_eta(i,j)/eta(i,j);
    end
end

%}
figure('units','normalized','position',[0.1,0.1,0.6,0.6])
color_5=1/255*[117 114 181;91 183 205; 197 86 89;203 180 123];
color_5=flip(color_5);

subplot(2,2,1)
histogram(eta(:,4),'normalization','probability', 'Facecolor',color_5(4,:),'EdgeColor', 'none');
xlabel('\eta');
ylabel('Frequency');
set(gca,'fontsize',16);
grid on;box on;grid minor;

subplot(2,2,2)
histogram(std_eta(:,4),'normalization','probability','Facecolor',color_5(4,:), 'EdgeColor', 'none');
xlabel('\sigma_\eta');
ylabel('Frequency');
set(gca,'fontsize',16);
grid on;box on;grid minor;

subplot(2,2,3)
histogram(CV(:,4),'normalization','probability','Facecolor',color_5(4,:), 'EdgeColor', 'none');
xlabel('$CV$', 'Interpreter', 'latex');
ylabel('Frequency');
set(gca,'fontsize',16);
grid on;box on;grid minor;

subplot(2,2,4)
for i = 1:4
    pd = fitdist(CV(:,i),'Normal'); % Fit a normal distribution to the data
    x_values = linspace(min(CV(:,i)), max(CV(:,i)), 100);
    y_values = cdf(pd, x_values);
    plt(i)=plot(x_values, y_values,'Color',color_5(i,:), 'LineWidth', 2);
    
    % Find the value corresponding to the 97.5th percentile
    value_97_5 = icdf(pd, 0.975);
    line([value_97_5 value_97_5], [0 0.975], 'Color', 'black', 'LineWidth', 2, 'LineStyle', '--');
    hold on;
end
xlabel('$CV$', 'Interpreter', 'latex');
ylabel('CDF');
set(gca,'fontsize',16);
forleg = arrayfun(@(x) ['$n=' num2str(x) '$'], L0, 'UniformOutput', false);
legend(forleg, 'NumColumns', 1, 'Location', 'Northeast', 'Interpreter', 'latex');
legend(plt,forleg)
grid on;box on;grid minor;
%}