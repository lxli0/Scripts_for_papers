%{
clc,clear
L0=[100 200 500 1e3];
b0=[0.7 1.0 1.3];
mmin_set=0;
mmin_true=mmin_set;

for j=1:length(b0)
    j
    b=b0(j);
    for k=1:length(L0)
        L=L0(k);
        parpool(8);
        parfor i=1:1e3
            Dm_GR=-1/b.*log10(rand(1,L))+mmin_true;  
            tic;
            [b_TGR_est_GS(i),mcorner_TGR_est_GS(i),loglikelihood_TGR_rand]=Estimation_TGR_gridsearch_continuous(Dm_GR,mmin_set);
            t_GS(i)=toc;
            tic;
            [b_TGR_est_NR(i),mcorner_TGR_est_NR(i),loglikelihood_TGR_rand]=Estimation_TGR_NR_continuous(Dm_GR,mmin_set);
            t_NR(i)=toc;
        end
        delete(gcp('nocreate'))
        mean_T_GS(j,k)=mean(t_GS);
        mean_T_NR(j,k)=mean(t_NR);
        mean_b_GS(j,k)=mean(b_TGR_est_GS);
        std_b_GS(j,k)=std(b_TGR_est_GS);
        mean_b_NR(j,k)=mean(b_TGR_est_NR);
        std_b_NR(j,k)=std(b_TGR_est_NR);
        mean_mcorner_GS(j,k)=mean(mcorner_TGR_est_GS);
        std_mcorner_GS(j,k)=std(mcorner_TGR_est_GS);
        mean_mcorner_NR(j,k)=mean(mcorner_TGR_est_NR);
        std_mcorner_NR(j,k)=std(mcorner_TGR_est_NR);
    end
end
%}
figure('units','normalized','position',[0.1,0.1,0.4,0.6])
subplot(3,2,1)
for j=1:length(b0)
    scatter(L0,mean_b_NR(j,:)./mean_b_GS(j,:),80,'filled');
    hold on;
end
forleg = arrayfun(@(x) ['$b=' num2str(x) '$'],b0', 'UniformOutput', false);
legend(forleg,'Numcolumns',2,'Location','South', 'Interpreter', 'latex');
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$b_{NR}/b_{GS}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.95 1.03])
subplot(3,2,2)
for j=1:length(b0)
    scatter(L0,std_b_NR(j,:)./std_b_GS(j,:),80,'filled');
    hold on;
end
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$\sigma_{b_{NR}}/\sigma_{b_{GS}}$', 'Interpreter', 'latex');set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.95 1.03])


subplot(3,2,3)
for j=1:length(b0)
    scatter(L0,mean_mcorner_NR(j,:)./mean_mcorner_GS(j,:),80,'filled');
    hold on;
end
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$m_{{corner}_{NR}}/m_{{corner}_{GS}}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.95 1.03])

subplot(3,2,4)
for j=1:length(b0)
    scatter(L0,std_mcorner_NR(j,:)./std_mcorner_GS(j,:),80,'filled');
    hold on;
end
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$\sigma_{m_{{corner}_{NR}}}/\sigma_{m_{{corner}_{GS}}}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.95 1.03])


subplot(3,2,5)
for j=1:length(b0)
    scatter(L0,mean_T_NR(j,:),80,'filled');
    hold on;
end
grid on;box on;
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$t_{NR}$ (s)', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;

subplot(3,2,6)
for j=1:length(b0)
    scatter(L0,mean_T_GS(j,:),80,'filled');
    hold on;
end
grid on;box on;
xlabel('$n$', 'Interpreter', 'latex');
ylabel('$t_{GR}$ (s)', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
