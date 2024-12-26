clc,clear
b=1;
L0=[100 500 1e3];
zeta0=[0.1 0.3 0.6 0.9];
mmin_set=0;
mmin_true=mmin_set;
Mmin_set=10^(1.5*(mmin_true+6.07));

for j=1:length(zeta0)
    j
    zeta=zeta0(j);
    for k=1:length(L0)
        L=L0(k);
        mcorner_TGR=zeta*log10(L)/b;
        M_corner=10^(1.5*(mcorner_TGR+6.07));
        parpool(8);
        parfor i=1:1e3
            Dm_TGR=gentgr(2/3*b,L,Mmin_set,M_corner);
            tic;
            [b_TGR_est_GS(i),mcorner_TGR_est_GS(i),loglikelihood_TGR_rand]=Estimation_TGR_gridsearch_continuous(Dm_TGR,mmin_set);
            t_GS(i)=toc;
            tic;
            [b_TGR_est_NR(i),mcorner_TGR_est_NR(i),loglikelihood_TGR_rand]=Estimation_TGR_NR_continuous(Dm_TGR,mmin_set);
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
fortitle = ['\textbf{\textit{b}=', num2str(b), '}'];
sgtitle(fortitle, 'Interpreter', 'latex');
subplot(3,2,1)
for k=1:length(L0)
    scatter(zeta0,mean_b_NR(:,k)./mean_b_GS(:,k),80,'filled');
    hold on;
end
forleg = arrayfun(@(x) ['$n=' num2str(x) '$'], L0', 'UniformOutput', false);
legend(forleg,'Numcolumns',2,'Location','North', 'Interpreter', 'latex');
xlabel('$\zeta$', 'Interpreter', 'latex');
ylabel('$b_{NR}/b_{GS}$', 'Interpreter', 'latex');set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.9 1.2])

subplot(3,2,2)
for k=1:length(L0)
    scatter(zeta0,std_b_NR(:,k)./std_b_GS(:,k),80,'filled');
    hold on;
end
grid on;box on;
xlabel('$\zeta$', 'Interpreter', 'latex');
ylabel('$\sigma_{b_{NR}}/\sigma_{b_{GS}}$', 'Interpreter', 'latex');set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.9 1.4])

subplot(3,2,3)
for k=1:length(L0)
    scatter(zeta0,mean_mcorner_NR(:,k)./mean_mcorner_GS(:,k),80,'filled');
    hold on;
end
xlabel('$\zeta$', 'Interpreter', 'latex');
ylabel('$m_{{corner}_{NR}}/m_{{corner}_{GS}}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.9 1.2])

subplot(3,2,4)
for k=1:length(L0)
    scatter(zeta0,std_mcorner_NR(:,k)./std_mcorner_GS(:,k),80,'filled');
    hold on;
end
grid on;box on;
xlabel('$\zeta$', 'Interpreter', 'latex');
ylabel('$\sigma_{m_{{corner}_{NR}}}/\sigma_{m_{{corner}_{GS}}}$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;
ylim([0.9 1.4])

subplot(3,2,5)
for k=1:length(L0)
    scatter(zeta0,mean_T_NR(:,k),80,'filled');
    hold on;
end
grid on;box on;
xlabel('$\zeta$', 'Interpreter', 'latex');
ylabel('$t_{NR}$ (s)', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;

subplot(3,2,6)
for k=1:length(L0)
    scatter(zeta0,mean_T_GS(:,k),80,'filled');
    hold on;
end
grid on;box on;
xlabel('$\xi$', 'Interpreter', 'latex');
ylabel('$t_{GR}$ (s)', 'Interpreter', 'latex');
set(gca, 'fontsize', 16);
grid on; box on; grid minor;

function [Dm]=gentgr(beta,L,Mmin,M_corner)
        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
end