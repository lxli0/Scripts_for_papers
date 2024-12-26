%{
clc,clear
b0=[0.6 1 1.4];
L0=[1000 500 200];
zeta0=[0.8 0.7 0.6];
correct_two=zeros(1,length(b0));
correct_LRT=zeros(1,length(b0));
correct_MDMC=zeros(1,length(b0));
correct_zero=zeros(1,length(b0));
times=1000;

for k=1:length(b0)
    b=b0(k);
    L=L0(k);
    zeta=zeta0(k);
    parpool(10);
    parfor i=1:times
        %% generate catalogs
        beta=2/3*b;
        mmin=0;
        Mmin=10^(1.5*(mmin+6.07));
        mcorner_TGR=zeta*log10(L)/b;
        M_corner=10^(1.5*(mcorner_TGR+6.07));

        U10=rand(1,L);
        U20=rand(1,L);
        U1=Mmin-M_corner*log(U10);
        U2=Mmin*U20.^(-1/beta);
        DM=min(U1,U2);
        Dm=2/3*log10(DM)-6.07;    
        
        dDm=diff(Dm);
        dDm_pos=dDm(dDm>=0);

        b_pos=1/log(10)/mean(dDm_pos);   

        %% Estimation
        [b_GR,loglikelihood_GR]=Estimation_GR(Dm,mmin);
        [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR(Dm,mmin);

        %% Judgement % Result: 1 means GR, 2 means TGR
        % Likelihood-ratio test
        R=loglikelihood_TGR-loglikelihood_GR;
        if 2*R>3.84
            Target_LRT(i)=2;
        else
            Target_LRT(i)=1;
        end

        % MDMC
        dmmax=max(dDm_pos);
        F_dm_pos=1-10.^(-b_pos*dmmax);
        P_dm_pos=F_dm_pos^(length(dDm_pos));
        
        if P_dm_pos>0.05
            Target_MDMC(i)=1;
        else
            Target_MDMC(i)=2;
        end
        
        if Target_LRT(i)==2&&Target_MDMC(i)==2
            Target_total(i)=1;
        elseif Target_LRT(i)==2&&Target_MDMC(i)==1
            Target_total(i)=2;
        elseif Target_LRT(i)==1&&Target_MDMC(i)==2
             Target_total(i)=3;
        elseif Target_LRT(i)==1&&Target_MDMC(i)==1
              Target_total(i)=4;
        end
    end
        delete(gcp('nocreate'))
        correct_two(k)=length(find(Target_total==1))/times;
        correct_LRT(k)=length(find(Target_total==2))/times;
        correct_MDMC(k)=length(find(Target_total==3))/times;
        correct_zero(k)=length(find(Target_total==4))/times;
end
save('compare_MDMC_LRT_TGR.mat')
%}

figure('units','normalized','position',[0.1,0.1,0.4,0.3])
color_map=flip(autumn(4));
Correction=[correct_two',correct_LRT',correct_MDMC',correct_zero'];
pl = bar(1:1:length(b0), Correction,1);
XTickLabel = cell(1, numel(b0));
for i = 1:numel(b0)
    XTickLabel{i} = sprintf('$b=%.1f,\\, n=%d,\\, \\zeta=%.1f$', b0(i), L0(i), zeta0(i));
end
set(gca, 'XTickLabel', XTickLabel, 'TickLabelInterpreter', 'latex'); % Setting the custom tick labels with LaTeX
ylabel('Frequency')
for i=1:4
    set(pl(i),'FaceColor',color_map(i,:));
end
ylim([0 1]);
set(gca,'fontsize',16);
grid on; box on;
legend('Both methods are correct','Only LRT is correct','Only MMDC is correct','Both methods are wrong','Numcolumns',2)