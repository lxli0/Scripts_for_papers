
clc,clear
b=1;
L=1e4;
mu=0.8;
sigma0=[0.2 0.4];
mc_select0=[-2 -1 0 1 2];
zeta=0.5;
correction_AIC=zeros(length(sigma0),length(mc_select0));
correction_BIC=zeros(length(sigma0),length(mc_select0));
correction_LRT=zeros(length(sigma0),length(mc_select0));
correction_CC=zeros(length(sigma0),length(mc_select0));
correction_KS=zeros(length(sigma0),length(mc_select0));
correction_LT=zeros(length(sigma0),length(mc_select0));
correction_MMC=zeros(length(sigma0),length(mc_select0));
times=1000;

for k=1:length(sigma0)
    sigma=sigma0(k);   
    for j=1:length(mc_select0)
        mc_select=mu+mc_select0(j)*sigma;
        parpool(20);
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

            cri=normcdf(Dm,mu,sigma);
            ran_cri=rand(1,L);
            jkf=ran_cri<=cri;
            Dm=Dm(jkf);
            Dm=Dm(Dm>=mc_select);
            Dm=Dm-mc_select;
            
            %% Estimation
            [b_GR,loglikelihood_GR]=Estimation_GR(Dm,mmin);
            [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR(Dm,mmin);

            %% Judgement % Result: 1 means GR, 2 means TGR
            Moment_max=10^(1.5*(max(Dm)+6.07));
            Moment_min=10^(1.5*6.07);
            Moment_corner=10^(1.5*(mcorner_TGR+6.07));
            n=length(Dm);

            % AIC
            AIC_GR=-2*loglikelihood_GR+2*1; % one parameter
            AIC_TGR=-2*loglikelihood_TGR+2*2; % two parameter
            dAIC=AIC_GR-AIC_TGR;
            if dAIC<0
                Target_AIC(i)=1;
            else
                Target_AIC(i)=2;
            end

            % BIC
            BIC_GR=-2*loglikelihood_GR+1*log(n); % one parameter
            BIC_TGR=-2*loglikelihood_TGR+2*log(n); % two parameter
            dBIC=BIC_GR-BIC_TGR;
            if dBIC<0
                Target_BIC(i)=1;
            else
                Target_BIC(i)=2;
            end

            % Likelihood-ratio test
            R=loglikelihood_TGR-loglikelihood_GR;
            if 2*R>3.84
                Target_LRT(i)=2;
            else
                Target_LRT(i)=1;
            end

            % Curve curvature 
            p_CC=curvature_judgement(Dm,b_GR);
            if p_CC>0.05
                Target_CC(i)=1;
            else
                Target_CC(i)=2;
            end

            % KS test
            p_KS=KStest(Dm,b_GR);
            if p_KS>0.05
                Target_KS(i)=1;
            else
                Target_KS(i)=2;
            end

            % Linearity test 
            Target_LT(i)=Linearity_test(Dm);

            % MMC
            p_MMC_GR=(1-(Moment_min/Moment_max)^(2/3*b_GR))^n;
            P_MMC_TGR=(1-(((Moment_min/Moment_max)^(2/3*b_TGR))*exp((Moment_min-Moment_max)/Moment_corner)))^n;
            if p_MMC_GR>0.05
                Target_MMC(i)=1;
            else
                Target_MMC(i)=2;
            end
        end
        delete(gcp('nocreate'))

        correction_AIC(k,j)=length(find(Target_AIC==2))/times;
        correction_BIC(k,j)=length(find(Target_BIC==2))/times;
        correction_LRT(k,j)=length(find(Target_LRT==2))/times;
        correction_CC(k,j)=length(find(Target_CC==2))/times;
        correction_KS(k,j)=length(find(Target_KS==2))/times;
        correction_LT(k,j)=length(find(Target_LT==2))/times;
        correction_MMC(k,j)=length(find(Target_MMC==2))/times;
    end
end
save('incompleteness_TGR.mat')
%}
%
%%
clc,clear
load('incompleteness_TGR.mat')
figure('units','normalized','position',[0.1,0.1,0.55,0.6])
color=1/255*[231 98 84;239 133 71;247 170 88;170 220 224;114 188 213;82 143 173;55 103 149];
color_MFD=[0.1, 0.1, 0.5;   
    0.58, 0.44, 0.86; 
    0.5, 0, 0; 
    0, 0.5, 0.5;
    0.72, 0.53, 0.04
];

for k=1:length(sigma0)
    sigma=sigma0(k);
    beta=2/3*b;
    subplot(2,length(sigma0),k) 
    m_range=0:0.01:log10(L)/b+1;
    M_range=10.^(1.5*(m_range+6.07));
    mmin_true=0;
    Mmin_set=10^(1.5*(mmin_true+6.07));
    mcorner_TGR=zeta*log10(L)/b;
    M_corner=10^(1.5*(mcorner_TGR+6.07));
            
    CDF_TGR=1-(Mmin_set./M_range).^beta.*exp((Mmin_set-M_range)/M_corner);
    pdf_TGR=diff(CDF_TGR)./diff(m_range);
    pdf_TGR=[pdf_TGR(1),pdf_TGR];
    pdf_TGR = L*pdf_TGR.*normcdf(m_range, mu, sigma);
    P_TGR=cumsum(pdf_TGR(1:end-1).*diff(m_range));
    P_TGR=[P_TGR(end),P_TGR(end)-P_TGR];
    semilogy(m_range,P_TGR,'-k','linewidth',2.5);hold on;
    for j=1:length(mc_select0)
        mc_select=mu+mc_select0(j)*sigma;
        jkf=find(abs(m_range-mc_select)<=1e-6);
        pl1(j)=semilogy([mc_select mc_select],[1 P_TGR(jkf)],'color',color_MFD(j,:),'linewidth',2.5);
        hold on;
    end
    ylim([1 L]);
    xlim([0 log10(L)/b]);
    set(gca,'fontsize',16);
    grid on;box on;grid minor;
    xlabel('$m$', 'Interpreter', 'latex');
    fortitle = ['\sigma=', num2str(sigma)];
    title(fortitle);
    if k==1
        ylabel('$N(\geq m)$', 'Interpreter', 'latex');
    else
        set(gca, 'YTickLabel', []) 
    end

    subplot(2,length(sigma0),k+length(sigma0)) 
    Correction=[
        correction_AIC(k,:)',correction_BIC(k,:)',correction_LRT(k,:)', ...
        correction_CC(k,:)', correction_KS(k,:)',correction_LT(k,:)',correction_MMC(k,:)'];
    pl = bar(1:1:length(mc_select0), Correction,1);
    for i=1:7
        set(pl(i),'FaceColor',color(i,:));
    end
    set(gca,'XTickLabel',num2str(mu+sigma*mc_select0'));
    grid on;box on;grid minor;
    xlabel('$m_C$', 'Interpreter', 'latex');
    set(gca,'fontsize',16);
    if k==1
            ylabel('$R$', 'Interpreter', 'latex');
    else
        set(gca, 'YTickLabel', []) 
    end
end
subplot(2,length(sigma0),1)
forleg = {'$m_C=\mu-2\sigma$','$m_C=\mu-\sigma$','$m_C=\mu$','$m_C=\mu+\sigma$','$m_C=\mu+2\sigma$'};
legend(pl1, forleg, 'NumColumns', 1, 'Location', 'Northeast', 'Interpreter', 'latex');


subplot(2,length(sigma0),2+length(sigma0))
hL=legend(["Akaike information criterion","Bayesian information criterion","likelihood-ratio test","curve curvature test", ...
    "Kolmogorov–Smirnov test", "linearity test","maximum magnitude criterion"],'Numcolumns',3,'Location','Southeast');
  


subplotPositions = findobj(gcf, 'Type', 'axes');
for i = 1:length(sigma0)*2
    pos = get(subplotPositions(i), 'OuterPosition');
    if mod(i,2)==0
        newOuterPos = [pos(1), pos(2)+0.03 , pos(3) + 0.05, pos(4)-0.03 ];
    else
        newOuterPos = [pos(1), pos(2)+0.09 , pos(3) + 0.05, pos(4)-0.03 ];
    end
    set(subplotPositions(i), 'OuterPosition', newOuterPos);
end

newPosition = [0.5 0.03 0.1 0.06]; % left, bottom, width, height
newUnits = 'normalized';
set(hL, 'Position', newPosition);