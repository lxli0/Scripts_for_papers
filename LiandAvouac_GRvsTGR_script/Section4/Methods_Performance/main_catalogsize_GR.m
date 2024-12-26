%{
clc,clear
b0=[0.6 1 1.5];
L0=[100 200 500 1000 2000];
correction_AIC=zeros(length(b0),length(L0));
correction_BIC=zeros(length(b0),length(L0));
correction_LRT=zeros(length(b0),length(L0));
correction_CC=zeros(length(b0),length(L0));
correction_KS=zeros(length(b0),length(L0));
correction_LT=zeros(length(b0),length(L0));
correction_MMC=zeros(length(b0),length(L0));
times=1000;

for k=1:length(b0)
    b=b0(k);   
    for j=1:length(L0)
        L=L0(j);
        parpool(20);
        parfor i=1:times
            %% generate catalogs
            Mmin=0;
            Dm=-1/b.*log10(rand(1,L))+Mmin; % generate random magnitude

            %% Estimation
            [b_GR,loglikelihood_GR]=Estimation_GR(Dm,Mmin);
            [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR(Dm,Mmin);

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
        
        correction_AIC(k,j)=length(find(Target_AIC==1))/times;
        correction_BIC(k,j)=length(find(Target_BIC==1))/times;
        correction_LRT(k,j)=length(find(Target_LRT==1))/times;
        correction_CC(k,j)=length(find(Target_CC==1))/times;
        correction_KS(k,j)=length(find(Target_KS==1))/times;
        correction_LT(k,j)=length(find(Target_LT==1))/times;
        correction_MMC(k,j)=length(find(Target_MMC==1))/times;
    end
end
save('catalogsize_GR.mat');
%}
%
%%
clc,clear
load('catalogsize_GR.mat')
figure('units','normalized','position',[0.1,0.1,0.8,0.6])
color=1/255*[231 98 84;239 133 71;247 170 88;170 220 224;114 188 213;82 143 173;55 103 149];
color_MFD=[0.1, 0.1, 0.5;   
    0.58, 0.44, 0.86; 
    0.5, 0, 0; 
    0, 0.5, 0.5;
    0.72, 0.53, 0.04
];

for k=1:length(b0)
    b=b0(k);
    beta=2/3*b;
    for j=1:length(L0)
        L=L0(j);
        subplot(2,length(b0),k) 
        m_range=0:0.02:log10(L)/b+1;
        M_range=10.^(1.5*(m_range+6.07));
        mmin_true=0;
        Mmin_set=10^(1.5*(mmin_true+6.07));

        P_GR=(Mmin_set./M_range).^beta;
        semilogy(m_range,L*P_GR,'color',color_MFD(j,:),'linewidth',2.5);hold on;
        ylim([1 max(L0)]);
        xlim([0 log10(L)/min(b0)]);
        set(gca,'fontsize',16);
        grid on;box on;grid minor;
        xlabel('$m$', 'Interpreter', 'latex');
       % fortitle = ['\textbf{\textit{b}=', num2str(b), '}'];
       % title(fortitle, 'Interpreter', 'latex');
        if k==1
            ylabel('$N(\geq m)$', 'Interpreter', 'latex');
        else
            set(gca, 'YTickLabel', []) 
        end
    end

    subplot(2,length(b0),k+length(b0))
    Correction=[
        correction_AIC(k,:)',correction_BIC(k,:)',correction_LRT(k,:)', ...
        correction_CC(k,:)', correction_KS(k,:)',correction_LT(k,:)',correction_MMC(k,:)'];
    pl = bar(1:1:length(L0), Correction,1);
    for i=1:7
        set(pl(i),'FaceColor',color(i,:));
    end
    set(gca,'XTickLabel',num2str(L0'));
    grid on;box on;grid minor;
    xlabel('$n$', 'Interpreter', 'latex');
    set(gca,'fontsize',16);
    if k==1
        ylabel('$R$', 'Interpreter', 'latex');
    else
        set(gca, 'YTickLabel', []) 
    end
end
subplot(2,length(b0),1)
forleg = arrayfun(@(x) ['$n=' num2str(x) '$'], L0, 'UniformOutput', false);
legend(forleg, 'NumColumns', 1, 'Location', 'Northeast', 'Interpreter', 'latex');


subplot(2,length(b0),3+length(b0))
hL=legend(["Akaike information criterion","Bayesian information criterion","likelihood-ratio test","curve curvature test", ...
    "Kolmogorov–Smirnov test", "linearity test","maximum magnitude criterion"],'Numcolumns',4,'Location','Southeast');


subplotPositions = findobj(gcf, 'Type', 'axes');
for i = 1:2*length(b0)
    pos = get(subplotPositions(i), 'OuterPosition');
    if mod(i,2)==1
        newOuterPos = [pos(1), pos(2)+0.05 , pos(3) + 0.05, pos(4)+0.01 ];
    else
        newOuterPos = [pos(1), pos(2)+0.01 , pos(3) + 0.05, pos(4)+0.01 ];
    end
    set(subplotPositions(i), 'OuterPosition', newOuterPos);
end

newPosition = [0.5 0.03 0.1 0.06]; % left, bottom, width, height
newUnits = 'normalized';
set(hL, 'Position', newPosition);