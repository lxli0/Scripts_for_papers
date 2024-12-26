%{
clc,clear
b0=[0.6 1 1.5];
L0=[100 200 500 1000 2000];
zeta0=[0.3 0.5 0.6 0.7 0.8 0.9];
correction_AIC=zeros(length(b0),length(L0),length(zeta0));
correction_BIC=zeros(length(b0),length(L0),length(zeta0));
correction_LRT=zeros(length(b0),length(L0),length(zeta0));
correction_CC=zeros(length(b0),length(L0),length(zeta0));
correction_KS=zeros(length(b0),length(L0),length(zeta0));
correction_LT=zeros(length(b0),length(L0),length(zeta0));
correction_MMC=zeros(length(b0),length(L0),length(zeta0));
times=1000;

for k=1:length(b0)
    b=b0(k);   
    for j=1:length(L0)
        L=L0(j);
        for jk=1:length(zeta0)
            zeta=zeta0(jk);
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

            correction_AIC(k,j,jk)=length(find(Target_AIC==2))/times;
            correction_BIC(k,j,jk)=length(find(Target_BIC==2))/times;
            correction_LRT(k,j,jk)=length(find(Target_LRT==2))/times;
            correction_CC(k,j,jk)=length(find(Target_CC==2))/times;
            correction_KS(k,j,jk)=length(find(Target_KS==2))/times;
            correction_LT(k,j,jk)=length(find(Target_LT==2))/times;
            correction_MMC(k,j,jk)=length(find(Target_MMC==2))/times;
        end
    end
end
save('catalogsize_TGR.mat');
%}
%
%%
clc,clear
load('catalogsize_TGR.mat')
figure('units','normalized','position',[0.1,0.1,0.8,0.9])
color=1/255*[231 98 84;239 133 71;247 170 88;170 220 224;114 188 213;82 143 173;55 103 149];
color_MFD=flip(summer(length(L0)));
mainColors = [
    0.1, 0.1, 0.5;   
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

        axesHandles(1, k) =subplot(1+length(zeta0),length(b0),k) ;
        for jk=1:length(zeta0)
            zeta=zeta0(jk);            
            blendFactor = 0.8*(1 - jk / length(zeta0));  % Blending towards white, less so for larger indices
            gradientColor = mainColors(j, :) + (1 - mainColors(j, :)) * blendFactor;

            m_range=0:0.1:log10(L)/b+1;
            M_range=10.^(1.5*(m_range+6.07));
            mmin_true=0;
            Mmin_set=10^(1.5*(mmin_true+6.07));

            mcorner_TGR=zeta*log10(L)/b;
            M_corner=10^(1.5*(mcorner_TGR+6.07));  
            P_TGR=(Mmin_set./M_range).^beta.*exp((Mmin_set-M_range)/M_corner);
        
            pl_MFD(j)=semilogy(m_range,L*P_TGR, 'Color', gradientColor,'linewidth',1);
            hold on;
        end
        ylim([1 max(L0)]);
        xlim([0 log10(L)/min(b0)]);
        set(gca,'fontsize',16);
        grid on;box on;grid minor;
        xlabel('$m$', 'Interpreter', 'latex');
        if k==1
            ylabel('$N(\geq m)$', 'Interpreter', 'latex');
        else
            set(gca, 'YTickLabel', []) 
        end
    end
end

for k=1:length(b0)
    b=b0(k);
    beta=2/3*b;
       
    for j=1:length(zeta0)
        zeta=zeta0(j);

      %  fortitle = ['\textbf{\textit{b}=', num2str(b), '}'];
      %  title(fortitle, 'Interpreter', 'latex');
            
        
        axesHandles(j+1, k) = subplot(1+length(zeta0),length(b0),length(b0)+3*(j-1)+k);
        %subplot(1+length(L0),length(b0),length(b0)+3*(j-1)+k)   
        Correction=[
            reshape(correction_AIC(k,:,j),[],1),reshape(correction_BIC(k,:,j),[],1),reshape(correction_LRT(k,:,j),[],1), ...
            reshape(correction_CC(k,:,j),[],1),reshape(correction_KS(k,:,j),[],1),reshape(correction_LT(k,:,j),[],1),reshape(correction_MMC(k,:,j),[],1)];
        pl = bar(1:1:length(L0), Correction,1);
        
        for i=1:7
            set(pl(i),'FaceColor',color(i,:));
        end
        ylim([0 1]);
        set(gca,'fontsize',16);
        
        if (length(b0)+3*(j-1)+k)>length(zeta0)*length(b0)
            set(gca,'XTickLabel',num2str(L0'));            
            xlabel('$n$', 'Interpreter', 'latex');
        else
            set(gca, 'XTickLabel', [])
        end
        
        if mod(length(b0)+3*(j-1)+k,length(b0))==1
            ylabel('$R$', 'Interpreter', 'latex');
        else
            set(gca, 'YTickLabel', [])
        end
       grid on;box on;grid minor;
    end
end
%{
subplot(1+length(L0),length(b0),1)
forleg = arrayfun(@(x) ['$n=' num2str(x) '$'], L0, 'UniformOutput', false);
legend(forleg, 'NumColumns', 1, 'Location', 'Northeast', 'Interpreter', 'latex');
%}
subplot(1+length(zeta0),length(b0),1)
forleg = arrayfun(@(x) ['$n=' num2str(x) '$'], L0, 'UniformOutput', false);
legend(pl_MFD,forleg, 'NumColumns', 2, 'Location', 'Southeast', 'Interpreter', 'latex');

subplot(1+length(zeta0),length(b0),1+length(b0))
hL=legend(["Akaike information criterion","Bayesian information criterion","likelihood-ratio test","curve curvature test", ...
    "Kolmogorov–Smirnov test", "linearity test","maximum magnitude criterion"],'Numcolumns',4,'Location','Southeast');


for i = 1:1 + length(zeta0)
    for j = 1:length(b0)
        ax = axesHandles(i, j);
        pos = get(ax, 'OuterPosition');
        if i >=2
            newOuterPos = [pos(1), pos(2)+0.02*(i-2)-0.06 , pos(3) + 0.05, pos(4) ];
        elseif i==1
            newOuterPos = [pos(1), pos(2)-0.04 , pos(3) + 0.05, pos(4)+0.06 ];
        end
        set(ax, 'OuterPosition', newOuterPos);
    end
end

newPosition = [0.5 0.01 0.1 0.06]; % left, bottom, width, height
newUnits = 'normalized';
set(hL, 'Position', newPosition);

%}