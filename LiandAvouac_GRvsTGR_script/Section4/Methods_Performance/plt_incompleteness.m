
figure('units','normalized','position',[0.1,0.1,0.7,0.6])
color=1/255*[231 98 84;239 133 71;247 170 88;170 220 224;114 188 213;82 143 173;55 103 149];
color_MFD=[0.1, 0.1, 0.5;   
    0.58, 0.44, 0.86; 
    0.5, 0, 0; 
    0, 0.5, 0.5;
    0.72, 0.53, 0.04
];

load('incompleteness_GR.mat')
for k=1:length(sigma0)
    sigma=sigma0(k);
    beta=2/3*b;
    subplot(2,2*length(sigma0),k) 
    m_range=0:0.01:log10(L)/b+1;
    M_range=10.^(1.5*(m_range+6.07));
    mmin_true=0;
    Mmin_set=10^(1.5*(mmin_true+6.07));
    pdf_GR=L*b*log(10).*10.^(-b.*(m_range-mmin_true)).*normcdf(m_range, mu, sigma);
    P_GR=cumsum(pdf_GR(1:end-1).*diff(m_range));%./sum(pdf_GR(1:end-1).*diff(m_range));
    P_GR=[P_GR(end),P_GR(end)-P_GR];
    semilogy(m_range,P_GR,'-k','linewidth',2.5);hold on;
    for j=1:length(mc_select0)
        mc_select=mu+mc_select0(j)*sigma;
        jkf=find(abs(m_range-mc_select)<=1e-6);
        pl1(j)=semilogy([mc_select mc_select],[1 P_GR(jkf)],'color',color_MFD(j,:),'linewidth',2.5);
        hold on;
    end
    ylim([1 L]);
    xlim([0 log10(L)/b]);
    set(gca,'fontsize',16);
    grid on;box on;grid minor;
    xlabel('$m$', 'Interpreter', 'latex');


    if k==1
        ylabel('$N(\geq m)$', 'Interpreter', 'latex');
    else
        set(gca, 'YTickLabel', []) 
    end

    subplot(2,2*length(sigma0),k+2*length(sigma0)) 
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

subplot(2,2*length(sigma0),2+2*length(sigma0))
hL=legend(["Akaike information criterion","Bayesian information criterion","likelihood-ratio test","curve curvature test", ...
    "Kolmogorov–Smirnov test", "linearity test","maximum magnitude criterion"],'Numcolumns',4,'Location','Southeast');

    %{
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
%}
newPosition = [0.5 0.03 0.1 0.06]; % left, bottom, width, height
newUnits = 'normalized';
set(hL, 'Position', newPosition);

clc,clear
load('incompleteness_TGR.mat')
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
    
    subplot(2,2*length(sigma0),k+length(sigma0)) 
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
    set(gca, 'YTickLabel', []) 

    subplot(2,2*length(sigma0),k+3*length(sigma0)) 
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
    set(gca, 'YTickLabel', []) 
end

subplot(2,2*length(sigma0),2*length(sigma0))
forleg = {'$m_C=\mu-2\sigma$','$m_C=\mu-\sigma$','$m_C=\mu$','$m_C=\mu+\sigma$','$m_C=\mu+2\sigma$'};
legend(pl1, forleg, 'NumColumns', 1, 'Location', 'Northeast', 'Interpreter', 'latex');

subplotPositions = findobj(gcf, 'Type', 'axes');
for i = 1:length(sigma0)*4
    pos = get(subplotPositions(i), 'OuterPosition');
    if mod(i,2)==0
        newOuterPos = [pos(1), pos(2)+0.02 , pos(3) + 0.04, pos(4)-0.03 ];
    else
        newOuterPos = [pos(1), pos(2)+0.07 , pos(3) + 0.04, pos(4)-0.03 ];
    end
    set(subplotPositions(i), 'OuterPosition', newOuterPos);
end
%}
