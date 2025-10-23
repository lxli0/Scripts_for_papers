clc; clear;
close all
load('OUTPUT_temporal.mat')
color = flip(turbo(Sample_number));  % Retain original color scheme

%% Categorization
cat_tran_potential = [];

for i = 1:Sample_number
    Temp_plr = allMedian_PLR{i};
    frac_below_005 = sum(Temp_plr < 0.05) / numel(Temp_plr);
    frac_below_010 = sum(Temp_plr < 0.10) / numel(Temp_plr);

    if frac_below_005 >0&&frac_below_010<1
        cat_tran_potential = [cat_tran_potential, i];
    end
end

[~, tran_order] = sort(cellfun(@length, sites(cat_tran_potential)));
cat_tran_sorted = cat_tran_potential(tran_order);
allMedian_PLR = allMedian_PLR(cat_tran_sorted);
allt_complete = allt_complete(cat_tran_sorted);
allm_complete = allm_complete(cat_tran_sorted);
mc = mc(cat_tran_sorted);
global_bGR = global_bGR(cat_tran_sorted);
global_bTGR = global_bTGR(cat_tran_sorted);
global_mcorner_TGR = global_mcorner_TGR(cat_tran_sorted);


legend_labels = sites(cat_tran_sorted);

%% Plot for first site
time_node0 = cell(1, Sample_number);
time_node0{1} = [4.2];
time_node0{2} = [11];
time_node0{3} = [2015.1, 2015.4];
time_node0{4} = [16 31];
time_node0{5} = [26.32];
time_node0{6} = [2022.33];
time_node0{7} = [2012.905, 2012.9135];
time_node0{8} = [0.55];
time_node0{9} = [0.66];
time_node0{10} = [2011.14 2011.15];
time_node0{11} = [10];
time_node0{12} = [3.5];

for i = 7:7%1:length(legend_labels)
    site = legend_labels{i};
    file = ['./Data/', site, '.txt'];
    D = load(file);
    T_all = D(:,1);
    M_all = D(:,5);
    % Call plotting function with current site's data
    plt_segmented(global_bGR(i), global_bTGR(i), global_mcorner_TGR(i), ...
                  T_all, M_all, allt_complete{i}, allm_complete{i}, ...
                  allMedian_PLR{i}, mc(i), time_node0{i});
    sgtitle(site);
end

%% Modified plt_segmented with color-aware titles and text
function [] = plt_segmented(global_bGR,global_bTGR,global_mcorner_TGR,T_all, M_all, T, M, pLRT, Mc, time_node0)
    figure('units','normalized','position',[0.1,0.1,0.3,0.7])
    subplot(length(time_node0)+3, 1, 1);

    yyaxis left
    scatter(T, M, 10, 'k', 'filled');
    hold on;
    xlabel('Time');
    ylabel('Magnitude');
    ax = gca; ax.YColor = 'k';
     for k = 1:length(time_node0)
        xline(time_node0(k), '--', 'Color', 'g', 'LineWidth', 3);
     end
    

    yyaxis right
    plot(T, pLRT, '-b', 'LineWidth', 2.5);
    hold on;
    ylabel('$p_{LRT}$', 'Interpreter', 'latex');
    ylim([0 1]);
    ax.YColor = 'k';
    set(gca, 'FontSize', 16);
    grid on; box on; grid minor;

    time_node = [min(T), time_node0, max(T)];
    diffm = diff(M);
    nonzero_elements = diffm(diffm ~= 0);
    delta = min(abs(nonzero_elements)) / 2;

    for i = 1:length(time_node)-1
        jkf = find(T >= time_node(i) & T < time_node(i+1));
        M_sub = M(jkf);
        [b_GR(i), b_TGR(i), mcorner_TGR(i)] = get_duration_MFD(M_sub, Mc, delta);
    end
    b_GR(i+1)=global_bGR;
    b_TGR(i+1)=global_bTGR;
    mcorner_TGR(i+1)=global_mcorner_TGR-Mc;
    

    for i = 1:length(time_node)-1
        jkf_all = find(T_all >= time_node(i) & T_all < time_node(i+1));
        M_all_sub = M_all(jkf_all);
        jkf = find(T >= time_node(i) & T < time_node(i+1));
        M_sub = M(jkf);
        Dm_sub = M_sub - Mc;

        for j = 1:length(time_node)
            loglikelihood_GR = loglike_GR_discrete(Dm_sub, 0, delta, b_GR(j));
            loglikelihood_TGR = loglike_TGR_discrete(Dm_sub, 0, delta, b_TGR(j), mcorner_TGR(j));
            R = loglikelihood_TGR - loglikelihood_GR;
            P_LR = 1 - chi2cdf(2 * R, 1);

            subplot_idx = (j-1)*(length(time_node0)+1) + length(time_node0)+1 + i;
            subplot(length(time_node0)+3, length(time_node0)+1, subplot_idx);

            % Title text
            if i == 1
                time_label = ['start to ', num2str(time_node(i+1))];               
            elseif i == length(time_node)-1
                time_label = [num2str(time_node(i)), ' to end'];
            else
                time_label = [num2str(time_node(i)), ' to ', num2str(time_node(i+1))];
            end
            
            if i==j
                colorflag=1;
            elseif j==length(time_node)
                colorflag=2;
            else
                colorflag=3;
            end

            plt_segmented_MFD(M_all_sub, Dm_sub, Mc, b_GR(j), b_TGR(j), mcorner_TGR(j), P_LR, colorflag);

            
            if i==j
                title(time_label, 'Color', 'magenta');
            elseif j==length(time_node)
                title(time_label, 'Color', 'cyan');
            else
                title(time_label);
            end
   
        end
    end
end

%% Helper functions (unchanged logic, added is_red param in plt_segmented_MFD)
function [b_GR,b_TGR,mcorner_TGR]=get_duration_MFD(Dm_sub,Mc,delta)
    Dm_sub=Dm_sub-Mc;
    Mmin=0;
    [b_GR,~]=Estimation_GR_discrete(Dm_sub,Mmin,delta);
    [b_TGR,mcorner_TGR,~]=Estimation_TGR_gridsearch_discrete(Dm_sub,Mmin,delta,b_GR);   
end

function [loglikelihood_GR]=loglike_GR_discrete(Dm,Mmin,delta,b_GR)
    DMoment=10.^(1.5*(Dm+6.07));
    Moment_min0=10^(1.5*(Mmin+6.07));
    Moment_min=10^(1.5*((Mmin-delta)+6.07));
    jkf=find((DMoment-Moment_min0)>=-1e-16);
    n=length(jkf);
    beta_GR=2/3*b_GR;
    loglikelihood_GR=n*(beta_GR*log(Moment_min)+log(beta_GR))-(1+beta_GR)*sum(log(DMoment(jkf)));
end

function [loglikelihood_TGR]=loglike_TGR_discrete(Dm,Mmin,delta,b_TGR,mcorner_TGR)
    DMoment=10.^(1.5*(Dm+6.07));
    Moment_min0=10^(1.5*(Mmin+6.07));
    Moment_min=10^(1.5*((Mmin-delta)+6.07));
    jkf=find((DMoment-Moment_min0)>=-1e-16);
    n=length(jkf);
    Moment_corner_TGR=10.^(1.5*(mcorner_TGR+6.07)); 
    beta_TGR=2/3*b_TGR;

    loglikelihood_TGR=n*beta_TGR*log(Moment_min)+1/Moment_corner_TGR*(n*Moment_min-sum(DMoment(jkf))) ...
        -beta_TGR*sum(log(DMoment(jkf)))+sum(log(beta_TGR./DMoment(jkf)+1/Moment_corner_TGR));  
end

function []=plt_segmented_MFD(M,Dm,mc,b_GR,b_TGR,Mcorner_TGR_minus_mc,P_LR,colorflag)
    L=length(M);
    dm=0.1;
    m=floor(min(M)/dm)*dm:dm:max(M);
    if max(m)<max(M)
        m=[m,max(M)];
    end
    n0=hist(M,m);   
    for i=1:length(m)
        cn0(i)=length(find((M-m(i))>=-1e-10));
    end
    n=log10(n0);
    cn=log10(cn0);
    color1=[0.75 0.75 0.75; 0.25 0.25 0.25];
    color2=1/255*[132 94 194;178 91 0;0 139 200];
    semilogy(m,n0,'o','color',0.4*[1 1 1],'markersize',6.5,'markerfacecolor',color1(1,:),'MarkerEdgeColor',color1(1,:)); hold on;
    semilogy(m,cn0,'o','color',0.4*[1 1 1],'markersize',6.5,'markerfacecolor',color1(2,:),'MarkerEdgeColor',color1(2,:)); hold on;

    if ~isnan(mc)
        m=0:dm:max(Dm);
        if max(m)<max(Dm)
            m=[m,max(Dm)];
        end
        Moment_min=10^(1.5*(min(m)+6.07));
        Moment_corner=10^(1.5*(Mcorner_TGR_minus_mc+6.07));
        F_GR=(Moment_min./10.^(3/2*(m+6.07))).^(2/3*b_GR);
        F_TGR=(Moment_min./10.^(3/2*(m+6.07))).^(2/3*b_TGR).*exp((Moment_min-10.^(3/2*(m+6.07)))/ Moment_corner);

        GR_con=GR_confidential_interval(b_GR,length(Dm),mc);
        semilogy(m+mc,length(Dm)*F_TGR,'--','color',color2(2,:),'linewidth',3); hold on;
        scatter(mc,length(Dm)*F_GR(1)*1.5,100,1/255*[255,165,0],'v','filled'); hold on;

        xMin = round(min(M)); xMax = ceil(max(M));
        xlim([xMin xMax]);
        xticks(xMin:xMax);
        xticklabels(arrayfun(@num2str, xMin:xMax, 'UniformOutput', false));
        yMin = 0; yMax = ceil(log10(L));
        ylim([10^yMin, 10^yMax]);
        yticks(10.^(yMin:yMax));
        yticklabels(arrayfun(@(y) sprintf('10^{%d}', round(log10(y))), yticks, 'UniformOutput', false));

        ax = gca;
        xTextPos = 0.95*ax.XLim(2);
        yTextPos = 0.9*10^(log10(ax.YLim(2)));

        textStr = sprintf('$b_{\\mathrm{GR}}$: %.2f\\\\\n$b_{\\mathrm{TGR}}$: %.2f\\\\\n$m_{\\mathrm{corner}}$: %.1f\\\\\n$p_{\\mathrm{LRT}}$: %.2f', ...
            b_GR, b_TGR, Mcorner_TGR_minus_mc+mc, P_LR);

        if colorflag==1
            text_color = 'magenta';
        elseif colorflag==2
            text_color = 'cyan';
        else
            text_color = 'black';
        end

        text(xTextPos, yTextPos, textStr, ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'right', ...
            'FontSize', 13, ...
            'Interpreter', 'latex', ...
            'Color', text_color);
    end

    grid on; box off;
    xlabel('\it m');
    ylabel('\it N');
    set(gca,'fontsize',13);
    hold on;
end

function [pl] = GR_confidential_interval(b, N, mc_set)
    color1 = 1/255 * [132 94 194; 178 91 0; 0 139 200];
    n0 = 1:N;
    dm = 0.01;
    m_range = 0:dm:3*log10(N)/b;
    P_large = 10.^(-b * m_range);
    alpha = 90;

    for j = 1:length(n0)
        n = n0(j);
        P = binopdf(n, N, P_large);
        [~, h0] = max(P);
        median(j) = m_range(h0);
        TT = cumsum(P) / sum(P);
        [~, h1] = min(abs(TT - (1 - alpha/100)/2));
        [~, h2] = min(abs(TT - (1 + alpha/100)/2));
        small(j) = m_range(h1);
        large(j) = m_range(h2);
    end

    x_fill = [small + mc_set, fliplr(large + mc_set)];
    y_fill = [n0, fliplr(n0)];
    fill(x_fill, y_fill, color1(3,:), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold on;

    pl = plot(median + mc_set, n0, '-', 'Color', color1(3,:), 'LineWidth', 3);
    hold on;
end



