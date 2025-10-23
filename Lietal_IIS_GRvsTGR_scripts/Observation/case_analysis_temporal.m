function [Dt,DmplusMc,median_PLR,median_bGR,median_bTGR,median_mcorner, ...
    median_EQrate_aGR,median_COV,median_LV,mmax_pre, mmax_obs]=case_analysis_temporal(D,Mc,global_b_GR, global_b_TGR, global_mcorner,time_flag)
    figure('units','pixels','position',[100, 100, 1000, 1600]);
    tiledlayout(4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');  % <--- use tight spacing
    D=sortrows(D,1);
    T0=D(:,1);
       
    Order0=1:1:length(T0);
    Order0=Order0';
    M0=D(:,5);
    
    diffm=diff(M0);
    nonzero_elements=diffm(diffm~=0);
    delta= min(abs(nonzero_elements))/2;

    jkf=find(M0-Mc>=-1e-11);
    Dm=M0(jkf);
    Dt=T0(jkf);
    Dorder=Order0(jkf);
    Dm=Dm-Mc;
    Mmin=0;
                         
    %% temporal evolution 
    Ntotal=length(Dm);
    min_size=100;
    max_size=1000;
    w0=ceil(Ntotal/max_size)+1:1:min([floor(Ntotal/min_size)-1,ceil(Ntotal/max_size)+1+2]); % number of windows
    
    [P_LR_median_eachwindow,bGR_median_eachnode,EQrate_aGR_median_eachnode, ...
        bTGR_median_eachnode,mcorner_median_eachnode, COV_median_eachnode, LV_median_eachnode]=get_temporal_variation(Dm,Dt,time_flag,delta,Mc,Mmin,min_size,max_size,w0);
    DmplusMc=Dm+Mc;
    
    nexttile
    median_PLR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,P_LR_median_eachwindow,1);
    nexttile
    median_bGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,bGR_median_eachnode,2);
    nexttile
    median_EQrate_aGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,EQrate_aGR_median_eachnode,3);
    nexttile
    median_bTGR=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,bTGR_median_eachnode,4);
    nexttile
    median_mcorner=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,mcorner_median_eachnode,5);
    nexttile
    median_COV=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,COV_median_eachnode,6);
    nexttile
    median_LV=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,LV_median_eachnode,7);
    
   %% Maximum magnitude evolution plot
    [mmax_pre, mmax_obs] = calc_mmax_obs_pred(DmplusMc, Mc, global_b_GR, global_b_TGR, global_mcorner, delta);
    nexttile

    % Axis limits based on observed values
    min_x = floor(min(mmax_obs));
    max_x = ceil(max(mmax_obs));
    plot([min_x, max_x], [min_x, max_x], '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]); 
    hold on;

    % Define fixed colors
    color_GR = [1 0 0];   % Red for GR
    color_TGR = [0 0 1];  % Blue for TGR

    % --- GR Stage 1 ---
    h(1) = plot(mmax_pre.GR_mode_st1, mmax_obs, '-', 'Color', color_GR, 'LineWidth', 1.5);
    h(2) = plot(mmax_pre.GR_5_st1, mmax_obs, '-*', 'Color', color_GR, 'LineWidth', 1.2);
    h(3) = plot(mmax_pre.GR_95_st1, mmax_obs, '-o', 'Color', color_GR, 'LineWidth', 1.2);

    % --- GR Stage 2 ---
    h(4) = plot(mmax_pre.GR_mode_st2, mmax_obs, '--', 'Color', color_GR, 'LineWidth', 1.5);
    h(5) = plot(mmax_pre.GR_5_st2, mmax_obs, '--*', 'Color', color_GR, 'LineWidth', 1.2);
    h(6) = plot(mmax_pre.GR_95_st2, mmax_obs, '--o', 'Color', color_GR, 'LineWidth', 1.2);

    % --- TGR Stage 1 ---
    h(7) = plot(mmax_pre.TGR_mode_st1, mmax_obs, '-', 'Color', color_TGR, 'LineWidth', 1.5);
    h(8) = plot(mmax_pre.TGR_5_st1, mmax_obs, '-*', 'Color', color_TGR, 'LineWidth', 1.2);
    h(9) = plot(mmax_pre.TGR_95_st1, mmax_obs, '-o', 'Color', color_TGR, 'LineWidth', 1.2);

    % --- TGR Stage 2 ---
    h(10) = plot(mmax_pre.TGR_mode_st2, mmax_obs, '--', 'Color', color_TGR, 'LineWidth', 1.5);
    h(11) = plot(mmax_pre.TGR_5_st2, mmax_obs, '--*', 'Color', color_TGR, 'LineWidth', 1.2);
    h(12) = plot(mmax_pre.TGR_95_st2, mmax_obs, '--o', 'Color', color_TGR, 'LineWidth', 1.2);

    % Add legend
    legend(h, {
        'GR mode st1', 'GR 5th st1', 'GR 95th st1', ...
        'GR mode st2', 'GR 5th st2', 'GR 95th st2', ...
        'TGR mode st1', 'TGR 5th st1', 'TGR 95th st1', ...
        'TGR mode st2', 'TGR 5th st2', 'TGR 95th st2' ...
    }, 'Location', 'bestoutside');

    % Labeling and formatting
    xlabel('Predicted M_{max}');
    ylabel('Observed M_{max}');
    set(gca, 'FontSize', 16);
    grid on; box on; grid minor;
    xlim([min_x, max_x]);
    ylim([min_x, max_x]);

end

function [mmax_pre, mmax_obs] = calc_mmax_obs_pred(M, Mc, global_b_GR, global_b_TGR, global_mcorner, delta)
    numpercal = 100;
    numcal = floor(length(M) / numpercal);

    mmax_obs = NaN(1, numcal);

    mmax_pre.GR_mode_st1 = NaN(1, numcal);
    mmax_pre.GR_5_st1 = NaN(1, numcal);
    mmax_pre.GR_95_st1 = NaN(1, numcal);
    mmax_pre.GR_mode_st2 = NaN(1, numcal);
    mmax_pre.GR_5_st2 = NaN(1, numcal);
    mmax_pre.GR_95_st2 = NaN(1, numcal);

    mmax_pre.TGR_mode_st1 = NaN(1, numcal);
    mmax_pre.TGR_5_st1 = NaN(1, numcal);
    mmax_pre.TGR_95_st1 = NaN(1, numcal);
    mmax_pre.TGR_mode_st2 = NaN(1, numcal);
    mmax_pre.TGR_5_st2 = NaN(1, numcal);
    mmax_pre.TGR_95_st2 = NaN(1, numcal);

    for j = 1:numcal
        jkf = 1:min(j * numpercal, length(M));
        mmax_obs(j) = max(M(jkf));

        % GR using global b-value
        mmax_pre.GR_mode_st1(j) = Mc + log10(length(jkf)) / global_b_GR;
        mmax_pre.GR_5_st1(j) = mmax_pre.GR_mode_st1(j) - (1 / global_b_GR) * log10(-log(0.05));
        mmax_pre.GR_95_st1(j) = mmax_pre.GR_mode_st1(j) - (1 / global_b_GR) * log10(-log(0.95));

        % TGR using global parameters
        [mmax_pre.TGR_mode_st1(j), mmax_pre.TGR_5_st1(j), mmax_pre.TGR_95_st1(j)] = ...
            get_mmax_tgr(global_b_TGR, global_mcorner - Mc, length(jkf), Mc);

        % GR with estimated b-value
        [b_GR, ~] = Estimation_GR_discrete(M(jkf) - Mc, 0, delta);
        mmax_pre.GR_mode_st2(j) = Mc + log10(length(jkf)) / b_GR;
        mmax_pre.GR_5_st2(j) = mmax_pre.GR_mode_st2(j) - (1 / b_GR) * log10(-log(0.05));
        mmax_pre.GR_95_st2(j) = mmax_pre.GR_mode_st2(j) - (1 / b_GR) * log10(-log(0.95));

        % TGR with estimated parameters
        [b_TGR, mcorner_TGR, ~] = Estimation_TGR_gridsearch_discrete(M(jkf) - Mc, 0, delta, b_GR);
        [mmax_pre.TGR_mode_st2(j), mmax_pre.TGR_5_st2(j), mmax_pre.TGR_95_st2(j)] = ...
            get_mmax_tgr(b_TGR, mcorner_TGR, length(jkf), Mc);
    end
end

function [mmax_mode, mmax_5, mmax_95] = get_mmax_tgr(b_TGR, mcorner_TGR, L, Mc)
    M_corner = 10^(1.5 * mcorner_TGR + 9.1);
    mmin_true = 0;
    Mmin = 10^(1.5 * mmin_true + 9.1);
    beta = (2 / 3) * b_TGR;

    mtest = 0:0.01:log10(L) / b_TGR;
    M = 10.^(1.5 * mtest + 9.1);
    misfit = abs((Mmin ./ M).^beta .* exp((Mmin - M) ./ M_corner) - 1 / L);
    [~, ind] = min(misfit);
    mmax_mode = mtest(ind) + Mc;

    max_TGR = NaN(1, 1e3);
    for k = 1:1e3
        Dm_TGR = gentgr(beta, L, Mmin, M_corner);
        max_TGR(k) = max(Dm_TGR) + Mc;
    end

    pct = prctile(max_TGR, [5, 95]);
    mmax_5 = pct(1);
    mmax_95 = pct(2);

end



function Dm = gentgr(beta, L, Mmin, M_corner)
    U10 = rand(1, L);
    U20 = rand(1, L);
    U1 = -M_corner * log(U10);
    U2 = Mmin * U20.^(-1/beta);
    DM = min(U1, U2);
    Dm = 2/3 * log10(DM) -  6.07;
end




function [P_LR_median_eachwindow,bGR_median_eachnode,EQrate_aGR_median_eachnode, ...
    bTGR_median_eachnode,mcorner_median_eachnode,COV_median_eachnode,LV_median_eachnode]=get_temporal_variation(Dm,Dt,time_flag,delta,Mc,mmin,min_size,max_size,w0)
    Ntotal=length(Dm);
    num_rand=200; % for each w0, number of random partition
    P_LR_median_eachwindow=zeros(length(w0),Ntotal);
    bGR_median_eachnode=zeros(length(w0),Ntotal);
    EQrate_aGR_median_eachnode=zeros(length(w0),Ntotal);
    bTGR_median_eachnode=zeros(length(w0),Ntotal);
    mcorner_median_eachnode=zeros(length(w0),Ntotal);
    COV_median_eachnode=zeros(length(w0),Ntotal);
    LV_median_eachnode=zeros(length(w0),Ntotal);

    if time_flag == 2
        Dt = (Dt - min(Dt)) * 365.25;
    end
    
    for k=1:length(w0)
        w=w0(k);
        P_LR=zeros(num_rand,Ntotal);
        bGR=zeros(num_rand,Ntotal);
        EQrate_aGR=zeros(num_rand,w+1);
        bTGR=zeros(num_rand,Ntotal);
        mcorner=zeros(num_rand,w+1);
        COV=zeros(num_rand,w+1);
        LV=zeros(num_rand,w+1);
        node=zeros(num_rand,w+1);
        
        for i=1:num_rand
            node(i,:)=generateNumbers(w,length(Dm),min_size,max_size);% w+1 nodes, including the first and last events
            for j=1:w
                Dm_win=Dm(node(i,j):node(i,j+1));
                Dt_win=Dt(node(i,j):node(i,j+1));
                
                [b_GR,loglikelihood_GR]=Estimation_GR_discrete(Dm_win,mmin,delta);
                [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_discrete(Dm_win,mmin,delta,b_GR);
                R=loglikelihood_TGR-loglikelihood_GR;
                
                P_LR(i,node(i,j):node(i,j+1))=1-chi2cdf(2*R,1);
                
                bGR(i,node(i,j):node(i,j+1))=b_GR;
                
                EQrate_aGR(i,node(i,j):node(i,j+1))=10^(b_GR*Mc+log10(length(Dt_win)))/(max(Dt_win)-min(Dt_win));
                EQrate_aGR(i,node(i,j):node(i,j+1))=log10(EQrate_aGR(i,node(i,j):node(i,j+1)));
                
                bTGR(i,node(i,j):node(i,j+1))=b_TGR;
                mcorner(i,node(i,j):node(i,j+1))=mcorner_TGR+Mc;    
                
                COV(i,node(i,j):node(i,j+1))=std(diff(Dt_win))/mean(diff(Dt_win));
                LV_dT = 0;
                diff_T=diff(Dt_win);
                nonzero_vals = diff_T(diff_T ~= 0);
                min_nonzero = min(abs(nonzero_vals));
                diff_T(diff_T == 0) = min_nonzero / 2;
                num_diff=length(diff_T);
                for jk = 1:num_diff-1
                    LV_dT = LV_dT + 3*(diff_T(jk) - diff_T(jk+1))^2 / (diff_T(jk) + diff_T(jk+1))^2;
                end

                LV_dT = LV_dT / (num_diff-1);
                LV(i,node(i,j):node(i,j+1))=LV_dT;
                
            end
        end
        P_LR_median_eachwindow(k,:)=median(P_LR,1); 
        bGR_median_eachnode(k,:)=median(bGR,1); 
        EQrate_aGR_median_eachnode(k,:)=median(EQrate_aGR,1);
        bTGR_median_eachnode(k,:)=median(bTGR,1); 
        mcorner_median_eachnode(k,:)=median(mcorner,1);    
        COV_median_eachnode(k,:)=median(COV,1);
        LV_median_eachnode(k,:)=median(LV,1);

    end
    
end

function nonde_order=generateNumbers(n,Ntotal,min_size,max_size) % n is the window number
    nonde_order=zeros(1,n+1);
    target_sum=Ntotal-n*min_size; % available number for arrangement
    scaled_numbers=inf*ones(1,n);
    while max(scaled_numbers)>=max_size-min_size
        random_numbers=rand(1,n);
        scaled_numbers=floor(random_numbers/sum(random_numbers)*target_sum); % number of events (-100) for each window
    end
    nonde_order(1)=1;
    for i=2:n
        nonde_order(i)=nonde_order(i-1)+min_size+scaled_numbers(i-1);
    end
    nonde_order(n+1)=Ntotal;
end

function [meadian_plt_variable]=plt_temporalMFD(T0,M0,Dt,DmplusMc,w0,plt_variable,plt_variable_flag)
    yyaxis left
    scatter(T0,M0,10,[0.75 0.75 0.75],'filled');
    hold on;
    scatter(Dt,DmplusMc,10,'k','filled');
    xlabel('Time')
    ylabel('Magntidue');
    ax = gca; 
    ax.YColor = 'k';

    yyaxis right
    col=autumn(length(w0));
    for k=1:length(w0)
        h(k)=plot(Dt,plt_variable(k,:),'-','color',col(k,:),'linewidth',2.5);hold on;
        leg{k}=['$w=$',num2str(w0(k))];
    end
    meadian_plt_variable=median(plt_variable,1);
    h(k+1)=plot(Dt,meadian_plt_variable,'-b','linewidth',2.5);hold on;
    leg{k+1}=['median'];
    if plt_variable_flag==1
        ylabel('$p_{LRT}$','Interpreter','latex')
        ylim([0 1]);
    elseif plt_variable_flag==2
        ylabel('$b_{GR}$','Interpreter','latex')
        legend(h,leg, 'Interpreter', 'latex','location','northwest','Numcolumns',3);
    elseif plt_variable_flag==3
        ylabel('$a_{GR} \, (\mathrm{d}^{-1})$', 'Interpreter', 'latex');
    elseif plt_variable_flag==4
        ylabel('$b_{TGR}$','Interpreter','latex')
    elseif plt_variable_flag==5
        ylabel('$m_{corner}$','Interpreter','latex')
    elseif plt_variable_flag==6
        ylabel('CV_{\Delta\itt}');
    elseif plt_variable_flag==7
        ylabel('LV_{\Delta\itt}');
    end  
    ax = gca; 
    ax.YColor = 'k';
    set(gca,'fontsize',16);
    grid on;box on;grid minor;
end