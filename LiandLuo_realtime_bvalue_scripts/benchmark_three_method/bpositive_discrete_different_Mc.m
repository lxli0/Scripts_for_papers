% In this script, we test the b-positive method for ideal catalogs under different magnitude difference threshold for discrete situation
clc,clear
Sim =1e4; % number of simulated catalogs
L = 1e4; % number of events in each simulated catalog
b=1; % set b-value, we should expect the methods provide estimation close to this value
mc=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0]; % Mc or Mc'
delta=0.05; % magnitudes are discretized by 2*delta
Magn_Start=-delta; % starting magnitude for simulations

for i = 1 : Sim 
    Pcm=rand(1,L);
    Magn=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitudes
    Magn=roundn(Magn,log10(2*delta)); % magnitude discretization
    diff_M=diff(Magn); % calculate magnitude difference
    for j=1:length(mc)
        [apo(i,j),bpo(i,j)]=cal_a_b(diff_M,2,50,delta,mc(j)); % calculate a- and b-values by b-positive method 
    end
    Sim-i
end

for j=1:length(mc)
    mean_b1(j)=mean(bpo(:,j));
    std_b1(j)=std(bpo(:,j));
end

% plot
figure;
errorbar(mc, mean_b1,std_b1,'-o','linewidth',2);
hold on;

% b-positive
function [a,b]=cal_a_b(Dm,pl,limit,delta,Mmin)
    L=length(Dm);
    dm=0.1;
    m=0:dm:floor(max(Dm)/dm)*dm;
    n0=hist(Dm,m);   
    cn0(1)=L;
    for i=2:floor(max(Dm)/dm)+1
        cn0(i)=L-sum(n0(1:i-1)); %number of samples no less than mM
    end
    n=log10(n0);
    cn=log10(cn0);
    % Calculate b value    
    ind=find(m-Mmin>=-1e-6);
    jkf=find(Dm-Mmin>=-1e-6);
    if length(jkf)>=limit
        Dm0=Dm(jkf);
        b=1/delta/log(10)*acoth(( mean(Dm0)-Mmin+delta )/delta);
        hhf=find(cn0>=10);
        for f=ind:hhf(end)
            a0(f)=log10(cn0(f))+b*m(f);
        end
        a0(a0==0)=NaN;
        z=find(~isnan(a0));
        a01=a0(z);
        mean_a=mean(a01);
        std_a=std(a01);
        if std_a/mean_a<0.1
            a=mean_a;
            if pl==1 %plot
                color=1/255*[132 94 194;178 91 0;0 139 200];
                plot_n0=semilogy(m,n0,'o', 'color', 0.4*[1 1 1], ...
                            'markersize', 6, 'markerfacecolor', color(3,:), ...
                            'MarkerEdgeColor', color(3,:));
                hold on;
                plot_cn0=semilogy(m,cn0,'o', 'color', 0.4*[1 1 1], ...
                            'markersize', 6, 'markerfacecolor', color(2,:), ...
                            'MarkerEdgeColor', color(2,:));
                grid on;
                box on;
                xlabel('Magnitude');
                ylabel('Number');
                set(gca,'fontsize',16);

                scatter(Mmin,cn0(ind),180,color(1,:),'v','filled');hold on;
                plot_range=Mmin:0.1:m(hhf(end));
                y=a-b*plot_range;
                plot(plot_range,10.^y,'-','color',color(1,:),'linewidth',2.5);
                legend([plot_n0,plot_cn0],'Non-cumulative','Cumulative','Location','Northeast');
                str1 = {'Mc=','a=','b='};
                str2 = {Mmin,a,b};
                text('string',str1, 'Units','normalized','position',[0.75,0.65],'fontsize',16); 
                text('string',str2, 'Units','normalized','position',[0.84,0.65],'fontsize',16); 
            end
        else
            b=nan;
            a=nan;
        end
    else 
        a=nan;
        b=nan;
    end
end