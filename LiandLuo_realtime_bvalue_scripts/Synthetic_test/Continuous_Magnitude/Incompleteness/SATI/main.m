%% The script evaluates the performance of the three methods under STAI (continuous magnitudes)
clc,clear
Sim=200;
tdelete=[0,0.2,1,2];
  
for i=1:Sim
    mp=7;
    [T,M]=gen_catalog(mp);
    for j=1:length(tdelete)
        jkf=find(T>=tdelete(j));
        T=T(jkf);
        Magn=M(jkf);
        num0(i)=length(Magn);
        [Mc_mle(i,j),b_mle(i,j),number_mle(i,j)]=MBS_MLE(Magn);
        [Mc_pos(i,j),b_pos(i,j),number_pos(i,j)]=MBS_pos(Magn);  
        Mc_kms(i,j)=Mc_mle(i,j);
        number_kms(i,j)=length(Magn(Magn>=Mc_kms(i,j)));
        b_kms(i,j)=KMS(Magn,Mc_kms(i,j));
    end   
    Sim-i
end


for j=1:length(tdelete)
    meanb_mle(j)=mean(b_mle(:,j));
    stdb_mle(j)=std(b_mle(:,j));
    meannum_mle(j)=mean(number_mle(:,j));
    stdnum_mle(j)=std(number_mle(:,j));
    meanMc_mle(j)=mean(Mc_mle(:,j));
    stdMc_mle(j)=std(Mc_mle(:,j));
    
    meanb_pos(j)=mean(b_pos(:,j));
    stdb_pos(j)=std(b_pos(:,j));
    meannum_pos(j)=mean(number_pos(:,j));
    stdnum_pos(j)=std(number_pos(:,j));
    meanMc_pos(j)=mean(Mc_pos(:,j));
    stdMc_pos(j)=std(Mc_pos(:,j));
    
    meanb_kms(j)=mean(b_kms(:,j));
    stdb_kms(j)=std(b_kms(:,j));
    meannum_kms(j)=mean(number_kms(:,j));
    stdnum_kms(j)=std(number_kms(:,j));
    meanMc_kms(j)=mean(Mc_kms(:,j));
    stdMc_kms(j)=std(Mc_kms(:,j));
end
%% plot
tdelete_plot=1:1:length(tdelete);
figure('units','normalized','position',[0.1,0.1,0.4,0.6])
color=1/255*[117 114 181; 197 86 89;203 180 123];
%color=1/255*[77 133 189;247 144 61;89 169 90];
subplot(3,1,1)
errorbar(tdelete_plot-0.07,meanb_mle,stdb_mle,'o','color',color(1,:),'linewidth',2)
hold on;
errorbar(tdelete_plot,meanb_pos,stdb_pos,'o','color',color(2,:),'linewidth',2)
hold on;
errorbar(tdelete_plot+0.07,meanb_kms,stdb_kms,'o','color',color(3,:),'linewidth',2)
legend('MLE','b-positive','KMS')
set(gca,'fontsize',16)
ylabel('b-value')
set(gca,'xticklabel',[]);
xlim([0.5,4.5])
ylim([0.6 1.2])
grid on;box on;
set(gca,'position',[0.15 0.69 0.75 0.25])

subplot(3,1,2)
yyaxis left
errorbar(tdelete_plot-0.07,meanMc_mle,stdMc_mle,'o','color',color(1,:),'linewidth',2)
hold on;
yyaxis right
errorbar(tdelete_plot,meanMc_pos,stdMc_pos,'o','color',color(2,:),'linewidth',2)
hold on;
set(gca,'ycolor',[0 0 0]);
ylabel('Mc''')
ylim([0 0.6])
yyaxis left
errorbar(tdelete_plot+0.07,meanMc_kms,stdMc_kms,'o','color',color(3,:),'linewidth',2)
set(gca,'fontsize',16)
xlim([0.5,4.5])
set(gca,'xticklabel',[]);
ylabel('Mc')
ylim([2 4])
grid on;box on;
set(gca,'ycolor',[0 0 0]);
set(gca,'position',[0.15 0.40 0.75 0.25])

subplot(3,1,3)
errorbar(tdelete_plot-0.07,meannum_mle,stdnum_mle,'o','color',color(1,:),'linewidth',2)
hold on;
errorbar(tdelete_plot,meannum_pos,stdnum_pos,'o','color',color(2,:),'linewidth',2)
hold on;
errorbar(tdelete_plot+0.07,meannum_kms,stdnum_kms,'o','color',color(3,:),'linewidth',2)
set(gca,'fontsize',16)
xlim([0.5,4.5])
ylim([0 1400])
xlabel('Days Delete')
ylabel('Sample Number')
set(gca,'xticklabel',{' ','0','','0.2','','1','','2',' '}); 
set(gca,'position',[0.15 0.11 0.75 0.25])
grid on;box on;
