clc,clear
b=1; % set b-value
beta=2/3*b;
mmin=0;
Mmin=10^(1.5*(mmin+6.07));

time_cut=10;

times=1000;
MFD_type=[1 1 2 2 2 2];% 1 GR; 2TGR
L=[1e3 2e3 1e3 1e3 2e3 2e3];
mcorner=[inf inf 2 2.5 2 2.5];
correct=zeros(length(MFD_type),times,time_cut);

for k=1:length(MFD_type)
    k
    parpool(8);
    parfor i=1:times 
        if MFD_type(k)==1
            Dm0=-1/b.*log10(rand(1,L(k)))+mmin; % generate random magnitude
        else
            M_corner=10^(1.5*(mcorner(k)+6.07));
            U10=rand(1,L(k));
            U20=rand(1,L(k));
            U1=Mmin-M_corner*log(U10);
            U2=Mmin*U20.^(-1/beta);
            DM=min(U1,U2);
            Dm0=2/3*log10(DM)-6.07;
        end
        for j=1:time_cut
            Dm=Dm0(1:j/time_cut*L);
            [b_GR,loglikelihood_GR]=Estimation_GR(Dm,mmin);
            [b_TGR,mcorner_TGR,loglikelihood_TGR]=Estimation_TGR_gridsearch_continuous(Dm,mmin)

            % Likelihood-ratio test
            R=loglikelihood_TGR-loglikelihood_GR;
            P_LR=1-chi2cdf(2*R,1);
            if P_LR>0.05
                Target_LRT=1; %GR
            else
                Target_LRT=2;
            end
            if Target_LRT==MFD_type(k)
                correct(k,i,j)=1;
            end
        end   
    end
    delete(gcp('nocreate'))
end


figure('units','normalized','position',[0.1,0.1,0.4,0.25])
color=1/255*[247 170 88;231 98 84;204 121 167;170 68 153;136 204 238;51 34 136];

for k=1:length(MFD_type)
    correct_tar=reshape(correct(k,:,:),times,[]);
    subplot(1,2,2)
    num_correct=sum(correct_tar,2);
    poten=0:1:10;
    for i=1:length(poten)
        num(i)=length(find(num_correct==poten(i)));
    end
    scatter(poten,num/times,50,color(k,:),'filled');hold on;
    xlabel('Number of subcatalogs');
    ylabel('Frequency');
    grid on;box on;grid minor;
    set(gca,'fontsize',16);
    ylim([0 1])
    if MFD_type(k)==1
        leg{k}=['GR,',' $n$=',num2str(L(k))];
    else
        leg{k}=['TGR,',' $n$=',num2str(L(k)),', $m_{corner}$=',num2str(mcorner(k))];
    end  

    subplot(1,2,1)
    order_correct=sum(correct_tar,1);
    scatter(1:10,order_correct/times,50,color(k,:),'filled');hold on;
    xlabel('# of the subcatalog')
    ylabel('$R$', 'Interpreter', 'latex');
    grid on;
    box on;
    set(gca,'fontsize',16);
    ylim([0 1])
    xlim([1,10])
end
subplot(1,2,2)
legend(leg,'location','north', 'Interpreter', 'latex');


