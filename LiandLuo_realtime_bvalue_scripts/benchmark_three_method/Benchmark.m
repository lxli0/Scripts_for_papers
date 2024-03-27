clc,clear
%% Initial setting
Sim=1e4; % number of random catalog
L=1e3; % catalog size
b=1; % the set b-value
%% continuous
figure('units','normalized','position',[0.1,0.1,0.35,0.8])

subplot(3,1,1)
[mean_MLE_con,std_MLE_con]=MLE_continuous(Sim,L,b);
hold on;
set(gca,'fontsize',16)
ylabel( 'Frequency' )
xlabel( 'b-value' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line

subplot(3,1,2)
[mean_bpos_con,std_bpos_con]=bpositive_continuous(Sim,L,b);
hold on;
set(gca,'fontsize',16)
ylabel( 'Frequency' )
xlabel( 'b-value' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line

subplot(3,1,3)
[mean_KMS_con,std_KMS_con]=KMS_continuous(Sim,L,b);
hold on;
xlabel( 'b-value' )
ylabel( 'Frequency' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
grid on;
set(gca,'fontsize',16)

%% discrete
figure('units','normalized','position',[0.1,0.1,0.35,0.8])

subplot(3,1,1)
[mean_MLE_dis,std_MLE_dis]=MLE_discrete(Sim,L,b);
hold on;
set(gca,'fontsize',16)
ylabel( 'Frequency' )
xlabel( 'b-value' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line

subplot(3,1,2)
[mean_bpos_dis,std_bpos_dis]=bpositive_discrete(Sim,L,b);
hold on;
set(gca,'fontsize',16)
ylabel( 'Frequency' )
xlabel( 'b-value' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line

subplot(3,1,3)
[mean_KMS_dis,std_KMS_dis]=KMS_discrete(Sim,L,b);
hold on;
xlabel( 'b-value' )
ylabel( 'Frequency' )
plot( [ 1 1 ] , [ 0 ,0.3 ] , '--k' ) % add the real b-value line
grid on;
set(gca,'fontsize',16)
