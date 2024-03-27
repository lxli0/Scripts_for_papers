clc,clear
Cat_Raw=load( 'Cat_16Nov_ZMAP.txt' ) ;
% select the events with the previous threshold criteria
Center   = [ 13.324 , 44.013 ] ;  % Long and Lat
DistMax  = 30;                   % Km
DepthMax = 30;                   % Km
Cat = Cat_Raw( Cat_Raw( : , 7 ) <= DepthMax & ...
      distance( Cat_Raw( : , 2 ) , Cat_Raw( : , 1 ) , Center(2) , Center(1) ) ...
      .* pi/180*6371 <= DistMax , : ) ;
% select events based on STAI
DeltaT = 1/6 ; %in days
Time = datenum( Cat( : , [ 3 : 5 , 8 : 10 ]) ) - ...
        datenum( Cat( 1 , [ 3 : 5 , 8 : 10 ]) ) ;
Cat_NoSTAI = Cat( Time >= DeltaT , : ) ;
Dm=Cat_NoSTAI(:,6);

[Mc_MBS,b_mle,number_mle]=MBS_MLE_discrete(Dm);
jkf=Dm>=(Mc_MBS-1e-6);
Dm_complete=Dm(jkf);
b_KMS=KMS(Dm_complete);
[Mc_pos,b_pos,number_pos]=MBS_pos_discrete(Dm);
