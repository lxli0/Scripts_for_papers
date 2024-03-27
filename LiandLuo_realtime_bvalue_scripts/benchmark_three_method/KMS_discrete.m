% The script tests the KMS method for the ideal catalog (discrete magnitude)
% Sim: number of random catalogs
% L: number of events in each simulated catalog
% b: set b-value, we should expect the methods provide estimation close to this value
function [mean_b,std_b]=KMS_discrete(Sim,L,b)   
    delta=0.05; % magnitudes are discretized by 2*delta
    Magn_Start = -delta; % starting magnitude for simulations
    b_es=zeros(1,Sim);
    for i=1:Sim
        Pcm=rand(1,L);
        Magn=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitude
        Magn=roundn(Magn,log10(2*delta)); % magnitude discretization
        b_es(i)=KMS(Magn);
        Sim-i
    end
    mean_b=mean(b_es);
    std_b=std(b_es);

    % plot
    edges = [0.85:0.01:1.15];
    h=histogram(b_es,edges,'Normalization','probability');
    h.FaceColor = 1/255*[250 192 61];
    h.EdgeColor=[1 1 1];
    xlim([0.85 1.15]);
    ylim([0 0.15])
    grid on;
end

