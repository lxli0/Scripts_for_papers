% The script tests the MLE method for the ideal catalog (continuous magnitude)
% Sim: number of random catalogs
% L: number of events in each simulated catalog
% b: set b-value, we should expect the methods provide estimation close to this value
function [mean_b,std_b]=MLE_continuous(Sim,L,b)
    Magn_Start = 0; % starting magnitude for simulations
    for i=1:Sim
        Pcm=rand(1,L);
        Magn=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitudes
        Mmin=0;
        b_es(i)=1/(log(10)*(mean(Magn)-Mmin)); % estimated b-value by MLE method
    end
    mean_b=mean(b_es);
    std_b=std(b_es);

    % plot
    edges = [0.85:0.01:1.15];
    h=histogram(b_es,edges,'Normalization','probability');
    h.FaceColor = 1/255*[46 89 167];
    h.EdgeColor=[1 1 1];
    xlim([0.85 1.15]);
    ylim([0 0.15])
    grid on;
end
