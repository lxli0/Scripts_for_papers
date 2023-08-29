% The script tests the b-positive method for the ideal catalog (continuous magnitude)
% Sim: number of random catalogs
% L: number of events in each simulated catalog
% b: set b-value, we should expect the methods provide estimation close to this value
function [mean_b,std_b]=bpositive_continuous(Sim,L,b)
    Magn_Start = 0; % starting magnitude for simulations
    for i=1:Sim
        Pcm=rand(1,L);
        Magn=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitude
        Dm0=diff(Magn); % calculate magnitude difference
        jkf=Dm0>=-1e-6;
        Dm=Dm0(jkf); % filter earthquakes by Mc'
        Mmin=0;
        b_es(i)=1/(log(10)*(mean(Dm)-Mmin)); % calculate b-value
    end
    mean_b=mean(b_es);
    std_b=std(b_es);

    % plot
    edges = [0.85:0.01:1.15];
    h=histogram(b_es,edges,'Normalization','probability');
    h.FaceColor = 1/255*[209 41 32];
    h.EdgeColor=[1 1 1];
    xlim([0.85 1.15]);
    ylim([0 0.15])
    grid on;
end
