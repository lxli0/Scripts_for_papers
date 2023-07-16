% The script tests the b-positive method for the ideal catalog (discrete magnitude)
clc,clear
Sim=2e3; % number of random catalogs
L=1e3; % number of events in each simulated catalog
b=1; % set b-value, we should expect the methods provide estimation close to this value
delta=0.05; % magnitudes are discretized by 2*delta
Magn_Start = -delta; % starting magnitude for simulations
Mc_prime=2*delta; % Mc'
for i=1:Sim
    Pcm=rand(1,L);
    Magn=-1/b.*log10(Pcm)+ Magn_Start; % generate random magnitude
    Magn=roundn(Magn,log10(2*delta)); % magnitude discretization
    Dm0=diff(Magn); % calculate magnitude difference
    jkf=find((Dm0-Mc_prime)>=-1e-6);
    Dm=Dm0(jkf); % filter earthquakes by Mc'
    b_es(i)=1/delta/log(10)*acoth((mean(Dm)-Mc_prime+delta )/delta); % calculate b-value
end
mean(b_es)
std(b_es)

% plot
figure('units','normalized','position',[0.1,0.1,0.2,0.2])
edges = [0.85:0.01:1.15];
h=histogram(b_es,edges,'Normalization','probability');
h.FaceColor = 1/255*[209 41 32];
h.EdgeColor=[1 1 1];
xlabel( 'b-value' )
ylabel( 'Frequency')
xlim([0.85 1.15]);
ylim([0 0.15])
set(gca,'fontsize',16);
grid on;


