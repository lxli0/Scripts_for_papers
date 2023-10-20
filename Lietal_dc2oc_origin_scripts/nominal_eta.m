clc,clear
sensitivity=[1.15e-7 1.65e-8 1.43e-8 4e-9 7.1e-8 1.4e-7 1.3e-7 1e-6 3.2e-8 ...
    1e-6 1e-7 1e-8 1e-6 3.5e-6 5e-7 5e-8 1.6e-7 5e-7 1.3e-9 ...
    5e-9 1.5e-9 6e-9 5.6e-9 2.65e-9 4e-7 5e-7 7e-9 1e-8 ];
frequency=[3 8;0.4 0.8;0.8 1.6;0.4 1.0;0.5 1;1 2;2 4;77 77;80 240;
    0.3 1;0.3 1;0.3 1;500 500;1 2;1 5;2 4;3 7;2 8;0.04 0.1;
    0.33 1;0.1 0.9;0.5 0.9;0.5 2;1 3;0.5 2;0.5 2;0.06 0.14;100 240];
number=[7 5 6 5 5];
color=1/255*[
    253,109,90
    254,180,11
    109,195,84
    153,68,135
    81,140,216
    ];
all_number=length(sensitivity);
k=0;
figure('units','normalized','position',[0.1,0.1,0.3,0.5])
for j=1:length(number)
    for i=1:number(j)
        k=k+1;
        if frequency(k,1)~=frequency(k,2)
            h(k)=loglog([sensitivity(k),sensitivity(k)],1./[frequency(k,1),frequency(k,2)], ...
            'linewidth',2.5,'color',color(j,:));
        else
            h(k)=loglog([sensitivity(k),sensitivity(k)],1./[frequency(k,1),frequency(k,2)],'o', ...
            'markersize',5, 'markerfacecolor',color(j,:),'MarkerEdgeColor', color(j,:));
        end   
        text(1.1*sensitivity(k),(1/frequency(k,1)+1/frequency(k,2))/2,num2cell(k));
        hold on;
    end
end
legend(h([1,9,15,20,25]),'Volcano','Atmosphere','Earth tide','Earthquake','Other','NumColumns',2,'Location','SouthEast')

xlabel('\eta_c (Pa^-^1)')
ylabel('Period (s)')
set(gca,'YDir','reverse') 
grid on;
box on;
set(gca,'fontsize',16);
