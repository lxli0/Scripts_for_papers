clc,clear
sensitivity=[1.15e-7 1.65e-8 1.43e-8 4e-9 7.1e-8 1.4e-7 1.3e-7 1e-6 3.2e-8 ...
    1e-6 1e-7 1e-8 1e-6 1.4e-6 5e-7 4.2e-7 1.6e-7 5e-7 1.3e-9 ...
    5e-9 1.5e-9 6e-9 5.6e-9 2.65e-9 4e-7 5e-7 7e-9 1e-8 ];
frequency=[3 8;0.4 0.8;0.8 1.6;0.4 1.0;0.5 1;1 2;2 4;77 77;80 240;
    0.3 1;0.3 1;0.3 1;500 500;1 2;1 5;2 4;3 7;2 8;0.04 0.1;
    0.33 1;0.1 0.9;0.5 0.9;0.5 2;1 3;0.5 2;0.5 2;0.06 0.14;100 240];
number=[7 5 6 5 5];
color=1/255*[
    31, 119, 180
    188 189 34
    255, 127, 14
    214 39 40
    148 103 189
    ];
all_number=length(sensitivity);
k=0;
figure('units','normalized','position',[0.1,0.1,0.4,0.6])
for j=1:length(number)
    for i=1:number(j)
        k=k+1;
        if frequency(k,1)~=frequency(k,2)
            h(k)=loglog([frequency(k,1),frequency(k,2)],[sensitivity(k),sensitivity(k)], ...
            'linewidth',2.5,'color',color(j,:));
        else
            h(k)=loglog([frequency(k,1),frequency(k,2)],[sensitivity(k),sensitivity(k)],'o', ...
            'markersize',5, 'markerfacecolor',color(j,:),'MarkerEdgeColor', color(j,:));
        end   
        text((frequency(k,1)+frequency(k,2))/2, 1.2*sensitivity(k), num2str(k), 'FontSize', 16);
        hold on;
    end
end
legend(h([1,9,15,20,25]),'Volcano','Atmosphere','Earth tide','Earthquake','Other','NumColumns',1,'Location','NorthWest')
ylim([1e-9 2e-6])
ylabel('$\eta_c$ (Pa$^{-1}$)', 'Interpreter', 'latex');
xlabel('Frequency (Hz)')
grid on;
box on;
set(gca,'fontsize',16);
