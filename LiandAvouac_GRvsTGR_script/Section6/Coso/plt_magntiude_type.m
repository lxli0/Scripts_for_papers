T=Dt;
for i=1:length(T)
    count_w(i)=length(find(MTYPE(1:i)==1));
    count_n(i)=length(find(MTYPE(1:i)==2));
    count_l(i)=length(find(MTYPE(1:i)==3));
    count_h(i)=length(find(MTYPE(1:i)==4));
    count_c(i)=length(find(MTYPE(1:i)==5));
    count_none(i)=length(find(MTYPE(1:i)==6));
end
color=1/255*[247 170 88;231 98 84;204 121 167;170 68 153;136 204 238;51 34 136];
plot(T,count_c,'-','color',color(1,:),'linewidth',2.5);hold on;
plot(T,count_l,'-','color',color(2,:),'linewidth',2.5);hold on;
plot(T,count_h,'-','color',color(3,:),'linewidth',2.5);hold on;
plot(T,count_n,'-','color',color(4,:),'linewidth',2.5);hold on;
plot(T,count_none,'-','color',color(5,:),'linewidth',2.5);hold on;
plot(T,count_w,'-','color',color(6,:),'linewidth',2.5);hold on;

grid on;box on;grid minor;
set(gca,'fontsize',16);
xlabel('Year');
ylabel('Cumulative Number')
xlim([1981 2022])
legend('Coda Amplitude','Local','Helicorder','No Magnitude','None','Moment','location','northwest', 'FontSize', 14);

%set(gca, 'Color', 'none');
%set(gcf, 'Color', 'none'); 