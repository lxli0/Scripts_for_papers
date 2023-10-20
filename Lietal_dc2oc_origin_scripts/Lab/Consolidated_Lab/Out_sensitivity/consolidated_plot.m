clc,clear
figure('units','normalized','position',[0.1,0.1,0.9,0.6])
subplot(1,2,1)
%{
namelist1 = dir('*_out.txt');
file_name1 = {namelist1.name};
l = length(file_name1); 
z=0;
for i=1:l
    D=load(file_name1{i});
    number=length(D(:,1));
    Depth(z+1:z+number)=D(:,1);
    eta(z+1:z+number)=D(:,2);
    z=z+number;
end
loglog(Depth/1e6,eta,'o','markersize',8, 'markerfacecolor',[0.75 0.75 0.75],'MarkerEdgeColor',[0.75 0.75 0.75]);
%}
namelist2 = dir('*_nihe.txt');
file_name2 = {namelist2.name};
l = length(file_name2); 
Sample_number=l;
%% create colormap
color=nan(Sample_number,3);
color_end=1/255*[220,20,60];
color_begin=1/255*[255,255,0];
for i=1:Sample_number
    color(i,:)=color_begin+(color_end-color_begin)*i/Sample_number;
end
rank_number=[9
18
25
1
24
15
5
11
3
12
23
13
10
14
17
7
20
16
21
22
19
4
6
2
8]; % Based on sample number

for j=1:l
    i=find(rank_number==j);
    D=load(file_name2{i});
    w1(i)=D(1);
    w2(i)=D(2);
    depth_min(i)=D(3);
    depth_max(i)=D(4);
    xi=linspace(depth_min(i),depth_max(i),20);
    yi=w1(i).*xi.^(w2(i));
    h(j)=loglog(xi/1e6,yi,'linewidth',2.5,'color',color(rank_number(i),:));
    hold on;
end

str=["Christensen and Wang, 1985";
"Yukutake et al., 1988";
"King, 2002";
"Winkler and McGowan, 2004";
"Hicks and Berry, 1956";
"Wyllie et al., 1958";
"Sano et al., 1992";
"Zinszner et al., 1997";
"Asef and Najibi, 2013";
"Njiekak and Schmitt, 2019";
"Khaksar et al., 1999";
"MacBeth and Schuett, 2007";
"Miller et al., 2021";
"Nur and Simmons, 1969";
"Grochau and Gurevich,2008";
"Smith et al., 2010";
"Peacock et al., 1994";
"Best et al., 2007";
"Wei et al., 2022";
"Simmons and Brace, 1965";
"Tsuji et al., 2006";
"Wang et al., 2005";
"Manghnani et al., 1974";
"Eberhart-Phillips et al., 1989";
"Birch, 1960"];
legend(h,str);
ylim([1e-11,1e-5])
xlabel('Confining Pressure (MPa)')
ylabel('\eta (Pa^-^1)')
grid on;
box on;
set(gca,'fontsize',16);
