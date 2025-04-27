clc,clear
data = load('cross_out.txt');
x = data(:, 1);
y = data(:, 2); 
d=sqrt((x+16).^2+(y+16).^2);
z=data(:,3);
stress = data(:, 6);


unique_d = unique(d);
unique_z = unique(z);
[D,Z] = meshgrid(unique_d, unique_z);
STRESS = reshape(stress, length(unique_z), length(unique_d));
figure;  % Create a new figure
pcolor(D,Z,STRESS);
shading interp;  % Smooths the color mapping
colorbar;  % Adds a color bar to indicate the scale of the depth-normal stress
xlabel('Distance (km)');
ylabel('Depth (km)');
load('vik.mat');
colormap(vik);
caxis([-5,5])
set(gca, 'YDir', 'reverse');
set(gca,'fontsize',16);
axis equal
xlim([0 54])
hold on;
plot([16*sqrt(2),16*sqrt(2)],[0 100],'k--','linewidth',2.5)
grid on;box on;grid minor;

