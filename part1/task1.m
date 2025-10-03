clear;
clc;
close all;

% Defining material properties
rho = [2810 8500 8500]; 
cp = [960 380 500]; 
k = [130 115 16.2];
names = ["Aluminum" "Brass" "Steel"];  
A = pi*0.0127^2; % m^2 --> cross-sectional area
x = linspace(0.034925,0.123825,8); % Position of thermocouples
xAxis = [0,x];

% Reading in fils to struct (credit: Jeff Glusman)

a=dir('../data/*mA');

for i=1:length(a)
    thermocouples = readmatrix(['../data/' a(i).name]);
    thermocouples = removeNaNs(thermocouples);

    b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3
    matTable(i).name = b{1};
    % {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
    v = strsplit(b{2},'V'); % volts are always in the second portion
    ampval= strsplit(b{3},'mA'); % amps are always in the third portion
    matTable(i).volts = str2num(v{1}); % convert string to number (vector)
    matTable(i).amps = str2num(ampval{1}) / 1000;
    matTable(i).rho = rho(names==matTable(i).name);
    matTable(i).cp = cp(names==matTable(i).name); 
    matTable(i).k = k(names==matTable(i).name); 
    matTable(i).ssValues = thermocouples(end,2:9);
    [p,S] = polyfit(x,matTable(i).ssValues,1);
    [matTable(i).y_fit,delta] = polyval(p,[0,x],S);
    matTable(i).slopeExp = p(1);
    matTable(i).T0Exp = matTable(i).y_fit(1);
    matTable(i).Q = matTable(i).volts * matTable(i).amps;
    matTable(i).slopeAna = matTable(i).Q / (matTable(i).k * A);
    matTable(i).y = matTable(i).T0Exp + matTable(i).slopeAna .* xAxis;
 
end

% Plotting Results

titles = ["Aluminum 25V" "Aluminum 30V" "Brass 25V" "Brass 30V" "Steel 22V"];
figure; 
for i=1:length(a) 
  plot(xAxis,matTable(i).y_fit,'Color','#FAA0A0','LineWidth',2);
  hold on; 
  grid on;
  plot(xAxis,matTable(i).y,'Color','#AF67DB','LineWidth',2);
  xlabel('Distance (m)');
  ylabel('Temperature (°C)');
  title(titles(i)); 
  legend('Experimental Solution','Analytical Solution','Location','northwest');
  hold off; 
  exportgraphics(gcf,['figures/', a(i).name, '.pdf']);
end

