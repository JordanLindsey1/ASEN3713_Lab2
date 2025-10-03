clear;
clc;
close all;

% Defining material properties
names = ["Aluminum" "Brass" "Steel"];  
rho   = [2810       8500    8500]; 
cp    = [960        380     500]; 
k     = [130        115     16.2];

A = pi*0.0127^2; % m^2 --> cross-sectional area
x = linspace(0.034925,0.123825,8); % Position of thermocouples
xAxis = linspace(0, x(end), 1000);
line = @(M, x, T0) M * x + T0;

% Reading in fils to struct (credit: Jeff Glusman)

a = dir('../data/*mA');

for i=1:length(a)
    thermocouples = readmatrix(['../data/' a(i).name]);
    thermocouples = removeNaNs(thermocouples);

    b = strsplit(a(i).name,'_'); % gives a cell array (b) that is 1x3
    matTable(i).name = b{1};
    % {'material','voltsV','ampsmA'} -- now split by 'V' and 'mA'
    v      = strsplit(b{2},'V'); % volts are always in the second portion
    ampval = strsplit(b{3},'mA'); % amps are always in the third portion

    matTable(i).volts          = str2num(v{1}); % convert string to number (vector)
    matTable(i).amps           = str2num(ampval{1}) / 1000;
    matTable(i).rho            = rho(names==matTable(i).name);
    matTable(i).cp             = cp (names==matTable(i).name); 
    matTable(i).k              = k  (names==matTable(i).name); 
    matTable(i).initial_values = thermocouples(1, 2:9);

    [p,S] = polyfit(x,matTable(i).initial_values,1);
    % [matTable(i).y_fit,delta] = polyval(p,[0,x],S);

    matTable(i).slopeExp = p(1);
    matTable(i).T0       = matTable(i).initial_values(1);
    matTable(i).Q        = matTable(i).volts * matTable(i).amps;
end

% Plotting Results

titles = ["Aluminum 25V" "Aluminum 30V" "Brass 25V" "Brass 30V" "Steel 22V"];

if ~isfolder('figures/task2') 
  mkdir figures task2;
end

for i=1:length(a) 
  figure(i); 
  hold on; 
  grid on;

  scatter(x, matTable(i).initial_values, 'Color', '#FAA0A0');
  plot(xAxis, line( matTable(i).slopeExp, xAxis, matTable(i).T0 ), 'Color', '#FAA0A0', 'LineWidth', 2);

  ylim([10 20]);
  xlabel('Distance (m)');
  ylabel('Temperature (°C)');
  title(titles(i)); 
  legend('Experimental Data', 'Linear Fit', 'Location', 'best');

  exportgraphics(gcf,['figures/task2/', a(i).name, '.pdf']);
end

