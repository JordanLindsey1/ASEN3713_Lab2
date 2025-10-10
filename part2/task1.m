clc;
clear;

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
    matTable(i).H_exp = p(1);
    matTable(i).T_0 = matTable(i).y_fit(1);
    matTable(i).Q = matTable(i).volts * matTable(i).amps;
    matTable(i).H_an = matTable(i).Q / (matTable(i).k * A);
    matTable(i).y = matTable(i).T_0 + matTable(i).H_an .* xAxis;
    matTable(i).init_vals = thermocouples(1,2:9);
end

t = 1:1000;
H = matTable(1).H_an;
L = .1492;
x = .1238;
alpha = 130 / (960*2810);

summation = 0;

% sin((2*n-1)*pi*x/(2*L))

for n=1:10
    summation = summation + ((-1)^n * 8*H / (((2*n-1)*pi)^2)) .* -1^(n+1) .* exp(-((2*n-1)*pi/(2*L))^2*alpha.*t);
end

u_Al25_T8 = matTable(1).T_0 + H*x + summation;

plot(t,u_Al25_T8,'LineWidth',2,'Color','red')
xlabel("Time (t)");
ylabel("u(T8,t)");
title("Plot 1")
