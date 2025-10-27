clc;
clear;
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

thermocouples.Al25 = readmatrix(['../data/' a(1).name]);
thermocouples.Al25 = removeNaNs(thermocouples.Al25);
thermocouples.Al30 = readmatrix(['../data/' a(2).name]);
thermocouples.Al30 = removeNaNs(thermocouples.Al30);
thermocouples.Br25 = readmatrix(['../data/' a(3).name]);
thermocouples.Br25 = removeNaNs(thermocouples.Br25);
thermocouples.Br30 = readmatrix(['../data/' a(4).name]);
thermocouples.Br30 = removeNaNs(thermocouples.Br30);
thermocouples.St22 = readmatrix(['../data/' a(5).name]);
thermocouples.St22 = removeNaNs(thermocouples.St22);

for i=1:length(a)

    thermo_temp = readmatrix(['../data/' a(i).name]);
    thermo_temp = removeNaNs(thermo_temp);

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
    matTable(i).ssValues = thermo_temp(end,2:9);
    [p,S] = polyfit(x,matTable(i).ssValues,1);
    [matTable(i).y_fit,delta] = polyval(p,[0,x],S);
    matTable(i).H_exp = p(1);
    matTable(i).T_0 = matTable(i).y_fit(1);
    matTable(i).Q = matTable(i).volts * matTable(i).amps;
    matTable(i).H_an = matTable(i).Q / (matTable(i).k * A);
    matTable(i).y = matTable(i).T_0 + matTable(i).H_an .* xAxis;
    matTable(i).init_vals = thermo_temp(1,2:9);

    [p, S] = polyfit(x, matTable(i).init_vals, 1);
    matTable(i).M_exp = p(1);

end

%% MATERIAL 1

H = matTable(1).H_exp;
L = .1905;
alpha = 130 / (960*2810);
t = (0:(length(thermocouples.Al25) - 1)) .* 10;
u_t = zeros(length(thermocouples.Al25),8);
M = matTable(1).M_exp;

for i = 0:7
    x = 0.0349 + 0.0127*i;

    summation = 0;
    for n=1:10
        lambda = (2*n-1)*pi/(2*L);
        bn = 2/L * (M-H) * (sin(lambda*L)-lambda*L*cos(lambda*L))/(lambda^2);
    
        summation = summation + bn * sin(lambda*x) .* exp(-lambda^2*alpha.*t);
    end
    
    u_t(:,i+1) = matTable(1).T_0 + H*x + summation;
end



figure(1);
hold on;
for i = 1:8
    plot((0:(length(u_t) - 1)).*10,u_t(:,i),linewidth=2,color="black")
    plot((0:(length(u_t) - 1)).*10,thermocouples.Al25(:,i+1),linewidth=2,color="red")
end
hold off;

xlabel("Time (s)");
ylabel("Temperature (C)");
title("Aluminum at 25 volts");
legend(["Modeled Data","Experimental Data"],Location="southeast");


%% MATERIAL 2

H = matTable(2).H_exp;
t = (0:(length(thermocouples.Al30) - 1)) .* 10;
u_t = zeros(length(thermocouples.Al30),8);
M = matTable(2).M_exp;

for i = 0:7
    x = 0.0349 + 0.0127*i;

    summation = 0;
    for n=1:10
        lambda = (2*n-1)*pi/(2*L);
        bn = 2/L * (M-H) * (sin(lambda*L)-lambda*L*cos(lambda*L))/(lambda^2);
    
        summation = summation + bn * sin(lambda*x) .* exp(-lambda^2*alpha.*t);
    end
    
    u_t(:,i+1) = matTable(2).T_0 + H*x + summation;
end

figure(2);
hold on;
for i = 1:8
    plot((0:(length(u_t) - 1)).*10,u_t(:,i),linewidth=2,color="black")
    plot((0:(length(u_t) - 1)).*10,thermocouples.Al30(:,i+1),linewidth=2,color="red")
end
hold off;

xlabel("Time (s)");
ylabel("Temperature (C)");
title("Aluminum at 30 volts");
legend(["Modeled Data","Experimental Data"],Location="southeast");


%% MATERIAL 3

H = matTable(3).H_exp;
t = (0:(length(thermocouples.Br25) - 1)) .* 10;
u_t = zeros(length(thermocouples.Br25),8);
M = matTable(3).M_exp;

for i = 0:7
    x = 0.0349 + 0.0127*i;

    summation = 0;
    for n=1:10
        lambda = (2*n-1)*pi/(2*L);
        bn = 2/L * (M-H) * (sin(lambda*L)-lambda*L*cos(lambda*L))/(lambda^2);
    
        summation = summation + bn * sin(lambda*x) .* exp(-lambda^2*alpha.*t);
    end
    
    u_t(:,i+1) = matTable(3).T_0 + H*x + summation;
end

figure(3);
hold on;
for i = 1:8
    plot((0:(length(u_t) - 1)).*10,u_t(:,i),linewidth=2,color="black")
    plot((0:(length(u_t) - 1)).*10,thermocouples.Br25(:,i+1),linewidth=2,color="red")
end
hold off;

xlabel("Time (s)");
ylabel("Temperature (C)");
title("Brass at 25 volts");
legend(["Modeled Data","Experimental Data"],Location="southeast");


%% MATERIAL 4

H = matTable(4).H_exp;
t = (0:(length(thermocouples.Br30) - 1)) .* 10;
u_t = zeros(length(thermocouples.Br30),8);
M = matTable(4).M_exp;

for i = 0:7
    x = 0.0349 + 0.0127*i;

    summation = 0;
    for n=1:10
        lambda = (2*n-1)*pi/(2*L);
        bn = 2/L * (M-H) * (sin(lambda*L)-lambda*L*cos(lambda*L))/(lambda^2);
    
        summation = summation + bn * sin(lambda*x) .* exp(-lambda^2*alpha.*t);
    end
    
    u_t(:,i+1) = matTable(4).T_0 + H*x + summation;
end

figure(4);
hold on;
for i = 1:8
    plot((0:(length(u_t) - 1)).*10,u_t(:,i),linewidth=2,color="black")
    plot((0:(length(u_t) - 1)).*10,thermocouples.Br30(:,i+1),linewidth=2,color="red")
end
hold off;

xlabel("Time (s)");
ylabel("Temperature (C)");
title("Brass at 30 volts");
legend(["Modeled Data","Experimental Data"],Location="southeast");


%% MATERIAL 5

H = matTable(5).H_exp;
t = (0:(length(thermocouples.St22) - 1)) .* 10;
u_t = zeros(length(thermocouples.St22),8);
M = matTable(5).M_exp;

for i = 0:7
    x = 0.0349 + 0.0127*i;

    summation = 0;
    for n=1:10
        lambda = (2*n-1)*pi/(2*L);
        bn = 2/L * (M-H) * (sin(lambda*L)-lambda*L*cos(lambda*L))/(lambda^2);
    
        summation = summation + bn * sin(lambda*x) .* exp(-lambda^2*alpha.*t);
    end
    
    u_t(:,i+1) = matTable(5).T_0 + H*x + summation;
end

figure(5);
hold on;
for i = 1:8
    plot((0:(length(u_t) - 1)).*10,u_t(:,i),linewidth=2,color="black")
    plot((0:(length(u_t) - 1)).*10,thermocouples.St22(:,i+1),linewidth=2,color="red")
end
hold off;

xlabel("Time (s)");
ylabel("Temperature (C)");
title("Steel at 22 volts");
legend(["Modeled Data","Experimental Data"],Location="southeast");


