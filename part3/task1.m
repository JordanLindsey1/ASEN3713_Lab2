clc;
clear;
close all;

range        = 2;
alpha_num    = 500;
error_matrix = @(actual, observed) sqrt( mean( rmse(actual, observed).^2 ) );
matTable     = readData('../data/');
len          = length(matTable);
u_t_before   = part2_task3(matTable);
error_before = zeros(len, 1);
min_errors   = zeros(len, 1);

for i = 1 : len
  error_before(i) = error_matrix(matTable(i).thermo_temp, u_t_before(i).data);
end

for i = 1 : len
  alpha      = linspace(matTable(i).alpha * (1 - range), matTable(i).alpha * (1 + range), alpha_num)';
  alpha_best = matTable(i).alpha;
  min_error  = Inf;

  for j = 1 : alpha_num
    matTable(i).alpha = alpha(j);

    u_t   = part2_task3(matTable(i)); 
    error = error_matrix(matTable(i).thermo_temp, u_t.data);

    if (error < min_error) 
      alpha_best = alpha(j);
      min_error = error;
    end
  end

  matTable(i).alpha = alpha_best;
  min_errors(i)     = min_error;
end

part2_task3(matTable, 1);
