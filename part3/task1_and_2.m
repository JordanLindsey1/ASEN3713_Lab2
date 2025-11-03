clc;
clear;
close all;


%% Task 1
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

u_t = part2_task3(matTable, 1);

%% Task 2
tolerance = 9e-4;
Fo        = zeros(len, 1);
t_ss      = zeros(len, 1);
L         = .1905;

for i = 1 : len
  [data_height, data_len] = size(u_t(i).data);

  t_ss_vec = zeros(data_len, 1);
  
  for j = 1 : data_len
    grad       = gradient( u_t(i).data(:, j) );
    small_grad = (-tolerance <= grad) & (grad <= tolerance);

    for k = data_height : -1 : 1
      if ( ~small_grad(k) )
        t_ss_vec(j) = (k + 1) * 10;
        break;
      end
    end
  end

  t_ss_avg = mean(t_ss_vec);
  t_ss(i)  = t_ss_avg;
  Fo(i)    = matTable(i).alpha * t_ss_avg / L^2;

  figure(i);
  hold on;
  xline(t_ss(i));
end


