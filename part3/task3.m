% /****************************************************************/
% /*  Author: ASEN3802 GROUP 27                                   */
% /*  Title:  Th8 vs Model                                         */
% /*  Purpose: Plot thermocouple (Th8) experimental data with      */
% /*           ±TC uncertainty bounds and compare to model curve.  */
% /*  Creation Date: 2025-11-03                                    */
% /*  Revisions: Yes, All the Revisions                            */
% /****************************************************************/

clc; clear; close all;

%% --- Settings ---
data_dir = '../data/';   % path 
Tc_err   = 1.0;          % ± thermocouple uncertainty

%% --- Load data & run model ---
matTable = readData(data_dir);
u_t      = part2_task3(matTable);

%% --- color ---
expColor   = [0.0 0.2 0.6]; 
boundColor = [0.5 0.5 0.5];   
modelColor = [0.8 0.2 0.0];   

%% --- Plot material/voltage ---
for i = 1:numel(matTable)
    t       = matTable(i).t(:);
    j8      = size(matTable(i).thermo_temp, 2);  % last thermocouple index
    Th8_exp = matTable(i).thermo_temp(:, j8); % experimental data
    Th8_mod = u_t(i).data(:, j8); % model prediciton data

    figure('Color','w'); hold on; grid on;

    % Experimental Th8 with error bounds
    hExp = plot(t, Th8_exp, '-', 'Color', expColor, 'LineWidth', 2,'DisplayName', 'Experimental Th_8');
    hBoundTop = plot(t, Th8_exp + Tc_err, '--', 'Color', boundColor,'LineWidth', 1, 'DisplayName', '±TC Error Bounds');
    plot(t, Th8_exp - Tc_err, '--', 'Color', boundColor,'LineWidth', 1, 'HandleVisibility', 'off');

    % Model prediction
    hModel = plot(t, Th8_mod, '-', 'Color', modelColor, 'LineWidth', 2, ...
        'DisplayName', 'Model Prediction');

    % Labels and stuff
    xlabel('Time (s)');
    ylabel('Temperature (°C)');
    title(sprintf('%s - %.0f V', matTable(i).name, matTable(i).volts));

    legend([hExp, hBoundTop, hModel], 'Location', 'southeast', 'Box', 'off');
    set(gca, 'FontSize', 12, 'Box', 'on', 'LineWidth', 1);
end
