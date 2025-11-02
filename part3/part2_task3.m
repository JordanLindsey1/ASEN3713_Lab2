function [u_t] = part2_task3(matTable, do_plot)
    arguments
        matTable struct;
        do_plot logical = 0;
    end

    L        = .1905;
    x        = linspace(0.034925,0.123825,8); % Position of thermocouples
    len      = length(matTable);
    u_t(len) = struct( ...
        'name', '',    ...
        'data', []     ...
    );

    for i = 1 : len
        H           = matTable(i).H_exp;
        u_t(i).name = strcat( matTable(i).name, " ", num2str(matTable(i).volts), "V" );

        for j = 1 : 8
            summation = 0;
            for n = 1 : 10
                lambda = (2 * n - 1) * pi / (2 * L);
                bn     = (-1)^n * (8 * H * L) / ((2 * n - 1) * pi)^2;
            
                summation = summation + bn * sin(lambda * x(j)) * exp(-lambda^2 * matTable(i).alpha * matTable(i).t);
            end
            
            u_t(i).data(:, j) = matTable(i).T_0 + (H * x(j)) + summation;
        end

        if (do_plot)
            figure(i);
            hold on;

            for j = 1 : 8
                plot(matTable(i).t, u_t(i).data(:, j), linewidth=2, color="black");
                plot(matTable(i).t, matTable(i).thermo_temp(:, j), linewidth=2, color="red")
            end

            xlabel("Time (s)");
            ylabel("Temperature (C)");
            title( strcat(matTable(i).name, " at ", num2str(matTable(i).volts ), " volts") );
            legend("Modeled Data", "Experimental Data" , "Location", "southeast");
        end
    end

end
