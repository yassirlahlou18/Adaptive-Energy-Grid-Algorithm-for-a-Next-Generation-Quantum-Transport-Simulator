function plotAdaptTransmission()
    % Parameters
    Nx_fine = 10000;
    Nx_coarse = 10;
    dE = 0.01;
    V0 = 0;
    m = 0.1;
    
    % Generate material structures
    x_fine = linspace(0, 1, Nx_fine);
    x_coarse = linspace(0, 1, Nx_coarse);
    mat_fine = struct('Nx', Nx_fine, 'dE', dE, 'V', V0*ones(1, Nx_fine), 'x', x_fine, 'm', m*ones(1, Nx_fine));
    % mat_coarse = struct('Nx', Nx_coarse, 'dE', dE, 'V', V0*ones(1, Nx_coarse), 'x', x_coarse, 'm', m*ones(1, Nx_coarse));
    
    % Obtain transmission using solve_transport for fine and coarse grids
    [E_fine, TEL_fine, ~, ~, ~] = solve_transport({{mat_fine}});

    % Interpolate fine grid transmission to get a function for true values
    trueFunction = @(energy) interp1(E_fine, TEL_fine, energy, 'spline', 'extrap');
    

    % Define the energy points for coarse grid
    energyPoints = x_coarse;

    % Apply adaptive trapezoidal rule to coarse grid
    NE_values3 = [50,60,70,80,90,100,110,120,200,250,300];

    plotErrorEstimateAccuracy(trueFunction, energyPoints, NE_values3);

end
