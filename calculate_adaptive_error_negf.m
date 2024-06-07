function calculate_adaptive_error_negf()
    % Parameters for NEGF
    Nx_fine = 1000;
    Nx_coarse = 10;
    dE = 0.01;
    V0 = 0;
    m = 0.1;  % Effective mass
    
    x_fine = linspace(0, 1, Nx_fine);
    x_coarse = linspace(0, 1, Nx_coarse);
    mat_fine = struct('Nx', Nx_fine, 'dE', dE, 'V', V0*ones(1, Nx_fine), 'x', x_fine, 'm', m*ones(1, Nx_fine), 'rgf', 0);

    [E_fine, TEL_fine, ~, ~, GR_fine] = solve_transport_negf({{mat_fine}});

    trueFunction = @(energy) interp1(E_fine, TEL_fine, energy, 'spline', 'extrap');

   
    energyPoints = x_coarse;

    NE_values3 = [50,60,70,80,90,100,110,120,200,250];
    maxPoints=200;

    [estimatedError2, estimatedErrorWithEstimation2, subdivisions, functionValues] = adaptiveTrapezoidalSecondDerivative2(trueFunction, energyPoints, maxPoints)
    Nx_coarse = length(subdivisions); % Number of grid points

    mat_adapt = struct('Nx', Nx_coarse, 'dE', dE,  'V', V0*ones(1, Nx_coarse), 'x', subdivisions, 'm', m*ones(1, Nx_coarse), 'rgf', 0);
   [E_adapt, TEL_adapt, ~, ~, GR_adapt] = solve_transport_negf({{mat_adapt}});


    % plotErrorEstimateAccuracy(trueFunction, energyPoints, NE_values3);

    % Extract and plot diagonal elements of the Green's function
    diag_GR_fine = diag(GR_fine);
    diag_GR_adapt = diag(GR_adapt);
    fprintf('Dimensions of data: %d x %d\n', size(GR_adapt, 1), size(GR_adapt, 2));



    figure;
    plot(E_fine, diag_GR_fine, 'r-');
    title('Diagonal Elements of Green''s Function');
    xlabel('Energy');
    ylabel('Diagonal Elements');
    legend('Optimized grid');

    figure;
    plot(E_adapt, diag_GR_adapt);
    title('Diagonal Elements of Green''s Function');
    xlabel('Energy');
    ylabel('Diagonal Elements');
    legend('Adaptive grid');

    % compare_fft_nufft(diag_GR_fine, diag_GR_adapt, E_fine, E_adapt)
end