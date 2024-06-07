function compare_fourier_transforms_single_column()
    % Params
    Nx_fine = 1000;
    dE = 0.01;
    V0 = 0;
    m = 0.1;  

    x_fine = linspace(0, 1, Nx_fine);
    mat_fine = struct('Nx', Nx_fine, 'dE', dE, 'V', V0*ones(1, Nx_fine), 'x', x_fine, 'm', m*ones(1, Nx_fine), 'rgf', 0);
    [E_fine, ~, ~, ~, GR_fine] = solve_transport_negf({{mat_fine}});
    
    % Select a specific column from the Green's function
    col_fine = floor(Nx_fine / 2); % Middle column for fine grid
    column_fine = GR_fine(:, col_fine);
    
    % Energy points for adaptive grid
    energyPoints = linspace(min(E_fine), max(E_fine), 10);  % Coarse initial grid
    maxPoints = 250;
    trueFunction = @(energy) interp1(E_fine, column_fine, energy, 'spline', 'extrap');
    [~, ~, subdivisions, ~] = adaptiveTrapezoidalSecondDerivative2(trueFunction, energyPoints, maxPoints);
    
    % Calculate non-uniform dE
    dE_adapt = diff(subdivisions);
    
    % Adaptive grid 
    Nx_adapt = length(subdivisions);
    mat_adapt = struct('Nx', Nx_adapt, 'dE', dE_adapt, 'V', V0*ones(1, Nx_adapt), 'x', subdivisions, 'm', m*ones(1, Nx_adapt), 'rgf', 0);
    [E_adapt, ~, ~, ~, GR_adapt] = solve_transport_negf({{mat_adapt}});
    
    col_adapt = floor(Nx_adapt / 2); % Middle column for adaptive grid
    column_adapt = GR_adapt(:, col_adapt);
    fprintf('Dimensions of column: %d x %d\n', size(GR_adapt, 1), size(GR_adapt, 2));

    % fft on the fine grid
    F_fine = fftshift(fft(column_fine));
    k_fine = linspace(-pi, pi, length(column_fine));
    
    % nufft on the adaptive grid
    [k_adapt, F_adapt] = perform_nufft(column_adapt, subdivisions, length(column_fine));
    
    % Plot
    figure;
    plot(k_fine, abs(F_fine), 'r-', 'DisplayName', 'Uniform FFT');
    hold on;
    plot(k_adapt, abs(F_adapt), 'b-', 'DisplayName', 'Adaptive NUFFT');
    title('Comparison of Fourier Transforms (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;
end

function [k, F] = perform_nufft(data, x, K)
    % Column vectors
    data = data(:);
    x = x(:);
    fprintf('Dimensions of x: %d x %d\n', size(x, 1), size(x, 2));
    fprintf('Dimensions of data: %d x %d\n', size(data, 1), size(data, 2));
    
    % Frequency range
    k = linspace(-pi, pi, K).';
    
    x = x * 2 * pi;  % Rescale x 
    
    F = nufft(data, x, k);
end
