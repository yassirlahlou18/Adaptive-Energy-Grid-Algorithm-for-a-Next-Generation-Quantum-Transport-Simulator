function compare_fourier_transforms2()
    
    % Params
    Nx_fine = 100;
    dE = 0.01;
    V0 = 0;
    m = 0.1;
    Temp=300;
    Ne_fine=250;
   



    x_fine = linspace(0, 0.3, Nx_fine);

    mat_fine = struct('Nx', Nx_fine, 'dE', dE*ones(1, Ne_fine-1), 'V', V0*ones(1, Nx_fine),'Temp', Temp,'x', x_fine, 'm', m*ones(1, Nx_fine), 'rgf', 1);
    [E_fine, ~, ~, ~, GR_fine] = solve_transport_negf({{mat_fine}});
        fprintf('Dimensions of E_fine: %d x %d\n', size(E_fine, 1), size(E_fine, 2));
                fprintf('Dimensions of GR_fine: %d x %d\n', size(GR_fine, 1), size(GR_fine, 2));


    
    % Select a specific column from the Green's function
    column_fine = GR_fine(:, 1);
                    fprintf('Dimensions of column_fine: %d x %d\n', size(column_fine, 1), size(column_fine, 2));

    % Energy points for adaptive grid
    energyPoints = linspace(min(x_fine), max(x_fine), 10)  % Coarse initial grid
    maxPoints = Nx_fine;
    trueFunction = @(energy) interp1(E_fine, column_fine, energy, 'spline', 'extrap');
    [~, ~, subdivisions, ~] = adaptiveTrapezoidalSecondDerivative2(trueFunction, energyPoints, maxPoints);
    
    % Calculate non-uniform dE
    dE_adapt = diff(subdivisions);
        
    % Adaptive grid 
    Nx_adapt = length(subdivisions)
    Nx_adaspt = length(dE_adapt)
    mat_adapt = struct('Nx', Nx_fine, 'dE', dE_adapt, 'V', V0*ones(1, Nx_adapt), 'Temp', Temp,'x', subdivisions, 'm', m*ones(1, Nx_adapt), 'rgf', 1);
    [E_adapt, ~, ~, ~, GR_adapt] = solve_transport_negf({{mat_adapt}});
    
    column_adapt = GR_adapt(:, 1);
    fprintf('Dimensions of column: %d x %d\n', size(GR_adapt, 1), size(GR_adapt, 2));

    % fft on the fine grid
  % nufft on the adaptive grid
    [k_fine, F_fine] = perform_nufft(column_fine, E_fine, length(column_fine));
    
    % nufft on the adaptive grid
    [k_adapt, F_adapt] = perform_nufft(column_adapt, subdivisions, length(column_adapt));
    N_F_fine=abs(F_fine);
    N_F_fine=N_F_fine-min(N_F_fine);
    



    N_F_Adapt=abs(F_adapt);
    N_F_Adapt=N_F_Adapt-min(N_F_Adapt);
    scale=max(N_F_Adapt)/max(N_F_fine);
    N_F_fine=N_F_fine*scale;

    
    % Plot
    figure;
    plot(k_fine, N_F_fine, 'r-', 'DisplayName', 'Uniform NUFFT');
    hold on;
    plot(k_adapt, N_F_Adapt, 'b-', 'DisplayName', 'Adaptive NUFFT');
    title('Comparison of Fourier Transforms (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;

    [GL_fine, GG_fine] = solve_GLGG(E_fine, mat_fine);
    energyPoints2 = linspace(min(E_fine), max(E_fine), 10);  % Coarse initial grid
    column_GL_fine = diag(GL_fine);
    column_GG_fine = diag(GG_fine);
    
    % Use adaptiveTrapezoidalSecondDerivative2 on column_GL_fine
    trueFunction2 = @(energy) interp1(E_fine, column_GL_fine, energy, 'spline', 'extrap');
    [~, ~, subdivisions2, ~] = adaptiveTrapezoidalSecondDerivative2(trueFunction2, energyPoints2, maxPoints);
    
    % Calculate non-uniform dE
    dE_adapt2 = diff(subdivisions2);
    
    % Adaptive grid
    Nx_adapt2 = length(subdivisions2);
    E_adapt2 = subdivisions2;  % Energy vector
    mat_adapt2 = struct('Nx', Nx_adapt2, 'dE', dE_adapt2, 'V', V0*ones(1, Nx_adapt2), 'Temp', Temp,'x', subdivisions2, 'm', m*ones(1, Nx_adapt2), 'rgf', 1);

    
    % Calculate Lesser and Greater Green's functions for the adaptive grid
    [GL_adapt, GG_adapt] = solve_GLGG(E_adapt2, mat_adapt2);
    column_GL_adapt = diag(GL_adapt);
    column_GG_adapt = diag(GG_adapt);
    
     % Perform NUFFT on the fine grid
    [k_fine, F_GL_fine] = perform_nufft2(GL_fine, x_fine, length(column_GL_fine));
    [~, F_GG_fine] = perform_nufft2(GG_fine, x_fine, length(column_GG_fine));
    
    % Perform NUFFT on the adaptive grid
    [k_adapt, F_GL_adapt] = perform_nufft2(GL_adapt, subdivisions, length(column_GL_adapt));
    [~, F_GG_adapt] = perform_nufft2(GG_adapt, subdivisions, length(column_GG_adapt));
    
    % Multiply the NUFFT results and perform inverse FFT on the product
    F_product_fine = ifft(F_GL_fine * F_GG_fine);
    fprintf('Dimensions of F_product_fine: %d x %d\n', size(F_product_fine, 1), size(F_product_fine, 2));

    F_product_adapt = ifft(F_GL_adapt * F_GG_adapt);
    fprintf('Dimensions of F_product_adapt: %d x %d\n', size(F_product_adapt, 1), size(F_product_adapt, 2));

    
    % Plot the result for both pairs on the same plot
    figure;
    plot(k_fine, abs(F_product_fine(:,1)), 'r-', 'DisplayName', 'Uniform Grid');
    hold on;
    plot(k_adapt, abs(F_product_adapt(:,1)), 'b-', 'DisplayName', 'Adaptive Grid');
    title('Comparison of Fourier Transforms of GL * GG');
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
    
    F = nufft(data, x,k);
end

function [k, F] = perform_nufft2(data, x, K)
    % Column vectors

    fprintf('Dimensions of x: %d x %d\n', size(x, 1), size(x, 2));
    fprintf('Dimensions of data: %d x %d\n', size(data, 1), size(data, 2));
    
    % Frequency range
    k = linspace(-pi, pi, K).';
    
    x = x * 2 * pi;  % Rescale x 
    
    F = nufft(data, x,k);
end

