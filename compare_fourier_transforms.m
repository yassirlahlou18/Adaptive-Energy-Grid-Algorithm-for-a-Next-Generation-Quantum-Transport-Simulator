function compare_fourier_transforms()
    
    structure = 'dqw';
    mat_fine = material(structure, 1);

    % Modify energy vector
    NE_fine = 1001;
    dE = 0.0015;
    %dE = mat_fine.dE;
    mat_fine.dE = dE*ones(1, NE_fine-1);
    
    [E_fine, ~, ~, ~, GR_fine] = solve_transport_negf({{mat_fine}});

    fprintf('Dimensions of E_fine: %d x %d\n', size(E_fine, 1), size(E_fine, 2));
    fprintf('Dimensions of GR_fine: %d x %d\n', size(GR_fine, 1), size(GR_fine, 2));


    % Select a specific column from the Green's function, here the first column
    column_fine = GR_fine(:, 1);

    fprintf('Dimensions of column_fine: %d x %d\n', size(column_fine, 1), size(column_fine, 2));

    % Energy points for adaptive grid
    energyPoints = linspace(min(E_fine), max(E_fine), 10);  % Coarse initial grid
    maxPoints = NE_fine;
    trueFunction = @(energy) interp1(E_fine, column_fine, energy, 'spline', 'extrap');
    [~, ~, subdivisions, ~] = adaptiveTrapezoidalSecondDerivative2(trueFunction, energyPoints, maxPoints-1);
    
    % Calculate non-uniform dE
    dE_adapt = diff(subdivisions);
        
    % Adaptive grid 
    mat_adapt = material(structure, 1);
    mat_adapt.dE = dE_adapt;
    [E_adapt, ~, ~, ~, GR_adapt] = solve_transport_negf({{mat_adapt}});
    
    % Important to select the same column as before for comparison
    column_adapt = GR_adapt(:, 1);
    fprintf('Dimensions of GR_adapt: %d x %d\n', size(GR_adapt, 1), size(GR_adapt, 2));

    % fft on the fine grid
    % nufft on the adaptive grid
    [k_fine, F_fine] = perform_nufft_selectgrid(column_fine, E_fine,  length(column_fine), dE);
    
    % nufft on the adaptive grid
    [k_adapt, F_adapt] = perform_nufft_selectgrid(column_adapt, E_adapt, 3 * length(column_adapt), min(dE_adapt));

    N_F_fine=abs(F_fine);
    N_F_fine=N_F_fine-min(N_F_fine);

    N_F_Adapt=abs(F_adapt);
    N_F_Adapt=N_F_Adapt-min(N_F_Adapt);

    %scale=max(N_F_Adapt)/max(N_F_fine);
    %N_F_fine=N_F_fine*scale;
    
    %scale_k = max(k_adapt)/max(k_fine);
    
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


    % Plot
    figure;
    h = plot(E_fine, column_fine, 'g-.', 'DisplayName', 'Original data (G^r)');
    set(h,{'LineWidth'},{5});
    hold on;
    % Whenever you see nufft(F, k, -E), that is the inverse transform from
    % k-space back to energy.
    g = plot(E_fine, nufft(F_fine, k_fine, -E_fine)/length(k_fine), 'r--', 'DisplayName', 'Inverse Uniform NUFFT');
    set(g, {'LineWidth'}, {3});
    plot(E_adapt, nufft(F_adapt, k_adapt, -E_adapt)/length(k_adapt), 'b-', 'DisplayName', 'Inverse Adaptive NUFFT');
    title('Comparison of Recovery (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;

    [GL_fine, GG_fine] = solve_GLGG(E_fine, mat_fine);
    energyPoints2 = linspace(min(E_fine), max(E_fine), 10);  % Coarse initial grid
    column_GL_fine = imag(GL_fine(:, 1));
    column_GG_fine = imag(GG_fine(:, 1));
    
    % Use adaptiveTrapezoidalSecondDerivative2 on column_GL_fine
    trueFunction2 = @(energy) interp1(E_fine, column_GL_fine, energy, 'spline', 'extrap');
    [~, ~, subdivisions2, ~] = adaptiveTrapezoidalSecondDerivative2(trueFunction2, energyPoints2, maxPoints-1);
    
    % Calculate non-uniform dE
    dE_adapt2 = diff(subdivisions2);
    
    % Adaptive grid
    E_adapt2 = subdivisions2;  % Energy vector
    mat_adapt2 = material(structure, 1);
    mat_adapt2.dE = dE_adapt2;

    % Calculate Lesser and Greater Green's functions for the adaptive grid
    [GL_adapt, GG_adapt] = solve_GLGG(E_adapt2, mat_adapt2);
    column_GL_adapt = imag(GL_adapt(:,1));
    column_GG_adapt = imag(GG_adapt(:,1));
    
     % Perform NUFFT on the fine grid
    [k_fine, F_GL_fine] = perform_nufft_selectgrid(column_GL_fine, E_fine, length(column_GL_fine), dE);
    [~, F_GG_fine] = perform_nufft_selectgrid(column_GG_fine, E_fine, length(column_GG_fine), dE);
    
    % Perform NUFFT on the adaptive grid
    % increase length of k space to avoid copies in energy domain, min(dE)
    % for good NUFFT representation.
    [k_adapt, F_GL_adapt] = perform_nufft_selectgrid(column_GL_adapt, E_adapt2, 6 * length(column_GL_adapt), min(dE_adapt2));
    [~, F_GG_adapt] = perform_nufft_selectgrid(column_GG_adapt, E_adapt2, 6 * length(column_GG_adapt), min(dE_adapt2));
    
    %Perform NUFFT for adaptive grid onto the same fine grid
    [~, F_GL_adapt_smallk] = perform_nufft_selectgrid(column_GL_adapt, E_adapt2, 1 * length(column_GL_adapt), dE);
    [~, F_GG_adapt_smallk] = perform_nufft_selectgrid(column_GG_adapt, E_adapt2, 1 * length(column_GG_adapt), dE);

    % Plot
    figure;
    plot(k_fine, F_GL_fine, 'r-', 'DisplayName', 'Uniform NUFFT');
    hold on;
    plot(k_adapt, F_GL_adapt, 'b-', 'DisplayName', 'Adaptive NUFFT');
    title('Comparison of Fourier Transforms Lesser (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;

    % Plot
    figure;
    h = plot(E_fine, column_GL_fine, 'g-.', 'DisplayName', 'Original data (G^L)');
    set(h,{'LineWidth'},{5});
    hold on;
    g = plot(E_fine, nufft(F_GL_fine, k_fine, -E_fine)/length(k_fine), 'r--', 'DisplayName', 'Inverse Uniform NUFFT');
    set(g, {'LineWidth'}, {3});
    plot(E_adapt2, nufft(F_GL_adapt, k_adapt, -E_adapt2)/length(k_adapt), 'b-', 'DisplayName', 'Inverse Adaptive NUFFT');
    title('Comparison of Recovery (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;

    % Plot
    figure;
    plot(k_fine, F_GG_fine, 'r-', 'DisplayName', 'Uniform NUFFT');
    hold on;
    plot(k_adapt, F_GG_adapt, 'b-', 'DisplayName', 'Adaptive NUFFT');
    title('Comparison of Fourier Transforms Lesser (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;

    % Plot
    figure;
    h = plot(E_fine, column_GG_fine, 'g-.', 'DisplayName', 'Original data (G^G)');
    set(h,{'LineWidth'},{5});
    hold on;
    g = plot(E_fine, nufft(F_GG_fine, k_fine, -E_fine)/length(k_fine), 'r--', 'DisplayName', 'Inverse Uniform NUFFT');
    set(g, {'LineWidth'}, {3});
    plot(E_adapt2, nufft(F_GG_adapt, k_adapt, -E_adapt2)/length(k_adapt), 'b-', 'DisplayName', 'Inverse Adaptive NUFFT');
    title('Comparison of Recovery (Single Column)');
    xlabel('Frequency');
    ylabel('Magnitude');
    legend;
    hold off;
    
    % Multiply the NUFFT results and perform inverse FFT on the product
    F_product_fine = nufft(F_GL_fine .* F_GG_fine, k_fine, -E_fine)/length(k_fine);
    fprintf('Dimensions of F_product_fine: %d x %d\n', size(F_product_fine, 1), size(F_product_fine, 2));

    F_product_adapt_smallk = nufft(F_GL_adapt_smallk .* F_GG_adapt_smallk, k_fine, -E_adapt2)/length(k_fine);
    fprintf('Dimensions of F_product_adapt: %d x %d\n', size(F_product_adapt_smallk, 1), size(F_product_adapt_smallk, 2));

    F_product_adapt = nufft(F_GL_adapt .* F_GG_adapt, k_adapt, -E_adapt2)/length(k_adapt);
    fprintf('Dimensions of F_product_adapt: %d x %d\n', size(F_product_adapt, 1), size(F_product_adapt, 2));

    %column_GL_interpolated = nufft(F_GL_adapt, k_adapt, -E_fine)/length(k_adapt);
    %column_GG_interpolated = nufft(F_GG_adapt, k_adapt, -E_fine)/length(k_adapt);

    column_GL_interpolated = interp1(E_adapt2, column_GL_adapt, E_fine);
    column_GG_interpolated = interp1(E_adapt2, column_GG_adapt, E_fine);


    [k_fine_int, F_GL_fine_int] = perform_nufft_selectgrid(column_GL_interpolated, E_fine, length(column_GL_fine), dE);
    [~, F_GG_fine_int] = perform_nufft_selectgrid(column_GG_interpolated, E_fine, length(column_GG_fine), dE);

    F_product_adapt_int = nufft(F_GL_fine_int .* F_GG_fine_int, k_fine_int, -E_fine)/length(k_fine_int);

    convolution_product = cconv(column_GG_fine, column_GL_fine, length(E_fine));
    
    
    
    % Plot the result for both pairs on the same plot
    figure;
    h = plot(E_fine, abs(convolution_product), 'g-.', 'DisplayName', 'Direct Convolution');
    set(h, {'LineWidth'}, {3});
    hold on;
    plot(E_fine, abs(F_product_fine(:,1)), 'r-', 'DisplayName', 'Uniform Grid FFT');
    g = plot(E_adapt2, abs(F_product_adapt_smallk(:,1)), 'b:', 'DisplayName', 'Adaptive Grid NUFFT A1');
    set(g, {'LineWidth'}, {3});
    n = plot(E_adapt2, abs(F_product_adapt(:,1)), 'k--', 'DisplayName', 'Adaptive Grid NUFFT A2');
    set(n, {'LineWidth'}, {3});
    m = plot(E_fine, abs(F_product_adapt_int(:,1)), 'm-.', 'DisplayName', 'Interpolation FFT');
    set(m, {'LineWidth'}, {3});
    title('Comparison of Fourier Transforms of F^-1(GL(t) * GG(t)) = P(E)');
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
    
    %F = nufft(data, x,k);
    F = nufft(data, x);
end


function [k, F] = perform_nufft_selectgrid(data, x, K, dE)
    % Column vectors
    data = data(:);
    x = x(:);
    fprintf('Dimensions of x: %d x %d\n', size(x, 1), size(x, 2));
    fprintf('Dimensions of data: %d x %d\n', size(data, 1), size(data, 2));
    
    % Frequency range
    k = linspace(0, 1/(dE) - 1/(dE)/K, K).';
    
    %x = x * 2 * pi;  % Rescale x 
    
    %F = nufft(data, x,k);
    F = nufft(data, x, k);
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

