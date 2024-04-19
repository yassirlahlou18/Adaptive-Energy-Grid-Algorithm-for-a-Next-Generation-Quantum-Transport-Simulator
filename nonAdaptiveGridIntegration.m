function [n_approx, rel_error] = nonAdaptiveGridIntegration(f, Emin, Emax, NE, n_true)
    DeltaE = (Emax - Emin) / NE; % Grid spacing
    E_grid = Emin:DeltaE:Emax;
    n_approx = 0;
    for i = 1:length(E_grid)-1
        E_mid = (E_grid(i) + E_grid(i+1))/2;
        n_approx = n_approx + f(E_mid) * DeltaE;
    end

    % Calculate relative error
    rel_error = abs(n_approx - n_true) / n_true;
end
