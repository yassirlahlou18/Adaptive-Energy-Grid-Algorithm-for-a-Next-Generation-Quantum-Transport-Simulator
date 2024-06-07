function plotGraphFixedNbPts(NE_values,f,Emin,Emax)
energyPoints=linspace(Emin,Emax,10);

relative_errors_adaptive = zeros(size(NE_values));

for i = 1:length(NE_values)

    % Adaptive Grid Integration
    % Assuming adaptive method returns relative error similar or less than non-adaptive
    rel_error_adaptive = adaptiveTrapezoidalFixedNbPts(f, energyPoints, NE_values(i)); % Compute relative error for adaptive
    relative_errors_adaptive(i) = rel_error_adaptive;
end

% Plotting
figure;
plot(NE_values, relative_errors_adaptive, 'r--*', 'LineWidth', 2);
xlabel('Number of Subintervals');
ylabel('Relative Error');
legend('Adaptive');
title('Relative Error vs. Number of Subintervals');
grid on;

end