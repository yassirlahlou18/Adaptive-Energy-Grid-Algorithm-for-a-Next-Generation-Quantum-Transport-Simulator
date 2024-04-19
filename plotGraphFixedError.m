function plotGraphFixedError(NE_values,f,Emin,Emax,n_true)
relative_errors_non_adaptive = zeros(size(NE_values));
relative_errors_adaptive = zeros(size(NE_values));
sub_intervals_adaptive = zeros(size(NE_values));

for i = 1:length(NE_values)
    % Non-Adaptive Grid Integration
    [n_approx, rel_error_non_adaptive] = nonAdaptiveGridIntegration(f, Emin, Emax, NE_values(i), n_true);
    relative_errors_non_adaptive(i) = rel_error_non_adaptive;

    % Adaptive Grid Integration
    % Assuming adaptive method returns relative error similar or less than non-adaptive
    maxError = rel_error_non_adaptive; % Adjust as needed
    [n_adaptive, subIntervals, subIntervalSizes] = adaptiveTrapezoidal(f, Emin, Emax, maxError);
    rel_error_adaptive = abs(n_adaptive - n_true) / n_true; % Compute relative error for adaptive
    relative_errors_adaptive(i) = rel_error_adaptive;
    sub_intervals_adaptive(i) = length(subIntervalSizes); % Store the number of subintervals used by the adaptive method
end

% figure;
% plot(sub_intervals_adaptive, relative_errors_adaptive, 'r--*', 'LineWidth', 2);
% xlabel('Number of Subintervals');
% ylabel('Relative Error');
% legend('Adaptive Trapezoidal');
% title('Relative Error vs. Number of Subintervals');
% grid on;

figure;
plot(NE_values, relative_errors_non_adaptive, 'b-o', 'LineWidth', 2);
xlabel('Number of Subintervals');
ylabel('Relative Error');
legend('Non-Adaptive');
title('Relative Error vs. Number of Subintervals');
grid on;

% figure;
% plot(NE_values, sub_intervals_adaptive, 'b-o', 'LineWidth', 2);
% xlabel('Number of Subintervals Non Adaptative');
% ylabel('Number of Subintervals Adaptative Trapezoidal');
% title('Number of Subintervals Adaptative Trapezoidal vs Non Adaptative'); % For at least the same error
% grid on;

% figure;
% plot(relative_errors_non_adaptive, relative_errors_adaptive, 'b-o', 'LineWidth', 2);
% xlabel('Number of Subintervals Non Adaptative');
% ylabel('Number of Subintervals Adaptative Trapezoidal');
% title('Number of Subintervals Adaptative Trapezoidal vs Non Adaptative'); % For at least the same error

grid on;
end