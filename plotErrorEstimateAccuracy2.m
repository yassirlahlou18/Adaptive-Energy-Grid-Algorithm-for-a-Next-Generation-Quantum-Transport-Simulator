function plotErrorEstimateAccuracy2(f, energyPoints, k)
    % Initialize arrays
    n_true = integral(f, energyPoints(1), energyPoints(end));
    richardsonRelativeErrors = zeros(size(k));
    secondDerivativeRelativeErrors = zeros(size(k));
    quadratureRelativeErrors = zeros(size(k));  % For the new method
    trueRelativeErrors = zeros(size(k));
    NonAdaptiveError = zeros(size(k));
    rel_error_non_adaptive = zeros(size(k));

    richardsonEstimatedErrors = zeros(size(k));
    secondDerivativeEstimatedErrors = zeros(size(k));
    quadratureEstimatedErrors = zeros(size(k));  % For the new method
    secondDerivativeEstimatedErrorWithEstimation = zeros(size(k));
    quadratureEstimatedErrorWithEstimation = zeros(size(k));
    RichardsonEstimatedErrorWithEstimation = zeros(size(k));
    n_approx = zeros(size(k));

    % Initialize arrays for timing
    timeNonAdaptive = zeros(size(k));
    timeRichardson = zeros(size(k));
    timeSecondDerivative = zeros(size(k));
    timeQuadrature = zeros(size(k));

    % Loop over different values of maximum points
    for i = 1:length(k)
        maxPoints = k(i);  % Current max points from array k

        % Timing non-adaptive method
        tic;
        [n_approx(i), rel_error_non_adaptive(i)] = nonAdaptiveGridIntegration(f, energyPoints(1), energyPoints(end), maxPoints, n_true);
        timeNonAdaptive(i) = toc;

        % Timing Richardson method
        tic;
        [richardsonEstimatedErrors(i), RichardsonEstimatedErrorWithEstimation(i)] = adaptiveTrapezoidalRichardson(f, energyPoints, maxPoints);
        timeRichardson(i) = toc;

        % Timing Second Derivative method
        tic;
        [secondDerivativeEstimatedErrors(i), secondDerivativeEstimatedErrorWithEstimation(i)] = adaptiveTrapezoidalSecondDerivative(f, energyPoints, maxPoints);
        timeSecondDerivative(i) = toc;

        % Timing Quadrature method
        tic;
        [quadratureEstimatedErrors(i), quadratureEstimatedErrorWithEstimation(i)] = adaptiveQuadratureSimpson(f, energyPoints, maxPoints);  % Assuming this function also returns estimated error
        timeQuadrature(i) = toc;

        % Assuming this function returns the true error
        trueRelativeErrors(i) = adaptiveTrapezoidalFixedNbPts(f, energyPoints, maxPoints);  

        richardsonRelativeErrors(i) = abs(richardsonEstimatedErrors(i) - trueRelativeErrors(i)) / trueRelativeErrors(i);
        secondDerivativeRelativeErrors(i) = abs(secondDerivativeEstimatedErrors(i) - trueRelativeErrors(i)) / trueRelativeErrors(i);
        quadratureRelativeErrors(i) = (quadratureEstimatedErrors(i) - trueRelativeErrors(i)) / trueRelativeErrors(i);
        NonAdaptiveError(i) = abs(n_approx(i) - n_true) / n_true;
    end
    
    % Plotting the results: Relative Errors of estimations
    figure;
    % plot(k, richardsonRelativeErrors, '-o', 'DisplayName', 'Richardson Estimation Relative Error to True error');
    hold on;
    plot(k, secondDerivativeRelativeErrors, '-x', 'DisplayName', 'Second Derivative Estimation Relative Error to True error');
    plot(k, quadratureRelativeErrors, '-s', 'DisplayName', 'Quadrature Estimation Relative Error to True error');
    %plot(k, trueRelativeErrors, '-^', 'DisplayName', 'True Relative Error');
    xlabel('Number of Points');
    ylabel('Relative Error between Estimated and True Error');
    title('Error Estimates Accuracy Comparison');
    legend show;
    grid on;
    
    % Plotting the results: Estimated Errors
    figure;
    semilogy(k, NonAdaptiveError, '^-', 'Color', '#0072BD', 'DisplayName', 'Non-Adaptive True Error'); % Blue solid
    hold on;
    % semilogy(k, richardsonEstimatedErrors, 'o-', 'Color', '#008000', 'DisplayName', 'Richardson Estimated Error with true integral');
    % semilogy(k, secondDerivativeEstimatedErrors, '^-', 'Color', '#EDB120', 'DisplayName', 'Second Derivative Estimated Error with true integral');
    % semilogy(k, quadratureEstimatedErrors, 'o-', 'Color', '#D95319', 'DisplayName', 'Quadrature Estimated Error with true integral');
    semilogy(k, trueRelativeErrors, '-^', 'DisplayName', 'Perfect Adaptive Error');
    % semilogy(k, RichardsonEstimatedErrorWithEstimation, '-s', 'DisplayName', 'Richardson Estimated Error');
    % semilogy(k, secondDerivativeEstimatedErrorWithEstimation, '^--', 'Color', '#EDB120','DisplayName', 'Second Derivative Estimated Error');
    % semilogy(k, quadratureEstimatedErrorWithEstimation, 'o--', 'Color', '#D95319' ,'DisplayName', 'Quadrature Estimated Error');
    xlabel('Number of Points');
    ylabel('Error');
    title('Error Comparison');
    legend show;
    grid on;

    % Plotting the timing results
    figure;
    % plot(k, timeNonAdaptive, '^-', 'DisplayName', 'Non-Adaptive Time');
    plot(k, timeRichardson, 'o-', 'Color', '#008000', 'DisplayName', 'Richardson Time');
    hold on;
    plot(k, timeSecondDerivative, '^-', 'Color', '#EDB120', 'DisplayName', 'Second Derivative Time');
    plot(k, timeQuadrature, 'o-', 'Color', '#D95319', 'DisplayName', 'Quadrature Time');
    xlabel('Number of Points');
    ylabel('Time (seconds)');
    title('Computation Time Comparison');
    legend show;
    grid on;


end
