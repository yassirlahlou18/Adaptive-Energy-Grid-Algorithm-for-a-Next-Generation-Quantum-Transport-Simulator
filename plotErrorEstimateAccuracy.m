function plotErrorEstimateAccuracy(f, Emin, Emax, k)
    % Initialize arrays
    richardsonRelativeErrors = zeros(size(k));
    secondDerivativeRelativeErrors = zeros(size(k));
    quadratureRelativeErrors = zeros(size(k));  % For the new method
    trueRelativeErrors = zeros(size(k));

    richardsonEstimatedErrors = zeros(size(k));
    secondDerivativeEstimatedErrors = zeros(size(k));
    quadratureEstimatedErrors = zeros(size(k));  % For the new method

    % Loop over different values of maximum points
    for i = 1:length(k)
        maxPoints = k(i);  % Current max points from array k

        % Using existing functions to calculate relative and estimated errors
        [richardsonEstimatedErrors(i), ~] = adaptiveTrapezoidalRichardson(f, Emin, Emax, maxPoints);
        [secondDerivativeEstimatedErrors(i), ~] = adaptiveTrapezoidalSecondDerivative(f, Emin, Emax, maxPoints);
        [quadratureEstimatedErrors(i), ~] = adaptiveQuadratureSimpson(f, Emin, Emax, maxPoints);  % Assuming this function also returns estimated error
        trueRelativeErrors(i) = adaptiveTrapezoidalFixedNbPts(f, Emin, Emax, maxPoints);  % Assuming this function returns the true error

        richardsonRelativeErrors(i) = (richardsonEstimatedErrors(i) - trueRelativeErrors(i)) / trueRelativeErrors(i);
        secondDerivativeRelativeErrors(i) = (secondDerivativeEstimatedErrors(i) - trueRelativeErrors(i)) / trueRelativeErrors(i);
        quadratureRelativeErrors(i) = (quadratureEstimatedErrors(i) - trueRelativeErrors(i)) / trueRelativeErrors(i);
    end
    
    % Plotting the results: Relative Errors of estimations
    figure;
    plot(k, richardsonRelativeErrors, '-o', 'DisplayName', 'Richardson Estimation Relative Error to True error');
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
    plot(k, richardsonEstimatedErrors, '-o', 'DisplayName', 'Richardson Estimated Error');
    hold on;
    plot(k, secondDerivativeEstimatedErrors, '-x', 'DisplayName', 'Second Derivative Estimated Error');
    plot(k, quadratureEstimatedErrors, '-s', 'DisplayName', 'Quadrature Estimated Error');
    plot(k, trueRelativeErrors, '-^', 'DisplayName', 'True Integral Error');
    xlabel('Number of Points');
    ylabel('Estimated Error');
    title('Estimated Errors Comparison');
    legend show;
    grid on;
end
