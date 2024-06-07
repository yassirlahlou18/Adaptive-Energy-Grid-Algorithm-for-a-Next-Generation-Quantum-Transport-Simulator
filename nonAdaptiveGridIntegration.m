function [n_approx, rel_error] = nonAdaptiveGridIntegration(f, Emin, Emax, maxpoints, n_true)
        
    subIntervalLength = (Emax - Emin) / maxpoints;

    % Initialize the integral
    subIntervalIntegral =  zeros(size(maxpoints));

    % Compute the integral over each subinterval
    for i = 1:maxpoints
        % Calculate the current subinterval endpoints
        x0 = Emin + (i-1) * subIntervalLength;
        x1 = Emin + i * subIntervalLength;
 

        % % Evaluate the function at the endpoints of the subinterval
        % f0 = f(x0);
        % f1 = f(x1);

        % Calculate the trapezoidal estimate for the subinterval
        subIntervalIntegral(i) = trapezoidalRule(f, x0, x1);

        % Sum up the estimates
        
    end
    n_approx=sum(subIntervalIntegral, 'omitnan');
   
    rel_error=abs(n_approx-n_true)/n_true;


end

