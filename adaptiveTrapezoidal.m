function [totalIntegral, numSubintervals, subintervalSizes] = adaptiveTrapezoidal(F, Emin, Emax, maxError)
    % Initialize the total integral to 0
    totalIntegral = 0;
    
    % Start with a single subinterval over the entire domain
    intervals = [Emin, Emax];
    
    % Initialize arrays to hold the number of subintervals and their sizes
    numSubintervals = 0;
    subintervalSizes = [];
    
    % Loop until all subintervals have been processed
    while ~isempty(intervals)
        newIntervals = [];
        
        for i = 1:2:length(intervals)-1
            % Extract the current subinterval limits
            a = intervals(i);
            b = intervals(i+1);
            
            % Compute the trapezoidal approximation over the current subinterval
            approxIntegral = trapz([a, b], [F(a), F(b)]);
            
            % Obtain the true integral over the current subinterval
            trueInt= integral(F, a, b);

            
            % Calculate the relative error
            relativeError = abs((trueInt - approxIntegral) / trueInt);
            
            % Check if the error exceeds the maximum allowed error
            if relativeError > maxError
                % If so, split the subinterval in half for further refinement
                midPoint = (a + b) / 2;
                newIntervals = [newIntervals, a, midPoint, midPoint, b];
            else
                % Otherwise, add the approximated integral to the total and record the subinterval size
                totalIntegral = totalIntegral + approxIntegral;
                numSubintervals = numSubintervals + 1;
                subintervalSizes = [subintervalSizes, b-a];
            end
        end
        
        % Prepare for the next iteration
        intervals = newIntervals;
    end
    
    %fprintf('Total Integral: %f\n', totalIntegral);
    %fprintf('Number of subintervals: %d\n', numSubintervals);
    %for i = 1:length(subintervalSizes)
        %fprintf('Subinterval %d size: %f\n', i, subintervalSizes(i));
    %end
end
