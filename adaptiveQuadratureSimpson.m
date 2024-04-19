function [estimatedError, matchCounter] = adaptiveQuadratureSimpson(f, Emin, Emax, maxPoints)
    subdivisions = [Emin, Emax];
    matchCounter = 0; % Counter for matches
    while length(subdivisions) - 1 < maxPoints
        [estimates, errors, trueErrors] = updateEstimatesSimpson(f, subdivisions);
        
        [~, maxErrorIdx] = max(errors);
        [~, maxTrueErrorIdx] = max(trueErrors);
        
        if maxErrorIdx == maxTrueErrorIdx
            matchCounter = matchCounter + 1; % Increment if the max error index matches
        end
        
        if length(subdivisions) - 1 < maxPoints
            newPoint = (subdivisions(maxErrorIdx) + subdivisions(maxErrorIdx + 1)) / 2;
            subdivisions = sort([subdivisions, newPoint]);
        else
            break;
        end
    end
    
    estimatedError = max(errors);
    adaptivePlot(f, subdivisions, Emin, Emax, 'Simpson Quadrature Estimation')
end

function [estimates, errors, trueErrors] = updateEstimatesSimpson(f, subdivisions)
    numSubdivisions = length(subdivisions) - 1;
    estimates = zeros(1, numSubdivisions);
    errors = zeros(1, numSubdivisions);
    trueErrors = zeros(1, numSubdivisions);
    
    for i = 1:numSubdivisions
        a = subdivisions(i);
        b = subdivisions(i+1);
        m = (a + b) / 2;
        simpsonEstimate = (b - a) / 6 * (f(a) + 4 * f(m) + f(b));
        trapezoidalEstimate = (f(a) + f(b)) * (b - a) / 2;
        
        trueIntegral = integral(f, a, b);
        
        estimate = simpsonEstimate;
        errorEstimate = abs(simpsonEstimate - trapezoidalEstimate) / 15; % Error estimate refinement
        
        estimates(i) = estimate;
        errors(i) = errorEstimate;
        trueErrors(i) = abs((estimate - trueIntegral) / trueIntegral);
    end
end


