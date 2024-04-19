function [estimatedError, matchCounter] = adaptiveTrapezoidalSecondDerivative(f, Emin, Emax, maxPoints)
    subdivisions = [Emin, Emax];
    matchCounter = 0; %counter for matches
    while length(subdivisions) - 1 < maxPoints
        [estimates, errors, trueErrors] = updateEstimatesSecondDerivative(f, subdivisions);
        
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
    
    totalEstimate = sum(estimates);
    estimatedError = max(errors);
    %trueTotal = integral(f, Emin, Emax);
    %totalRelativeError = abs((totalEstimate - trueTotal) / trueTotal);
    %adaptivePlot(subdivisions, f, totalEstimate, 'Second Derivative Estimation')

    adaptivePlot(f, subdivisions, Emin, Emax, 'Second Derivative Estimation')
end

function [estimates, errors, trueErrors] = updateEstimatesSecondDerivative(f, subdivisions)
    numSubdivisions = length(subdivisions) - 1;
    estimates = zeros(1, numSubdivisions);
    errors = zeros(1, numSubdivisions);
    trueErrors = zeros(1, numSubdivisions); % true errors
    
    for i = 1:numSubdivisions
        a = subdivisions(i);
        b = subdivisions(i+1);
        estimate = trapezoidalRule(f, a, b);
        trueIntegral = integral(f, a, b); 
        
        secondDerivativeEstimate = estimateSecondDerivative(f, (a + b) / 2);
        errorEstimate = ((b - a)^3 / 12) * abs(secondDerivativeEstimate); 
        
        estimates(i) = estimate;
        errors(i) = errorEstimate; 
        trueErrors(i) = abs((estimate - trueIntegral) / trueIntegral); 
    end
    
end

function secondDeriv = estimateSecondDerivative(f, x)
    h = 1e-4; % A small step for the second derivative approximation
    secondDeriv = (f(x + h) - 2 * f(x) + f(x - h)) / h^2; % Second derivative estimation
end

function estimate = trapezoidalRule(f, a, b)
    estimate = (f(a) + f(b)) * (b - a) / 2; 
end
