function [estimatedError, matchCounter] = adaptiveTrapezoidalRichardson(f, Emin, Emax, maxPoints)
    subdivisions = [Emin, Emax];
    matchCounter = 0; 
    while length(subdivisions) - 1 < maxPoints
        [estimates, errors, trueErrors] = updateEstimatesRichardson(f, subdivisions);
        
        [~, maxErrorIdx] = max(errors);
        [~, maxTrueErrorIdx] = max(trueErrors);
        
        if maxErrorIdx == maxTrueErrorIdx
            matchCounter = matchCounter + 1; 
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
    %adaptivePlot(subdivisions, f, totalEstimate, 'Richardson Estimation')
    adaptivePlot(f, subdivisions, Emin, Emax, 'Richardson estimatiom')
end

function [estimates, errors, trueErrors] = updateEstimatesRichardson(f, subdivisions)
    numSubdivisions = length(subdivisions) - 1;
    estimates = zeros(1, numSubdivisions);
    errors = zeros(1, numSubdivisions);
    trueErrors = zeros(1, numSubdivisions); 
    
    for i = 1:numSubdivisions
        a = subdivisions(i);
        b = subdivisions(i+1);
        T1 = trapezoidalRule(f, a, b);
        T2 = (trapezoidalRule(f, a, (a+b)/2) + trapezoidalRule(f, (a+b)/2, b)) / 2;
        trueIntegral = integral(f, a, b); 
        
        estimates(i) = T2;
        errors(i) = abs(T2 - T1); % Richardson error estimate
        trueErrors(i) = abs((T2 - trueIntegral) / trueIntegral); 
    end
end
