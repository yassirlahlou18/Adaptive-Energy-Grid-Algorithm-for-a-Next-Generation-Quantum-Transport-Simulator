function [estimatedError, estimatedErrorWithEstimation] = adaptiveTrapezoidalRichardson(f, energyPoints, maxPoints)
    subdivisions = energyPoints;
    matchCounter = 0; 
    validSubdivisionFound = true; % To check if a valid subdivision was found

    while (length(subdivisions) - 1 < maxPoints) && validSubdivisionFound
        [estimates, errors, trueErrors] = updateEstimatesRichardson(f, subdivisions);
        
        [~, maxErrorIdx] = max(errors);
        [~, maxTrueErrorIdx] = max(trueErrors);
        
        if maxErrorIdx == maxTrueErrorIdx
            matchCounter = matchCounter + 1; 
        end
        
    
        intervalSizes = diff(subdivisions);
        maxIntervalSize = max(intervalSizes);
        
        validSubdivisionFound = false;
        [~, errorSortedIndices] = sort(errors, 'descend'); % Sorting error indices in descending order

        for idx = errorSortedIndices
            newPoint = (subdivisions(idx) + subdivisions(idx + 1)) / 2;
            newSubdivisions = sort([subdivisions, newPoint]);
            
            % Check the size of the new sub-intervals
            newIntervalSizes = diff(newSubdivisions);
            minNewIntervalSize = min(newIntervalSizes);

            % Only accept the new subdivision if it is not more than 5 times smaller than the largest interval
            if maxIntervalSize / minNewIntervalSize <= 0.1*maxPoints
                subdivisions = newSubdivisions;
                validSubdivisionFound = true;
                break; % Exit loop after valid subdivision is found
            end
        end

        if ~validSubdivisionFound
            fprintf('No valid subdivisions found that satisfy the maximum size ratio constraint.\n');
            break; % Exit while loop if no valid subdivisions can be found
        end
    end
    
    
    totalEstimate = sum(estimates);
    trueTotalIntegral=integral(f, energyPoints(1), energyPoints(length(energyPoints)));
    estimatedError = abs(totalEstimate-trueTotalIntegral)/trueTotalIntegral;
    estimatedErrorWithEstimation=max(errors);
    %trueTotal = integral(f, Emin, Emax);
    %totalRelativeError = abs((totalEstimate - trueTotal) / trueTotal);
    %adaptivePlot(subdivisions, f, totalEstimate, 'Richardson Estimation')
    % adaptivePlot(f, subdivisions, Emin, Emax, 'Richardson estimatiom')
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
