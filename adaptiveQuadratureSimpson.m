function [estimatedError, estimatedErrorWithEstimation] = adaptiveQuadratureSimpson(f, energyPoints, maxPoints)
    subdivisions = energyPoints;
    matchCounter = 0; % counter for matches
    validSubdivisionFound = true; % To check if a valid subdivision was found

    while (length(subdivisions) - 1 < maxPoints) && validSubdivisionFound
        [estimates, errors, trueErrors] = updateEstimatesSimpson(f, subdivisions);
        
        [~, maxErrorIdx] = max(errors);
        [~, maxTrueErrorIdx] = max(trueErrors);
        
        if maxErrorIdx == maxTrueErrorIdx
            matchCounter = matchCounter + 1; % Increment if the max error index matches
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
    trueTotalIntegral = integral(f, energyPoints(1), energyPoints(length(energyPoints)));
    estimatedError = abs(totalEstimate - trueTotalIntegral) / trueTotalIntegral;
    estimatedErrorWithEstimation = max(errors);
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


