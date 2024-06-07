function [estimatedError, estimatedErrorWithEstimation, subdivisions, functionValues] = adaptiveTrapezoidalSecondDerivative2(f, energyPoints, maxPoints)
    subdivisions = energyPoints;
    functionValues = f(energyPoints); % Initialize function values at initial points
    matchCounter = 0; % Counter for matches
    validSubdivisionFound = true; % To check if a valid subdivision was found

    while (length(subdivisions) - 1 < maxPoints) && validSubdivisionFound
        [estimates, errors, trueErrors] = updateEstimatesSecondDerivative(f, subdivisions);
        
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
            newFunctionValues = f(newSubdivisions); % Get function values at new subdivisions
            
            % Check the size of the new sub-intervals
            newIntervalSizes = diff(newSubdivisions);
            minNewIntervalSize = min(newIntervalSizes);

            % Only accept the new subdivision if it is not more than 5 times smaller than the largest interval
            if maxIntervalSize / minNewIntervalSize <= 0.1 * maxPoints
                subdivisions = newSubdivisions;
                functionValues = newFunctionValues; % Update function values with new points
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
    trueErrors = zeros(1, numSubdivisions); % True errors
    
    for i = 1:numSubdivisions
        a = subdivisions(i);
        b = subdivisions(i + 1);
        estimate = trapezoidalRule(f, a, b);
        trueIntegral = integral(f, a, b);
        
        secondDerivativeEstimate = estimateSecondDerivative(f, (a + b) / 2);
        errorEstimate = ((b - a) ^ 3 / 12) * abs(secondDerivativeEstimate);
        
        estimates(i) = estimate;
        errors(i) = errorEstimate;
        trueErrors(i) = abs((estimate - trueIntegral) / trueIntegral);
    end
end

function secondDeriv = estimateSecondDerivative(f, x)
    h = 1e-4; % A small step for the second derivative approximation
    secondDeriv = (f(x + h) - 2 * f(x) + f(x - h)) / h ^ 2; % Second derivative estimation
end

function estimate = trapezoidalRule(f, a, b)
    estimate = (f(a) + f(b)) * (b - a) / 2;
end
