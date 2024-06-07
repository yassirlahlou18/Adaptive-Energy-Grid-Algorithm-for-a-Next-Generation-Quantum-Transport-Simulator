function totalRelativeError = adaptiveTrapezoidalFixedNbPts(f, energyPoints, maxPoints)

    subdivisions = energyPoints;
    trapezoidalEstimates = [];
    trueIntegrals = [];
    errors = [];
    validSubdivisionFound = true; % To check if a valid subdivision was found


    while (length(subdivisions) - 1 < maxPoints) && validSubdivisionFound
       
        [trapezoidalEstimates, trueIntegrals, errors] = updateEstimatesAndErrors(f, subdivisions);
      
        [~, maxErrorIdx] = max(errors);
        

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
            if maxIntervalSize / minNewIntervalSize <= 3.5
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
    

    
    % total relative error over whole domain
    totalEstimate = sum(trapezoidalEstimates);
    trueTotal = sum(trueIntegrals);
    totalRelativeError = abs((totalEstimate - trueTotal) / trueTotal);





    % adaptivePlot(f, subdivisions, Emin, Emax, 'true integral')
    
end



function [trapezoidalEstimates, trueIntegrals, errors] = updateEstimatesAndErrors(f, subdivisions)
    numSubdivisions = length(subdivisions) - 1;
    trapezoidalEstimates = zeros(1, numSubdivisions);
    trueIntegrals = zeros(1, numSubdivisions);
    errors = zeros(1, numSubdivisions);
    
    for i = 1:numSubdivisions
        trapezoidalEstimates(i) = trapezoidalRule(f, subdivisions(i), subdivisions(i + 1));   
        trueIntegrals(i) = integral(f, subdivisions(i), subdivisions(i + 1));
        % relative error
        errors(i) = abs((trapezoidalEstimates(i) - trueIntegrals(i)) / trueIntegrals(i));
    end
end

function estimate = trapezoidalRule(f, a, b)
    estimate = (f(a) + f(b)) * (b - a) / 2;
end
