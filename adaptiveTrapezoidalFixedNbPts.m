function totalRelativeError = adaptiveTrapezoidalFixedNbPts(f, Emin, Emax, maxPoints)

    subdivisions = [Emin, Emax];
    trapezoidalEstimates = [];
    trueIntegrals = [];
    errors = [];
    

    while length(subdivisions) - 1 < maxPoints
       
        [trapezoidalEstimates, trueIntegrals, errors] = updateEstimatesAndErrors(f, subdivisions);
      
        [~, maxErrorIdx] = max(errors);
        

        if length(subdivisions) - 1 < maxPoints
            newPoint = (subdivisions(maxErrorIdx) + subdivisions(maxErrorIdx + 1)) / 2;
            subdivisions = sort([subdivisions, newPoint]);
        else
            break;
        end
    end
    
    % total relative error over whole domain
    totalEstimate = sum(trapezoidalEstimates);
    trueTotal = sum(trueIntegrals);
    totalRelativeError = abs((totalEstimate - trueTotal) / trueTotal);





    adaptivePlot(f, subdivisions, Emin, Emax, 'true integral')
    
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
