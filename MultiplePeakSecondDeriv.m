function [totalIntegral, subdivisions] = MultiplePeakSecondDeriv(f, Emin, Emax, maxSubIntervals)
    
    numSamples = 1000; 
    thresholdRatio = 0.1;  

  
    E = linspace(Emin, Emax, numSamples);
    values = arrayfun(f, E);
    [~, locs] = findpeaks(values, E, 'MinPeakProminence', 0.1*max(values));  
    

    peakRegions = [];
    for i = 1:length(locs)
        peakCenter = locs(i);
        peakValue = f(peakCenter);
        threshold = peakValue * thresholdRatio;
        
       
        leftIndex = find(values(1:find(E == peakCenter)) < threshold, 1, 'last');
        if isempty(leftIndex)
            leftBound = Emin;
        else
            leftBound = E(leftIndex);
        end
        
        
        rightIndex = find(values(find(E == peakCenter):end) < threshold, 1, 'first');
        if isempty(rightIndex)
            rightBound = Emax;
        else
            rightBound = E(find(E == peakCenter) + rightIndex - 1);
        end
        
        peakRegions = [peakRegions; leftBound, peakCenter, rightBound];
    end

    
    peakRegions = sortrows(peakRegions, 1);  
    mergedPeakRegions = [];
    currentRegion = peakRegions(1, :);
    for i = 2:size(peakRegions, 1)
        if peakRegions(i, 1) <= currentRegion(3)
           
            currentRegion(3) = max(currentRegion(3), peakRegions(i, 3));
        else
           
            mergedPeakRegions = [mergedPeakRegions; currentRegion];
            currentRegion = peakRegions(i, :);
        end
    end
    mergedPeakRegions = [mergedPeakRegions; currentRegion];  

   
    totalIntegral = 0;
    subdivisions = [];

  
    for i = 1:size(mergedPeakRegions, 1)
        [peakIntegral, subDivs] = adaptiveTrapezoidalSecondDerivative(f, mergedPeakRegions(i, 1), mergedPeakRegions(i, 3), maxSubIntervals);
        totalIntegral = totalIntegral + peakIntegral;
        subdivisions = [subdivisions, subDivs];  
    end

   
    figure;
    plot(E, values, 'b-'); hold on;
    for i = 1:size(mergedPeakRegions, 1)
        fill([mergedPeakRegions(i, 1), mergedPeakRegions(i, 3), mergedPeakRegions(i, 3), mergedPeakRegions(i, 1)], ...
            [0, 0, max(values), max(values)], 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    hold off;
    title('Function with Merged Peak Integration Regions Highlighted');
    xlabel('E');
    ylabel('f(E)');
    legend('Function', 'Adaptive Integration Regions');

    % Output
    fprintf('Total Estimated Integral: %.4f\n', totalIntegral);
    fprintf('Subdivisions Used: %s\n', mat2str(unique(subdivisions)));
end