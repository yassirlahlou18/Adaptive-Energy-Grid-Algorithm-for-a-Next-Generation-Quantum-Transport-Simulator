function adaptivePlot(f, subdivisions, Emin, Emax, methodName)

    fig = figure; 
    set(fig, 'NumberTitle', 'off', 'Name', ['Adaptive Trapezoidal Rule with ', methodName ]);

    % Loop through each interval to plot trapezoids
    for i = 1:length(subdivisions)-1
        xLeft = subdivisions(i);
        xRight = subdivisions(i+1);
        yLeft = f(xLeft);
        yRight = f(xRight);

        % Plotting the trapezoids (approximation)
        patch([xLeft xRight xRight xLeft], [0 0 yRight yLeft], [0.9, 0.9, 0.9], 'EdgeColor', 'r', 'LineWidth', 1.5);
    end

    % Now plot the original function over the interval [Emin, Emax] on top of the approximation
    hold on;
    fplot(f, [Emin, Emax], 'b', 'LineWidth', 2);

    % Finalize the plot
    hold off;
    grid on;
    legend('Trapezoidal Approximation', 'Location', 'Best');
    xlabel('x');
    ylabel('f(x)');
    title(['Adaptive Trapezoidal Rule with ',methodName]);
end
