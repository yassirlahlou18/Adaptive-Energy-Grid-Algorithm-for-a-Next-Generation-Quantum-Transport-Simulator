function estimate = trapezoidalRule(f, a, b)
    estimate = (f(a) + f(b)) * (b - a) / 2;
end