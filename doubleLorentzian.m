function y = doubleLorentzian(E)
    E01 = -5e-3;  % Center of the first peak
    Gamma1 = 1e-3;  % Width of the first peak

    E02 = 5e-3;  % Center of the second peak
    Gamma2 = 2*1e-3;  % Width of the second peak

    y = (Gamma1/2) ./ (pi * ((E - E01).^2 + (Gamma1/2)^2)) + (Gamma2/2) ./ (pi * ((E - E02).^2 + (Gamma2/2)^2));
end

