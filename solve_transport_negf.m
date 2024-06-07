function [E, TEL, TER, gEL, gER, GR] = solve_transport_negf(varargin)
    tollim = 1e-6;  % criterion for numerical tolerance limit
    mat = varargin{1}{1}{1};  % structure containing the material parameters
    rgf = mat.rgf;  % variable that determines whether rgf is turned on or off
    Nx = mat.Nx;  % number of discretization points
    dE = mat.dE;  % discretization length of energy vector, can be a vector
    V = mat.V;  % equilibrium band diagram
    x = mat.x;  % position vector
    dx = x(2) - x(1);  % discretization length of position vector

    H = calc_discrete_hamiltonian2(Nx, mat, V);  % Hamiltonian matrix
    Emin = min(V);  % minimum energy
    Emax = max(V) + 0.3;  % maximum energy
    
    % Construct energy vector E based on dE
    if isscalar(dE)
        E = Emin - 5 * dE:dE:Emax;  % homogeneous energy vector
    else
        E = Emin + [0, cumsum(dE)]; % non-uniform energy vector
    end
    NE = length(E);  % number of energy points

    TEL = zeros(1, NE);  % transmission from left to right
    TER = zeros(1, NE);  % transmission from right to left
    gEL = zeros(NE, Nx);  % DOS induced from the left contact
    gER = zeros(NE, Nx);  % DOS induced from the right contact

    switch rgf
        case 0
            % Loop over all the energy points
            for IE = 1:NE
                % Initialize the retarded self-energy SigmaR
                SigmaR = sparse(Nx, Nx);
                % Prepare left contact
                TL = -H(1, 2);
                DL = E(IE) - H(1, 1);
                kL = 1 / dx * acos(-DL / (2 * TL));
                if real(kL) < 0 || (imag(kL) < 0 && abs(imag(kL)) > tollim)
                    kL = -kL;
                end
                SigmaR11 = -TL * exp(1i * kL * dx);
                Gamma11 = 1i * (SigmaR11 - SigmaR11');
                % Prepare right contact
                TR = -H(Nx-1, Nx);
                DR = E(IE) - H(Nx, Nx);
                kR = 1 / dx * acos(-DR / (2 * TR));
                if real(kR) < 0 || (imag(kR) < 0 && abs(imag(kR)) > tollim)
                    kR = -kR;
                end
                SigmaRNN = -TR * exp(1i * kR * dx);
                GammaNN = 1i * (SigmaRNN - SigmaRNN');
                SigmaR(1, 1) = SigmaR11;
                SigmaR(Nx, Nx) = SigmaRNN;
                M = E(IE) * spdiags(ones(Nx, 1), 0, Nx, Nx) - H - SigmaR;
                GR = inv(M);
                gEL(IE, :) = (1 / (2 * pi * dx * 1e-9) * abs(GR(:, 1)).^2 * Gamma11)';
                gER(IE, :) = (1 / (2 * pi * dx * 1e-9) * abs(GR(:, Nx)).^2 * GammaNN)';
                TEL(IE) = real(GR(1, Nx) * GammaNN * GR(1, Nx)' * Gamma11);
                TER(IE) = real(GR(Nx, 1) * Gamma11 * GR(Nx, 1)' * GammaNN);
            end
        case 1
            eta = 1e-8;
            gR = zeros(NE, Nx);
            GR = zeros(NE, Nx);
            GRn1 = zeros(NE, Nx);
            TL = -H(1, 2);
            DL = E' + 1i * eta - H(1, 1);
            kL = 1 / dx * acos(-DL / (2 * TL));
            SigmaR11 = -TL * exp(1i * kL * dx);
            Gamma11 = 1i * (SigmaR11 - conj(SigmaR11));
            TR = -H(Nx-1, Nx);
            DR = E' + 1i * eta - H(Nx, Nx);
            kR = 1 / dx * acos(-DR / (2 * TR));
            SigmaRNN = -TR * exp(1i * kR * dx);
            GammaNN = 1i * (SigmaRNN - conj(SigmaRNN));
            gR(:, Nx) = 1 ./ (E' - H(Nx, Nx) - SigmaRNN);
            for I = Nx-1:-1:2
                gR(:, I) = 1 ./ (E' - H(I, I) - H(I, I+1) * gR(:, I+1) * H(I+1, I));
            end
            GR(:, 1) = 1 ./ (E' - H(1, 1) - H(1, 2) * gR(:, 2) * H(2, 1) - SigmaR11);
            GRn1(:, 1) = GR(:, 1);
            for I = 2:Nx
                GR(:, I) = gR(:, I) + gR(:, I) .* H(I, I-1) .* GR(:, I-1) .* H(I-1, I) .* gR(:, I);
                GRn1(:, I) = gR(:, I) .* H(I, I-1) .* GRn1(:, I-1);
            end
            gEL = 1 / (2 * pi * dx * 1e-9) * abs(GRn1).^2 .* (Gamma11 * ones(1, Nx));
            gER = real(1i * (GR - conj(GR))) / (2 * pi * dx * 1e-9) - gEL;
            TER = (abs(GRn1(:, Nx)).^2 .* Gamma11 .* GammaNN)';
            TEL = TER;
    end
end
