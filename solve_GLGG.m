function [GL, GG] = solve_GLGG(E, varargin)

hbar = 1.0546e-34;  % reduced Planck's constant
e = 1.6022e-19;  % elementary charge
kB = 1.38e-23;  % Boltzmann constant

mat = varargin{1};  % structure containing the material parameters

Nx = mat.Nx;  % number of discretization points
V = mat.V;  % equilibrium band diagram
x = mat.x;  % position vector

dx = x(2) - x(1);  % discretization length of position vector

Temp = mat.Temp;  % lattice temperature

H = calc_discrete_hamiltonian2(Nx, mat, V);  % Hamiltonian matrix

NE = length(E);  % number of energy points

EfL = 0.2;  % Fermi level of the left contact
EfR = -0.4;  % Fermi level of the right contact

fEL = 1 ./ (1 + exp(e * (E - EfL) / (kB * Temp)));  % distribution function of the left contact
fER = 1 ./ (1 + exp(e * (E - EfR) / (kB * Temp)));  % distribution function of the right contact

ndE = zeros(NE, Nx);  % energy-resolved charge density
IdE = zeros(NE, Nx);  % energy-resolved current density

eta = 1e-8;  % imaginary part of the energy

gR = zeros(NE, Nx);  % small retarded Green's Function
gL = zeros(NE, Nx);  % small lesser Green's Function
gG = zeros(NE, Nx);  % small greater Green's Function
GR = zeros(NE, Nx);  % retarded Green's Function
GL = zeros(NE, Nx);  % lesser Green's Function 
GG = zeros(NE, Nx);  % greater Green's Function 
GLNN1 = zeros(NE, Nx);  % first off-diagonal element of lesser Green's Function 
SigmaL = zeros(NE, Nx);  % lesser Self-energy
SigmaG = zeros(NE, Nx);  % greater Self-energy

% Calculation of the retarded Green's Function
SigmaR = 1/2 * (SigmaG - SigmaL);

% Recursive algorithm for the calculation of all the required Green's Functions
TL = -H(1, 2);
DL = E' + 1i * eta - H(1, 1) - SigmaR(:, 1);

kL = 1 / dx * acos(-DL / (2 * TL));
SigmaR11 = -TL * exp(1i * kL * dx);
Gamma11 = 1i * (SigmaR11 - conj(SigmaR11));

SigL11 = 1i * Gamma11 .* fEL';
SigG11 = 1i * Gamma11 .* (fEL' - 1);

TR = -H(Nx-1, Nx);
DR = E' + 1i * eta - H(Nx, Nx) - SigmaR(:, Nx);

kR = 1 / dx * acos(-DR / (2 * TR));
SigmaRNN = -TR * exp(1i * kR * dx);
GammaNN = 1i * (SigmaRNN - conj(SigmaRNN));

SigLNN = 1i * GammaNN .* fER';
SigGNN = 1i * GammaNN .* (fER' - 1);

gR(:, Nx) = 1 ./ (E' - H(Nx, Nx) - SigmaRNN - SigmaR(:, Nx));
gL(:, Nx) = gR(:, Nx) .* (SigLNN + SigmaL(:, Nx)) .* conj(gR(:, Nx));
gG(:, Nx) = gR(:, Nx) .* (SigGNN + SigmaG(:, Nx)) .* conj(gR(:, Nx));

for I = Nx-1:-1:2
    gR(:, I) = 1 ./ (E' - H(I, I) - H(I, I+1) * gR(:, I+1) * H(I+1, I) - SigmaR(:, I));
    gL(:, I) = gR(:, I) .* (H(I, I+1) .* gL(:, I+1) * H(I+1, I) + SigmaL(:, I)) .* conj(gR(:, I));
    gG(:, I) = gR(:, I) .* (H(I, I+1) .* gG(:, I+1) * H(I+1, I) + SigmaG(:, I)) .* conj(gR(:, I));
end

GR(:, 1) = 1 ./ (E' - H(1, 1) - H(1, 2) * gR(:, 2) * H(2, 1) - SigmaR11 - SigmaR(:, 1));
GL(:, 1) = GR(:, 1) .* (SigL11 + H(1, 2) * gL(:, 2) .* H(2, 1) + SigmaL(:, 1)) .* conj(GR(:, 1));
GG(:, 1) = GR(:, 1) .* (SigG11 + H(1, 2) * gG(:, 2) .* H(2, 1) + SigmaG(:, 1)) .* conj(GR(:, 1));

for I = 2:Nx
    GR(:, I) = gR(:, I) + gR(:, I) .* H(I, I-1) .* GR(:, I-1) .* H(I-1, I) .* gR(:, I);
    ML = gR(:, I) * H(I, I-1) .* GR(:, I-1) * H(I-1, I) .* gL(:, I);
    GL(:, I) = gL(:, I) + gR(:, I) .* H(I, I-1) .* GL(:, I-1) .* H(I-1, I) .* conj(gR(:, I)) + (ML - conj(ML));
    MG = gR(:, I) * H(I, I-1) .* GR(:, I-1) * H(I-1, I) .* gG(:, I);
    GG(:, I) = gG(:, I) + gR(:, I) .* H(I, I-1) .* GG(:, I-1) .* H(I-1, I) .* conj(gR(:, I)) + (MG - conj(MG));
end
end

