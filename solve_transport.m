%function to calculate the transmission and density-of-states in the
%structures defined in material.m
function [E TEL TER gEL gER]=solve_transport(varargin)

tollim=1e-6;                            %criterion for numerical tolerance limit

hbar=1.0546e-34;                        %reduced Plank's constant
m0=9.1095e-31;                          %rest effective mass
e=1.6022e-19;                           %elementary charge

mat=varargin{1}{1}{1};                  %structure containing the material parameters

Nx=mat.Nx;                              %number of discretization points
dE=mat.dE;                              %discretization length of energy vector
V=mat.V;                                %equilibrium band diagram
x=mat.x;                                %position vector
m=mat.m;                                %position-dependent effective mass

dx=x(2)-x(1);                           %discretization length of position vector

mL=m(1);                                %effective mass of the left contact
mR=m(Nx);                               %effective mass of the right contact

H=calc_discrete_hamiltonian(Nx,mat,V);  %Hamiltonian matrix

Emin=min(V);                            %minimum energy
Emax=max(V)+0.3;                        %maximum energy

E=Emin-5*dE:dE:Emax;                    %homogeneous energy vector
NE=length(E);                           %number of energy points

TEL=zeros(1,NE);                        %transmission from left to right
TER=zeros(1,NE);                        %transmission from right to left
gEL=zeros(NE,Nx);                       %DOS induced from the left contact
gER=zeros(NE,Nx);                       %DOS induced from the right contact

%loop over all the energy points
for IE=1:NE,
    
    %Initialize the possibility to inject from right and left contact
    injectionL=0;
    injectionR=0;
    
    %Self energy Sigma and injection vector S
    Sigma=sparse(Nx,Nx);
    S=zeros(Nx,2);
    
    %Prepare left contact
    TL=-H(1,2);
    DL=E(IE)-H(1,1);

    kL=1/dx*acos(-DL/(2*TL));
    
    %Make sure that wave vector kL corresponds to a propagating state ->
    condition1=real(kL)<0;
    condition2=(imag(kL)<0)&&abs(imag(kL)>tollim);
    
    if condition1||condition2,
        kL=-kL;
    end
    
    %Convert wave vector kL to self-energy and injection vector
    Sigma11=-TL*exp(i*kL*dx);
    S11=-TL*(1-exp(2*i*kL*dx));
    
    %calculate electron velocity at x(1)
    if (-DL/(2*TL))<1,
        vL=hbar*kL*1e9/(mL*m0);
        dEkL_dk=hbar^2*kL*1e9/(mL*m0)/e;
        injectionL=1;
    else
        vL=0;
    end
    
    %Prepare right contact
    TR=-H(Nx-1,Nx);
    DR=E(IE)-H(Nx,Nx);
    
    kR=1/dx*acos(-DR/(2*TR));
    
    %Make sure that wave vector kR corresponds to a propagating state ->
    condition1=real(kR)<0;
    condition2=(imag(kR)<0)&&abs(imag(kR)>tollim);
    
    %calculate electron velocity at x(Nx)
    if condition1||condition2,
        kR=-kR;
    end
    
    %Convert wave vector kR to self-energy and injection vector
    SigmaNN=-TR*exp(i*kR*dx);
    SNN=-TR*(1-exp(2*i*kR*dx));
    
    if (-DR/(2*TR))<1,
        vR=hbar*kR*1e9/(mR*m0);
        dEkR_dk=hbar^2*kR*1e9/(mR*m0)/e;
        injectionR=1;
    else
        vR=0;
    end
    
    %Update self-energy matrix and injection vector
    Sigma(1,1)=Sigma11;
    Sigma(Nx,Nx)=SigmaNN;
    
    S(1,1)=S11;
    S(Nx,2)=SNN;
    
    %Prepare matrix M=(E-H-Sigma)
    M=E(IE)*spdiags(ones(Nx,1),0,Nx,Nx)-H-Sigma;
    
    %Solve M*phi=S (1) LU factorization (2) Solve
    [L U P Q]=lu(M);
    
    phi=Q*(U\(L\(P*S)));
    phiL=phi(:,1);
    phiR=phi(:,2);
    
    %Extract bR when state injected from the left
    bR_L=-1/TR*(DR*phiL(Nx)+TR*phiL(Nx-1));
    
    %Extract bL when state injected from the right
    bL_R=-1/TL*(DL*phiR(1)+TL*phiR(2));
    
    %Calculate transmission from left to right TEL and from right to left TER
    if injectionL&&injectionR
        TEL(IE)=abs(bR_L).^2*vR/vL;
        TER(IE)=abs(bL_R).^2*vL/vR;
    end
    
    %Calculate DOS from left contact (do not include the 2 for spin)
    if injectionL,
        gEL(IE,:)=1/(2*pi)*(abs(phiL).^2/dEkL_dk)';
    end
    
    %Calculate DOS from right contact (do not include the 2 for spin)
    if injectionR,
        gER(IE,:)=1/(2*pi)*(abs(phiR).^2/dEkR_dk)';
    end
    
end