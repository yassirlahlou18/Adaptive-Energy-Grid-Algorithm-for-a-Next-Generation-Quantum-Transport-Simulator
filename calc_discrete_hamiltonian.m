function MG=calc_discrete_hamiltonian(Nx,mat,V_pot)

hbar=1.0546e-34;
e=1.6022e-19;
mo=9.1095e-31;

T0=hbar.^2/(mo*e*1e-18);

m=mat.m;
x=mat.x;

MG=sparse(Nx,Nx);

for I=2:Nx-1,
    
    dxp=x(I+1)-x(I);
    dxm=x(I)-x(I-1);
    
    MG(I,I-1)=-2*T0/(dxm*(m(I)+m(I-1))*(dxm+dxp));
    MG(I,I+1)=-2*T0/(dxp*(m(I)+m(I+1))*(dxm+dxp));
    MG(I,I)=-MG(I,I-1)-MG(I,I+1)+V_pot(I);
    
end

dx=x(2)-x(1);
MG(1,1)=T0/(dx.^2*m(1))+V_pot(1);
MG(1,2)=-T0/(2*dx.^2*m(1));

dx=x(Nx)-x(Nx-1);
MG(Nx,Nx)=T0/(dx.^2*m(Nx))+V_pot(Nx);
MG(Nx,Nx-1)=-T0/(2*dx.^2*m(Nx));