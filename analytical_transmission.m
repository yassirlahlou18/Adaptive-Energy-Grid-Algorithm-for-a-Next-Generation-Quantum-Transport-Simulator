function [E,TE]=analytical_transmission()

hbar=1.0546e-34;
e=1.6022e-19;
m0=9.1095e-31;

mL=0.065;
mC=0.1;
mR=0.065;

L=5e-9;

Vbarrier=0.3;

E=0:1e-3:1.0;

kL=sqrt(2*m0*mL*e*E/hbar.^2);
kR=sqrt(2*m0*mR*e*E/hbar.^2);
kappaC=sqrt(2*m0*mC*e*(Vbarrier-E)/hbar.^2);

F1=(kappaC-i*kR*mC/mR)./(2*kappaC).*exp(kappaC*L);
F2=(kappaC+i*kR*mC/mR)./(2*kappaC).*exp(-kappaC*L);
F3=(i*kL-kappaC*mL/mC).*(kappaC+i*kR*mC/mR)./(4*i*kL.*kappaC).*exp(-kappaC*L);
F4=(i*kL+kappaC*mL/mC).*(kappaC-i*kR*mC/mR)./(4*i*kL.*kappaC).*exp(kappaC*L);

F=F1+F2-F3-F4;

TE=1./abs(F).^2;