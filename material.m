function mat=material(structure,rgf)

%Units
%distance: nm
%potential: eV
%doping: cm^{-3}
 % Initialize structure fields
    mat.dx = [];
    mat.dE = [];
    mat.Temp = [];
    mat.L = [];
    mat.m = [];
    mat.V = [];
    mat.Ndop = [];
    mat.Epsr = [];
    mat.rgf = rgf; % Assign the input rgf to the output structure
switch structure,
    
    case 'constant_potential',
        
        dx=0.1;
        dE=1e-4;
        Temp=300;

        L=30;
        m=0.1;
        V=0;
        Ndop=1e19;
        Epsr=12;
    
    case 'quantum_well',
        
        dx=0.1;
        dE=1e-4;
        Temp=300;

        Lbarrier_left=20;
        Lqw=5;
        Lbarrier_right=20;

        mbarrier=0.1;
        mqw=0.065;

        Vbarrier=0.3;
        Vqw=0;
        
        Nbarrier=0;
        Nqw=0;
        
        Epsbarrier=10;
        Epsqw=12;

        L=[Lbarrier_left Lqw Lbarrier_right];
        m=[mbarrier mqw mbarrier];
        V=[Vbarrier Vqw Vbarrier];
        Ndop=[Nbarrier Nqw Nbarrier];
        Epsr=[Epsbarrier Epsqw Epsbarrier];
        
    case 'barrier',
        
        dx=0.2;
        dE=5e-4;
        Temp=300;

        Lside=30;
        Lqw_left=15;
        Lbarrier=4;
        Lqw_right=15;

        mbarrier=0.1;
        mqw=0.065;

        Vbarrier=0.3;
        Vqw=0;

        Nside=1e19;
        Nbarrier=0;
        Nqw=0;
        
        Epsbarrier=6;
        Epsqw=12;
        
        L=[Lside Lqw_left Lbarrier Lqw_right Lside];
        m=[mqw mqw mbarrier mqw mqw];
        V=[Vqw Vqw Vbarrier Vqw Vqw];
        Ndop=[Nside Nqw Nbarrier Nqw Nside];
        Epsr=[Epsqw Epsqw Epsbarrier Epsqw Epsqw];
        
    case 'rtd',
        
        dx=0.2;
        dE=2e-4;
        Temp=300;
        
        Lext_left=30;
        Lside_left=15;
        Lbarrier_left=2.6;
        Lqw=4;
        Lbarrier_right=2.6;
        Lside_right=15;
        Lext_right=30;
        
        mside=0.065;
        mbarrier=0.1;
        mqw=0.065;
        
        Vside=0;
        Vbarrier=0.3;
        Vqw=0;
        
        Next=1e19;
        Nside=0;
        Nbarrier=0;
        Nqw=1e18;
        
        Epsside=12;
        Epsbarrier=11;
        Epsqw=12;
        
        L=[Lext_left Lside_left Lbarrier_left Lqw Lbarrier_right Lside_right Lext_right];
        m=[mside mside mbarrier mqw mbarrier mside mside];
        V=[Vside Vside Vbarrier Vqw Vbarrier Vside Vside];
        Ndop=[Next Nside Nbarrier Nqw Nbarrier Nside Next];
        Epsr=[Epsside Epsside Epsbarrier Epsqw Epsbarrier Epsside Epsside];
        
    case 'dqw'
        
        dx=0.1;
        dE=1e-4;
        Temp=300;
        
        Lside_left=20;
        Lqw1=5;
        Lbarrier=4;
        Lqw2=6;
        Lside_right=20;
        
        mside=0.1;
        mbarrier=0.1;
        mqw1=0.065;
        mqw2=0.065;
        
        Vside=0.3;
        Vbarrier=0.3;
        Vqw1=0;
        Vqw2=0;
        
        Nside=0;
        Nbarrier=0;
        Nqw1=0;
        Nqw2=0;
        
        Epsside=10;
        Epsbarrier=10;
        Epsqw1=12;
        Epsqw2=12;
        
        L=[Lside_left Lqw1 Lbarrier Lqw2 Lside_right];
        m=[mside mqw1 mbarrier mqw2 mside];
        V=[Vside Vqw1 Vbarrier Vqw2 Vside];
        Ndop=[Nside Nqw1 Nbarrier Nqw2 Nside];
        Epsr=[Epsside Epsqw1 Epsbarrier Epsqw2 Epsside];
        
end

%dx (discretization length), L (segment length), m (segment effective mass), V (segment potential), Ndop (doping profile), Temp (temperature), and dE (distance between two energy discretization points) must absolutely be defined%
%No need to change the rest if you want a new structure%

xmin=0;
xmax=sum(L);

Nx=round((xmax-xmin)/dx)+1;

mat.Nx=Nx;
mat.x=zeros(1,Nx);
mat.V=zeros(1,Nx);
mat.Ndop=zeros(1,Nx);
mat.Epsr=zeros(1,Nx);

mat.x(1)=xmin;
mat.m(1)=m(1);
mat.V(1)=V(1);
mat.Ndop(1)=Ndop(1);
mat.Epsr(1)=Epsr(1);

mat.dE=dE;
mat.Temp=Temp;

mat.rgf=rgf;

index=2;
Ltot=L(1);

for I=1:length(L),
    
    while (mat.x(index-1)+dx)<=(Ltot+dx/10),
        
        mat.x(index)=mat.x(index-1)+dx;
        mat.m(index)=m(I);
        mat.V(index)=V(I);
        mat.Ndop(index)=Ndop(I);
        mat.Epsr(index)=Epsr(I);
        index=index+1;
        
    end
    
    if I<length(L),
        Ltot=Ltot+L(I+1);
    end
    
end