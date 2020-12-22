%% Fresnel Coefficients 
% This script provies Fresnel reflection and transmission coefficients
% for light incident on a planar slab of thickness 'd' in arbitrary
% directions characterized by ($\theta,\phi$). Here $\theta$ denotes
% the angle made by the incident light with $z$-axis while $\phi$
% denotes the angle made by the in-plane wavevector of the incident
% light with $x$-axis. The thickness 'd' is in units of $c/\omega$
% where 'c' is vacuum speed of light and $\omega$ is the frequency.
%
% For each example, a material matrix 'MM' is defined which contains
% dimensionless permeability tensor ($\epsilon$), permeability tensor
% ($\mu$), and magneto-electric coupling tensors ($\xi,\zeta$). For
% each material example, these response tensors are not arbitrary but
% also satisfy passivity constaint. We ensure that this passivity
% is satisfied since our theory is not applicable for optically active
% gain media. If the medium is not passive, an error message will be
% displayed. 
%
% As discussed in our manuscript, the spin-resolved emissivity along
% the incidence direction ($\theta,\phi$) requires the calculation
% of Fresnel reflection and transmission coefficients for light incident
% along the reflection direction ($\theta,\phi+\pi$) and the transmission
% direction ($\pi-\theta, \phi+\pi$) on the other side of the slab. 
% We calculate all these coefficients for each example. For simplicity,
% we focus on $\theta \in [0,\pi/2]$ and $\phi \in =[0,2\pi]$ which
% denotes light incident in the top hemisphere.
%
% The user can specify the parameters $\theta,\phi$, thickness 'd',
% and the material matrix MM in this script.
% The spin-resolved Kirchhoff's laws are derived based on certain
% relations between the Fresnel coefficients. These can be verified
% for various material (MM) and thickness (d) parameters,
% and indcidence directions ($\theta,\phi$). The user only needs to
% choose the parameters and run the command 'publish('make_report.m','pdf')'
% to produce a report on Fresnel coefficients for several examples. 
% On a standard computer with MATLAB 2018, it takes approximately 10
% minutes to put together all results and provide this universal
% perspective of Fresnel coefficients for several material classes. 
% Below, we fix some of the parameters.

theta=pi/3;  % angle between incident light and z-axis (normal to slab)
phi=pi/6; % angle between parallel wavevector and x-axis 
d=0.4; % thickness of slab in units of c/w. 

%%
% Below, we first consider reciprocal media that satisfy
% SKL-1 and then nonreciprocal media which satisfy SKL-2 or SKL-3.
% We also provide these parameters for a double layered thin film. 
% 

%% Uniaxial anisotropic (reciprocal SKL-1) medium
nf=1  % example 1 and figure 1
ep=[2+0.1i 0 0; 0 -2+0.1i 0; 0 0 2+0.1i]; % ep_xx not equal to ep_yy
mu=(1+1e-6*1i)*eye(3);   % small loss added for passivity 
xi=zeros(3); zeta=zeros(3);
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  
		%

%% Biaxial anisotropic (reciprocal SKL-1) medium
nf=nf+1
ep=[2+0.1i 0 0; 0 -2+0.1i 0; 0 0 3+0.1i];
mu=(1+1e-6*1i)*eye(3);
xi=zeros(3); zeta=zeros(3);
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

%
%% Isotropic magnetoelectric (reciprocal SKL-1) medium (Pasteur medium)
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=0.2i*eye(3);   % nonzero magnetoelectric coupling 
zeta=-transpose(xi);  % reciprocity condition 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

%
%% Anisotropic magnetoelectric (reciprocal SKL-1) medium [xy-coupling]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=[0 0.2i 0; 0 0 0; 0 0 0;];  % Ex-Hy, Ey-Hx coupling 
zeta=-transpose(xi);  % reciprocity condition 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

				%
%% Anisotropic magnetoelectric (reciprocal SKL-1) medium [xz-coupling]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=[0 0 0.2i; 0 0 0; 0 0 0;];  % Ex-Hz, Ez-Hx coupling 
zeta=-transpose(xi);  % reciprocity condition 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

				%
%% Anisotropic magnetoelectric (reciprocal SKL-1) medium [yz-coupling]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=[0 0 0; 0 0 0.2i; 0 0.1i 0;];  % Ey-Hz, Ez-Hy coupling 
zeta=-transpose(xi);  % reciprocity condition 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

				%
%% Composite anisotropic magnetoelectric (reciprocal SKL-1) medium 
nf=nf+1 
ep=[2+0.1i 0 0; 0 -2+0.1i 0; 0 0 2+0.1i]; % uniaxial anisotropy 
mu=(1+1e-6*1i)*eye(3);
xi=[0 0 0; 0 0 0.2i; 0 0.1i 0;];  % magnetoelectric anisotropy 
zeta=-transpose(xi);  % reciprocity condition 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

				%
%% Multilayered anisotropic magnetoelectric (reciprocal SKL-1) medium
nf=nf+1
% Material 1 
ep1=[4+0.1i 0 0; 0 -3+0.1i 0; 0 0 4+0.1i];  % uniaxial anisotropic 
mu1=(1+1e-6*1i)*eye(3);
xi1=zeros(3); zeta1=zeros(3);
MM1=[ep1 xi1; zeta1 mu1];
% Material 2
ep2=(6+0.1i)*eye(3);
mu2=(1+1e-6*1i)*eye(3);
xi2=0.2i*eye(3);   % isotropic magnetoelectric coupling 
zeta2=-transpose(xi2);
MM2=[ep2 xi2; zeta2 mu2];
% Geometry
td=[0.4 0.5];  % thicknesses of layers.
make_2layered_table(theta,phi,MM1,MM2,td,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incidence along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=-t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-1.  

 %
%% Gyroelectric (nonreciprocal SKL-2) medium [gyrotropy axis along z]
nf=nf+1 
ep=[4+0.1i 0.2i 0; -0.2i 4+0.1i 0; 0 0 4+0.1i;]; % must be anti-symmetric 
mu=(1+1e-6*1i)*eye(3);
xi=zeros(3); zeta=zeros(3); 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-2.  

				%
%% Gyroelectric (nonreciprocal SKL-3) medium [gyrotropy axis along x]
nf=nf+1 
ep=[4+0.1i 0 0; 0 4+0.1i 0.2i; 0 -0.2i 4+0.1i;]; % must be anti-symmetric 
mu=(1+1e-6*1i)*eye(3);
xi=zeros(3); zeta=zeros(3); 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-2.  

				%
%% Gyroelectric (nonreciprocal SKL-3) medium [gyrotropy axis along y]
nf=nf+1 
ep=[4+0.1i 0 0.2i; 0 4+0.1i 0; -0.2i 0 4+0.1i;]; % must be anti-symmetric 
mu=(1+1e-6*1i)*eye(3);
xi=zeros(3); zeta=zeros(3); 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-3.  
 
				%
%% Gyromagnetic (nonreciprocal SKL-2) medium [gyrotropy axis along z]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=[-2+0.1i 0.2i 0; -0.2i -2+0.1i 0; 0 0 -2+0.1i]; % must be anti-symmetric
xi=zeros(3); zeta=zeros(3); 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-2.  

				%
%% Gyromagnetic (nonreciprocal SKL-3) medium [gyrotropy axis along x]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=[2+0.1i 0 0; 0 2+0.1i 0.2i; 0 -0.2i 2+0.1i]; % must be anti-symmetric
xi=zeros(3); zeta=zeros(3); 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-3.  

				%
%% Isotropic magnetoelectric (nonreciprocal SKL-2) medium (Tellegen medium)
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=0.2*eye(3);  %  nonzero magnetoelectric coupling 
zeta=transpose(xi);  % non-reciprocity  
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-2.  

				%
%% Anisotropic magnetoelectric (nonreciprocal SKL-3) medium [xz-coupling]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=[0 0 0.2; 0 0 0; 0 0 0;]; % Ex-Hz, Ez-Hx coupling 
zeta=transpose(xi);  % non-reciprocity
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-3.  

      %
%% Anisotropic magnetoelectric (nonreciprocal SKL-3) medium [yz-coupling]
nf=nf+1 
ep=(4+0.1i)*eye(3); 
mu=(1+1e-6*1i)*eye(3);
xi=[0 0 0; 0 0 0.2; 0 0.1 0;]; % Ey-Hz, Ez-Hy coupling  
zeta=transpose(xi);  % non-reciprocity  
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-3.  

				%
%% Composite gyroelectric magnetoelectric (nonreciprocal SKL-2) medium
nf=nf+1 
ep=[4+0.1i 0.2i 0; -0.2i 4+0.1i 0; 0 0 4+0.1i];  % gyrotropy axis along z  
mu=(1+1e-6*1i)*eye(3);
xi=0.2*eye(3);  % Isotropic magnetoelectric 
zeta=transpose(xi);  % rotational symmetry preserved 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-2.  

				%
%% Composite gyroelectric magnetoelectric (nonreciprocal SKL-3) medium
nf=nf+1 
ep=[4+0.1i 0 0; 0 4+0.1i 0.2i; 0 -0.2i 4+0.1i];  % gyrotropy axis alon x  
mu=(1+1e-6*1i)*eye(3);
xi=[0 0 0.2; 0 0 0; 0.1 0 0;];  % Ex-Hz, Ez-Hx coupling
zeta=transpose(xi);
% rotational symmetry not preserved for this nonreciprocal medium 
MM=[ep xi; zeta mu];
make_table(theta,phi,MM,d,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-3.  

%
%% Multilayered gyroelectric magnetoelectric (nonreciprocal SKL-2) medium
nf=nf+1
% Material 1 
ep1=[4+0.1i 0.2i 0; -0.2i 4+0.1i 0; 0 0 4+0.1i];  % gyrotropy axis along z  
mu1=(1+1e-6*1i)*eye(3);
xi1=zeros(3); zeta1=zeros(3);
MM1=[ep1 xi1; zeta1 mu1];
% Material 2
ep2=(6+0.1i)*eye(3);
mu2=(1+1e-6*1i)*eye(3);
xi2=0.2*eye(3);  % isotropic magnetoelectric 
zeta2=transpose(xi2); % nonreciprocity 
MM2=[ep2 xi2; zeta2 mu2];
% Geometry
td=[0.4 0.5];  % thicknesses of layers.
make_2layered_table(theta,phi,MM1,MM2,td,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{ss,pp}(\theta,\phi)=r_{ss,pp}(\theta,\phi+\pi)$, 
% $r_{sp}(\theta,\phi)=r_{ps}(\theta,\phi+\pi)$,
% $t_{ss,pp}(\theta,\phi)=t_{ss,pp}(\pi-\theta,\phi+\pi)$,
% $t_{sp}(\theta,\phi)=t_{ps}(\pi-\theta,\phi+\pi)$.
% These relations lead to SKL-2.  

				%
%% Multilayered gyroelectric magnetoelectric (nonreciprocal SKL-3) medium
nf=nf+1
% Material 1 
ep1=[4+0.1i 0 0.2i; 0 4+0.1i 0; -0.2i 0 4+0.1i];  % gyrotropy axis along y  
mu1=(1+1e-6*1i)*eye(3);
xi1=zeros(3); zeta1=zeros(3);
MM1=[ep1 xi1; zeta1 mu1];
% Material 2
ep2=(6+0.1i)*eye(3);
mu2=(1+1e-6*1i)*eye(3);
xi2=[0 0 0; 0 0 0.2; 0 0.1 0;];  % Ey-Hz, Ez-Hy coupling 
zeta2=transpose(xi2);
MM2=[ep2 xi2; zeta2 mu2];
% Geometry
td=[0.4 0.5];  % thicknesses of layers.
make_2layered_table(theta,phi,MM1,MM2,td,nf);

%%
% In the table, incidence direction is $(\theta,\phi)$, reflection direction
% is $(\theta,\phi+\pi)$, and transmission direction is $(\pi-\theta,\phi+\pi)$.
% The reflection and transmission coefficients for light incident along these directions 
% satisfy the relations: 
% $r_{sp}(\theta,\phi)=-r_{ps}(\theta,\phi)$,
% $t_{ss,pp}(\pi-\theta,\phi+\pi)=t_{ss,pp}(\theta,\phi+\pi)$,
% $t_{sp}(\pi-\theta,\phi+\pi)=-t_{ps}(\theta,\phi+\pi)$.
% These relations lead to SKL-3.  


