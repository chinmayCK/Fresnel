function[rss, rps, rsp, rpp, tss, tps, tsp, tpp]=multilayer_top(theta,phi)

  kp=sin(theta);  % parallel component of wavevector.
  
  %% define materials to be used here.
  exx=3+0.1i; 
  ep=[exx 0 0.1i; 0 exx 0; -0.1i 0 exx]; 
  %mu=(1+1e-6*1i)*eye(3);
  mu=[2+0.1i 0 0; 0 2+0.1i 0.1i; 0 -0.1i 2+0.1i;];  
  %xi=[0 0.1 0; 0.1 0 0; 0 0 0;]; zeta=1*transpose(xi); 
  %xi=0.2*eye(3);  zeta=xi; 
  xi=zeros(3); zeta=zeros(3);
  MM1=[ep xi; zeta mu;];
  Mc=-1i*MM1; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);   % is pas=1, then this material is passive 

  %%
  exx=4+0.1i; 
  ep=[exx 0 0; 0 exx 0; 0 0 exx]; 
  mu=(1+1e-6*1i)*eye(3);
  %mu=[2+0.1i 0 0.1i; 0 2+0.1i 0; -0.1i 0 2+0.1i;];  
  xi=[0 0 0.1; 0 0 0.2; 0 0 0;]; zeta=1*transpose(xi); 
  %xi=0.2*eye(3);  zeta=xi; 
  %xi=zeros(3); zeta=zeros(3);
  MM2=[ep xi; zeta mu;];
  Mc=-1i*MM2; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);   % is pas=1, then this material is passive 

  %% Define the geometry (layers and materials etc)
  td=[0.3, 0.4, 0.5];  % define the number of layers here and thicknesses
  N=length(td);

%define the material in each layer. for now we will do it for three layers.
% solve for kz and field profiles in each layer.
% there are either 2 or 4 solutions of kz. and 4 solutions for fields inside. 
  [kz,ad,bd,au,bu,s]=getkzEH(kp,phi,MM1);  % For material 1.
  if s==0  % if we cannot find the solutions for the given material then s=0. 
    disp('special attention required for this material'); return;
  else
  end;
  %%% the folowwing part puts material 1 in layer 1 and layer 3. 
  kzd1(1)=kz(1); kzd2(1)=kz(2); kzu1(1)=kz(3); kzu2(1)=kz(4);  % kz(downward going)d(first or second wave)1--(layer no.)
  afd(:,1)=ad; bfd(:,1)=bd; afu(:,1)=au; bfu(:,1)=bu;  % the associated fields.
  % Because the same material is there in the third layer, we have
  kzd1(3)=kz(1); kzd2(3)=kz(2); kzu1(3)=kz(3); kzu2(3)=kz(4);  % kz(downward going)d(first or second wave)1--(layer no.)
  afd(:,3)=ad; bfd(:,3)=bd; afu(:,3)=au; bfu(:,3)=bu;  % the associated fields.

  [kz,ad,bd,au,bu,s]=getkzEH(kp,phi,MM2);  % For material 2.
  if s==0
    disp('special attention required for this material'); return;
  else
  end;
  %%%% the following part puts material 2 in layer 2. 
  kzd1(2)=kz(1); kzd2(2)=kz(2); kzu1(2)=kz(3); kzu2(2)=kz(4);  % kz(downward going)d(first or second wave)1--(layer no.)
  afd(:,2)=ad; bfd(:,2)=bd; afu(:,2)=au; bfu(:,2)=bu;  % the associated fields.

  
  % Solve the boundary conditions.
  % The unknowns are rss,rsp,rpp,rps, and 4 coefficients of up and down waves in each layer.
  tN=4*N+4;  % there are 4N+4 unknowns to solve for incidence of s and p polarization. 
  As=zeros(tN);  
  
  MMv=eye(6);   % material matrix 
  kz0=sqrt(1-kp^2);
  [esi, epi]=findspvectors(kp,-kz0,phi,MMv);  % interface at z=0 so -kz0 is the incident wave
  [esr, epr]=findspvectors(kp,kz0,phi,MMv);   % kz0 is the reflected wave   

	      % s-polarized incidence: unknowns are rss, rps, tss, tps
  % First 4 B.C.s 
  As(1,:)=[-esr(1) -epr(1) afd(1,1) bfd(1,1) afu(1,1) bfu(1,1) zeros(1,tN-6)];   % Ex continuity
  As(2,:)=[-esr(2) -epr(2) afd(2,1) bfd(2,1) afu(2,1) bfu(2,1) zeros(1,tN-6)];   % Ey continuity 
  As(3,:)=[-esr(4) -epr(4) afd(4,1) bfd(4,1) afu(4,1) bfu(4,1) zeros(1,tN-6)];   % Hx continuity 
  As(4,:)=[-esr(5) -epr(5) afd(5,1) bfd(5,1) afu(5,1) bfu(5,1) zeros(1,tN-6)];   % Hy
  % B.C. at intermediate interfaces
  for j=1:(N-1)
    % compute phases 
    pad=exp(-1i*kzd1(j)*td(j)); pbd=exp(-1i*kzd2(j)*td(j));  % minus sign because we are going in -z direction 
    pau=exp(-1i*kzu1(j)*td(j)); pbu=exp(-1i*kzu2(j)*td(j));  x=tN-4*j-6; 
    As(4*j+1,:)=[zeros(1,4*j-2) pad*afd(1,j) pbd*bfd(1,j) pau*afu(1,j) pbu*bfu(1,j) -afd(1,j+1) -bfd(1,j+1) -afu(1,j+1) -bfu(1,j+1) zeros(1,x)];
    As(4*j+2,:)=[zeros(1,4*j-2) pad*afd(2,j) pbd*bfd(2,j) pau*afu(2,j) pbu*bfu(2,j) -afd(2,j+1) -bfd(2,j+1) -afu(2,j+1) -bfu(2,j+1) zeros(1,x)];
    As(4*j+3,:)=[zeros(1,4*j-2) pad*afd(4,j) pbd*bfd(4,j) pau*afu(4,j) pbu*bfu(4,j) -afd(4,j+1) -bfd(4,j+1) -afu(4,j+1) -bfu(4,j+1) zeros(1,x)];
    As(4*j+4,:)=[zeros(1,4*j-2) pad*afd(5,j) pbd*bfd(5,j) pau*afu(5,j) pbu*bfu(5,j) -afd(5,j+1) -bfd(5,j+1) -afu(5,j+1) -bfu(5,j+1) zeros(1,x)];
  end
				% B.C. at the last interface
  pad=exp(-1i*kzd1(N)*td(N)); pbd=exp(-1i*kzd2(N)*td(N));  % minus sign because we are going in -z direction 
  pau=exp(-1i*kzu1(N)*td(N)); pbu=exp(-1i*kzu2(N)*td(N));  
  As(tN-3,:)=[zeros(1,tN-6) pad*afd(1,N) pbd*bfd(1,N) pau*afu(1,N) pbu*bfu(1,N) -esi(1) -epi(1)];
  As(tN-2,:)=[zeros(1,tN-6) pad*afd(2,N) pbd*bfd(2,N) pau*afu(2,N) pbu*bfu(2,N) -esi(2) -epi(2)];
  As(tN-1,:)=[zeros(1,tN-6) pad*afd(4,N) pbd*bfd(4,N) pau*afu(4,N) pbu*bfu(4,N) -esi(4) -epi(4)];
  As(tN,:)=[zeros(1,tN-6) pad*afd(5,N) pbd*bfd(5,N) pau*afu(5,N) pbu*bfu(5,N) -esi(5) -epi(5)];

  Bs=[esi(1); esi(2); esi(4); esi(5); zeros(tN-4,1)];  % s-polarized incidence 
  rs=linsolve(As,Bs);
  rps=rs(2);  rss=rs(1); tss=rs(tN-1); tps=rs(tN); 
  Bp=[epi(1); epi(2); epi(4); epi(5); zeros(tN-4,1)];  % p-polarized incidence 
  rp=linsolve(As,Bp);
  rsp=rp(1);  rpp=rp(2); tsp=rp(tN-1); tpp=rp(tN); 

  if isnan(rss)||isnan(rsp)||isnan(rps)||isnan(rpp)
    disp('Some matrices were ill-conditioned, NaN discovered, removing this solution'); 
    rss=0; rsp=0; rps=0; rpp=0; tss=0; tsp=0; tps=0; tpp=0;
    return; 
  else
  end
			   % round off-extremely small values to zero 
  rss=fresnelround(rss);    tss=fresnelround(tss);
  rsp=fresnelround(rsp);    tsp=fresnelround(tsp);
  rps=fresnelround(rps);    tps=fresnelround(tps);
  rpp=fresnelround(rpp);    tpp=fresnelround(tpp);

  return;


function[kz,ad,bd,au,bu,s]=getkzEH(kp,phi,MM)

  ep=MM(1:3,1:3);
  if (kp==sqrt(ep(1,1)))||(kp==1)   % whenever kz=0 is analytically known, the computer makes errors in finding solutions very close to that value of kp
      kp=1.01*kp;
  else
  end  

  % first find kz inside the material 
  kx=kp*cos(phi); ky=kp*sin(phi);
  % if the dispersions are known don't waste your time. Use them.
  if isdiag(MM(1:3,1:3)) && isdiag(MM(4:6,4:6)) && ~sum(any(MM(1:3,4:6))) && 0
    epm=MM(1,1); mum=MM(4,4);
    disp('Isotropic');
    ks=sqrt(epm*mum-kx^2-ky^2);
    if imag(ks)>0   % so that inside the medium the solution doesnt grow with |z-z0|
       ks=-ks;
    else
    end
    kz=ks;    
  else
    % search for solution around randomly distributed initial points 
        m=1;
    Nk=25;   % number of random points  
    km=max([1,1.5*kp]);   % range of initial points from 0 to +km 
    % Increase Nk and km if solution is not found easily
    %kzr=-km+2*km*rand(1,Nk);
    kzr=-km+2*km*rand(1,Nk);
    kzi=-km+2*km*rand(1,Nk); % Im(kz) < 0 to avoid exponential growth inside the medium 
    % random distribution guarantees solution with less number of points than
    % uniform distribution of these initial points
    % fminsearch options
    opt=optimset('Display','off');
    ks=[];
    
    for j=1:Nk
        for l=1:Nk
            x0=[kzr(j) kzi(l)];
            [x, fval]=fminsearch(@(x)matcond(MM,kx,ky,x),x0,opt);  % this works better than others e.g. fmincon,
            z=x(1)+1i*x(2);
            if (fval<=-1e15) && (abs(z)<50) %  && (x(2)<=0) 
                % condition number singularity check 
                % reasonable solution (not 1e43 eg) check
                % only exponentially decaying solutions are allowed inside
                % the medium (z<0) hence Im(kz)<0 
                   
                % the following part is to not repeat solutions
                if m>1
                    if compare(ks,z)>0  % if the solution z is not obtained previously
                        ks(m)=z;
                        m=m+1;
                    else 
                        % drop the solution (uniquetol doesn't work well. we want very precise solutions here)
                    end
                else
                    ks(m)=z; m=m+1;  % first ever solution 
                end
            
            else
                % the solution obtained does not satisfy the conditions 
            end   
        end
    end
  end
			       %kz=sort(ks,'ComparisonMethod','real');
  [~, idx]=sort(real(ks)); kz=ks(idx);	
  % waves going in negative z are listed first. 
  
  if (isempty(kz))
    kz=0; au=0; bu=0; ad=0; bd=0; s=0; % no solutions found. Special attention required. (s==0)  
    return; 
  else
  end
  
  % Find the fields associated with these waves
  k=1;
  E=[];
  for j=1:length(kz)
    f=findfields(kx,ky,kz(j),MM);
    %pause; % null space : dim depends on kz and other factors
    for l=1:size(f,2)
        % Make sure that the solution is divergence free (validate the solution)
        D=MM*f(:,l); kv=[kx ky kz(j)]; d1=kv*D(1:3); d2=kv*D(4:6);  % last two are div.D=0 and div.B=0
	%d1=1e-9; d2=1e-9; 
	if (abs(d1)<1e-8) && (abs(d2)<1e-8)  % Divergence values are extremely small 
            E(:,k)=f(:,l);  
            k=k+1;
        else
        end
    end
  end

  if size(E,2)==4
    % Two eigensolutions inside the material 
    ad=E(:,1); bd=E(:,2);   % waves going downward in the slab 
    au=E(:,3); bu=E(:,4);   % waves going upward in the slab 
  else
    ad=0; bd=0; au=0; bu=0; s=0; return; 
  end

  if length(kz)==2
    kz=[kz(1), kz(1), kz(2), kz(2)];  
  else
  end
  s=1; 
  return;


function[A]=fresnelround(a)
    if abs(a)<=1e-7
        A=0;
    else
        A=a;
    end
    return;

function[E]=findfields(kx,ky,kz,MM)
% for a give kx,ky,kz solution find the eigenvectors that are
% solutions to maxwell's equations.
  kd=[0 -kz ky; kz 0 -kx; -ky kx 0;];
  MMk=[zeros(3) kd; -kd zeros(3)];
  M=MM+MMk;
  E=null(M);
  if isempty(E)
    disp('null command did not work. Using svd to compute the null basis');
    [V, D]=eigs(M);
    [m idx]=min(abs(diag(D)));
    E=V(:,idx);
  else
  end
  return;
    
function[y]=compare(A,z)
% uniquetol functionality for complex numbers.
% compare given complex number z with previously obtained numbers
% in A
  h=abs(A-z)/max(abs(A));
  if ~isempty(find(h<1e-6, 1))
      y=0;  % z is a solution obtained previously
  else
      y=1;  % z is not obtained previously
  end

  return;
  
function[C]=matcond(MM,kx,ky,x0)
    kz=x0(1)+1i*x0(2);
    kd=[0 -kz ky; kz 0 -kx; -ky kx 0;];
    MMk=[zeros(3) kd; -kd zeros(3)];
    M=MM+MMk;
    C=-1*abs(cond(M));
    % det is not a good test for testing singularity of the matrix
    % C=abs(det(M));
    return;
    
function[esm epm]=findspvectors(kp,kz,phi,MM)
    % For a given kp,kz,phi wavevector, this finds s-p polarized waves 
    % THESE MAY NOT BE SOLUTIONS their linear combinations are. 
    kx=kp*cos(phi); ky=kp*sin(phi);
    kd=[0 -kz ky; kz 0 -kx; -ky kx 0;];
    ep=MM(1:3,1:3); xi=MM(1:3,4:6); zeta=MM(4:6,1:3); mu=MM(4:6,4:6);
    % s or TE polarization 
    ev=[sin(phi) -cos(phi) 0]';
    hv=linsolve(mu,(kd-zeta)*ev);
    esm=[ev; hv;]/norm([ev; hv;]); 
    % p or TM polarization 
    hv=[sin(phi) -cos(phi) 0]'; 
    ev=linsolve(ep,-(xi+kd)*hv);
    epm=[ev; hv;]/norm([ev; hv;]);      
return;   
  
function[t]=isIllConditioned(A)
    if abs(rcond(A))<1e-15
        t=1;
    else
        t=0;
    end
return;
