function[rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_2layered_top(theta,phi,MM1,MM2,td)

  kp=sin(theta);
% please define materials MM1, MM2, and geometry td for this two layered slab
  
  N=length(td);
% solve for kz and field profiles in each layer.
% there are either 2 or 4 solutions of kz. and 4 solutions for fields inside. 
  [kz,ad,bd,au,bu,s]=getkzEH(kp,phi,MM1);  % For material 1.
  if s==0
    disp('special attention required for this material');
    rss=0; rsp=0; rps=0; rpp=0; tss=0; tsp=0; tps=0; tpp=0;
    return;
  else
  end;
  % this puts material 1 in layer 1.
  kzd1(1)=kz(1); kzd2(1)=kz(2); kzu1(1)=kz(3); kzu2(1)=kz(4);  % kz(downward going)d(first or second wave)1--(layer no.)
  afd(:,1)=ad; bfd(:,1)=bd; afu(:,1)=au; bfu(:,1)=bu;  % the associated fields.
  
  [kz,ad,bd,au,bu,s]=getkzEH(kp,phi,MM2);  % For material 2.
  if s==0
    disp('special attention required for this material');
    rss=0; rsp=0; rps=0; rpp=0; tss=0; tsp=0; tps=0; tpp=0;
    return;
  else
  end;
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

    % first find kz inside the material 
  kx=kp*cos(phi); ky=kp*sin(phi);

  syms kzs;
  kd=[0 -kzs ky; kzs 0 -kx; -ky kx 0;];
  MMk=[zeros(3) kd; -kd zeros(3)];
  M=MM+MMk;

  expr=det(M);
  p=collect(expr,kzs);
  r=solve(p,kzs);
  ks=double(r); Nk=length(ks);
  [~, idx]=sort(real(ks)); kz=ks(idx); % downward and upward going waves
			       % kz=sort(ks,'ComparisonMethod','real')

  % the above solution may contain a repeated solution. Each solution has two null states like vacuum. 
  if abs(kz(1)-kz(2)) <= 1e-6
    kz=[kz(1) kz(3)];
  else
  end
  
  k=1;
  E=[];
  for j=1:length(kz)
    f=findfields(kx,ky,kz(j),MM);  
    for l=1:size(f,2)
      % Make sure that the solution is divergence free                        
      D=MM*f(:,l); kv=[kx ky kz(j)]; d1=kv*D(1:3); d2=kv*D(4:6);
      %d1=1e-10; d2=1e-10; 
      % last two are div.D=0 and div.B=0                                                         
      if (abs(d1)<1e-4) && (abs(d2)<1e-4)  % Divergence values are zero
        E(:,k)=f(:,l);
        k=k+1;
      else
      end
    end
  end
  s=1; % solution is reliable or not 
  if size(E,2)==4
    ad=E(:,1); bd=E(:,2);
    au=E(:,3); bu=E(:,4); 
  else
    ad=zeros(6,1); bd=zeros(6,1); au=zeros(6,1); bu=zeros(6,1);
    s=0; disp('All solutions not found'); 
  end 
  %disp(E)

  if length(kz)==2
    kz=[kz(1), kz(1), kz(2), kz(2)];  
  else
  end
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

function[A]=abs2(a)
  A=a.*conj(a);
  return; 
