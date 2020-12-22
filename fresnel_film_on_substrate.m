function[rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi,MM,d)

  % d is in dimensionless units: (wd/c) for a given frequency
  % phi is angle with x-axis
  % kp is parallel component of wavevector
  % MM: 6 by 6 material matrix
  % output : reflection and transmission coefficients
  
  Mc=-1i*MM; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);  % check if all eigenvalues are positive 
  if pas~=1
    disp('This is not a passive medium');
  else
  end; 
  
  % pas should be 1   
  kp=sin(theta); 

  % substrate permittivity and permeability are assumed to be constant
  ep_slab=1.5^2+1e-6*1i;
  mu_slab=1+1e-6*1i;  
  kz_slab=sqrt(ep_slab-kp^2);
  
  % first find kz inside the material and then the corresponding eigenstates
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
    E=zeros(6,4); s=0; disp('All solutions not found'); 
  end 
  %disp(E)
 
  % s-p polarized modes in vacuum 
  MMv=eye(6);   % material matrix 
  kz0=sqrt(1-kp^2);
  [esi, epi]=findspvectors(kp,-kz0,phi,MMv);  % interface at z=0 so -kz0 is the incident wave
  [esr, epr]=findspvectors(kp,kz0,phi,MMv);   % kz0 is the reflected wave 
% the transmitted waves on the other side of interface at z=d have the following eigenvectors.
  MM_slab=[ep_slab*eye(3) zeros(3); zeros(3) mu_slab*eye(3)];
  [est, ept]=findspvectors(kp,-kz_slab,phi,MM_slab);
  
  if length(kz)==2
    pa=exp(1i*kz(1)*(-d)); pb=pa;  % phase for downward going waves 
    par=exp(1i*kz(2)*(-d)); pbr=par; % phase for upward going waves
  else
    pa=exp(1i*kz(1)*(-d)); pb=exp(1i*kz(2)*(-d));  % phase for downward going waves
    par=exp(1i*kz(3)*(-d)); pbr=exp(1i*kz(4)*(-d)); % phase for upward going waves
  end
    
				% incidence of s-polarized light 
  As=[esr(1) epr(1) -ad(1) -bd(1) -au(1) -bu(1) 0 0;
      esr(2) epr(2) -ad(2) -bd(2) -au(2) -bu(2) 0 0;
      esr(4) epr(4) -ad(4) -bd(4) -au(4) -bu(4) 0 0;
      esr(5) epr(5) -ad(5) -bd(5) -au(5) -bu(5) 0 0;
      0 0 pa*ad(1) pb*bd(1) par*au(1) pbr*bu(1) -est(1) -ept(1);
      0 0 pa*ad(2) pb*bd(2) par*au(2) pbr*bu(2) -est(2) -ept(2);
      0 0 pa*ad(4) pb*bd(4) par*au(4) pbr*bu(4) -est(4) -ept(4);
      0 0 pa*ad(5) pb*bd(5) par*au(5) pbr*bu(5) -est(5) -ept(5);
     ];
  Bs=-1*[esi(1); esi(2); esi(4); esi(5); 0; 0; 0; 0;];
  rs=linsolve(As,Bs);
  rps=rs(2);  rss=rs(1); tss=rs(7); tps=rs(8); 
			     % incidence of p-polaried incident wave; 
  Bp=-1*[epi(1); epi(2); epi(4); epi(5); 0; 0; 0; 0;];
  rp=linsolve(As,Bp);
  rsp=rp(1); rpp=rp(2); tsp=rp(7); tpp=rp(8);

  if isnan(rss)||isnan(rsp)||isnan(rps)||isnan(rpp)
    disp('Some matrices were ill-conditioned, NaN discovered, removing this solution'); 
    rss=0; rsp=0; rps=0; rpp=0; tss=0; tsp=0; tps=0; tpp=0;
  else
  end
  
  % round off-extremely small values to zero 
  rss=fresnelround(rss);    tss=fresnelround(tss);
  rsp=fresnelround(rsp);    tsp=fresnelround(tsp);
  rps=fresnelround(rps);    tps=fresnelround(tps);
  rpp=fresnelround(rpp);    tpp=fresnelround(tpp);
  
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
  % We might want to use svd here 
  % [U S V]=svd(M,'econ')
  E=null(M);
  if isempty(E)
    %disp('null command did not work. Using svd to compute the null basis');
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
    
