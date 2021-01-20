function[rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_bottom(theta,phi,MM,d)

  %{
  % uncomment for the example material below  
  %phi=0;
  % single frequency calculation. Collects data for various kp-values. 
  %d=0.5;    % thickness dimensionless units.
  % gyroelectric material
  exx=2+0.1i; 
  exy=0; 
  ep=[exx 0 exy; 0 exx 0; -exy 0 exx]; 
  mu=(1+1e-6*1i)*eye(3);
  %mu=[2+0.1i 0 0.1i; 0 2+0.1i 0; -0.1i 0 2+0.1i;];  
  %xi=[0 0.1 0; 0.1 0 0; 0 0 0;]; zeta=1*transpose(xi); 
  %xi=0.2*eye(3);  zeta=xi; 
  xi=zeros(3); zeta=zeros(3);
  MM=[ep xi; zeta mu;];
  %}

  Mc=-1i*MM; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);  % check if all eigenvalues are positive 
  if pas~=1
    disp('This is not a passive medium');
  else
  end;
  
  kp=sin(theta); 
  phi=phi+pi;   % incidence along (theta, phi) means the parallel wavevector makes an angle phi+pi with x-axis
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
    E=zeros(6,4); s=0; disp('All solutions not found'); 
  end 
  %disp(E)
  
  % s-p polarized modes in vacuum 
  MMv=eye(6);   % material matrix 
  kz0=sqrt(1-kp^2);
  [esi, epi]=findspvectors(kp,kz0,phi,MMv);  % interface at z=-d so kz0 is the incident wave from the bottom 
  [esr, epr]=findspvectors(kp,-kz0,phi,MMv);   % -kz0 is the reflected wave for incidence from negativ inf.  
  % the transmitted waves on the other side of interface at z=0 also has same eigenvectors. 
  
  if size(E,2)==4
    % Two eigensolutions inside the material 
    a=E(:,1); b=E(:,2);   % waves going downward in the slab 
    ar=E(:,3); br=E(:,4);   % waves going upward in the slab 
			  % phases till the waves reach z=0 interface
    if length(kz)==2
        pa=exp(1i*kz(1)*(d)); pb=pa; 
	par=exp(1i*kz(2)*(d)); pbr=par;
    else
        pa=exp(1i*kz(1)*(d)); pb=exp(1i*kz(2)*(d)); 
	par=exp(1i*kz(3)*(d)); pbr=exp(1i*kz(4)*(d));
    end

    %{
    if length(kz)==1
        pa=exp(1i*kz(1)*(d)); pb=pa; 
        par=exp(-1i*kz(1)*(d)); pbr=par;
    else
        pa=exp(1i*kz(1)*(d)); pb=exp(1i*kz(2)*(d)); 
        par=exp(-1i*kz(1)*(d)); pbr=exp(-1i*kz(2)*(d));
    end
    %}
    
    % incidence of s-polarized light 
    As=[esr(1) epr(1) -a(1) -b(1) -ar(1) -br(1) 0 0;
        esr(2) epr(2) -a(2) -b(2) -ar(2) -br(2) 0 0;
        esr(4) epr(4) -a(4) -b(4) -ar(4) -br(4) 0 0;
        esr(5) epr(5) -a(5) -b(5) -ar(5) -br(5) 0 0;
        0 0 pa*a(1) pb*b(1) par*ar(1) pbr*br(1) -esi(1) -epi(1);
        0 0 pa*a(2) pb*b(2) par*ar(2) pbr*br(2) -esi(2) -epi(2);
        0 0 pa*a(4) pb*b(4) par*ar(4) pbr*br(4) -esi(4) -epi(4);
        0 0 pa*a(5) pb*b(5) par*ar(5) pbr*br(5) -esi(5) -epi(5);
        ];
    Bs=-1*[esi(1); esi(2); esi(4); esi(5); 0; 0; 0; 0;];
    rs=linsolve(As,Bs);
    rps=rs(2);  rss=rs(1); tss=rs(7); tps=rs(8); 
    %rsd=(-kz0-kz)/(-kz0+kz)  % vac-eps2 interface check for validity -kz0
    % is the incident wave. Note the sign for the given convention where
    % z<0 is the semi-infinite medium 
    
    % incidence of p-polaried incident wave; 
    Bp=-1*[epi(1); epi(2); epi(4); epi(5); 0; 0; 0; 0;];
    rp=linsolve(As,Bp);
    rsp=rp(1); rpp=rp(2); tsp=rp(7); tpp=rp(8);
    %rpd=(-em*kz0-kz)/(-em*kz0+kz)  % vac-eps2 interface check for validity 
  else
    disp('Special attention required for this material.');
    th=asin(kp)/pi; disp(th); disp(phi); disp('kz='); disp(kz);
    rss=0; rsp=0; rps=0; rpp=0; tss=0; tsp=0; tps=0; tpp=0;
  end

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
    
