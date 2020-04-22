function[rss rps rsp rpp]=fresnel(theta,phi,MM)

  %{
  exx=2+0.1*1i;
  exy=0.1i;
  ep=[exx exy 0; -exy exx 0; 0 0 exx;];
  mu=(1+1e-6*1i)*eye(3);
  xi=zeros(3); zeta=zeros(3);
  MM=[ep xi; zeta mu;]; 
  %}
  
  % Check the passivity constraint of the material
  Mc=-1i*MM; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);  % check if all eigenvalues are positive 
			  % pas should be 1
  if pas~=1
    disp('This is not a passive medium');
  else
  end; 

  kp=sin(theta);   % the code works for kp>1 (evanescent waves)
  % first find kz inside the material 
  kx=kp*cos(phi); ky=kp*sin(phi);

  epx=sqrt(real(MM(1,1)));  % this is needed for finding kz-solutions inside the medium. 
  ep=MM(1:3,1:3);
  if (kp==sqrt(ep(1,1)))||(kp==1)   % whenever kz=0 is analytically known, the computer makes errors in finding solutions very close to that value of kp
      kp=1.01*kp;
  else
  end
  
  if isdiag(MM(1:3,1:3)) && isdiag(MM(4:6,4:6)) && 0
    epm=MM(1,1); mum=MM(4,4);   
    ks=sqrt(epm*mum-kx^2-ky^2);
    if imag(ks)>0   % so that inside the medium the solution doesnt grow with |z-z0|
       ks=-ks;
    else
    end
   
  else
    % search for solution around randomly distributed initial points 
    m=1;
    Nk=20;   % number of random points  
    km=max([epx,1.5*kp]);   % range of initial points to search for solutions of kz 
	 % Increase Nk and km if all required solutions are not found.
    % for hyperbolic media, finding all solutions can be tricky. change the parameters here 
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
            if (fval<=-1e15) && (abs(z)<100) % && (x(2)<=0) 
                % condition number singularity check (we want precise solution else nullspace is difficult to find)
                % reasonable solution check 
                % only exponentially decaying solutions are allowed inside
                % the medium (z<0) hence Im(kz)<0 
                   
                % the following part is to omit repeated solutions
                if m>1
                    if compare(ks,z)>0  % if the solution z is not obtained previously
                        ks(m)=z;
                        m=m+1;
                    else 
                        % drop the solution (uniquetol doesn't work well.)
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
  kz=sort(ks,'ComparisonMethod','real');
  kz=kz(1:length(kz)/2);  % select the first two -z-propagating solutions
  
  % if there is no solution of kz inside the medium then all light is
  % essentially reflected (Some crazy material might do this)
  if (isempty(kz))
    rss=-1; rsp=0; rps=0; rpp=-1; disp('No solution found inside the material'); return; 
  else
  end 
  
  % It is better to find the eigenvectors inside the medium rather than s p
  % polarized vectors which are not necessarily solutions. But we can
  % nontheless decompose the eigenvectors into s and p polarizations later on if required.  
  k=1;
  E=[];
  for j=1:length(kz)
    f=findfields(kx,ky,kz(j),MM);  % null space : dim depends on kz and other factors
    for l=1:size(f,2)
        % Make sure that the solution is divergence free 
        D=MM*f(:,l); kv=[kx ky kz(j)]; d1=kv*D(1:3); d2=kv*D(4:6);  % last two are div.D=0 and div.B=0
        if (abs(d1)<1e-8) && (abs(d2)<1e-8)  % ensuring that divergence values are zero.  
            E(:,k)=f(:,l);  
            k=k+1;
        else
        end
    end
  end
  %E
  
  % s-p polarized modes in vacuum 
  MMv=eye(6);   % material matrix 
  kz0=sqrt(1-kp^2);
  % fields on the vacuum side of the geometry 
  [esi epi]=findspvectors(kp,-kz0,phi,MMv);  % interface at z=0 so -kz0 is the incident wave
  [esr epr]=findspvectors(kp,kz0,phi,MMv);   % kz0 is the reflected wave 

  if size(E,2)==2
    % Two eigensolutions inside the material 
    a=E(:,1); b=E(:,2);
    As=[esr(1) epr(1) -a(1) -b(1);
        esr(2) epr(2) -a(2) -b(2);
        esr(4) epr(4) -a(4) -b(4);
        esr(5) epr(5) -a(5) -b(5);];
    Bs=-1*[esi(1); esi(2); esi(4); esi(5);];
    rs=linsolve(As,Bs);
    rps=rs(2);  rss=rs(1);
    %rsd=(-kz0-kz)/(-kz0+kz)  % vac-eps2 interface check for validity -kz0
    % is the incident wave. Note the sign for the given convention where
    % z<0 is the semi-infinite medium 
    % p-polaried incident wave; 
    Bp=-1*[epi(1); epi(2); epi(4); epi(5);];
    rp=linsolve(As,Bp);
    rsp=rp(1); rpp=rp(2);
    %rpd=(-em*kz0-kz)/(-em*kz0+kz)  % vac-eps2 interface check for validity 
  elseif size(E,2)==1  
    rss=0; rsp=0; rps=0; rpp=0; 
    disp('Only one kz solution found. Please try to find all solutions');
  else
    disp('Special attention required for this material. There are multiple solutions');
    rss=0; rsp=0; rps=0; rpp=0; 
  end

  if isnan(rss)||isnan(rsp)||isnan(rps)||isnan(rpp)
    disp('Some matrices were ill-conditioned, NaN discovered, removing this solution'); 
    rss=0; rsp=0; rps=0; rpp=0; 
  else
  end
  % round off-extremely small values to zero 
  rss=fresnelround(rss);
  rsp=fresnelround(rsp);
  rps=fresnelround(rps);
  rpp=fresnelround(rpp);  
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
