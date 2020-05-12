function[rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_top(theta,phi,MM,d)

  %{
  %phi=0;
  % single frequency calculation. Collects data for various kp-values. 
  %d=0.5;    % thickness dimensionless units.
  % gyroelectric material
  exx=2+0.1i; 
  exy=0; 
  ep=[exx 0 exy; 0 exx 0; -exy 0 exx]; 
  mu=(1+1e-6*1i)*eye(3);
  %mu=[2+0.1i 0 0.1i; 0 2+0.1i 0; -0.1i 0 2+0.1i;];  
  xi=[0 0.1 0; 0.1 0 0; 0 0 0;]; zeta=1*transpose(xi); 
  %xi=0.2*eye(3);  zeta=xi; 
  %xi=zeros(3); zeta=zeros(3);
  MM=[ep xi; zeta mu;];
  %}
  
  Mc=-1i*MM; Mc=(Mc+Mc');
  pas=all(eig(Mc)>1e-8);  % check if all eigenvalues are positive 
  if pas~=1
    disp('This is not a passive medium');
  else
  end; 
  
  % pas should be 1   
  kp=sin(theta); 
  epx=sqrt(real(MM(1,1)));
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
    %km=max([kp,3.5]);   % Use this when kz is expected to be large like 2.5
    km=max([epx,1.5*kp]);   % range of initial points from 0 to +km 
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
  
  % if there is no solution of kz inside the medium then all light is
  % essentially reflected (Some crazy unrealistic material might do this)
  if (isempty(kz))
    rss=0; rsp=0; rps=0; rpp=0; tsp=0; tps=0; tss=0; tpp=0;  
    disp('Please check material'); return; 
  else
  end 
  kzn=-1*kz;   % waves going in upward direction. 
  %kzt=horzcat(kz,kzn)  % all possible kz-wavevectors inside the medium 
  kzt=kz;
  k=1;
  E=[];
  for j=1:length(kzt)
    f=findfields(kx,ky,kzt(j),MM);
    %pause; % null space : dim depends on kz and other factors
    for l=1:size(f,2)
        % Make sure that the solution is divergence free (validate the solution)
        D=MM*f(:,l); kv=[kx ky kzt(j)]; d1=kv*D(1:3); d2=kv*D(4:6);  % last two are div.D=0 and div.B=0
	%d1=1e-9; d2=1e-9; 
	if (abs(d1)<1e-8) && (abs(d2)<1e-8)  % Divergence values are extremely small 
            E(:,k)=f(:,l);  
            k=k+1;
        else
        end
    end
  end
  %E
  % return; 
  % s-p polarized modes in vacuum 
  MMv=eye(6);   % material matrix 
  kz0=sqrt(1-kp^2);
  [esi, epi]=findspvectors(kp,-kz0,phi,MMv);  % interface at z=0 so -kz0 is the incident wave
  [esr, epr]=findspvectors(kp,kz0,phi,MMv);   % kz0 is the reflected wave 
  % the transmitted waves on the other side of interface at z=d also have
  % the same eigenvectors. 
  
  if size(E,2)==4
    % Two eigensolutions inside the material 
    a=E(:,1); b=E(:,2);   % waves going downward in the slab 
    ar=E(:,3); br=E(:,4);   % waves going upward in the slab 
    % phases till the waves reach z=-d interface 
    if length(kz)==2
        pa=exp(1i*kz(1)*(-d)); pb=pa; 
				%par=exp(-1i*kz(1)*(-d)); pbr=par;
	par=exp(1i*kz(2)*(-d)); pbr=par;
    else
        pa=exp(1i*kz(1)*(-d)); pb=exp(1i*kz(2)*(-d)); 
		    %par=exp(-1i*kz(1)*(-d)); pbr=exp(-1i*kz(2)*(-d));
	par=exp(1i*kz(3)*(-d)); pbr=exp(1i*kz(4)*(-d));
    end
    
    % incidence of s-polarized light 
    As=[esr(1) epr(1) -a(1) -b(1) -ar(1) -br(1) 0 0;
        esr(2) epr(2) -a(2) -b(2) -ar(2) -br(2) 0 0;
        esr(4) epr( 4) -a(4) -b(4) -ar(4) -br(4) 0 0;
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
    
function[ep]=epsInSb(w,Bx,By,Bz)

% B : magnetic field in Tesla 
% B=1;  % applied magnetic field in Tesla 
epsinf=15.7;
wl=3.62e13;
wt=3.39e13;
G=5.65e11;    %phonon decay rate
gp=3.39e12;    % plasma decay 
ne=1.07e17*1e6;   % doping density
me=9.1094e-31;  % electron mass 
mf=0.022*me;   % effective mass this is tunable 
wp=3.14e13;   % plasma frequency 
ec=1.6022e-19;  % electron charge
% wc=ec*B/mf;  % cyclotron frequency 
wcz=ec*Bz/mf;
wcx=ec*Bx/mf; 
wcy=ec*By/mf;

%{
wc=wcz;
exx=epsinf*(1+(wl^2-wt^2)./(wt^2-w.^2-1i*G*w)+wp^2*(w+1i*gp)./(w.*(wc^2-(w+1i*gp).^2)));
exy=1i*epsinf*wp^2*wc./(w.*((w+1i*gp).^2-wc^2)); 
ezz=epsinf*(1+(wl^2-wt^2)./(wt^2-w.^2-1i*G*w)-wp^2./(w.*(w+1i*gp))); 
ep=[exx exy 0; -exy exx 0; 0 0 ezz;];
return;
%}

w0=0;
Dw=w0^2-w^2-1i*gp*w;
M=[Dw -1i*w*wcz 1i*w*wcy; 1i*w*wcz Dw -1i*w*wcx; -1i*w*wcy 1i*w*wcx Dw;]; 
epd=epsinf*wp^2*inv(M);
elor=epsinf*(1+(wl^2-wt^2)./(wt^2-w.^2-1i*G*w));
ep=epd+elor*eye(3);

return;
