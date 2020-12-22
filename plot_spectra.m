clear all;
c=2.997928e8;  % speed of light
td=1e-6;    % thickness of InSb film on top of glass of const. n=2.25

theta=pi/4;  % angle with z-axis
phi=pi/4;  % angle with x-axis


%%%%% SKL-1 %%%%%%
Nw=100;
w=linspace(0.1,10,Nw)*1e13;    % frequency in rad/s
Bx=0; By=0; Bz=0;  % InSb placed in magnetic field. Zero magnetic field --- SKL1

for j=1:Nw
	disp(j);
  ep=epsInSb(w(j),Bx,By,Bz);
  mu=(1+1e-6*1i)*eye(3);
  xi=zeros(3); zeta=zeros(3);
  MM=[ep xi; zeta mu;];
  d=w(j)/c*td;
  kp=sin(theta);
  % emissivity calculation in the direction (theta,phi)
  [rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi+pi,MM,d);
  Rpp=1/4*abs(rss+rpp+1i*(rsp-rps)).^2;
  Rmp=1/4*abs(rss-rpp+1i*(rsp+rps)).^2; 
  Rmm=1/4*abs(rss+rpp-1i*(rsp-rps)).^2;
  Rpm=1/4*abs(rss-rpp-1i*(rsp+rps)).^2;
  ap(j)=1-Rpm-Rpp;  % RCP emissivity  
  am(j)=1-Rmp-Rmm;  % LCP emissivity 
  % absorptivity calculation in the direction (theta,phi)
  [rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi,MM,d);
  Rpp=1/4*abs(rss+rpp+1i*(rsp-rps)).^2;
  Rmp=1/4*abs(rss-rpp+1i*(rsp+rps)).^2; 
  Rmm=1/4*abs(rss+rpp-1i*(rsp-rps)).^2;
  Rpm=1/4*abs(rss-rpp-1i*(rsp+rps)).^2;
  abp(j)=1-Rmp-Rpp; % RCP absorptivity 
  abm(j)=1-Rpm-Rmm; % LCP absorptivity 
end
figure(1); plot(w,ap,'r',w,am,'b','LineWidth',1.5);
legend('\eta_{(+)}','\eta_{(-)}'); xlabel('Frequency (rad/s)'); ylabel('Emissivity (\eta)');

figure(2); plot(w,abp,'r',w,abm,'b','LineWidth',1.5);
legend('\alpha_{(+)}','\alpha_{(-)}'); xlabel('Frequency (rad/s)'); ylabel('Absorptivity (\eta)');

%return;

%%%%%% SKL-2 %%%%%%%

Bx=0; By=0; Bz=1;  % InSb placed in magnetic field along z-axis (perp. to slab)
		  % it is better to calculate this part using parpool

for j=1:Nw
  ep=epsInSb(w(j),Bx,By,Bz);
  mu=(1+1e-6*1i)*eye(3);
  xi=zeros(3); zeta=zeros(3);
  MM=[ep xi; zeta mu;];
  d=w(j)/c*td;
  kp=sin(theta);
  % emissivity calculation in the direction (theta,phi)
  [rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi+pi,MM,d);
  Rpp=1/4*abs(rss+rpp+1i*(rsp-rps)).^2;
  Rmp=1/4*abs(rss-rpp+1i*(rsp+rps)).^2; 
  Rmm=1/4*abs(rss+rpp-1i*(rsp-rps)).^2;
  Rpm=1/4*abs(rss-rpp-1i*(rsp+rps)).^2;
  ap(j)=1-Rpm-Rpp;  % RCP emissivity  
  am(j)=1-Rmp-Rmm;  % LCP emissivity 
  % absorptivity calculation in the direction (theta,phi)
  [rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi,MM,d);
  Rpp=1/4*abs(rss+rpp+1i*(rsp-rps)).^2;
  Rmp=1/4*abs(rss-rpp+1i*(rsp+rps)).^2; 
  Rmm=1/4*abs(rss+rpp-1i*(rsp-rps)).^2;
  Rpm=1/4*abs(rss-rpp-1i*(rsp+rps)).^2;
  abp(j)=1-Rmp-Rpp; % RCP absorptivity 
  abm(j)=1-Rpm-Rmm; % LCP absorptivity 
end
figure(3); plot(w,ap,'r',w,am,'b','LineWidth',1.5);
legend('\eta_{(+)}','\eta_{(-)}'); xlabel('Frequency (rad/s)'); ylabel('Emissivity (\eta)');

figure(4); plot(w,abp,'r',w,abm,'b','LineWidth',1.5);
legend('\alpha_{(+)}','\alpha_{(-)}'); xlabel('Frequency (rad/s)'); ylabel('Absorptivity (\eta)');

%return; 
%%%%%%%% SKL-3  %%%%%%%%%%%

Bx=2; By=0; Bz=0;  % InSb placed in magnetic field along x axis: SKL-3 

Nw=100;
w=linspace(2,5,Nw)*1e13; % we zoom in over smaller freq. range because change wrt B-field is small

for j=1:Nw
  ep=epsInSb(w(j),Bx,By,Bz);
  mu=(1+1e-6*1i)*eye(3);
  xi=zeros(3); zeta=zeros(3);
  MM=[ep xi; zeta mu;];
  d=w(j)/c*td;
  kp=sin(theta);
  % emissivity, absorptivity calculation 
  [rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi+pi,MM,d);
  Rpp=1/4*abs(rss+rpp+1i*(rsp-rps)).^2;
  Rmp=1/4*abs(rss-rpp+1i*(rsp+rps)).^2; 
  Rmm=1/4*abs(rss+rpp-1i*(rsp-rps)).^2;
  Rpm=1/4*abs(rss-rpp-1i*(rsp+rps)).^2;
  ap(j)=1-Rpm-Rpp;  % RCP emissivity in (theta,phi) direction 
  am(j)=1-Rmp-Rmm;  % LCP emissivity in (theta,phi) direction
  abp_c(j)=1-Rmp-Rpp;  % RCP absorptivity in (theta,phi+pi) conjugate direction
  abm_c(j)=1-Rpm-Rmm;  % LCP absorptivity in (theta,phi+pi) conjugate direction
  % absorptivity, absorptivity calculation  
  [rss, rps, rsp, rpp, tss, tps, tsp, tpp]=fresnel_film_on_substrate(theta,phi,MM,d);
  Rpp=1/4*abs(rss+rpp+1i*(rsp-rps)).^2;
  Rmp=1/4*abs(rss-rpp+1i*(rsp+rps)).^2; 
  Rmm=1/4*abs(rss+rpp-1i*(rsp-rps)).^2;
  Rpm=1/4*abs(rss-rpp-1i*(rsp+rps)).^2;
  abp(j)=1-Rmp-Rpp; % RCP absorptivity in (theta,phi) direction
  abm(j)=1-Rpm-Rmm; % LCP absorptivity in (theta,phi) direction
  ap_c(j)=1-Rpm-Rpp; % RCP emissivity in (theta,phi+pi) conjugate direction 
  am_c(j)=1-Rmp-Rmm; % LCP emissivity in (theta,phi+pi) conjugate direction  
end
figure(5); plot(w,ap,'r',w,am,'b','LineWidth',1.5);
hold on;   plot(w,ap_c,'r--',w,am_c,'b--','LineWidth',1.5);
legend('\eta_{(+)}','\eta_{(-)}','\eta_{(+)} conj direction','\eta_{(-)} conj direction');
xlabel('Frequency (rad/s)'); ylabel('Emissivity (\eta)');

figure(6); plot(w,abp,'r',w,abm,'b','LineWidth',1.5);
hold on;  plot(w,abp_c,'r--',w,abm_c,'b--','LineWidth',1.5);
legend('\alpha_{(+)}','\alpha_{(-)}','\alpha_{(+)} conj direction','\alpha_{(-)} conj direction');
xlabel('Frequency (rad/s)'); ylabel('Absorptivity (\eta)'); 

return; 




