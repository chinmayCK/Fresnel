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
