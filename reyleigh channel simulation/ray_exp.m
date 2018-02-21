function y=ray_exp(x,taumaxf,Trms,nprofile)

%     this subroutine computes the channel correlation matrix 
%     (freq. domain);
%	nprofil: channel profile: 
%		  1 - 	path delay distribution: exponential - tau(n) ~ exp(-tau(n)/Trms);
%			delay-power profile    : uniform - 1/sqrt(N) - Hoeher's suggestion
%			Doppler power spectrum : Clarke's model
%
%		  2 - 	path delay distribution: uniform
%			delay-power profile    : exponential
%			Doppler power spectrum : flat
%	


pi2=2*pi;

%size(x)
taumax=5*Trms;
f_car = taumaxf/taumax;

% nprofile=1 or 2 gives the same freq. correlation

y = (1- exp(-taumax*(1/Trms + 1i * pi2 * f_car * x)))./ ...
	(Trms*(1-exp(-5))*(1/Trms + 1i*pi2*f_car*x));
%y = 1./(Trms*(1/Trms + 1i*pi2*f_car*x));



