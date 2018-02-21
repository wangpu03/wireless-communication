function y=ray_jakes(x,fmaxt,nprofile)
%
%     this subroutine computes the normalized channel-correlation 
%	matrix in time 
%	nprofile: channel profile: 
%		  1 - 	path delay distribution: exponential
%			delay-power profile    : uniform
%			Doppler power spectrum : Clarke's model
%
%		  2 - 	path delay distribution: uniform
%			delay-power profile    : exponential
%			Doppler power spectrum : flat
%

pi2=2*pi;

if (nprofile==1)
	y=besselj(0,pi2*fmaxt*x);
else 
	y=pi*sinc(2*fmaxt*x);
end

