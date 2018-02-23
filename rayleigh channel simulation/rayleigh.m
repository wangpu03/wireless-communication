function fk=rayleigh(fdts,fdmax,fcar,Trms,nprofile,K,L)

% This function generates a frequency-selective (2D)Rayleigh fading process. 
% Specs as follows:
%	1) Doppler power spectrum - Clarke's model;
%	2) N paths, path delay tau(n) has exponential distribution C*exp(-tau(n)/Trms);
%	3) Count only tau(n)'s in [0, 5Trms];
%	4) Equally-powered paths, according to Hoeher (Trans. VT, 1992)
%	5) Composite transmitting/receiving filter: sinc function
%
%	fk	    : complex output vector (channel frequency responses)
%	fdts	: fd_max*T_symbol
%	fdmax	: fd_max
%	fcar	: carrier separation (OFDM: fcar=1/t. t=ts-delta)
%	Trms	: rms delay
%	nprofile: channel profile: 
%		  1 - path delay distribution: exponential
%			  delay-power profile    : uniform
%			  Doppler power spectrum : Clarke's model
%
%		  2 - path delay distribution: uniform
%			  delay-power profile    : exponential
%			  Doppler power spectrum : flat
%	K 	    : time index, k>= 0, k=0; initialization
%	L 	    : number of carriers
%
% Prepared by Chengyang Li, 11/27/01

pi2 = 2*pi;
im=sqrt(-1);
taumax = 5*Trms;

if (nprofile==1)
	N = 50;
	tau = zeros(1,N);
	np = 1;	

	% collect N paths with delay in [0, taumax]
	while (np<=N)
		%temp = exprnd(Trms);
		temp = -Trms * log(1-rand(1,1));
		if (temp<=taumax)
			tau(1,np) = temp;
			np = np+1;
		end
	end

	% generate N fd's (normalized) and N random phases
	fd = fdts * cos(rand(1,N) * pi2);
	Phi = rand(1,N);	% normalized by 2pi

	h = zeros(L,K);
	for k=1:K
		for l=1:L
			%temp = exp(1i*pi2*(Phi+fd*(k-1)-fcar*(l-1)*tau));
                      % typo corrected on 4/5/02                        
			temp = exp(im*pi2*(Phi+fd*(k-1)-fcar*(l-1)*tau));
			h(l,k) = sum(temp);
		end
	end

	h = h./sqrt(N);

    % 4/5/02 comment out next two lines
	%fk = h./norm(h(:,1));	
	%fk = fk.';
    fk=h.';


% 4/5/02 Ufuk: the following accomplish pulse shaping   
else
	N = 5;
	tau_n = rand(N,1) * taumax;	% path delay uniformly distributed in [0,taumax]
	nsample = K;
	Fs = fdmax/fdts;		% Fs = 1/T_symbol - sampling rate
	[b,a] = butter(8,fdmax/(Fs/2));
		
	chn = 9000+nsample;
	alpha=diag(sqrt(exp(-tau_n/Trms)/2)) * (randn(N,chn) + 1i*randn(N,chn));
	temp = zeros(size(alpha));
	for n = 1:N
		temp(n,:) = filter(b,a,alpha(n,:));
	end
	alpha = temp(:,9001:9000+nsample);
	
	temp = zeros(L,N);
	for n = 1:N
		temp(:,n) = exp(-1i*pi2*tau_n(n)*fcar*(0:L-1)');
	end
	
	fek = temp * alpha;
	fk = fek.'/norm(fek(:,1));
end		



     


