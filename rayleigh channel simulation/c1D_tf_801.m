% 2x1D Wiener CE for OFDM systems 
% first 1D in time domain, followed by 1D in frequency domain
%
% i)  Channel: exponentially decaying profile, Jakes' Doppler spectrum
% ii) 802.11a parameters
% iii)802.11a pilot spacing (short + long + pilot-tones, Nf=12,Nt=1), 
%     with rectangular (as the standard) or diagonal pattern.
%
% i)  After time domain interpolation, more interploated channel coefficients should 
%     be used in freq. domain interpolation. (changed from Ufuk's code.)
% ii) For the pilots pattern of 802.11a, freq. first is not an appropriate order
%     since N_f = 12 is larger than the rule-of-thumb Nyquist rate. So OneD_ft.m
%     is NOT implemented FOR 802.11a parameters.
%
% Note: coded for OFDM trained CE survey project

clear all;
 
im=sqrt(-1);
randn('seed',sum(100*clock));
rand('seed',sum(100*clock));

L=64;   % OFDM block length 
L2=52;  % modulated OFDM carriers
C=16;   % length of cyclic prefix 

% system sampling frequency 
fs = 20e+6;	

% fcar: carrier separation in Hz
fcar=fs/L;  % subcarrier spacing in Hz

% Assume portable receivers with moderate doppler 
% random Doppler frequency (Hz)
% Doppler from the following values: 0Hz, 5Hz, 25Hz, 50Hz, or 100Hz
fDmax = 0; 	

% pilot pattern selection:
%   0: rect. grid (as the standard)
%   1: diagonal distri
%           xoooooo
%           oxooooo
%           ooxoooo
%           xooxooo
%           oxooxoo
%           ooxooxo
pilotpattern = 1;

% let tau_max = 5 * tau_rms;
tau_rms = 50e-9;
tau_max = 5 * tau_rms;
nprofile = 1;

Ts=4e-6;                             %ofdm block duration
% It seems the above Ts includes CP,(64+16)/20e6, since C=16 the length of cyclic prefix 

% Decision on number of taps used in Wiener filter:
%   Notice 10 taps used in wiener2d.m, i.e., mtap = 10. We use 5 taps first 
%   in time domain (for this case, assumer higher correlation in the time 
%   domain; frequency domain pilots are insufficient anyway) and then another 
%   5 in frequency domain. 
%   The matrices involved will be smaller, and thus the code will run 
%   approximately 400 % faster. 

%total number of taps in frequency domain
ftap=10;
% 10 points will participate in the calculation
%total number of taps in time domain
ttap=10;

iswitch=3;
% Assume we use weighted distances such as 
% |k-k'|f_{dmax}T_{s}+ |l-l'| tau_{max} \Delta F


% OFDM symbol/block duration
T_symbol=(L+C)/fs;
fmaxt=fDmax*T_symbol;
% normalization
taumaxf=tau_max*fcar;

SNR=0:10:40;

% pilot tones
stones=[-21 -7 7 21];
ftones=stones+32;

% short training symbols
sshort=[ -24 -20 -16 -12 -8 -4 4 8 12 16 20 24 ];  
fshort=sshort+32;

% long training symbols
slong=[-25:1:26];
flong=slong+32;

Nshort=length(fshort);
Ntones=length(ftones);
Nlong=length(flong);

% OFDM symbols in a packet (K)
Kvec=[10 20 50 100 200];
Kvec = [59];
lKvec=length(Kvec);

% MC runs, let's set it to 1 for speed. 
mcrun=10;

%mcsw = 0;   %   switch used in displaying progress

for Krun=1:lKvec
    K=Kvec(Krun);
  
    npilot=2*12+2*52+(K-4)*4; 
    coord_pilot=zeros(npilot,2);
    
    % compute pilot symbols ("rect grid")
    % 802.11a pilot pattern can be calculated in the following way
    j=0;
    for k=1:K
        if k < 3 
            for l=1:Nshort
                j=j+1; 
                coord_pilot(j,:)=[fshort(l) k];
            end
        elseif k < 5 
            for l=1:Nlong
                j=j+1;
                coord_pilot(j,:)=[flong(l) k];
            end                 
        else 
            for l=1:Ntones
                j=j+1;   
                if (pilotpattern == 0) 
                    % rect. grid (as the standard)
                    coord_pilot(j,:)=[ftones(l) k];
                elseif (pilotpattern == 1)
                    % diag. pattern 1
                    coord_pilot(j,:)=[mod(k-1+ftones(l)-7,L-12)+7 k];
                end               
            end             
        end
    end

 
    % Plot pilot positions to verify 
    % plot(coord_pilot(:,1),coord_pilot(:,2),'x');
    % xlabel('OFDM symbol (frequency) index'); ylabel('subcarrier (time) index');
    
    % the vector of all pilots
    Ntap = j;   % Ques?
    % Ntap=348=6*58, 6 is along frequency, 58 is along time
    
    % 4/5/02 "Ufuk's" comments: j from previous stage is the total number of pilot tones 

    mse=zeros(length(SNR),1); % MSE vector
    %mse_th= zeros(length(SNR),1);
  
    for drun=1:mcrun    
        % Channel H(L,K)
        fprintf('*');
        H = rayleigh(fmaxt,fDmax,fcar,tau_rms,nprofile,K,L).';
        %	H: complex output vector (channel frequency responses)
        % L: the number of carriers, K: the number of symbols
  
        for SNRcount=1:length(SNR)            
            N0_Eb=10^(-SNR(SNRcount)/20);
            sigma=sqrt(N0_Eb/2);
            noise=(randn(L,K)+im*randn(L,K));
            noises=sigma*noise;
            R=H+noises;
    
            % Channel est
            H_hat=zeros(L,K);
            H_t=zeros(L,K);
            
            for l=7:58
                % the valid frequency pins are between 7--58, the other 12
                % are zero
                for k=1:K     
                    % Search ttap nearest pilot positions in TIME domain                               
                    % In time domain filtering, there may only be two pilots 
                    % at a given frequency in this means
                    numlarge=1e6;     
                    rat=1e6*ones(Ntap,1);
         
                    for j=1:Ntap        
                        templ = (l-coord_pilot(j,1));
                        % l is freq direction, k is time direction,
                        % coord_pilot 1st column is frequency, 2nd column
                        % is time
                        %only time direction first                
                        tempk=(k-coord_pilot(j,2)); 
            
                        if (templ==0)   % count only the pilots @ the same sub carrier              
                            if (iswitch == 0)       %compute Euclidean distances
                                % iswitch is the norm definition of
                                % Euclidean distances
                                dist= tempk^2;                  
                            elseif (iswitch == 1)   %absolute distance 
                                dist=abs(tempk);          
                            elseif (iswitch == 2)   %compute weighted distances I 
                	            dist_t=abs(tempk);               
                  	            rho_t= rho_jakes(dist_t,fmaxt,nprofile);                     
                                dist=1-abs(rho_t);	     
                            elseif (iswitch == 3)   %compute weighted distances II
                                dist=abs(tempk)*fmaxt;
                            else
	                            error( 'switch > 3 not implemented')
                            end
           
	                        rat(j)=dist;                           
                        end          
                    end

                    [asortt, rbt]=sort(rat);    
                    rbval=(asortt< 1e5);   
                    t2tap=min(sum(rbval),ttap);
                    yt=zeros(t2tap,1);  
                    coord_yt=zeros(t2tap,2);  

                    for j=1:t2tap   
                        temp = [coord_pilot(rbt(j),1) coord_pilot(rbt(j),2)];    
                        coord_yt(j,:)=temp;
                        yt(j,:)=R(temp(1),temp(2)); 
                    end 
                    
                    % Find autocorrelation of pilots
                    thetat=zeros(t2tap,1); 
                    Phit=zeros(t2tap,t2tap);
                    
                    for i=1:t2tap 
                        for j=1:t2tap   % find correlation of pilots
                            distemp=coord_yt(i,:)- coord_yt(j,:);                            
                            dist_t=distemp(:,2);                                            
                            Phit(i,j) = ray_jakes(dist_t,fmaxt,nprofile);; 
                        end
                        
                        %find correlation of pilots with data
                        distth=coord_yt(i,2)-k;                                      
                        thetat(i) = ray_jakes(distth,fmaxt,nprofile);
                    end     
                    
                    % The pilots will have noise on them so include. 
                    Phit=Phit+N0_Eb*eye(t2tap);             
                    Wot=thetat'*inv(Phit);   
                    H_t(l,k)=Wot*yt; 
                end
            end     
            
            % now start freq. domain interpolation
            
          % When doing freq. domain interp., pilots include both the real pilots and those interp. from 
            % time domain. Coord_newpilot keeps coordinates of all those references (same for all k index - 
            % projection to l)
            temp = coord_pilot(:,1); 
            %tempnot = [1:6 59:64];
            coord_newpilot = [];
            for tmpind=1:length(temp)              
                TF = ismember(temp(tmpind),coord_newpilot);
                if (TF == 0)
                    coord_newpilot = [coord_newpilot temp(tmpind)];
                end
            end
            newpil_len = length(coord_newpilot);
            
            for l=7:58
                for k=1:K    
                    raf=1e6*ones(newpil_len,1);
                    for j=1:newpil_len              
                	    templ= l- coord_newpilot(j);                        
                                    
                        if (iswitch == 0)   %compute Euclidean distances,	
                            dist= templ^2;    
                        elseif (iswitch == 1) 
                            dist=abs(templ);          
                        elseif (iswitch == 2) 
                            dist_f=abs(templ);          
                            rho_f = rho_exp(dist_f,taumaxf,tau_rms,nprofile);
                            dist=1-abs(rho_f);	     
                        elseif (iswitch == 3) 	                                  
                            dist= abs(templ)*taumaxf;
                        else
                            error( 'switch > 3 not implemented');
                        end
                            
                        raf(j)=dist;                                                         
                        %end for iswitch loop
                    end

                    [asortf, rbf]=sort(raf);
                    rbval=(asortf< 1e5);   
                    f2tap=min(sum(rbval),ftap);   
                    yk=zeros(f2tap,1);  
                    coord_yk=zeros(f2tap,2);  
    
                    for j=1:f2tap   
                        temp   = coord_newpilot(rbf(j));
                        coord_yk(j,:)=[temp k];
                        %Instead of using the data use the time filtered estimates
                        yk(j,:)=H_t(temp,k);
                        %end for yk loop 
                    end       

                    thetaf=zeros(f2tap,1); 
                    Phif=zeros(f2tap,f2tap); 

                    for i=1:f2tap        
                        for j=1:f2tap     
                            distemp=coord_yk(i,:)- coord_yk(j,:);
                            dist_f=distemp(:,1);          
                            Phif(i,j) = ray_exp(dist_f,taumaxf,tau_rms,nprofile);
                        end                        

                        disfth=coord_yk(i,1)-l;
                        thetaf(i)=ray_exp(disfth,taumaxf,tau_rms,nprofile);           
                    end     

                    % The pilots will have noise on them so include. 
                    Phif=Phif+N0_Eb*eye(f2tap);
     
                    Wof=thetaf'*inv(Phif);   
                    H_hat(l,k)=Wof*yk;         
                end 
            end 
            
            for k=1:K
                err=H_hat(7:58,k)-H(7:58,k);
                temper=(norm(err,'fro')/norm(H(:,k),'fro'))^2;  % normalized MSE
                mse(SNRcount)=mse(SNRcount)+temper; 
            end           
            
            %end for SNR loop 
        end
        
        %end for Drunloop (mcrun times)
        
    end     
    
    fprintf('\n');
    mse=mse/(L2*K*mcrun);
    %mse_th=mse_th/(L2*K*mcrun);    
    
    %end for K loop
end

figure(1);
%semilogy(SNR,mse,'-*', SNR, mse_th, '-.x');
semilogy(SNR,mse,'-*');
%legend('MSE','MSE theory'); 

xlabel('SNR (dB)');
ylabel('MSE'); 
grid;

%title('Wiener 2D Estimator: MSE vs. SNR'); 
