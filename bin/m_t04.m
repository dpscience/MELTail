%function [T0]=m_t0('namedat',FWHM,FWHMint,T0shift)
% m_t0.m
% function [T0]=m_t0('namedat',FWHM,FWHMint,T0shift,psperchannel)
% T0: 		scalar or row vector, temporal origin of spectrum
% namedat:	string, data file name.
% FWHM:		scalar or row vector, FWHM of the different gaussian functions 
%		composing the resolution function.
% FWHMint:	scalar or row vector, the relative intensities of the above 
% 		functions. the sum should be =1.
% T0shift:	T0shift is defined in the file m_input which m_t0 reads. 
%		it is assumed the first element is the t0 to be determined.
%		
%		if this first element is negative, m_t0 is automatically called
%		by melt and the calculated t0 is then used. if it is zero or
%		negative, m_t0 is not called but the value given is used.
%
%		the programme determines the point where the first gaussian is 
%		centered. the other elements should give the shift, in 
%		fractions of channel widths, with respect to the first gaussian.
% 		
% example: the resolution function is a sum of two gaussians, the first with 
% FWHM=250 ps, relative intensity 0.7, and the second with FWHM=150 ps, 
% relative intensity 0.3, centered 0.2 channel widths to 'the right'(increasing
% time): FWHM=[250 150], FWHMint=[0.7 0.3], T0rel=[0 0.2]. if the programme
% finds that the first gaussian is centered at 0.7, T0shift=[0.7 0.9];
%
% T0shift contains the points inside a channel where the different gaussians
% composing the resolution function are centered.
% this is just the temporal origin of the spectrum expressed as a 
% fraction of the channel width. it is a scalar if the resolution function
% is a single gaussian and a vector if it is a sum of many gaussians.
%
% the linear filter is used to calculate
% the intensity vector for different values of t0shift, and the intensity 
% vector with the least oscillations then corresponds to the best value
% of t0shift. this is an approximate value, and fine-tuning of the order
% of +- 0.1 channel widths (at most) might be needed to get to the optimal
% t0shift. how can you know that the t0shift is not yet optimal? one tell-tale
% sign is a very short lifetime component, in the first two or three channels
% of the intensity vector F. this means that you have to move t0shift towards
% 0. on the other hand oscillations in the residual spectrum given by m_plot
% in the first few channels of data indicate that t0shift must be moved
% towards 1.
%
% version 3.2, june 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

if exist('Dt0')==0
%m_input;
m_data;
%yesfigure=1;
timetaken=cputime;
end

Ntausmall=round(Ntau/1); 		% neglect very long lifetime components
firstchannel=20;
TC=zeros(NDt0+firstchannel+8,Ntausmall);		% initialize TC matrix

%Fsave=zeros(Ntausmall,11);
%fsave=zeros(Ntausmall,11);
gridtime=psperchannel*[-firstchannel:NDt0+10];
gridlambda=1./(exp(const+[1:Ntausmall]/increment));
lambdat=gridtime'*gridlambda;
lambda=ones(NDt0+firstchannel+1+10,1)*gridlambda;
t=gridtime'*ones(1,Ntausmall);
ij=0;
PP=eye(Ntausmall)/Ntausmall;
Cn=diag(Dt0);
Ci=sum(Dt0-bg)^2*eye(Ntausmall);

for t0shiftin=.1:.1:1
ij=ij+1;
t0shiftsave(ij)=t0shiftin;
TC=zeros(NDt0+firstchannel+8,Ntausmall);		% initialize TC matrix

% make analysis matrix for application of linear filter
for m=1:max(size(FWHM))
fwhm=FWHM(m);
if m==1
t0shift=t0shiftin*psperchannel;
elseif m>1
t0shift=(t0shiftin+T0shift(m))*psperchannel;
end
sigma=fwhm/(2*sqrt(log(2)));tdivlambda=(erf((t-t0shift)/sigma))./(2*lambda);
Ydivlambda=(exp(-lambdat + (t0shift*lambda) + (sigma^2/4*lambda.^2)) .*erfc(sigma/2*lambda-t/sigma+t0shift/sigma))./(2*lambda);
TC=TC+FWHMint(m)*(Ydivlambda(1:NDt0+firstchannel+8,:)-Ydivlambda(2:NDt0+firstchannel+1+8,:)-tdivlambda(1:NDt0+firstchannel+8,:)+tdivlambda(2:NDt0+firstchannel+1+8,:));
end			% m=1:max(size(FWHM))

% calculate startana and stopana as in m_mat1 and prepare analysis matrix
if ij==1
[maxmat,indmaxmat]=max(TC);

	for itau=1:Ntausmall
for ishift=1:3
starttc=indmaxmat(itau)-left_of_maxt0+ishift-2;
stoptc=starttc+NDt0-1;
tc=TC(starttc:stoptc,itau)/sum(TC(starttc:stoptc,itau));
difsqshift(ishift)=sum((tc-(Dt0-bg)/volDt0).^2);	
end
	[difsqtau(itau),indshift]=min(difsqshift);
	tshiftvec(itau)=indshift-2;
	end
		[difsqt0,indt0]=min(difsqtau);
tshift=tshiftvec(indt0);

startana=indmaxmat(indt0)-left_of_maxt0+tshift;	% first bin of analysis matrix
stopana=startana+NDt0-1;			% last bin of analysis matrix
end

TCC=TC(startana:stopana,1:Ntausmall);	% define the portion of the analysis
for ijk=1:Ntausmall			% matrix to be used.
TCC(:,ijk)=TCC(:,ijk)/sum(TCC(:,ijk));	% normalize matrix vectors
end
%%%%%

% apply m_itert0 on a coarse grid

m_itert0;
chi2t0one(ij)=chi2;	% detect oscillations in the intensity vector
clear chi2save
%%%%%

end					% end of coarse calculation
clear F

% finer calculation by interpolation

[minchi2one,iminchi2one]=min(chi2t0one);
ghi=0;

for t0shiftin=t0shiftsave(iminchi2one)-0.1:0.02:t0shiftsave(iminchi2one)+0.1
ghi=ghi+1;
t0shiftin=rem(t0shiftin+1,1);
t0shiftsavebis(ghi)=t0shiftin;
TC=zeros(NDt0+firstchannel+8,Ntausmall);		% initialize TC matrix

% make analysis matrix for application of linear filter
for m=1:max(size(FWHM))
fwhm=FWHM(m);
if m==1
t0shift=t0shiftin*psperchannel;
elseif m>1
t0shift=(t0shiftin+T0shift(m))*psperchannel;
end
sigma=fwhm/(2*sqrt(log(2)));tdivlambda=(erf((t-t0shift)/sigma))./(2*lambda);
Ydivlambda=(exp(-lambdat + (t0shift*lambda) + (sigma^2/4*lambda.^2)) .*erfc(sigma/2*lambda-t/sigma+t0shift/sigma))./(2*lambda);
TC=TC+FWHMint(m)*(Ydivlambda(1:NDt0+firstchannel+8,:)-Ydivlambda(2:NDt0+firstchannel+1+8,:)-tdivlambda(1:NDt0+firstchannel+8,:)+tdivlambda(2:NDt0+firstchannel+1+8,:));
end			% m=1:max(size(FWHM))

TCC=TC(startana:stopana,1:Ntausmall);	% define the portion of the analysis
for ijk=1:Ntausmall			% matrix to be used.
TCC(:,ijk)=TCC(:,ijk)/sum(TCC(:,ijk));	% normalize matrix vectors
end

% apply m_iter
m_itert0;
chi2t0two(ghi)=chi2;
clear chi2save

end					% end of fine calculation
%splt0shiftsavebis=t0shiftsavebis(1):0.005:t0shiftsavebis(length(t0shiftsavebis));
%splchi2t0two=spline(t0shiftsavebis,chi2t0two,splt0shiftsavebis);
%[minchi2two,iminchi2two]=min(splchi2t0two);
%t0=splt0shiftsavebis(iminchi2two);
[minchi2two,iminchi2two]=min(chi2t0two);
t0=t0shiftsavebis(iminchi2two);

fprintf('- m_t0 took')
fprintf('%5.3g',cputime-timetaken)
fprintf(' seconds of cpu time \n')
fprintf('- for the given data t0 is approximately = ')
fprintf('%4.2g \n',t0)
plot(t0shiftsavebis,chi2t0two,'o')%,splt0shiftsavebis,splchi2t0two),grid

clear Ydivlambda sigma gridtime gridlambda lambdat lambda t tdivlambda
clear firstchannel PP Ci Cn Cs t0shiftsave