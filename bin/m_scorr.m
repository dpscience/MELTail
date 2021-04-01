% m_scorr.m
% subtracts the given source components from the data before analysis.
% makes the source components using the same formulae as for the 
% exponentials in the analysis matrix TC, and uses in particular the values
% for t0shift and fwhm defined earlier. no restrictions on the number of source
% components but this is always a risky operation in that an uncertainty
% about the source component will propagate through to the final result.
% if you know them, subtract them, if you don't, don't.
%
% version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

if exist('sourcetime')>0,
ndat=ND+10;			% number of channels at psperchannel ps/channel

firstchannel=round(psperchannel/2);

TCs=zeros(ndat+firstchannel,max(size(sourcetime)));
gridtime=psperchannel*[-firstchannel:ndat];
gridlambdas=1./sourcetime;
lambdast=gridtime'*gridlambdas;
lambdas=ones(ndat+firstchannel+1,1)*gridlambdas;
t=gridtime'*ones(1,max(size(sourcetime)));

for m=1:max(size(FWHM))
fwhm=FWHM(m);
t0shift=T0shift(m);
sigma=fwhm/(2*sqrt(log(2)));
tdivlambdas=(erf((t-t0shift)/sigma))./(2*lambdas);
Ydivlambdas=(exp(-lambdast + (t0shift*lambdas) + (sigma^2/4*lambdas.^2)) .*erfc(sigma/2*lambdas-t/sigma+t0shift/sigma))./(2*lambdas);
TCs=TCs+FWHMint(m)*(Ydivlambdas(1:ndat+firstchannel,:)-Ydivlambdas(2:ndat+firstchannel+1,:)-tdivlambdas(1:ndat+firstchannel,:)+tdivlambdas(2:ndat+firstchannel+1,:));
end				% m=1:max(size(FWHM))

for i=1:min(size(TCs))
TCs(:,i)=TCs(:,i)/sum(TCs(:,i))*sourceint(i)*volD;
end

TCs=round(TCs(startana:stopana,:));
sourcecorr=zeros(max(size(TCs)),1);
for i=1:min(size(TCs))
sourcecorr=sourcecorr+TCs(:,i);
end
D=D-sourcecorr;
D_bg=round((D-bg)+abs(D-bg))/2;		% fix negative values to 0


clear ndat Ydivlambdas sigma gridtime gridlambdas lambdast lambdas t tdivlambdas 
end				% exist('sourcetime')>0
