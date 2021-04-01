% m_tcmat.m
% makes the exponential matrix TC using the analytical expression 
% for the convolution between the exponential and the resolution curve
% (see for example: kirkegaard and eldrup, comp. phys. comm. 3 (1972) 240-255)
% which may be a sum of more than one gaussians, not centered at the same point.
% if you wish to save a particular analysis matrix (if you use the same one
% constantly for example while analysing simulated data) you must uncomment 
% the line corresponding to the name of the saved matrix in m_input.
% if this line is uncommented while running melt, the analysis matrix  will be
% loaded and used in the analysis instead of a new one being calculated.
% for T0shift=0.5 the resolution function is centered in the middle of a 
% channel, for T0shift =0 on the left edge, and for T0shift=1 on the right edge
% so T0shift=0 and T0shift=1 are equivalent.
% 
% version 3.2, june 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

if exist('Dt0')==0
m_input;
m_data;
end

ndat=ND+10;	% number of channels at psperchannel ps/channel
%firstchannel=round(psperchannel/2);
firstchannel=20;

if use_tail_fit_model > 0
    firstchannel=1;
end

TC=zeros(ndat+firstchannel,Ntau);		% initialize TC matrix
gridtime=psperchannel*[-firstchannel:ndat];
gridlambda=1./(ll_tau+exp(const+[1:Ntau]/increment));
lambdat=gridtime'*gridlambda;
lambda=ones(ndat+firstchannel+1,1)*gridlambda;
t=gridtime'*ones(1,Ntau);

if use_tail_fit_model  < 1
    for m=1:max(size(FWHM))
        fwhm=FWHM(m);
        
        %if m>1
        %    t0shift=(T0shift(1)+T0shift(m))*psperchannel;
        %else
        %    t0shift=T0shift(1)*psperchannel;
        %end
    
        t0shift=T0shift(m)*psperchannel;
        
        sigma=fwhm/(2*sqrt(log(2)));		% !sigma here includes a sqrt(2) factor! 
        tdivlambda=(erf((t-t0shift)/sigma))./(2*lambda);
        Ydivlambda=(exp(-lambdat + (t0shift*lambda) + (sigma^2/4*lambda.^2)) .*erfc(sigma/2*lambda-t/sigma+t0shift/sigma))./(2*lambda);
        TC=TC+FWHMint(m)*(Ydivlambda(1:ndat+firstchannel,:)-Ydivlambda(2:ndat+firstchannel+1,:)-tdivlambda(1:ndat+firstchannel,:)+tdivlambda(2:ndat+firstchannel+1,:));
    end
else % tail fit model ...
    t0shift=T0shift(1)*psperchannel;
   
    sigma=1.;%0.5*psperchannel;
    tdivlambda=(erf((t-t0shift)/sigma))./(2*lambda);
    Ydivlambda=(exp(-lambdat + (t0shift*lambda) + (sigma^2/4*lambda.^2)) .*erfc(sigma/2*lambda-t/sigma+t0shift/sigma))./(2*lambda);
    TC=TC+1.*(Ydivlambda(1:ndat+firstchannel,:)-Ydivlambda(2:ndat+firstchannel+1,:)-tdivlambda(1:ndat+firstchannel,:)+tdivlambda(2:ndat+firstchannel+1,:));
end

if savetcmat==1,
savemat=['save ',nameana,' ndat Ntau fwhm psperchannel t0shift const increment TC '];
eval(savemat);
end
clear tdivlambda Ydivlambda t lambdat lambda gridlambda gridtime
