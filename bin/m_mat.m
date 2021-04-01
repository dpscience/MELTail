% m_mat.m
% loads analysis matrix or calculates it using m_tcmat.
% performs preliminary analysis for the determination of startana and stopana
%
% version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

if exist('nameana')>0,
eval(['load ',nameana]);
else
m_tcmat;
end	
[maxmat,indmaxmat]=max(TC);

Ntauhalf=round(Ntau/2);		% neglect longer lifetimes in pre-treatment

	for itau=1:Ntauhalf
for ishift=1:5
starttc=indmaxmat(itau)-left_of_maxt0+ishift-3;

if use_tail_fit_model > 0
    starttc=indmaxmat(itau)+ishift;
end
    
stoptc=starttc+NDt0-1;

tc=TC(starttc:stoptc,itau)/sum(TC(starttc:stoptc,itau));

difsqshift(ishift)=sum((tc-(Dt0-bg)/volDt0).^2);	

end
	[difsqtau(itau),indshift]=min(difsqshift);
	tshiftvec(itau)=indshift-3;
	end
		[difsqt0,indt0]=min(difsqtau);
tshift=tshiftvec(indt0);

% first bin of analysis matrix
startana=indmaxmat(indt0)-left_of_max+tshift+chanshift;

if use_tail_fit_model > 0
    startana=1;
end

stopana=startana+ND-1;				% last bin of analysis matrix
[tau]=m_tgrid(Ntau,const,increment,ll_tau);
TCC=TC(startana:stopana,1:Ntau);	% define the portion of the analysis
for ijk=1:Ntau				% matrix to be used.
TCC(:,ijk)=TCC(:,ijk)/tau(ijk);	% normalize matrix vectors
end