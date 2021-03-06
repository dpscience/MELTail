% m_iter.m
% see: shukla, peter and hoffmann NIM A 335 (1993) 310.
% iteration procedure for the maximisation of Q,the difference between the 
% entropy and likelihood functions using the quantified maximum entropy 
% principle. (see Bryan in Maximum Entropy and Bayesian Methods, Kluwer,
% Dordrecht, 1990, p221) 
% calculates the lifetime intensity vector F and the model
% reconstruction MOD. a pre-filtering stage applies a linear filter
% (see hoffmann, shukla, peter, barbiellinni and manuel NIM A 335 (1993) 276.)
% to the data, filters noise and provides an a-priori solution for F.
% this pre-filtering stage stabilizes the iteration in the maxent algorithm
% which is sometimes unstable due to the non-linearity introduced by the
% exponential term.
% the iteration loop is automatically exited on convergence.
%
% version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, W?rzburg Germany
% danny.petschke@uni-wuerzburg.de

if exist('entwght')==0

% the input parameters of the fitting procedure
tolsing=cutoff;			% singular value cutoff
a=entwghtstart;

% apply inverse linear filter on data. m_iter then uses the filtered data
% Df and the intensity vector f for further calculations. this increases
% both stability and precision.
P=eye(Ntau)/Ntau;
Cn=diag(D);
Ci=sum(D-bg)^2*eye(Ntau);

Cs=TCC*Ci*P*TCC';
f=Ci*P*TCC'*inv(Cs+Cn)*(D-bg);
Df=TCC*f;
Df=round(Df+bg);
Ds=Df(1:NDs);
clear Cn Cs Ci P
%%%%%%

f=(f+abs(f))/2;				% fix negative values to 0 and create
maxf=max(f);				% a better a priori F than the flat 
for sortfi=1:Ntau			% prior with the help of the linear 
if f(sortfi)<maxf/5;			% filter
f(sortfi)=maxf/5;
end
end
F=f/sum(f);				% normalized kick-off solution.

% singular value decomposition of TCC.
% delimit the singular space of the problem (s), calculate matrices which remain
% constant during the iteration. f index signifies the full matrix, without this
% index the matrix with only the first s columns is considered.
[Vf,Sf,Uf]=svd(TCC);
s=rank(TCC,tolsing);
U=Uf(:,1:s);
S=Sf(1:s,1:s);
Vb=Vf(:,1:s);
V=Vf(1:NDs,1:s);

clear TC TCC Uf Sf Vf griddat gridres i ii iii iiii j k;
%%%%%

nit=0;					% initialize iteration count

%INIMAP=ones(Ntau,1)./m_tgrid(Ntau,const,increment)';
%INIMAP=INIMAP/sum(INIMAP);
INIMAP=ones(Ntau,1)/Ntau;		% reference for calculation of entropy
DDL=diag((Ds).\ones(NDs,1));		% 2nd derivative of likelihood function

M=S*V'*DDL*V*S;				% M constant as MOD=TC*F is linear
NEWUU=ones(Ntau,1);			% initialisation

else
a=entwght;				% entropy weight
end 					% if exist('entwght')==0

while (nit<=enditer);			% fix number of iteration steps

%%%%%%%%%%%%%%%%%%%%%%%%%% BEGIN ITERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nit=nit+1;				% iteration count

MOD=Vb*S*U'*F;				% construct model with initial F

if nit==1,				% preliminary adjustements
vol=sum(D-bg)/sum(MOD);			% adjust volumes:data/'normalized model'
F=vol*F;				% adjust F accordingly
INIMAP=vol*INIMAP;			% adjust initial map :entropy reference
UU=U\(log(INIMAP.\F));			% UU for 1st iter. UU=U\zeros(size(F));
MOD=vol*MOD;
end

MOD=round(MOD+bg);			% whole numbers
MODs=MOD(1:NDs);

% diagonalise K, K is now of dimension s*s.
K=U'*diag(F)*U;

if sum(sum(isnan(K)))+sum(sum(isinf(K)))~0;
if exist('nit2')==1
fprintf('melt aborted before reaching entwghtstop but results \n')
fprintf('can be checked using m_plot and can probably be used ! \n')
warnflag=1;
else
fprintf('convergence problems! to deal with this you can :\n')
fprintf('1) increase entwghtstart or entwghtstop so that you stop\n')
fprintf('   iterations before reaching the problem point\n')
fprintf('   (convergence is more difficult for smaller entwght)\n')
fprintf('2) or else, increase cutoff, as long as this does\n')
fprintf('   not affect the fit\n\n')
warnflag=1;
end
break
end

K=triu(K)-diag(diag(K))+triu(K)';	% correct slight assymetries
[P,DK]=eig(K,'nobalance');		% eigendecomposition of K
DK=diag(diag(DK));
for i=1:s				% correct slightly negative eigenvalues 
	if DK(i,i)< 0,			% to zero and calculate the square roots
	DK(i,i)=0;			% of the other eigenvalues
	else DK(i,i)=(DK(i,i))^0.5;
	end
end
J=DK*P'*M*P*DK;
J=triu(J)-diag(diag(J))+triu(J)';	% correct slight assymetries
[R,DJ]=eig(J,'nobalance');		% eigendecomposition of J
DJ=diag(diag(DJ));
Yi=R'*DK*P'; 					

%%%%%%%%%%%%%%%%%%%%%%%%% maximise Q
G=S*V'*(Ds.\(MODs-Ds));	% singular space gradient

% calculate X, used to calculate iterative correction on UU. xx=step length.
% a=alpha(entropy weight),b=mu(step length restriction), to be adjusted such
% that sum(INIMAP)~xx.
X=((a+b)*eye(s)+DJ)\(-a*Yi*UU-Yi*G);		
xx=X'*X;
%xxsave(nit)=xx;

DELUU=(-a*UU-G-M*Yi'*X)/(a+b);		% iterative correction on UU
UUBUF=UU;
UU=UU+DELUU;

NEWUUBUF=NEWUU;				% store the exponent of the exponential
NEWUU=U*UU;				% useful for troubleshooting

F=((exp(NEWUU)).*INIMAP);		% reconstruct F
fsave(:,nit)=sum(F);

if exist('b')==1,
stepb=xx/sum(INIMAP);		% adjust iteration step b:  xx <= sum(INIMAP)
if b<a
b=b/stepb;
else
b=b*stepb;
end
end

% test for convergence
UUG=-a*UU-G;
conv=2*((K*UUG)'*UUG)^2/(a*((K*UU)'*UU) + ((K*G)'*G))^2;

if use_tail_fit_model < 1
    if exist('sourcetime')==1,
        chi2=(sum(((MOD-D).^2)./(MOD+sourcecorr)))/ND;	% normalised chi square value
    else
        chi2=(sum(((MOD-D).^2)./MOD))/ND;		% normalised chi square value
    end
else
    chi2=(sum(((MOD-D).^2)./MOD))/ND;		% normalised chi square value
end
%if exist('sourcetime')==1,
%chi2=(sum(((MOD-D).^2)./(D+sourcecorr)))/ND;	% normalised chi square value
%else
%chi2=(sum(((MOD-D).^2)./D))/ND;		% normalised chi square value
%end
% store parameters indicating convergence history. useful for troubleshooting
chi2save(nit)=chi2;
convsave(nit)=conv;

if nit>5
if chi2<cut_off_chi2
if (mean(chi2save(nit-1:nit))<cut_off_chi2 && std(chi2save(nit-1:nit))<.02)
break
end
end
end

% if exist('nit2')==0
% if nit>25
% if (mean(chi2save(nit-5:nit))>1.25 & std(chi2save(nit-5:nit))<.02)
% fprintf('cannot reduce chi2 below ');
% fprintf('%5.3f\n',chi2save(nit))
% fprintf('try and change parameters: \n')
% fprintf('1.entwghtstart may need to be decreased \n')
% fprintf('2.FWHM or T0shift may not be optimal \n')
% no_stop=0;				% stop melt
% end
% error('exiting melt!')
% end
% end

end