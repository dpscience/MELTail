% m_res.m: prepare results to be visualised on line or to be saved using m_save
%
% version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

% diagonalise K, K is now of dimension s*s.
K=U'*diag(F)*U;

if sum(sum(isnan(K)))+sum(sum(isinf(K)))~0;
if warnflag==0;
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
end
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

L=diag(Yi*M*Yi');					% store eigenvalues
ent=sum(F-INIMAP-F.*log((F./INIMAP)+1e-50));		% entropy(avoid log(0)!)
ng1=-2*a*ent;						% shannon number
ng2=sum(L./(a+L));
q=a*ent-chi2/2;	
pa=exp(q)*prod((a./(a+L)).^0.5);
proda=prod((a./(a+L)).^0.5);
pf=prod((a+L).^0.5);				% posterior probability of f		% posterior probability of a
expq=exp(q);

pfa=pa*pf;		% pointwise joint posterior probability of f and a

% calculate errors on the intensities
if rank(Yi)>s-1,
no_errors=0;
Y=inv(Yi);
COVa=(1/a)*diag(F,0);
COVb=diag(F,0)*U*Y*diag(L./(a*(a+L)),0)*Y'*U'*diag(F,0);
COVARIANCE=COVa-COVb;				% covariance matrix
volF=sum(F);
end