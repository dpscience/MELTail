function [tau]=m_tgrid(nend,const,increment,ll_tau,print)
%function [tau]=m_tgrid(nend,const,increment)
% THIS FILE TO BE SAVED AS m_tgrid.m
% makes the exponentially increasing grid of lifetimes as used in m_tcmat.
%
% version 3.2, june 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

tau=ll_tau+exp(const+[1:nend]/increment);
n1=round(1*nend/6);tau1=tau(n1);
n2=round(2*nend/6);tau2=tau(n2);
n3=round(3*nend/6);tau3=tau(n3);
n4=round(4*nend/6);tau4=tau(n4);
n5=round(5*nend/6);tau5=tau(n5);
taumin=tau(1);
taumax=tau(nend);
if nargin>4
fprintf('the lowest lifetime in the grid is');
fprintf('%6.1f ',taumin);
fprintf('ps ');
fprintf(' and the highest is ');
fprintf('%6.1f',taumax);
fprintf('ps\n');
fprintf('at the following lifetime values on the grid (ps):\n')
fprintf('%6.1f\t %6.1f\t %6.1f\t %6.1f\t %6.1f\t \n',tau(n1),tau(n2),tau(n3),tau(n4),tau(n5))
fprintf('the interval between neighbouring lifetimes is (ps):\n')
fprintf('%6.1f\t %6.1f\t %6.1f\t %6.1f\t %6.1f\t \n',tau(n1)-tau(n1-1),tau(n2)-tau(n2-1),tau(n3)-tau(n3-1),tau(n4)-tau(n4-1),tau(n5)-tau(n5-1))
end
