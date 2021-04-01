% m_save.m 
% saves the analysis results in a .mat file if save_results=1 in m_input.
%
% version 3.2, june 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

if exist('COVARIANCE')==1,
dt = datestr(now,'mmmm dd yyyy HH MM SS FFF AM.txt')
nameres=['R',namedat];

%nameres=['R',namedat,'s']
%source='sourcetime sourceint sourcecorr';
saveres=['save ',nameres,' namedat stopdat startbg stopbg FWHM FWHMint', ...
' T0shift FWHMint psperchannel binfact cutoff entwghtstart entwghtstop', ...
' enditer left_of_max left_of_maxt0 Ntau const increment convsave bg bgstd', ...
' F f2 D MOD chi2save asave pa2saven pf2saven pfa2saven f2saven ng1 ng2', ...
' f2save COVARIANCE no_stop no_errors '];
eval(saveres);
fprintf('- saving results \n');
else
fprintf('nothing to save \n')
end


