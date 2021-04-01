% m_input.m
% user defined input. this file serves as input for melt, for m_to and for 
% m_tcmat at the same time and is in principle the only file you need to edit
% and modify for running these programs.  more details in the file itself.
% each input is accompanied with a brief description. the second half of
% this input file contains either variables that do not need to be altered
% ordinarily (enditer, b, left_of_max, left_of_maxt0), or those which might
% be altered from time to time (Ntau, const, increment) according to the kind
% of data being analysed (metallic lifetimes, polymers, etc).
%
% version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

%%%%%%%%%%%%%%%%%%%%%%%% INPUT TO BE INSERTED BY USER 

namedat='filename';	% datafile+variable name
stopdat=7000; % last bin of data
startbg=6500; % first bin for background calculation
stopbg=7000; % last bin for background calculation

% Setting 'use_tail_fit_model' = 1 restricts the optimization procedure only to the
% tail-region of the data, i.e. the region where the longest components are predominant.
% The tail-fitting optimization procedure ignores the instrumental response, meaning
% that the convolution with the resolution function is treated as neglegible
% and, therefore, the model reduces to exclusively a distribution of lifetimes without convoluted
% kernel (params 'FWHM', 'T0shift', 'FWHMint', 'left_of_max' and
% 'left_of_maxt0' will be ignored).
% The beginning of the tail-region of interest with respect to the maximum is set by 'right_of_max_tail_fit_analysis' in units of channels. 

use_tail_fit_model = 1; % 0 = false, 1 = true
right_of_max_tail_fit_analysis = 250; % the channel to the right of the maximum where the tail-region of interest starts

% Setting 'use_simulation' = 1 applies the optimization procedure on
% simulated data ('namedat' is ignored). The simulation generates lifetime
% spectra based on Gaussian distributions of characteristic lifetimes defined by the
% params: 'tau_mean_sim', 'tau_stddev_sim', 'Int_mean_sim'
% incorporating the instrumental reponse as given by 'FWHM', 'FWHMint',
% 'T0shift'.

use_simulation = 1; % 0 = false, 1 = true
sim_smoothness = 1000; % defines the smoothness of the lifetime grid in the simulation mode

countsInSpectrum = 12000000; % number of counts in the spectrum
constBkgrd = 15; % constant background of the spectrum

tau_mean_sim = [170. 380. 1500.]; % mean of the Gaussian distributed charact. lifetimes
tau_stddev_sim = [5. 5. 200.]; % standard deviations of the Gaussian distributed charact. lifetimes
Int_mean_sim = [0.15 0.15 0.70]; % contributions of the Gaussian distributed charact. lifetimes (must equal to 1.0)

% input values used by m_tcmat for the analysis matrix TC. The resolution
% function used may be a sum of any number of gaussians, not centered at
% the same point.

FWHM=[200.]; % full width at half max for each gaussian in ps
T0shift=[0.0];	% t0 for each component. all values are in channel
			% units (ps/psperchannel). the first element by
			% convention determines the temporal origin of the
			% spectrum  and should be between 0 and 1. the others
			% give the shift RELATIVE to this first component.			 

FWHMint=[1]; % intensity of each component. the sum should = 1.
psperchannel=5.0;	
binfact=1; % groups every binfact channels together before
			% analysis. useful for long data especially
			% when FWHM > 5-7 * psperchannel. If you change this
			% T0shift ALSO CHANGES. For example if T0shift=[0.5 3]
			% when binfact is 1, it can be either [0.25 1.5] or
			% [0.75 1.5] for binfact=2, according to the binning.
			
% if no source corrections are to be applied, CHANGE the next TWO lines to
% comments by typing in a % sign at the beginning.
% source correction will be completely ignored if the tail-fitting is
% enabled as optimization procedure ('use_tail_fit_model' == 1).
%sourcetime=[400];	% source correction lifetimes in ps
%sourceint=[.03];	% source correction intensities

cutoff=1e-5; % singular value cutoff
entwghtstart=1e-4; % entropy weight for initial convergence
entwghtstop=1e-8; % lower limit for entropy weight

cut_off_chi2 = 1.5; % minimum chi-square value used to stop the optimization

savetcmat=0; % 1/0 : save/do not save TC matrix calculated by m_tcmat.
% nameana=0;		% analysis matrix. when feeding input
% to melt, UNCOMMENT only if a previously calculated analysis matrix is to be
% loaded and used. when feeding input to m_tcmat UNCOMMENT always as
% nameana specifies the name of the matrix in case it is to be saved.

log_results = 0; % 1/0 : saves the results and input of each single run.

%%%%%%%%%%%%%% the following can be left at their default values %%%%%%%%%%%%%%%

enditer=100; % maximum number of iterations. 

left_of_max=0;	% first bin of data. start left_of_max bins before the maximum. default value=3.
left_of_maxt0=9;	% first bin of data for determination of t0. default value=9.
% to be changed according to the range of lifetimes investigated. the grid can
% be generated and examined using m_tgrid.
Ntau=130; % define analysis lifetime range
			    % should be <= ntau in the matrix TC!
const=4; % determines the first lifetime on the analysis grid
increment=29; % determines the spacing of the lifetimes on this grid % 29

ll_tau = 100.; % [ps] % offset of the lifetime to be added to 'const', i.e. the lifetimes where the tc_mat actually starts

% initial value for iteration step. change ONLY if serious convergence problems!
b=0; % default=0 and in this case b=0 all through the iteration.
chanshift=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 