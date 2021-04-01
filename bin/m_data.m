% m_data.m
% loads and edits data. calculates and subtracts background (bg).
% if the datafile is ascii, it should be called variable_name.dat so
% that the data can be identified in m_input by namedat='variable_name'. 
% an ascii datafile may contain any number of lines and columns with
% data reading from top to bottom and from left to right.
% the number of columns in each line must be the same!
% the columns must be aligned at left and separated by spaces.
% 
% version 3.2, june 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

if use_simulation < 1
    eval(['load ',namedat,'.dat']);		% load datafile: ascii format 
    DATFULL=round(eval(namedat)); 
else
    fprintf('- simulating lifetime spectrum ... \n')
    m_simulate;
    DATFULL = spectrum;
end

if use_tail_fit_model > 0
    left_of_max = 0;
    left_of_maxt0 = 0;
    
    % cut first channels ...
    [max_value,max_index]=max(DATFULL);
    start_index = max_index+right_of_max_tail_fit_analysis;
    stopdat=stopdat-right_of_max_tail_fit_analysis-max_index-1;
    DATFULL=DATFULL(start_index:size(DATFULL));
    
    startbg=startbg-start_index;
    stopbg=stopbg-start_index;
end

if binfact> 1
for ibin=1:floor(max(size(DATFULL))/binfact)
DATFUL(ibin)=sum(DATFULL((ibin-1)*binfact+1:ibin*binfact));
end
psperchannel=psperchannel*binfact;
startbg=round(startbg/binfact);
stopbg=round(stopbg/binfact)-binfact;
stopdat=round(stopdat/binfact);
left_of_maxt0=floor(left_of_maxt0/binfact);
left_of_max=max(floor(left_of_max/binfact),1);
%T0shift=T0shift/binfact;
clear DATFULL
DATFULL=DATFUL;
clear DATFUL
end

DATFULL=DATFULL(:);
bg=mean(DATFULL(startbg:stopbg));
bgstd=std(DATFULL(startbg:stopbg));
[maxdat,indmaxdat]=max(DATFULL);

% Define D, the data vector to be analysed
startD=indmaxdat-left_of_max;		% first bin of data. start 3 bins before maximum.

if use_tail_fit_model > 0
    startD=1;
end

stopD=stopdat;

D=round(DATFULL(startD:stopD));
D_bg=round((D-bg)+abs(D-bg))/2;		% fix negative values to 0
volD=sum(D-bg);
ND=length(D);

% Define Ds, the range within D used by m_iter for calculating MOD. neglects
% the parts which are essentially background

for zeroindex=26:max(size(D))
    if abs(mean(D(zeroindex-25:zeroindex))-bg) < bgstd
        break
    end
end

if zeroindex<ND-16
    zeroindex=zeroindex+15;
end

if use_tail_fit_model > 0
    zeroindex=stopD;
end

Ds=D(1:zeroindex);
volDs=sum(Ds-bg);
NDs=length(Ds);

% Define Dt0, the range within D used by m_t0 for calculating t0 and also
% by m_mat for calculating startana and stopana. shorter than Ds to save
% calculating time and also because t0 is more sensitive to variations in
% the shorter lifetimes
%startDt0=indmaxdat-left_of_maxt0;	% first bin of data.
Dt0=round(DATFULL(indmaxdat-left_of_maxt0:indmaxdat+left_of_maxt0));
volDt0=sum(Dt0-bg);
NDt0=length(Dt0);

if use_tail_fit_model > 0
    Dt0=round(DATFULL(startD:startD+1));
    volDt0=sum(Dt0-bg);
    NDt0=length(Dt0);
end

no_errors=1; 			% see m_res.m
no_stop=1;				% see m_iter.m
warnflag=0;				% see m_iter.m and m_res.m