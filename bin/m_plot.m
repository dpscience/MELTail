% m_plot.m
% plots and prints out results from melt on screen, including :
%-- lifetimes, intensities and errors (+-) on the intensities. an error of NaN
% signifies that the program was unable to calculate the error (a matrix
% inversion is involved and the matrix may sometimes be singular).
%-- chi2: the normalised chi square value for the fit. 
%-- pa: prior probablility of entweight.the 'quantified' part of maximum entropy.
% for a given cutoff (= tolsing) and a given data vector.
% MELT calculates the whole range of possible solutions, each
% with a corresponding probability and finally gives an average solution which
% is the weighted sum of each solution and the corresponding probability
%-- conv: a measure of the convergence. 2=no convergence, 0=full convergence
%-- ent :the solution 'entropy'(as compared to a reference INIMAP,see m_iter)
%-- q: entropy - likelihood.
% the graphs contain:
%-- residuals: they should be structureless, reflecting the noise. structure
% in the first few channels indicates unoptimal t0shift or FWHM.
% convergence can be checked by looking at the variations of various parameters
% as the iteration proceeds.
%-- chi square should attain a constant value ~ 1.
%-- conv should decrease rapidly after convergence ( <1e-3 is good enough)

% version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

% calculate the number distinct maxima in the spectrum, delimit the range 
% for determining lifetime center of mass and intensity.
% make top hat masks to single out each component according to the
% determined range and calculate the standard deviation(error bar) on the
% intensity of each component

FC=f2/sum(f2);
FHEC=F/sum(F);			% normalize F
tau=m_tgrid(Ntau,const,increment,ll_tau);
nit=length(chi2save);
ntau=length(F);
intau=[1:ntau];
ndat=length(D);
in=[1:ndat];
TAU=tau';

%TAUinterp=[min(TAU):5:max(TAU)]';
%[FCinterp]=interp1(TAU,FC,TAUinterp);
%Ntauinterp=max(size(TAUinterp));
%FCinterp=FCinterp./TAUinterp;
%FCinterp=FCinterp/sum(FCinterp);

% graphics
figure(1)
clf
subplot(2,2,1)
%tau,max(FC)-abs(log10(FC+1e-25))/max(abs(log10(FC+1e-25)))*max(FC+1e-25),'.y',
semilogx(tau,FC,'-m',tau,FHEC,'r');
axis([tau(1),tau(ntau),0,max(FHEC)]);
set(gca,'xtick',[1e2,1e3,1e4])
xlabel('lifetime [ps]');
ylabel('normalized intensity');
title('avg. and highest entropy solutions ')

if no_stop==1
subplot(2,2,3)
semilogx(asave,pfa2saven,'.')
set(get(gca,'Children'),'markersize',10)
xlabel('entropy weight')
ylabel('probability');
title('corresp. solutions: see figure-->')

axes('position',[0.55 0.1 0.4 0.8])
[if2save,jf2save]=find(f2save>50*max(f2save(:,1)));
if jf2save~[]
%waterfall(f2save(:,[1:min(jf2save)-1 max(jf2save)+1:max(size(pf2save))])'),view(7,7),colormap(cool),brighten(-0.8);
pcolor(flipud(f2save')),shading('interp'),colormap('jet');xlabel('log(TAU)'),ylabel('log(entweight)');
else
%waterfall(f2save'),view(7,7),colormap(cool),brighten(-0.8);
pcolor(flipud(f2save')),shading('interp'),colormap('jet');xlabel('log(TAU)'),ylabel('log(entweight)');
end
%plot(f2save+ones(max(size(f2save)),1)*(0:0.03:(min(size(f2save))-1)*0.03),'-')
%axis([tau(1),tau(ntau),0,inf]);
set(gca,'xtick',[]),set(gca,'ytick',[])
title('solution variation with ent. weight')
end

figure(2)
subplot(2,2,1)
residuals=((D-MOD).*abs(D-MOD))./MOD;
plot(in*psperchannel,residuals,in*psperchannel,zeros(ndat,1));
axis([0,ndat*psperchannel,-max(abs(residuals)),max(abs(residuals))]);
xlabel('ps');
ylabel('weighted residuals');

if log_results > 0
    clear MODDannyPrintOut_X  MODDannyPrintOut_Y MODDannyFullPrintOut;
    MODDannyPrintOut_Y=MOD';
    MODDannyPrintOut_X=in*psperchannel;
    MODDannyFullPrintOut(:,1)=MODDannyPrintOut_X;
    MODDannyFullPrintOut(:,2)=MODDannyPrintOut_Y;

    filetext = fileread('m_input.m');
    dt = datestr(now,'mmmm dd yyyy HH MM SS FFF AM.txt');

    fileID = fopen(dt,'w');
    fprintf(fileID,'%s\n\ntau [ps]\tfit [#]\n',filetext);
    for i=1:size(in,2)
        fprintf(fileID,'%f\t%f\n', MODDannyFullPrintOut(i,1), MODDannyFullPrintOut(i,2));
    end
    fclose(fileID);
end

subplot(2,2,2)
semilogy(chi2save);
xlabel('iteration number');
ylabel('chi square');

subplot(2,2,3)
semilogy(in*psperchannel,D,'.',in*psperchannel,MOD);
axis([0,ndat*psperchannel,0,max(D)]);
xlabel('ps');
ylabel('counts');
legend('data','fit')

subplot(2,2,4)
semilogy(convsave);
xlabel('iteration number');
ylabel('conv');

drawnow

datime=clock;
fprintf('%c',namedat),fprintf(' analysed on ')
fprintf('%2.0f.%2.0f.%3.0f ',datime(3),datime(2),datime(1));
fprintf('at ')
fprintf('%2.0f:%2.0f:%2.0f',datime(4),datime(5),datime(6));
fprintf('\n\n')
fprintf('average solution \n')
[av_int,av_error,av_taumean,av_width]=m_ltint(f2,tau,no_errors,COVARIANCE);
fprintf('\n\n')
fprintf('highest entropy solution \n')
[he_int,he_error,he_taumean,he_width]=m_ltint(F'+1e-10*(f2/sum(f2)),tau,no_errors,COVARIANCE);

if f2(1)/max(f2)>0.005;		% warning for unoptimal t0
fprintf('THERE IS A VERY SHORT LIFETIME COMPONENT \n')
fprintf('IF SPURIOUS, CHECK RESIDUALS AND REDUCE t0shift\n')
fprintf('\n')
end

if exist('sourcetime')==1
fprintf('\n')
fprintf('your source corrections are:\n');
fprintf('%6.1f\t\t %6.1f\n',[sourcetime;sourceint*100]);
end
fprintf('\n')

fprintf('FWHM (ps)\t intensity\t t0shift(channels):\n');
fprintf('%6.1f\t\t %6.1f\t\t %6.1f\n',[FWHM;FWHMint*100;T0shift]);
fprintf('\n')

fprintf('chi2\t\t sum(data-bg)\t sum(data-model)\t\n');
fprintf('%5.3f\t\t %d\t  %d\n\n',chi2save(nit),sum(D),sum(D)-sum(MOD));

if no_stop==1
fprintf('~ no. of components (information content/2)= ')
fprintf('%5.3f\n',ng2/2);
end
fprintf('total no. of iterations                    = ')
fprintf('%d\n',nit);