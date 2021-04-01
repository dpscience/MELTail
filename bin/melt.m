% melt.m
% a program for extracting lifetimes and intensities from a multiexponential
% decay curve as encountered in a positron lifetime experiment.
% see: shukla, peter and hoffmann NIM A 335 (1993) 310.
% uses the following matlab script files:
% m_input; m_data; m_mat1; m_mat2; m_scorr; m_iter; m_res m_save;
% the results may be visualised using m_plot
%
% <= version 4.0, december 1996
% Abhay SHUKLA, High Energy Group, European Synchrotron Radiation Facility
%               BP 220 F-38043, Grenoble France
% shukla@esrf.fr

% version 5.0, april 2021
% Danny Petschke, Department of Chemistry and Pharmacy, University Wuerzburg
%               Roentgenring 11, Würzburg Germany
% danny.petschke@uni-wuerzburg.de

clear 

m_input;
fprintf('- loading and preparing data and analysis matrix \n')
m_data;

if use_tail_fit_model < 1 % no tail-fitting ...
    if T0shift(1)<0
        fprintf('- finding t0 ...\n')
        m_t04;
        fprintf('- for this data t0shift is approximately = ')
        fprintf('%4.2g',t0')
        fprintf('\n')
        T0shift(1)=t0;
    end
end

m_mat;

if use_tail_fit_model < 1 % source values are far below the tail-cut-off region ...
    m_scorr;
end

fprintf('- pre-filtering and iterating... \n')
m_iter;

if (no_stop == 1 && chi2>=cut_off_chi2)
    fprintf('convergence problems!\n')
    fprintf('try and change parameters: \n')
    fprintf('1.entwghtstart may need to be decreased \n')
    if use_tail_fit_model < 1
        error('2.FWHM or T0shift may not be optimal')
    end
end

m_res;
fprintf('- initial convergence reached after')
fprintf(' %d ', nit)
fprintf('iterations with chi2=')
fprintf('%5.3f\n',chi2save(nit))

nit2=0;
entwght=a;

fprintf('- varying entropy weight...\n')
while a>entwghtstop
entwght=10^(-0.1)*entwght;		% entwght=10^(log10(entwght)-0.1);
m_iter;
m_res;
if warnflag==1;
break
end
nit2=nit2+1;
asave(nit2)=entwght;
pa2save(nit2)=pa;
pf2save(nit2)=pf;
pfa2save(nit2)=pa*pf;
f2save(:,nit2)=F;
end
if warnflag==1
sizedown=max(size(asave))-3;
asave=asave(1:sizedown);
pa2save=pa2save(1:sizedown);
pf2save=pf2save(1:sizedown);
pfa2save=pfa2save(1:sizedown);
f2save=f2save(:,1:sizedown);
F=f2save(:,sizedown);
chi2save=chi2save(1:nit-3);
convsave=convsave(1:nit-3);
chi2=chi2save(nit-3);
conv=convsave(nit-3);
nit2=sizedown;
MOD=Vb*S*U'*F;
MOD=round(MOD+bg);			% whole numbers
MODs=MOD(1:NDs);
end
pa2saven=pa2save/sum(pa2save);
pf2saven=pf2save/sum(pf2save);
pfa2saven=pfa2save/sum(pfa2save);
f2saven=f2save.*(pfa2saven'*ones(1,length(F)))';
f2=sum(f2saven');

m_plot