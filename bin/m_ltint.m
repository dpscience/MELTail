function [intensity,error,taumean,width]=m_ltint(F,tau,no_errors,COVARIANCE)
% interpolates the solution F from an exponential time scale onto a linear
% one and corrects it accordingly so that the summed intensities of each peak
% are the same.
% calculates these intensities as also the associated lifetimes, errors
% and peak-widths (FWHM).

sizeF=size(F);
if sizeF(1)<sizeF(2)
F=F';
end
Ntau=length(F);
TAU=tau';

FC=F/sum(F);			% normalize F

clear compnumber comparezero df dlogic lim compmax endj 
compnumber=0;
df=diff(FC);
df=[df(1);diff(FC)];
comparezero=zeros(Ntau,1);
logic=df>comparezero;
dlogic=[diff(logic);0];
lim(1)=5;				% fix limit on left
for y=1:Ntau,
if dlogic(y)==-1,
compnumber=compnumber+1;		% number of components= number of maxima
compmax(compnumber)=y;
elseif dlogic(y)==1 & compnumber>0,
lim(2*compnumber)=y;
lim(2*compnumber+1)=y+1;
end
end

if length(lim)>(2*compnumber),
lim=lim(1:2*compnumber);
elseif length(lim)<2*compnumber,
lim(2*compnumber)=Ntau-5;		% fix limit on right
end

endj=length(lim)-1;

% make top hat masks and calculate standard deviations
if no_errors==0,
indstandev=0;
for j=1:2:endj,
indstandev=indstandev+1;
mask=zeros(1,Ntau);
range=lim(j+1)-lim(j)+1;
mask(lim(j):lim(j+1))=ones(1,range);
standevf(indstandev)=(mask*COVARIANCE*mask')^.5;
end
else
standevf=ones(endj/2)*NaN;
end

%clear intensity error taumean
% calculates and prints out centre of masses of peaks (lifetime in ps),
% intensities ( % ) and corresponding errors ( +-% )
i=0;
for j=1:2:endj
i=i+1;
thispeak=FC(lim(j):lim(j+1)-1)./TAU(lim(j):lim(j+1)-1);
[thispeakmax,thispeakmaxind]=max(thispeak);
thispeakleft=thispeak(1:thispeakmaxind);
thispeakright=thispeak(thispeakmaxind:length(thispeak));
midindexleft2=min(find(thispeakleft>=thispeakmax/2));
midindexleft1=midindexleft2-1;
midindexright1=max(find(thispeakright>=thispeakmax/2))+thispeakmaxind-1;
midindexright2=midindexright1+1;
if midindexright2>max(size(thispeak)) | midindexleft1<1
width(i)=NaN;
fprintf('badly resolved component! may be an artifact \n\n')
else
widthright=interp1([thispeak(midindexright1) thispeak(midindexright2)],[TAU(midindexright1+lim(j)-1) TAU(midindexright2+lim(j)-1)],thispeakmax/2);
widthleft=interp1([thispeak(midindexleft1) thispeak(midindexleft2)],[TAU(midindexleft1+lim(j)-1) TAU(midindexleft2+lim(j)-1)],thispeakmax/2);
width(i)=abs(widthright-widthleft);
end
intensity(i)=sum(FC(lim(j):lim(j+1)));
taumean(i)=sum(FC(lim(j):lim(j+1)).*TAU(lim(j):lim(j+1)))/intensity(i);
error(i)=standevf(i)/sum(F)*100;
end
fprintf('lifetime(ps)\t intensity\t error on intensity\t FWHM(ps)\n')
fprintf('%6.1f\t\t %6.1f\t\t %6.1f\t\t\t %6.0f\n',[taumean;intensity*100;error;width])
