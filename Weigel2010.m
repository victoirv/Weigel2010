close all;clear all;
load('OMNI_OMNI2_merged')
%getyears=Year>=2000;
%getyears=getyears+(Year<=2011);
%getyears=(getyears==2);
%KP=Kp_index(getyears);
VBS=1/2*Plasma_bulk_speed.*(abs(Bz_GSM)-Bz_GSM);
%VBS=VBS(getyears);
%DST=Dst_index(getyears);

N=50;
Na=0;
Nb=36;


[xnew,xnew2,corr,corr2,ca,ca2,cb,cb2,cc,cc2,casd,ca2sd,cbsd,cb2sd]=IRsortboot(Dst_index,VBS,Ion_density,Na,Nb);


figure; plot(flipud(cb)); hold on; plot(flipud(cb2),'r'); legend('Low Density','High Density','Location','SouthEast')
plot(cb+2*cbsd,'b:');plot(cb-2*cbsd,'b:');
plot(cb2+2*cb2sd,'r:');plot(cb2-2*cb2sd,'r:');
ylabel('DST response to VBs')
xlabel('Time Lags')
title(sprintf('N:%d Nx:%d Nf:%d',N,Na,Nb))
print -dpng dstvbs.png