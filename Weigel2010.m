close all;clear all;
load('OMNI_OMNI2_merged')
getyears=Year>=2000;
getyears=getyears+(Year<=2011);
getyears=(getyears==2);
%KP=Kp_index(getyears);
VBS=1/2*Plasma_bulk_speed.*(abs(Bz_GSM)-Bz_GSM);


N=50;
Na=0;
Nb=48;
lag=0;
advance=6;

[xnew,xnew2,corr,corr2,ca,ca2,cb,cb2,cc,cc2,casd,ca2sd,cbsd,cb2sd]=IRsortboot(Dst_index,VBS,Ion_density,N,Na,Nb,lag,advance);


figure; plot(fliplr(cb)); hold on; plot(fliplr(cb2),'r'); legend('Low Density','High Density','Location','SouthEast')
plot(fliplr(cb+2*cbsd),'b:');plot(fliplr(cb-2*cbsd),'b:');
plot(fliplr(cb2+2*cb2sd),'r:');plot(fliplr(cb2-2*cb2sd),'r:');
ylabel('DST response to VBs')
xlabel('Time Lags')
title(sprintf('N:%d Nx:%d Nf:%d',N,Na,Nb))
print -dpng dstvbs.png

VBS=VBS(getyears);
DST=Dst_index(getyears);
ION=Ion_density(getyears);

deltas=Weigel2010Sortmap(DST,VBS,ION,N,Na,Nb,lag,advance);

%imagesc(deltas)
surf(repmat([1:10]',1,Nb),deltas(:,:,1),deltas(:,:,2))
colorbar
ylabel('Window Width')
xlabel('Window Shift')
print -dpng map.png