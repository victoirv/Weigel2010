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

if(0)
[xnew,xnew2,corr,corr2,ca,ca2,cb,cb2,cc,cc2,casd,ca2sd,cbsd,cb2sd]=IRsortboot(Dst_index,VBS,Ion_density,N,Na,Nb,lag,advance);

figure; plot(fliplr(cb)); hold on; plot(fliplr(cb2),'r'); legend('Low Density','High Density','Location','SouthEast')
plot(fliplr(cb+2*cbsd),'b:');plot(fliplr(cb-2*cbsd),'b:');
plot(fliplr(cb2+2*cb2sd),'r:');plot(fliplr(cb2-2*cb2sd),'r:');
ylabel('DST response to VBs')
xlabel('Time Lags')
title(sprintf('N:%d Nx:%d Nf:%d',N,Na,Nb))
print -dpng dstvbs.png
end

VBS=VBS(getyears);
DST=Dst_index(getyears);
ION=Ion_density(getyears);
HOUR=Hour(getyears);

Sorter=HOUR;

%profile on
N=1;
tic
deltas=Weigel2010Sortmap(DST,VBS,Sorter,N,Na,Nb,lag,advance);
toc
tic
deltas=Weigel2010SortmapGPU(DST,VBS,Sorter,N,Na,Nb,lag,advance);
toc

%profile viewer

%imagesc(deltas)
%In this case, the maximum window is 10, and it goes from -advance to
%Nb-advance, but best to check because of differing window widths
%figure; imagesc([max(max(deltas(2:2:end,:,1))) min(min(deltas(2:2:end,:,1)))],[2 10],deltas(2:2:end,:,3))
%surf(repmat([1:10]',1,Nb),deltas(:,:,1),deltas(:,:,2))
[xi,yi,zi]=griddata(deltas(:,:,1),deltas(:,:,2),deltas(:,:,3),repmat(min(min(deltas(:,:,1))):0.5:max(max(deltas(:,:,1))),10,1),repmat(1:10,Nb*2-1,1)');
imagesc([min(min(xi)) max(max(xi))],[min(min(yi)) max(max(yi))], zi)
colorbar
ylabel('Window Width')
xlabel('Window Center')
print -dpng map.png

