%function testsortmap

N=1000;
x=rand(1,N);
f=rand(1,N);
sorter=mod(1:N,24);
dt=0.1;

for i=7:N
    if(sorter(i)<=12)
        x(i)=(1-dt/10)*x(i-1)+f(i-1)*dt;
    else
        x(i)=(1-dt/1)*x(i-1)+f(i-1)*dt;
    end
end

Nb=40;

deltas=Weigel2010Sortmap(x,f,sorter,50,0,Nb,0,4);

[xi,yi,zi]=griddata(deltas(:,:,1),deltas(:,:,2),deltas(:,:,3),repmat(min(min(deltas(:,:,1))):0.5:max(max(deltas(:,:,1))),10,1),repmat(1:10,Nb*2-1,1)');
imagesc([min(min(xi)) max(max(xi))],[min(min(yi)) max(max(yi))], abs(zi))

colorbar
ylabel('Window Width')
xlabel('Window Center')
print -dpng testmap.png