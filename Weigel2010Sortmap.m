function deltas = Weigel2010Sortmap(x,f,sorter,N,numxcoef,numfcoef,lag, advance)
%Usage: [xnew, xnew2, corr, corr2, ca, ca2, cb, cb2, cc, cc2, casd,ca2sd,cbsd,cb2sd] = Weigel2010Sortmap(x,f,sorter,N,numxcoef,numfcoef,lag, advance)
%Where ca are the x coefficients, cb the f coefficients
%Allows for a matrix of impulses
%***Important: Assumes more data points than impulses***%
 
if (nargin < 7) || (nargin > 8)
    disp('Usage: [xnew, xnew2, corr, corr2, ca, ca2, cb, cb2, cc, cc2] = IRsort(x,f,sorter,numxcoef,numfcoef,lag)');
    disp('Where ca are the x coefficients, cb the f coefficients');
    disp('***Important: Assumes more data points than impulses***');
    error('');
end
if nargin == 6
    lag=0;
    advance=0;
end
if nargin == 7
    advance=0;
end


%Make x and f row vectors for standardization purposes
if(length(x)~=size(x,2))
    x=x';
end
if(length(f)~=size(f,2))
    f=f';
end
if(length(sorter)~=size(sorter,2))
    sorter=sorter';
end

%Index for first predicted timestep
predstart=max(numxcoef,numfcoef)+1+lag-advance;

xstart=predstart-numxcoef-lag+advance;
fstart=predstart-numfcoef-lag+advance;




numimpulses=min(size(f));
    
Nsort=10;
    len=floor(length(x)-predstart-advance);
    A=zeros(len,numxcoef+numfcoef*numimpulses+1); %+1 for mean-subtraction,
    for i=1:numxcoef
        A(1:len,i)=x(xstart+i-1:xstart+i-1+len-1);
    end
    for i=1:numfcoef
        for j=1:numimpulses
            A(1:len,i+numxcoef+(j-1)*numfcoef)=f(j,fstart+i-1:fstart+i-1+len-1);
        end
    end
    
    A(:,end)=1;
    b=x(predstart:predstart+len-1);
    A=[A(1:end,:) b'];
    
    goodrows=find(~isnan(A*zeros(size(A,2),1)));

    A=A(goodrows,:);
    b=A(:,end);
    A=A(:,1:end-1);
    
deltas=zeros(Nsort,numfcoef-1+1,3);
    cas=zeros(N,numxcoef);
    ca2s=cas;
    cbs=zeros(N,numfcoef);
    cb2s=cbs;
    ccs=zeros(N,1);
    cc2s=zeros(N,1);
    
    sortern=zeros(1,len);
    
    sorterm=zeros(Nsort,len+numfcoef);
    for i=1:Nsort
       sorterm(i,1:end-i)=sorter(fstart+i-1:fstart+i-1+len+numfcoef-1-i);
    end
    
disp('Going into loops')


for sortwindow=1:Nsort
for sortstart=0:numfcoef-sortwindow

    %Add sorter. A=[xs... fs... 1 sort]
    %for i=1:len
        %sortern(i)=mean(sorter(fstart+i+sortstart:fstart+sortstart+sortwindow+i));
    %end
    sortern=mean(sorterm(1:sortwindow,1+sortstart:len+sortstart),1);
    

    %time 0 is at numfcoef-advance, but sortstart starts at point furthest past
    %so numfcoef-advance-(sortstart+sortwindow/2) is lags from 0
    
    
    %A=[xs... fs... 1 sort b]
    %{
    for a=1:(numxcoef+numfcoef+1+1+1) %+1 for mean-normalization coef (cc, column of 1s), +1 for sorter, and +1 for column of 'b'
        A1(isnan(A1(:,a)),:)=[];
    end
    %}
    
    %Now separate into two parts via sortern (n for current iteration)
    %Make sorter only use good rows, same as A matrix
    sortern=sortern(goodrows);
    badsortrows=isnan(sortern);
    split=mean(sortern(~badsortrows));
    fprintf(' %2.1f ',split);
    lowrows=sortern<=split;
    highrows=sortern>split;
    lowrows(badsortrows)=0;
    highrows(badsortrows)=0;
    %A1=A(sortern>split,:);
    %A2=A(sortern<=split,:);
    
    
    %b=A1(:,end);
    %A1=A1(:,1:end-2);
    %b2=A2(:,end);
    %A2=A2(:,1:end-2);
    
    

    
    %Time to bootstrap

    
  
    [~,lowrowx]=find(lowrows,len);
    [~,highrowx]=find(highrows,len);
    len1=length(lowrowx);
    len2=length(highrowx);
    
    for n=1:N
        %rs=randsample(lowrowx,floor(len1/2));
        %rs2=randsample(highrowx,floor(len2/2));
        rs=randsample(lowrowx,floor(len1/2));
        rs2=randsample(highrowx,floor(len2/2));
        %Asamp=A1(rs,:);
        %A2samp=A2(rs2,:);
        
        %coef=Asamp(1:end,:)\b(rs);
        %coef2=A2samp(1:end,:)\b2(rs2);
        coef=A(rs,:)\b(rs);
        coef2=A(rs2,:)\b(rs2);
        
        if(numxcoef>0)
            cas(n,:)=coef(1:numxcoef);
            ca2s(n,:)=coef2(1:numxcoef);
        end
        cbs(n,:)=coef(numxcoef+1:end-1);
        ccs(n)=coef(end);
        
        cb2s(n,:)=coef2(numxcoef+1:end-1);
        cc2s(n)=coef2(end);
        
    end
    
    
    %ca=mean(cas);
    cb=mean(cbs);
    %ca2=mean(ca2s);
    cb2=mean(cb2s);

    %Sortstart - sortstart+window, middle should be 
    %Left: numcoef-advance-sortstart 
    %Right: numcoef-advance-sortstart-(sortwindow-1) because the width of
    %the window is actually sortwindow-1
    deltas(sortwindow,sortstart+1,:)=[numfcoef-advance-sortstart-(sortwindow-1)/2, sortwindow, sum(cb2)-sum(cb)];
    %deltas(sortwindow,sortstart+1,:)=[numfcoef-advance-sortstart-(sortwindow-1)/2, sortwindow, max(abs(cb2-cb))];
    
    %Also return both IR coeffs
    %Plot cb and cb2 every time
end
fprintf('\n%d - %d',sortwindow,split);
end

%{

    xnew=zeros(1,length(x));
    xnew(1:predstart)=xtemp(1:predstart);
    xnew2=xnew;
    
    
    %Anywhere f is nan, don't predict, just copy data
    iter=1:(length(f));
    iter=iter(iter>=predstart+advance); %Don't use copied variables
    iter=iter(iter<=length(f)-lag); %Allow space to predict
    
    for i=iter
        
        xnew(i+lag)=(xtemp(i-numxcoef:1:i-1)*ca')+(reshape(ftemp(:,i-numfcoef:1:i-1),1,[])*cb')+cc';
        xnew2(i+lag)=(xtemp(i-numxcoef:1:i-1)*ca2')+(reshape(ftemp(:,i-numfcoef:1:i-1),1,[])*cb2')+cc2';
        
    end



%Calculate correlation here to save program from needing to strip NaNs
skip=(isnan(xnew) | isnan(xtemp));
skip(1:predstart+lag)=1;
corr=corrcoef(xnew(~skip),xtemp(~skip)); %Ignore first added bit
corr=corr(1,2);
corr2=corrcoef(xnew2(~skip),xtemp(~skip)); %Ignore first added bit
corr2=corr2(1,2);

%}
