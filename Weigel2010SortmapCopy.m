function deltas = Weigel2010SortmapCopy(x,f,sorter,N,numxcoef,numfcoef,lag, advance)
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
    A=zeros(len,numxcoef+numfcoef*numimpulses+1+1); %+1 for mean-subtraction, +1 for sorter
        for i=1:numxcoef
        A(1:len,i)=x(xstart+i-1:xstart+i-1+len-1);
    end
    for i=1:numfcoef
        for j=1:numimpulses
            A(1:len,i+numxcoef+(j-1)*numfcoef)=f(j,fstart+i-1:fstart+i-1+len-1);
        end
    end
    
    A(:,end-1)=1;

deltas=zeros(Nsort,numfcoef-1+1);
for sortwindow=1:Nsort
for sortstart=0:numfcoef-sortwindow
tic
    %Add sorter. A=[xs... fs... 1 sort]
    for i=1:len
        A(i,end)=mean(sorter(fstart+i+sortstart:fstart+sortstart+sortwindow+i));
    end

    b=x(predstart:predstart+len-1);
    A1=[A(1:end,:) b'];
    
    %A=[xs... fs... 1 sort b]
    
    for a=1:(numxcoef+numfcoef+1+1+1) %+1 for mean-normalization coef (cc, column of 1s), +1 for sorter, and +1 for column of 'b'
        A1(isnan(A1(:,a)),:)=[];
    end
    
    %Now separate into two parts via sortern (n for current iteration)
    sortern=A1(:,end-1);
    A2=A1;
    split=mean(sortern);
    A1(sortern>split,:)=[];
    A2(sortern<=split,:)=[];
    
    
    b=A1(:,end);
    A1=A1(:,1:end-2);
    b2=A2(:,end);
    A2=A2(:,1:end-2);
    
    

    
    %Time to bootstrap
    cas=zeros(N,numxcoef);
    ca2s=cas;
    cbs=zeros(N,numfcoef);
    cb2s=cbs;
    ccs=zeros(N,1);
    cc2s=zeros(N,1);
    
    len1=size(A1,1);
    len2=size(A2,1);
    toc
    for n=1:N
        rs=randsample(len1,floor(len1/2));
        rs2=randsample(len2,floor(len2/2));
        Asamp=A1(rs,:);
        A2samp=A2(rs2,:);
        
        coef=Asamp(1:end,:)\b(rs);
        coef2=A2samp(1:end,:)\b2(rs2);
        
        if(numxcoef>0)
            cas(n,:)=coef(1:numxcoef)';
            ca2s(n,:)=coef2(1:numxcoef)';
        end
        cbs(n,:)=coef(numxcoef+1:end-1)';
        ccs(n)=coef(end);
        
        cb2s(n,:)=coef2(numxcoef+1:end-1)';
        cc2s(n)=coef2(end);
        
    end
    
    
    
    ca=mean(cas);
    cb=mean(cbs);
    ca2=mean(ca2s);
    cb2=mean(cb2s);

    
    deltas(sortwindow,sortstart+1)=max(abs(cb2-cb));
    
end
fprintf('%d ',sortwindow);
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
