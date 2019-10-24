function x=coffee_house_design1(Dset,n,tDomain,dist)

%Copyright (c) 2019-   Jia Liu

 

if nargin<4
     dist = 'euc';
 
elseif ischar(dist)
    methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
            'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
            'spearman'; 'hamming'; 'jaccard'; 'squaredeuclidean'};
    i = find(strncmpi(dist,methods,length(dist)));
    if length(i) > 1
        error(message('stats:AmbiguousDistance', dist));
    elseif isempty(i)
        error(message('stats:UnrecognizedDistance', dist));
    else
        if i == 12 %'squaredeuclidean'
            dist = 'sqe'; % represent squared Euclidean
        else
            dist = methods{i}(1:3);
        end
    end
end

 
x=[];
xbar=Dset;

x(1,:)=min(Dset,[],1);
x_tmp=max(Dset,[],1);
 
[x_tmp,n1]=reject_design(tDomain,x_tmp);
 
if n1
     x(2,:) = x_tmp; %find initial two points with maximal distane among all pairs.
 else
     ik=0;
     distX = pdist2(x,Dset,dist);
     sort_Mat=unique(distX);
     while n1<1        
         x_tmp= Dset(distX==sort_Mat(end-ik),:);
         [x_tmp,n1]=reject_design(tDomain,x_tmp);
         ik=ik+1;
     end
     x(2,:) = datasample(x_tmp,1,1,'Replace',false); 
 end
 
 
 
 
 
 i=3;
 temp=zeros(n-2,1);
 xbar = setdiff(xbar,x,'rows');
 while i<=n             
        distS = pdist2(x,xbar,dist);         
        min_distS=unique(min(distS,[],1));
        n1=0;
        i1=0;
        while n1<1           
            temp(i-2)=min_distS(end-i1); 
            [~,J1]=find(distS==temp(i-2));
            x_tmp=xbar(J1,:);             
            [x_tmp,n1]=reject_design(tDomain,x_tmp);
            i1=i1+1;
        end
        
        x(i,:) = datasample(x_tmp,1,1,'Replace',false);
        xbar = setdiff(xbar,x(i,:),'rows'); 
        i=i+1 ;
  end
    


end
        
    

