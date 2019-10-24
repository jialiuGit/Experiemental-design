function x=coffee_house_design(Dset,n,dist)

%Copyright (c) 2019-   Jia Liu
if nargin<3
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
 
x(1,:)=max(Dset,[],1);

x(2,:)=min(Dset,[],1);

i=3;
temp=zeros(n-2,1);
xbar = setdiff(xbar,x,'rows');
while i<=n      

    distS = pdist2(x,xbar,dist);
    temp(i-2)=max(min(distS,[],1));  %find the minimal distance between all poins in x and the candidates in xbar 
    [~,J1]=find(distS==temp(i-2));
    x(i,:) = datasample(xbar(J1,:),1,1,'Replace',false);
    xbar = setdiff(xbar,x(i,:),'rows'); 
    i=i+1 ;
end
    


end
        
    


