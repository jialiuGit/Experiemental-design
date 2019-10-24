function Min_dis_design=minDis_reject(n,dim,tDomain,delta,varargin)
%Copyright (c) 2019-   Jia Liu

Min_dis_design.fh.mindis_short=@mindis_short;

nargs = nargin;

if nargs < 4
    error(message('stats:minDis_reject:TooFewInputs'));
end
    
if nargin<5
     dist = 'euc';
else 
    dist=varargin{1};
end

if nargin <6
    weights =[];     
else
    weights=varargin{2};
    weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
end




x=rand(n,dim);                %dim(x)=n

Min_dis_design.xnoincPr=Min_dis_design.fh.mindis_short(n,dim,delta,x); %without incPr
  
 
[x,k]=reject_design(tDomain,x); %dim(x)=k

if isstruct(x)
     strX=x.x;
     Xorg=x.xorg;
     
     clear x
     
     x=strX; 
     
     clear strX
end
    

 
n11=n-k;
ik=0;
while ik<n11 
    x_one=rand(1,dim);
    
    [x_one,n111]=reject_design(tDomain,x_one);
    if n111
        ik=ik+1;
        if isstruct(x_one)
            strX=x_one.x;
            
            Xorg=[Xorg; x_one.xorg];
            clear x_one
            x_one=strX; 
            
            clear strX
        end      
        x=[x;x_one];  
    end
    
end
 



i=1;
nreplace=0;

% tic
while i<=n
    
     
    x2=x;
    x2(i,:)=[];
    if nargin <6
        distS = pdist2(x(i,:),x2,dist);
    else        
        distS = pdist2(x(i,:),x2,@(Xi,Xj) weuc(Xi,Xj,weights));
    end
        
        
        
    
    
    if min(distS) >= delta         
        i=i+1;         
    else
       
        n22=0;
         
        while n22<1  
            x_one=rand(1,dim);
            [xtmp,n22]=reject_design(tDomain,x_one);             
        end
        
        if isstruct(xtmp)
            strX=xtmp.x;
            
            Xorg=[Xorg; xtmp.xorg];
            clear xtmp
            xtmp=strX; 
            
            clear strX
        end  
            
 
        x(i,:)=xtmp;
        nreplace=nreplace+1;
        if nreplace == 1e5
            error('The threshold of the distance is not proper, pls reduce the threshold')
        end
            
            
    end    
end
 
Min_dis_design.delta=delta;
Min_dis_design.x=x;

if exist('Xorg','var')
    Min_dis_design.xorg=Xorg(1:n);
end
if  isempty(nreplace < n)
  fprintf(['Find ' num2str(nreplace) ' optimal new sampling locations in all the candidates.'] )
else 
  fprintf(['Find ' num2str(n) ' optimal new sampling locations in the candidates.'] )
end
  
end
 

function x=mindis_short(n,dim,delta,x)

 
 
i=1;
nreplace=0;
while i<=n
    x2=x;
    x2(i,:)=[];
    distS = pdist2(x(i,:),x2,'euc');   
    
    if min(distS) >= delta
        i=i+1;  
         
    else
        x(i,:)=rand(1,dim);
        nreplace=nreplace+1;
        if nreplace == 1e4
            error('The threshold of the distance is not proper, pls reduce the threshold')
        end
            
            
    end
end
 
        
 
end
    
 
