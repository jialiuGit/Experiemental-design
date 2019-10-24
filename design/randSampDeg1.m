function Min_dis_design=randSampDeg1(Dset,n,tDomain,delta,dist,weights)
%Copyright (c) 2019-   Jia Liu

%weights----is a n*1 vector

%grid_reg---A Degsin set with pre-defined resolution

%dist--methods of distance, see pdist, pdist2

%delta--threthold distance between two points in sampling

%n--size of samples

%a distance that weights each coordinate contribution differently
%only work with Euclidean distance.
 
Min_dis_design.fh.reject_design=@reject_design;
Min_dis_design.fh.randshort=@randSamp_short;

x=datasample(Dset,n,1,'Replace',false); 
Min_dis_design.xnoincPr=x;
xbar = setdiff(Dset,x,'rows');

[x,n111]=reject_design(tDomain,x);
    

 
n11=n-n111;
ik=0;
while ik<n11  
    
    x_one=datasample(xbar,1,1,'Replace',false);
    xbar = setdiff(xbar,x_one,'rows');
    
    
    if isempty(xbar)
        message('The study region is too small to find suitable starting points, change the initial design points')
        x=datasample(Dset,n,1,'Replace',false); 
        Min_dis_design.xnoincPr=x;
        [x,n111]=reject_design(tDomain,x);
        xbar = setdiff(Dset,x,'rows');
        n11=n-n111;
    else
        [x_one,n111]=reject_design(tDomain,x_one);
        if n111
            ik=ik+1;
            x=[x;x_one];     
        
        end
    end
end

x00=x;
Min_dis_design.xnoincPr=Min_dis_design.fh.randshort(Dset,n,delta,Min_dis_design.xnoincPr);
xbar0=xbar;



if nargin<5
     dist = 'euc';
end

if nargin ==6
    weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
end

i=1;
nreplace=0;


while i<=n
    
     
    x2=x;
    x2(i,:)=[];
    if nargin ==5
        distS = pdist2(x(i,:),x2,@(Xi,Xj) weuc(Xi,Xj,weights));
    else
        distS = pdist2(x(i,:),x2,dist);
    end
    
    
    if min(distS) >= delta
        xbar=setdiff(xbar0,x(i,:),'rows'); 
        i=i+1;  
      
    else
       
        n22=0;       
        while n22<1 && ~isempty(xbar)
            x_one=datasample(xbar,1,1,'Replace',false);
            [xtmp,n22]=reject_design(tDomain,x_one);     
            xbar = setdiff(xbar,x_one,'rows');

            
        end
        if isempty(xbar)
            delta=delta*0.95;          
            x=x00;            
            i=1;
            xbar = setdiff(Dset,x,'rows');
         
         elseif n22==1
            x(i,:)=xtmp;        
            nreplace=nreplace+1;
            
        end
            
        
    end
    
    
    
end

Min_dis_design.delta=delta;
Min_dis_design.x=x;

    
    
%   figure,
%   plot3(x(:,1),x(:,2),x(:,3),'ro');
%   zlabel('time')
%   title(['min\_dist with incPr, design size = ' num2str(n) 'and delta = ', num2str(delta)]);

 if  isempty(nreplace < n)
      fprintf(['Find ' num2str(nreplace) ' optimal new sampling locations in all the candidates.'] )
 else 
      fprintf(['Find ' num2str(n) ' optimal new sampling locations in the candidates.'] )
 end
  
end




function [x,n]=reject_design(tDomain,design)

%x is a design candidate point generated from [0,1]

% note that with mean function the time domain of t must be same as the study region on time


 
 
x=[];
for ii1=1:size(design,1)
    
    [~,ind_t]=min(pdist2(design(ii1,3), tDomain.t'));
       
    u=rand(numel(ind_t),1);
    ind_t=ind_t( tDomain.aptRate(ind_t) >=u);
        
    if ~isempty(ind_t) 
         
        x=[x;design(ii1,:)];
    else
        x=x;
    end
    
end
n=size(x,1) ;

  
    
end
 

function x=randSamp_short(Dset,n,delta,x)

 
xbar = setdiff(Dset,x,'rows');
 
xbar0=xbar;
i=1;
nreplace=0;
while i<=n
    x2=x;
    x2(i,:)=[];
    distS = pdist2(x(i,:),x2,'euc');
        
    
    if min(distS) >= delta
        i=i+1;  
        xbar=setdiff(xbar0,x,'rows');  %correction on 24.10.17!
    else
        x(i,:)=datasample(xbar,1,1,'Replace',false);
        xbar = setdiff(xbar,x(i,:),'rows'); 
        nreplace=nreplace+1;
       if isempty(xbar) 
            x=datasample(Dset,n,1,'Replace',false);
            i=1;
            xbar = setdiff(Dset,x,'rows');
            xbar0=xbar;
       end
     
        
    end
end
      
   
end
    
%   figure,
%   plot3(x(:,1),x(:,2),x(:,3),'ro');
%   zlabel('time')
%   title(['min\_dist w/o incPr, design size = ' num2str(n) 'and delta = ', num2str(delta)]);
%  if  isempty(nreplace < n)
%       fprintf(['Find ' num2str(nreplace) ' optimal new sampling locations in all the candidates.'] )
%  else 
%       fprintf(['Find ' num2str(n) ' optimal new sampling locations in the candidates.'] )
%  end

%end
