function x=minDis(n,dim,delta,dist,weights)
%Copyright (c) 2019-   Jia Liu
 
x=rand(n,dim);
 
if nargin<4
     dist = 'euc';
end

if nargin ==5
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
        i=i+1;  
         
    else
        x(i,:)=rand(1,dim);
        nreplace=nreplace+1;
        if nreplace == 1e5
            error('The threshold of the distance is not proper, pls reduce the threshold')
        end
            
            
    end
end
        
  %figure,
  %plot3(x(:,1),x(:,2),x(:,3),'ro');
  %zlabel('time')
  %title(['min\_dist w/o incPr, design size = ' num2str(n) 'and delta = ', num2str(delta)]);
 if  isempty(nreplace < n)
      fprintf(['Find ' num2str(nreplace) ' optimal new sampling locations in all the candidates.'] )
 else 
      fprintf(['Find ' num2str(n) ' optimal new sampling locations in the candidates.'] )
 end
 

end
