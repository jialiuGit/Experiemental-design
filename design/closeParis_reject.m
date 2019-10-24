function x_design= closeParis_reject(dim,n, k, delta,tDomain,dist)

%Copyright (c) 2017-   Jia Liu    

rand.fh.circ = @rand_circ;

rand.fh.ball=@randsphere;

if k>0.5*n
     error('k must be less or equal to n/2')
end

deltak=delta*sqrt(n/(n-k));
x_design.delta=deltak;
if nargin<6
    x=minDis_reject(n-k,dim,tDomain,deltak);
else
    x=minDis_reject(n-k,dim,tDomain,deltak,dist);
end
   
 
x2= datasample(x.x,k,1,'Replace',false); 

 
 
r=0.5*deltak;
    for j1=1:k
        
        if dim==2
            n11=0;
             
            while n11<1
                x_tmp=rand.fh.circ(1,x2(j1,:),r); 
                 
                if (any(x_tmp<0) || any(x_tmp>1))
                    continue
                end                 
                [x_tmp,n11]=reject_design(tDomain,x_tmp);  
                
                if isfield(x,'org')
                    x.org=[x.org;x_tmp.org];
                    
                    strX=x_tmp.x;
                     
                    clear x_tmp
                    x_tmp=strX; 
                    clear strX
                end 
            end
            x.x(n-k+j1,:)=x_tmp;
            
         else
            n11=0;
            while n11<1                
                x_tmp=rand.fh.ball(1,x2(j1,:),r); 
                if (any(x_tmp<0) || any(x_tmp>1))
                    continue
                end%                 
                [x_tmp,n11]=reject_design(tDomain,x_tmp);   
                
                 if isfield(x,'org')
                    x.org=[x.org;x_tmp.org];
                    
                    strX=x_tmp.x;
                     
                    clear x_tmp
                    x_tmp=strX; 
                    clear strX
                end 
                
                
                
            end
            x.x(n-k+j1,:)=x_tmp;          
              
        end
    end
    
    x_design.x=x.x;
    
    
    if isfield(x,'org')
        
        x_design.x_orig= x.org(1:n);
        
    else
        
        
    
    
        x2_orig= datasample(x.xnoincPr,k,1,'Replace',false);  


        for i1=1:k

            if dim==2
                n1=0;
                while n1<1
                    x_tmp=rand.fh.circ(1,x2_orig(i1,:),r);
                     if (any(x_tmp<0) || any(x_tmp>1))
                        continue
                     else
                         n1=1;
                     end
                end

                x.xnoincPr(n-k+i1,:)=x_tmp;

            else
                n1=0;
                while n1<1 
                    x_tmp=rand.fh.ball(1,x2_orig(i1,:),r);
                    if (any(x_tmp<0) || any(x_tmp>1))
                        continue
                     else
                         n1=1;
                     end
                end

                x.xnoincPr(n-k+i1,:)=x_tmp;

            end
        end
        x_design.x_orig= x.xnoincPr;
        
    end
     
%     figure,
%     plot3(x_design.x(:,1),x_design.x(:,2),x_design.x(:,3),'ro');
%    
%     zlabel('time')
%     title(['close pair with incPr, design size = ' num2str(n) ' and delta = ', num2str(delta)]);

end

function X = rand_circ(N,x,r)
% Generates N random points in a circle.
% RAND_CIRC(N) generates N random points in the unit circle at (0,0).
% RAND_CIRC(N,x,r) generates N random points in a circle with radius r 
% and center at point x.

    if nargin<2
       x(1) = 0;
       x(2) = 0;
       r = 1;
    end    
    
    theta = rand(1,N)*(2*pi);
    r = r*sqrt(rand(1,N));
    xc = x(1) + r.*cos(theta);
    yc = x(2) + r.*sin(theta);

    X=[xc' yc'];
end

function X = randsphere(N,x,r)
 

theta = rand(1,N)*(2*pi);
phi = pi*rand(1,N);

xc=x(1)+r.*sin(phi).*cos(theta);
yc=x(2)+r.*sin(phi).*sin(theta);
zc= x(3)+r.*cos(phi);
 
X=[xc' yc' zc'];
 
end