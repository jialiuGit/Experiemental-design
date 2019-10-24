function x= closeParisSamp(Dset,n, k, delta,dist) %Dset is dim


%Copyright (c) 2019-   Jia Liu
 

rand.fh.circ = @rand_circ;

rand.fh.ball=@randsphere;

if k>0.5*n
     error('k must be less or equal to n/2')
end

deltak=delta*sqrt(n/(n-k));
 
x=minDis(n,Dset,deltak,dist);   

 
x2= datasample(x,k,1,'Replace',false);  

 
r=0.5*deltak;
    for j1=1:k
        
        if Dset==2
            n1=0;
            while n1<1
                x_tmp=rand.fh.circ(1,x2(j1,:),r);
                if (any(x_tmp<0) || any(x_tmp>1))
                        continue
                else
                    n1=1;
                end
            end
            x(n-k+j1,:)=x_tmp;
        else
            
            n1=0;
            while n1<1 
                x_tmp=rand.fh.ball(1,x2(j1,:),r);
                if (any(x_tmp<0) || any(x_tmp>1))
                    continue
                 else
                     n1=1;
                 end
            end
            x(n-k+j1,:)=x_tmp;
        end
    end
    
    %x=[x;xk];
%   figure,
%   plot3(x(:,1),x(:,2),x(:,3),'ro');
%    
%   zlabel('time')
%   title(['close pair w/o incPr, design size = ' num2str(n) ' and delta = ', num2str(delta)]);


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

    X=[xc yc];
end

function X = randsphere(N,x,r)
 

theta = rand(1,N)*(2*pi);
r = r*sqrt(rand(1,N));
phi = pi*rand(1,N);

xc=x(1)+r.*sin(phi)*cos(theta);
yc=x(2)+r.*sin(phi)*sin(theta);
zc= x(3)+r.*cos(phi);


X=[xc' yc' zc'];
 
end