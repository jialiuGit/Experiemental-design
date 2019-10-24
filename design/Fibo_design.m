

function x = Fibo_design(n) 
%Copyright (c) 2019-   Jia Liu
% The Fibonacci lattice design 2D

knumus1= @k_nu;

k=knumus1(n);

Dset = zeros(n,2);

for i = 1:n    
     Dset(i,1) = mod((i+0.5 )/ n, 1);    
     Dset(i,2) = mod((k*(i-1)+0.5 )/ n, 1);
    
    
end

x=Dset;

if n>1
    plot(x(:,1),x(:,2),'.');
else
    plot(x(1),x(2),'.');
end

end


function k_numus1 = k_nu(n) 

k_nv2=1;
 
k_nv= 2;
k_nv1 = k_nv;

while k_nv < n   
                  
           k_nv= k_nv1+k_nv2;
         
           k_nv2=k_nv1;
           k_nv1=k_nv;
end

k_numus1 = k_nv2;
end