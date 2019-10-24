function x = Fibo_design_3D(n,index) 
%Copyright (c) 2019-   Jia Liu
 
delta=0.5;
tau1=sqrt(5)-1;
tau2=sqrt(2);

switch index
    case 'vec'    
        ind=linspace(1,n,n);

        cmplx=@(n)n(:).^delta.*exp(1j*2*pi*tau1*n(:));

        tmp_vec= 2*pi*tau2*ind';

        cmp_tem= cmplx(ind);

        re_osa=0.5+real(cmp_tem);
        im_osa_cos=0.5+imag(cmp_tem).*cos(tmp_vec);
        im_osa_sin=0.5+imag(cmp_tem).*sin(tmp_vec);

        x=zeros(n,3);
        x(:,1)=re_osa; 
        x(:,2)=im_osa_cos;
        x(:,3)=im_osa_sin;
    
    case 'scalar'

        cmplx=@(n)n.^delta.*exp(1j*2*pi*tau1*n);

        tmp_vec= 2*pi*tau2*n; %here n is a index augment

        cmp_tem= cmplx(n);

        x(1)=0.5+real(cmp_tem);
        x(2)=0.5+imag(cmp_tem)*cos(tmp_vec);
        x(3)=0.5+imag(cmp_tem)*sin(tmp_vec);

end




 
    
    



 

    
 

 

% if n>1
%     plot3(x(:,1),x(:,2),x(:,3),'.');
% else
%     plot(x(1),x(2),x(3),'.');
% end

end

 