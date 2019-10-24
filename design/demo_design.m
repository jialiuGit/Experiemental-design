
n=10;
delta=0.46; deltak=0.35;
% 
%  
% n=30;% 
% delta=0.25;deltak=0.2;

tDomain=priori_time(n);


z1=design_fix_3D_incPr(n,tDomain,'samp_method','Halton');
hold on
z2=design_fix(n,'samp_method','Halton','dim',3);

legend('with incPr','without incPr','Location','NorthEast')
hold off


zz1=design_fix_3D_incPr(n,tDomain,'samp_method','Sobol');
hold on
zz2=design_fix(n,'samp_method','Sobol','dim',3);
legend('with incPr','without incPr','Location','NorthEast')
hold off

r1=design_fix_3D_incPr(n,tDomain,'samp_method','Random');
hold on
r2=design_fix(n,'samp_method','Random','dim',3);
legend('with incPr','without incPr','Location','NorthEast')
hold off





rr1=design_fix_3D_incPr(n,tDomain,'D_N',100,'samp_method','Random','name','min_dist','delta',delta);
hold on
rr2=design_fix(100,'n',n,'samp_method','Random','name','min_dist','dim',3,'delta',delta);
legend('with incPr','without incPr','Location','NorthEast')
hold off




q1=design_fix_3D_incPr(n,tDomain,'D_N',100,'samp_method','Random','name','space_fill');
hold on
qq2=design_fix(100,'n',n,'samp_method','Random','name','space_fill','dim',3);
legend('with incPr','without incPr','Location','NorthEast')
hold off 



g1=design_fix_3D_incPr(n,tDomain,'D_N',100,'samp_method','Random','name','close_pair','delta',deltak);
hold on
gg2=design_fix(100,'n',n,'samp_method','Random','name','close_pair','dim',3,'delta',deltak);
legend('with incPr','without incPr','Location','NorthEast')
hold off

%design without rejection

n=10;
delta=0.4; deltak=0.2;


figure
subplot(2,3,1) 

z2=design_fix(n,'samp_method','Halton','dim',2);
subplot(2,3,2) 
zz2=design_fix(n,'samp_method','Sobol','dim',2);


subplot(2,3,3)
r2=design_fix(n,'samp_method','Random','dim',2);
subplot(2,3,4)
rr2=design_fix(100,'n',n,'samp_method','Random','name','min_dist','dim',2,'delta',delta);
subplot(2,3,5)
qq2=design_fix(100,'n',n,'samp_method','Random','name','space_fill','dim',2);
subplot(2,3,6)
gg2=design_fix(100,'n',n,'samp_method','Random','name','close_pair','dim',2,'delta',deltak);

 


n=30;
delta=0.25; deltak=0.1;


figure
subplot(2,3,1) 

z2=design_fix(n,'samp_method','Halton','dim',2);
subplot(2,3,2) 
zz2=design_fix(n,'samp_method','Sobol','dim',2);


subplot(2,3,3)
r2=design_fix(n,'samp_method','Random','dim',2);
subplot(2,3,4)
rr2=design_fix(100,'n',n,'samp_method','Random','name','min_dist','dim',2,'delta',delta);
subplot(2,3,5)
qq2=design_fix(100,'n',n,'samp_method','Random','name','space_fill','dim',2);
subplot(2,3,6)
gg2=design_fix(100,'n',n,'samp_method','Random','name','close_pair','dim',2,'delta',deltak);




 
 




%same designs include the rejection algorithm.
rr1=design_fix_3D_incPr(n,tDomain,'D_N',27000,'samp_method','Halton','name','min_dist','delta',delta);  
rr2=rr1;
rr2.x=rr1.x_orig;


f1=design_fix_3D_incPr(n,tDomain,'name','min_dist_rand','delta',delta);
f2=f1;
f2.x=f1.x_orig;

q1=design_fix_3D_incPr(n,tDomain,'D_N',64000,'samp_method','Halton','name','space_fill');

qq2=design_fix(64000,'n',n,'samp_method','Halton','name','space_fill','dim',3);


g1=design_fix_3D_incPr(n,tDomain,'name','close_pair','delta',delta);

gg2=g1;
gg2.x=g1.x_orig;


k1=design_fix_3D_incPr(n,tDomain,'D_N',27000,'samp_method','Fibonacci lattice');
kk1=k1;
kk1.x=k1.x_orig;

 
 
