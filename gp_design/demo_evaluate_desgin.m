


Dset=rand(1000,3);

tDomain=priori_time(n);
k=5;
lik='Poisson';

model ='h+g';

nMC=50;
mfpar.a = 2;
mfpar.b = 0.5;
mfpar.c = 30;    


delta=0.4; deltak=0.3;

n=10;

utility='aEPV';
%utility= 'EKL';
%utility= 'EKL_D';  see lemma 1.
%----------when hyperparameter are fixed!--------------
len_t=0.85;
s2_t=1;
len_s= 1;
z1=design_fix_3D_incPr(n,tDomain,'samp_method','Halton');        
design1=design_evaluate(z1,len_s,lik,model,nMC,mfpar,utility,len_t,s2_t);

 
%----------go through hyperparameter space!--------------
%utility= 'EKL_D'; 
design11=design_mcm_hyperpara_new(z1,lik,model,nMC,mfpar,utility);