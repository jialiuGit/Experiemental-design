
function design= design_evaluate(design,len_s,likhood,model,nMC,mfpar,utility,len_t,s2_t)

%we evaluate parameter tDomain.tdelta=max(m(t))-min(m(t)), controlling by the meanfunction parameters a and b, which
%should be much less than var(f)=3,a fixed value; Otherwise the effect from the concavity of the mean function will be vanish
%in the posterior prediction.,


%a convex fixed mean function

%load tDomain
%n--size of the design
%a,b---parameters of the mean function.
%len_s---length scale parameter
%likhood--likelihood
%model--model chosen to be analyzed


%demo
% n=5;
% a=0.2; b=0.5;  %depdending on the region of the design space! rescaling 100
%
% %set a convex mean function which has a linear trend of the frequency
% %fun, f(t), considering uneven distributed status of the time in the design space
%
% len_s=1;
%
% model='h+g';
%
% %likhood = 'Poisson';
%
%
% likhood='Gaussian';
%
%
% N=2000;
 
m = @(t) mfpar.a-mfpar.c*(t(:)-mfpar.b).^2;


design.fh.utility_KL=@utility_KL;
design.fh.utility_lambda=@utility_lambda;
 
time= design.x(:,3);

% Create the covariance functions
pl = prior_t('s2',10);
pm = prior_sqrtunif();
gpcf2 = gpcf_sexp('lengthScale',len_t,'magnSigma2',s2_t,'magnSigma2_prior',prior_t,'SelectedVariables',3);
 


switch likhood
    case 'Poisson'
        lik = lik_poisson();
        
    case 'Gaussian'
        lik = lik_gaussian();
end
Dset=design_fix(1000,'samp_method','lattice','dim',3,'grid_reg',linspace(0,1,10),'gridZreg',linspace(0,1,10)); %A test grid cell with 1000 grids
xpred =Dset.x; %evaluate on a [0,1] cube.

switch likhood
    case 'Poisson'
        
        
        switch utility
            
            case 'aEPV'
                
                % Generate samples from the prior predictive distribution
                meanfun = m(time);
                meanfun2=m(xpred(:,3));
                grid_PV=zeros(numel(len_s),1);
                var_grid_PV=grid_PV;
                
                grid_PE=grid_PV;
                var_grid_lambda=var_grid_PV;
                for i1 = 1:numel(len_s)
                    gpcf1 = gpcf_matern32('lengthScale', len_s(i1), 'magnSigma2', 2,'lengthScale_prior', pl, 'magnSigma2_prior', pm,'SelectedVariables',[1 2]);
                    
                    switch model
                        case 'h+g'
                            gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                        case 'h*g'
                            gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                            gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                    end
                    K= gp_trcov(gp,design.x); %calculate covariance matrix of the coordinates.
                    f = repmat(meanfun,1,nMC) + chol(K,'lower')*randn(numel(time),nMC); %more safe by using chol!
                    sampyd = poissrnd(exp(f));% generate future data from the prior predictive distribution!
                    
                    tPV=zeros(size(sampyd,2),1);
                    tPE=tPV;
                    for j1= 1: size(sampyd,2)
                        
                        [Ef, Varf] = gp_pred(gp, design.x, sampyd(:,j1), xpred, 'z', exp(meanfun),'zt', exp(meanfun2)); %double check here 03.10.17.
                       
                        tPV(j1) = 1/numel(~isnan(Varf))*sum(Varf(~isnan(Varf)));
                        
                        tPE(j1) = 1/numel(~isnan(Varf))*sum(Ef(~isnan(Varf)));
                        
                    end
                    grid_PV(i1)=mean(tPV);
                    var_grid_PV(i1)= var(tPV);
                    
                    grid_PE(i1)=mean(tPE);
                    var_grid_lambda(i1)= var(design.fh.utility_lambda(tPE,tPV));
                end
                design.utility =mean(grid_PV);
                design.u_lens=grid_PV;
                design.utility_var=var_grid_PV;
                
                
                design.u_lens_lambda=design.fh.utility_lambda(grid_PE,grid_PV);
                design.utility_var_lamda=var_grid_lambda;
                design.utility_lambda =mean(design.fh.utility_lambda(grid_PE,grid_PV));
                
            case 'EKL'      
                
                
                meanfun = m(time);           
                
                meanfun2=m(xpred(:,3));
                grid_EKL=zeros(numel(len_s),1);
                var_grid_EKL=grid_EKL;
                for i1 = 1:numel(len_s)
                    gpcf1 = gpcf_matern32('lengthScale', len_s(i1), 'magnSigma2', 2,'lengthScale_prior', pl, 'magnSigma2_prior', pm,'SelectedVariables',[1 2]);
                    
                    switch model
                        case 'h+g'
                            gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                        case 'h*g'
                            gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                            gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                    end
                    
                    Kpred= gp_trcov(gp,Dset.x);
                    K= gp_trcov(gp,design.x); %calculate covariance matrix of the coordinates.
                    f = repmat(meanfun,1,nMC) + chol(K,'lower')*randn(numel(time),nMC); %more safe by using chol!
                    sampyd = poissrnd(exp(f));
                    
                    mEKL=zeros(size(sampyd,2),1);                    
                    
                    for j1= 1: size(sampyd,2)
                        [mu2, K2] = gpla_jpred(gp, design.x, sampyd(:,j1), xpred,'z', exp(meanfun), 'zt', exp(meanfun2));                          
                        EKL=utility_KL(Kpred,K2,meanfun2,mu2);                       
                        mEKL(j1) = EKL;
                       
                    end
                    
                    grid_EKL(i1)=1/sum(~isnan(mEKL))*sum(mEKL(~isnan(mEKL)));
                    var_grid_EKL(i1)=var(mEKL(~isnan(mEKL)));
                    
                end
                
                design.utility =mean(grid_EKL);
                design.u_lens=grid_EKL;
                design.utility_var=var_grid_EKL;
                
                
                case 'EKL_D'
                
                
                
                    meanfun = m(time);  


                    grid_EKL=zeros(numel(len_s),1);
                    var_grid_EKL=grid_EKL;
                    for i1 = 1:numel(len_s)
                        gpcf1 = gpcf_matern32('lengthScale', len_s(i1), 'magnSigma2', 2,'lengthScale_prior', pl, 'magnSigma2_prior', pm,'SelectedVariables',[1 2]);

                        switch model
                            case 'h+g'
                                gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                                    'infer_params', 'covariance');
                            case 'h*g'
                                gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                                gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                                    'infer_params', 'covariance');
                        end

                    %Kpred= gp_trcov(gp,Dset.x);
                    K= gp_trcov(gp,design.x); %calculate covariance matrix of the coordinates.
                    f = repmat(meanfun,1,nMC) + chol(K,'lower')*randn(numel(time),nMC); %more safe by using chol!
                    sampyd = poissrnd(exp(f));

                    mEKL=zeros(size(sampyd,2),1);


                    for j1= 1: size(sampyd,2)

                        [mu2, K2] = gpla_jpred(gp, design.x, sampyd(:,j1), design.x,'z', exp(meanfun));  
                         
                        
                        EKL=utility_KL(K,K2,meanfun,mu2+meanfun);                       
                        mEKL(j1) = EKL;                        
                    end
                     
                    grid_EKL(i1)=1/sum(~isnan(mEKL))*sum(mEKL(~isnan(mEKL)));                    
                    var_grid_EKL(i1)=var(mEKL(~isnan(mEKL)));

                    end                     
                    design.utility =mean(grid_EKL);                     
                    design.u_lens=grid_EKL;
                    design.utility_var=var_grid_EKL;
                
        end
        
        
        
    case 'Gaussian'
        switch utility
            case 'aEPV'
                
                grid_PV=zeros(numel(len_s),1);
                var_grid_PV=grid_PV;
                for i1 = 1:numel(len_s)
                    gpcf1 = gpcf_matern32('lengthScale', len_s(i1), 'magnSigma2', 2,'lengthScale_prior', pl, 'magnSigma2_prior', pm,'SelectedVariables',[1 2]);
                    
                    switch model
                        case 'h+g'
                            gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                        case 'h*g'
                            gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                            gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                    end
                   
                    K= gp_trcov(gp,design.x); %calculate covariance matrix of the coordinates.
                    
                    f = chol(K,'lower')*randn(numel(time),nMC); %correction on 6.10.17, generate f'
                    sampyd=f+ 0.01*randn(size(design.x,1),nMC); % with Gaussian likelihood, the inference is doing with yt = y-m (generate data y')
                   
                    PV=zeros(size(sampyd,2),1);
                    for j1= 1: size(sampyd,2)                       
                        [~, PV1] = gp_pred(gp, design.x,sampyd(:,j1), xpred); %correction on 6.10.17
                        PV(j1)=1/numel(~isnan(PV1))*sum(PV1(~isnan(PV1)));
                    end
                    grid_PV(i1)=mean(PV);
                    var_grid_PV(i1)= var(PV);
                end
                design.utility =mean(grid_PV);
                design.u_lens=grid_PV;
                design.utility_var=var_grid_PV;
                
                
                
                
            case 'EKL'  % the general case
                
                grid_EKL=zeros(numel(len_s),1);
                var_grid_EKL=grid_EKL;
                for i1 = 1:numel(len_s)
                    gpcf1 = gpcf_matern32('lengthScale', len_s(i1), 'magnSigma2', 2,'lengthScale_prior', pl, 'magnSigma2_prior', pm,'SelectedVariables',[1 2]);
                    
                    switch model
                        case 'h+g'
                            gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                        case 'h*g'
                            gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                            gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                                'infer_params', 'covariance');
                    end
                    
                    Kpred= gp_trcov(gp,Dset.x);
                    K= gp_trcov(gp,design.x); %calculate covariance matrix of the coordinates.
                    sampyd = exp(0.5*diag(K));
                    mEKL=zeros(size(sampyd,2),1);
                    for j1= 1: size(sampyd,2)                       
                        [mu2, K2] = gp_jpred(gp, design.x, sampyd(:,j1), xpred); %correction on 6.10.17                     
                        EKL=utility_KL(Kpred,K2,0,mu2);
                        
                        mEKL(j1) = EKL;                       
                    end
                    
                    grid_EKL(i1)=1/sum(~isnan(mEKL))*sum(mEKL(~isnan(mEKL)));
                    var_grid_EKL(i1)=var(mEKL(~isnan(mEKL)));
                    
                end
                                
                design.utility =mean(grid_EKL);
                design.u_lens=grid_EKL;
                design.utility_var=var_grid_EKL;        
                
                
                
                 case 'EKL_D'
                 
                    grid_EKL=zeros(numel(len_s),1);
                    var_grid_EKL=grid_EKL;
                    for i1 = 1:numel(len_s)
                        gpcf1 = gpcf_matern32('lengthScale', len_s(i1), 'magnSigma2', 2,'lengthScale_prior', pl, 'magnSigma2_prior', pm,'SelectedVariables',[1 2]);

                        switch model
                            case 'h+g'
                                gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                                    'infer_params', 'covariance');
                            case 'h*g'
                                gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                                gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                                    'infer_params', 'covariance');
                        end

                        K= gp_trcov(gp,design.x); %calculate covariance matrix of the coordinates.
                         
                        sampyd = exp(0.5*diag(K));
                        mEKL=zeros(size(sampyd,2),1);

                        
                        for j1= 1: size(sampyd,2)                            
                            [mu2, K2] = gp_jpred(gp, design.x, sampyd(:,j1), design.x);  
                            EKL=utility_KL(K,K2,0,mu2);                          
                            mEKL(j1) = EKL;
                             
                        end
                        
                        grid_EKL(i1)=1/sum(~isnan(mEKL))*sum(mEKL(~isnan(mEKL)));                         
                        var_grid_EKL(i1)=var(mEKL(~isnan(mEKL)));

                    end
                    design.utility =mean(grid_EKL);                    
                    design.u_lens=grid_EKL;
                    design.utility_var=var_grid_EKL;
        end
        
        
        
end

end



function utility=utility_KL(A,B,muA,muB)


try
    cholB=chol(B, 'lower');
catch
    cholB=[];
end


if ~isempty(cholB)
    MU_dif=(muA-muB);
    factor= sum(log(diag(chol(A,'lower'))./diag(cholB)));
    utility=factor+ 0.5*( -size(B,1)+trace(A\B) + MU_dif'/A*MU_dif);    
else
     utility = nan;
end



end

function utility=utility_lambda(Ef,Varf)

utility=exp(2*Ef + Varf).*(exp(Varf) -1);


end

