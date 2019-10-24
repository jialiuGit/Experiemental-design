function design= design_mcm_hyperpara_new(design,likhood,model,nMC,mfpar,utility)

%we evaluate parameter tDomain.tdelta=max(m(t))-min(m(t)), controlling by the meanfunction parameters a and b, which
%should be much less than var(f)=3,a fixed value; Otherwise the effect from the concavity of the mean function will be vanish
%in the posterior prediction.,


m = @(t) mfpar.a-mfpar.c*(t(:)-mfpar.b).^2;


design.fh.utility_KL=@utility_KL;
design.fh.utility_lambda=@utility_lambda;

time= design.x(:,3);


% Create the covariance functions
gpcf1 = gpcf_matern32('lengthScale', 1, 'magnSigma2', 2);
gpcf2 = gpcf_sexp('lengthScale', 0.85, 'magnSigma2', 1); %l_t and sigma2_t initial values

% --- informative priors  ----------
pl_t = prior_gaussian('mu', 0.85, 's2', 0.05);    % temporal length-scale
ps_t = prior_gamma('sh',20, 'is', 20);              % variance
pl_s = prior_gaussian('mu', 0.85, 's2', 0.05);    % spatial length-scale
gpcf2 = gpcf_sexp(gpcf2,'magnSigma2_prior',ps_t,'SelectedVariables',3 , 'lengthScale_prior', pl_t);
gpcf1 = gpcf_matern32(gpcf1, 'lengthScale_prior', pl_s, 'magnSigma2_prior', prior_fixed, 'SelectedVariables',[1 2]);

switch likhood
    case 'Poisson'
        lik = lik_poisson();
        
    case 'Gaussian'
        lik = lik_gaussian();
end


opt=optimset('TolX',1e-6,'TolFun',1e-6);

% --- MAP estimate with Laplace approximation ---

% Set the approximate inference method to Laplace approximation


Dset=design_fix(1000,'samp_method','lattice','dim',3,'grid_reg',linspace(0,1,10),'gridZreg',linspace(0,1,10)); %A test grid cell with 1000 grids

xpred =Dset.x; %evaluate on a [0,1] cube.

switch likhood
    case 'Poisson'

        % Generate samples from the prior predictive distribution
        meanfun = m(time);
        meanfun2=m(xpred(:,3));
        
        switch model
            case 'h+g'
                gp = gp_set('lik', lik, 'cf', {gpcf1 gpcf2}, 'jitterSigma2', 1e-4,...
                    'infer_params', 'covariance');
            case 'h*g'
                gpcf=gpcf_prod('cf', {gpcf1,gpcf2});
                gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-4,...
                    'infer_params', 'covariance');
        end
       
        lt = pl_t.mu + sqrt(pl_t.s2)*randn(nMC,1);  % sample the temporal length-scale
        ls = pl_s.mu + sqrt(pl_s.s2)*randn(nMC,1);  % sample the spatial length-scale
        st = gamrnd(ps_t.sh,1./ps_t.is,nMC,1);  % sample the temporal variance
        for s1 = 1:nMC
            gps1 = gp;
            switch model
                case 'h+g'
                    gps1.cf{1}.lengthScale = ls(s1);
                    gps1.cf{2}.magnSigma2 = st(s1);
                    gps1.cf{2}.lengthScale = lt(s1);
                case 'h*g'
                    gps1.cf{1}.cf{1}.lengthScale = ls(s1);
                    gps1.cf{1}.cf{2}.magnSigma2 = st(s1);
                    gps1.cf{1}.cf{2}.lengthScale = lt(s1);
            end
            Kpred(:,:,s1)= gp_trcov(gps1,Dset.x);  % samples from the prior predictive covariance matrix
            K= gp_trcov(gps1,design.x); %calculate covariance matrix of the coordinates.
            f(:,s1) = meanfun + chol(K,'lower')*randn(numel(time),1); %more safe by using chol!
        end
        sampyd = poissrnd(exp(f));% generate future data from the prior predictive distribution!
       
        switch utility
            
            case 'aEPV'
                
                
                tPV=zeros(size(sampyd,2),1);
                tPE=tPV;
                for j1= 1: size(sampyd,2)
                 
                    gp=gp_optim(gp,design.x,sampyd(:,j1),'z',exp(meanfun),'opt',opt);  %optimize
                     
                    % Set the approximate inference method
                    
                    gp = gp_set(gp, 'latent_method', 'MCMC', 'jitterSigma2', 1e-6);
                   
                    gp_rec=gp_mc(gp, design.x, sampyd(:,j1), 'nsamples', 200, 'display', 20);
                    % Remove burn-in and thin
                    gp_rec=thin(gp_rec,21);
                    
                   
                   [Ef, Varf] = gp_pred(gp_rec, design.x, sampyd(:,j1), xpred, 'z', exp(meanfun),'zt', exp(meanfun2)); %double check here 03.10.17.
                    
                    tPV(j1) = 1/numel(~isnan(Varf))*sum(Varf(~isnan(Varf)));
                    
                    tPE(j1) = 1/numel(~isnan(Varf))*sum(Ef(~isnan(Varf)));
                    
                end
               
                design.utility =mean(tPV);
                
                design.utility_var=var(tPV);
                
                lambda_E=design.fh.utility_lambda(tPE,tPV);
                
                design.utility_lambda = lambda_E;
                design.utility_var_lamda = var(lambda_E);
            
            case 'EKL_D'                
              
                mEKL=zeros(size(sampyd,2),1);               
                gpMAP = gp;
                for j1= 1: size(sampyd,2)                    
                    
                    
                    gpMAP=gp_optim(gpMAP,design.x,sampyd(:,j1),'z',exp(meanfun),'opt',opt);  %optimize
                    gp = gpMAP;
                   
                    gp = gp_set(gp, 'latent_method', 'MCMC', 'jitterSigma2', 1e-6);
                   
                    [gp_rec,g]=gp_mc(gp, design.x, sampyd(:,j1), 'nsamples', 200, 'display', 20);
                    % Remove burn-in and thin
                    gp_rec=thin(gp_rec,21);
                  
                    w = gp_pak(gp_rec);
                  
                    LOGlikXprior = - gpla_e(median(w), gpMAP, design.x, sampyd(:,j1), 'z', exp(meanfun) ) ;
                    EKL = mean(-(gp_rec.etr-gp_rec.e)) - ( 0.5*log(det(cov(w))) + (size(w,2)/2)*log(2*pi) + LOGlikXprior );
                    mEKL(j1) = EKL;
                   
                end
                
                design.mEKL = mEKL;
                design.utility=1 /sum(~isnan(mEKL))*sum(mEKL(~isnan(mEKL)));
                
                design.utility_var=var(mEKL(~isnan(mEKL)));
                
        end
        
        
        
end
end


function utility=utility_lambda(Ef,Varf)
utility=exp(2*Ef + Varf).*(exp(Varf) -1);
end


