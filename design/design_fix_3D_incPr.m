function design=design_fix_3D_incPr(n,tDomain,varargin)

%Edit 02.10.2017

%Copyright (c) 2017-   Jia Liu
%              2019-   Jia Liu
%tDomain ---is a structure arrary

%weights----is a n*1 vector

%grid_reg---A Degsin set with pre-defined resolution, the default
%resolution for 2D is 40x40, for 3D is 40x40x40.

%dist--methods of distance, see pdist, pdist2, the default dist is euc.

%delta--threthold distance between two points in sampling, the default
%value for min_dist design is 0.06.

%n--size of design

%nclosepair--number of close pair when calling close pair min_dist design,
%the default value is n/2.

%D_N--number of all the points in the design space;  default value is
%equivent to the resolution 40x40x40 lattice

%dim-- dimensions of the design space

%x-- is the design based minimun distance

%nCandidate -- number of candiate points of design which satisfy basic
%designs but did not consider the uneven probablility in time domain.

%name--the design, the default design is sampling design equvilent to the
%sampling method.

%a distance that weights each coordinate contribution differently
%only work with Euclidean distance.


ip=inputParser;
ip.FunctionName = 'RANDSAMPDEG';
ip.addRequired('n', @(x) ~isempty(x) &&isscalar(x) && isfinite(x));
ip.addRequired('tDomain',@(x) isstruct(x)|| isempty(x));
ip.addOptional('design', [], @isstruct);
ip.addParameter('D_N',27000, @(x) isscalar(x) && isfinite(x)); %
ip.addParameter('delta',0.06, @(x) isscalar(x) && isfinite(x));
ip.addParameter('nCandidate',400, @(x) isscalar(x) && isfinite(x));
ip.addParameter('nclosepair', @(x) isscalar(x) && isfinite(x));
ip.addParameter('dim',3, @(x) isscalar(x) && isfinite(x)); %dim
ip.addParameter('grid_reg',linspace(0,1,30), @(x) isvector(x) && isreal(x));
ip.addParameter('gridZreg','grid_reg', @(x) isvector(x) && isreal(x));
ip.addParameter('weights',[], @(x) isvector(x) && isreal(x));
ip.addParameter('samp_method','lattice', @(x) ischar(x) || iscell(x));
ip.addParameter('name','Random',@(x) ischar(x) || iscell(x));
ip.addParameter('dist','euclidean', @(x) ischar(x) || iscell(x));
ip.addParameter('Pc',5, @(x) isscalar(x) && isfinite(x) ||ismatrix(x||isvector(x) && isreal(x) ));
ip.addParameter('Num',500, @(x) isscalar(x) && isfinite(x))
ip.addParameter('data',[], @(x) ismatrix(x))


ip.parse(n,tDomain,varargin{:});
design=ip.Results.design;

if isempty(tDomain)
    %tDomain=tDomain_set();
    error('Set uneven distributed design dimension!')
end

if isempty(design)
    % Initialize a design
    init=true;
end

% sampling method
if init || ~ismember('samp_method',ip.UsingDefaults)
    design.samp_method=ip.Results.samp_method;
end
%   if ~ismember(design.name,{'min_dist' 'rand_lattice' 'close_pair' 'space_fill'})
%       design.name=design.samp_method;
%   end
% name of the design
if  init ||~ismember('name',ip.UsingDefaults)
    if ~strcmp(design.samp_method,'lattice')&&~ismember(ip.Results.name,{'min_dist'  'min_dist_rand' 'rand_lattice' 'close_pair' 'space_fill'})
        design.name=design.samp_method;
    else
        design.name=ip.Results.name;
    end
end

if ismember(ip.Results.name,{'min_dist_rand' 'close_pair'})
      design.samp_method='sample method';
end

% type of distance (minimun distance)
if  init ||~ismember('dist',ip.UsingDefaults)%&& ...
    %ismember(design.name,{'min_dist' 'space_fill'})
    design.dist=ip.Results.dist;
end
%  number of candidates points in design set
if  init ||~ismember('D_N',ip.UsingDefaults)
    D_N=ip.Results.D_N;
end

%  number of sampling points w/0 consider uneven probability on time domain
if  init ||~ismember('nCandidate',ip.UsingDefaults)
    c_N=ip.Results.nCandidate;
end

% dimension of the design space
if init || ~ismember('dim',ip.UsingDefaults)
    design.dim=ip.Results.dim;
end
% threshold (e.g., minimum distance in min_dist design)
if init || ~ismember('delta',ip.UsingDefaults)
     if ismember(ip.Results.name,{'min_dist' 'close_pair' 'min_dist_rand'})
         design.threshold=ip.Results.delta;
     end
end

%    % time domain
if strcmp(design.samp_method,'lattice')
     if  init ||~ismember('grid_reg',ip.UsingDefaults)
          grid_reg=ip.Results.grid_reg; 
          if isempty(grid_reg)
              grid_reg=linspace(0,1,nthroot(D_N,design.dim));
          end
          design.grid_reg=grid_reg;
      end
      if  init ||~ismember('gridZreg',ip.UsingDefaults)
          gridZreg=ip.Results.gridZreg; 
          if isempty(gridZreg)
              gridZreg=grid_reg;     
          end
          design.gridZreg=gridZreg;
      end
end

if strcmp(design.name,'close_pair')
    if init ||~ismember('nclospair',ip.UsingDefaults)
        design.nclosepair=floor(0.5*ip.Results.n);
    end
end

if ismember(design.dist,{'Minkoski' 'mahalanobis' 'seuclidean'})
    if  init ||~ismember('Pc',ip.UsingDefaults)
        Pc=ip.Results.Pc;
    end
end
%  data
if  init ||~ismember('data',ip.UsingDefaults)
    data=ip.Results.data;
end

% dimension of the design space
if  ismember('weights',ip.UsingDefaults)
    design.weights=false;
else
    design.weights=true;
end

design.size=n;
design.fh.reject_design = @reject_design;

design.fh.rescal_Fibo= @rescal_to_unitcube;
desing.fh.vainwekscale=@vain_data_scale_week;

switch  design.samp_method
    
    case 'Sobol'
        %sobol random sequence in 3D
        
        p = sobolset(design.dim,'Skip',1e3,'Leap',1e2); %d---dimensions of design space
        
       
        if strcmp(design.name,'Sobol')
           
           Dset0 = net(p,n);
           Dset= design.fh.reject_design(tDomain,Dset0);  
           n1=1;
            
           while size(Dset,1)<n  
               tmp = net(p,n+n1);
               Dset0=tmp(n+n1,:);                
               [Dset_x,n11] = design.fh.reject_design(tDomain,Dset0);  
               
               n1=n1+1 ;
                
               if n11==1
                   Dset=[Dset;Dset_x];
               end
           end
        else
            Dset = net(p,D_N);
        end
        
    case 'rand_sobol'
        
         p = sobolset(design.dim+1,'Skip',1e3,'Leap',1e2); %d---dimensions of design space
        
         Dset=[];
         Dset_org=Dset;
         k=1;
         while size(Dset,1)<n%
             u = randi([1 D_N-n], design.dim+1,1);

             Dset0 = net(p,max(u)+k);             

             xtmp=[Dset0(u(1)+k,1) Dset0(u(2)+k,2) Dset0(u(3)+k,3)];              
             
             if size(Dset_org,1)<n
                 Dset_org=[Dset_org;xtmp];
             end
             
            
              %--------Reject design
              
              [xtmp,n11] = design.fh.reject_design(tDomain,xtmp);  
               
               
                
               if n11
                   Dset=[Dset;xtmp];
               end
              
              
         end
         design.x_orig=Dset_org;
     

    case 'Halton'
        %sobol random sequence in 2 D
        p = haltonset(design.dim,'Skip',1e3,'Leap',1e2); %d---dimensions of design space
        if strcmp(design.name,'Halton')
            Dset0 = net(p,n);
            
            Dset = design.fh.reject_design(tDomain,Dset0);  
            n1=1;
            while size(Dset,1)<n  
               tmp = net(p,n+n1);
               Dset0=tmp(n+n1,:);                
               [Dset_x,n11] = design.fh.reject_design(tDomain,Dset0);  
               
               n1=n1+1 ;
                
               if n11==1
                   Dset=[Dset;Dset_x];
               end
           end
        else
            Dset = net(p,D_N);
        end
        
        
    case 'lattice'
        if design.dim == 3
            %             if ~exist('gridZreg') %Improve this command line later
            %                 gridZreg = grid_reg;
            %             end
            [X, Y,Z] = meshgrid(grid_reg,grid_reg,gridZreg);
            p=[X(:) Y(:) Z(:)];
            disp('The default resolution is 30x30x30 3D grids.')
        else
            [X, Y] = meshgrid(grid_reg,grid_reg); %here d is resoluation of the z axis.
            p=[X(:) Y(:)];
            disp('The default resolution (grid_reg) is 30x30 2D grids.')
        end
        
        if strcmp(design.name,'rand_lattice')
            Dset0=datasample(p,n,1,'Replace',false);
            [Dset,n11] = design.fh.reject_design(tDomain,Dset0);  
            
            while n11<n
                 p=setdiff(p,Dset0,'rows');
                 Dset0=datasample(p,n-n11,1,'Replace',false);
                 [Dset_x,n111] = design.fh.reject_design(tDomain,Dset0);   
                 n11=n11+n111;
                 Dset=[Dset;Dset_x];
%                  
                 if isempty(Dset0)
                     error('The study region is too small to find suitable starting points!')
                 end

            end
        else
           Dset=p;
        end
        
    case 'Fibonacci lattice'
     
            Dset0=Fibo_design_3D(n,'vec');
            Dset=design.fh.rescal_Fibo(Dset0,Dset0);
            design.x_orig = Dset;
            [Dset,~,apt_ind] = design.fh.reject_design(tDomain,Dset); 
            if size(Dset,2)>3               
                Dset0=Dset0(apt_ind,:);                  
            else
                Dset0=Dset0(apt_ind,:);
            end
            size(Dset,1)
            if size(Dset,1)>=n
                Dset=Dset(1:n,:);
            else


                i1=1;

                while size(Dset,1)<(n)                  
                    xtmp=Fibo_design_3D(n+i1,'scalar');  

                    Dset0 = [Dset0;xtmp];

                    Dset=design.fh.rescal_Fibo(Dset0,Dset0);%first rescaling then apply reject rule

                    xtmp=Dset(end,:); 

                    [xtmp,n11] = design.fh.reject_design(tDomain,xtmp); 
                    %for new reject design               

                    if ~n11
                        Dset(end,:)=[];
                        Dset0(end,:)=[];
                    elseif size(xtmp,2)>3                          
                        Dset=[Dset;xtmp(1:3)];                        
                    else
                        Dset=[Dset;xtmp];
                    end 

                    i1=i1+1;

                end
            end
            

    case 'Random'   
        if strcmp(design.name,'Random')            
            if ~isempty(data)
                Dset0= datasample(data,n, 'Replace',false);
            else
                Dset0=rand(n,design.dim);
            end
            design.x_orig=Dset0;
            [Dset,n11] = design.fh.reject_design(tDomain,Dset0);  
            n1=0;
            k=n-n11;  
             
            while n1<k
                x_tmp=rand(0.5*n,design.dim);

                x_tmp= x_tmp(~ismember(x_tmp,Dset,'rows'),:);
                if isempty(x_tmp)
                    continue
                else
                    [Dset_x,n11] = design.fh.reject_design(tDomain,x_tmp);
                     if n11>0
                         n1=n1+n11;
                         Dset=[Dset;Dset_x];
                     end
                end
            end
            Dset=Dset(1:n,:);
        else            
            Dset=rand(D_N,design.dim);
        end       
end
 

if strcmp(design.name,'min_dist')
    delta=design.threshold;
    if ~isempty(data)
        Dset=data;
    end    
    min_dis_deg=randSampDeg1(Dset,n,tDomain,delta);
    x=min_dis_deg.x;
    design.x_orig=min_dis_deg.xnoincPr;
    design.threshold=min_dis_deg.delta;
elseif strcmp(design.name,'close_pair')
    delta=design.threshold;
    k=design.nclosepair;
    cPari_deg=closeParis_reject(design.dim,n, k, delta,tDomain);
    x=cPari_deg.x;
    design.x_orig=cPari_deg.x_orig;
    design.deltak=cPari_deg.delta;
    
elseif   strcmp(design.name,'min_dist_rand')
    design.samp_method='random';
    delta=design.threshold;
    min_dis_rand=minDis_reject(n,design.dim,tDomain,delta);
    x=min_dis_rand.x;
    design.x_orig = min_dis_rand.xnoincPr;
elseif strcmp(design.name,'space_fill')    
    x=coffee_house_design1(Dset,n,tDomain);
else
    x=Dset;
end


if ~isempty(x)
    design.x=x;
    figure
    if design.dim==3
        if n>1
            plot3(x(:,1),x(:,2),x(:,3),'ro');
        else
            plot3(x(1),x(2),x(3),'ro');
        end
        xlim([0 1]);ylim([0 1]);zlim([0 1]);
       zlabel('time');
    end
    if strcmp(design.name,'rand_lattice')
        title([design.name ' design with incPr, design size= ' num2str(n) ],'Interpreter','none')
    elseif ~strcmp(design.samp_method,design.name)         
        title([design.name ' design with incPr, design size= ' num2str(n)],'Interpreter','none')
    else
        title([design.name ' design with incPr, design size= ' num2str(n)],'Interpreter','none')
    end
end

end




function x=rescal_to_unitcube(x2,x)

if isequal(x,x2)
    min_x= min(x(:,1));
    min_y= min(x(:,2));
    min_z= min(x(:,3));
    x(:,1)=x(:,1)-min_x;
    x(:,1)=x(:,1)./max(x(:,1));
    x(:,2)=x(:,2)-min_y;
    x(:,2)=x(:,2)./max(x(:,2));
    x(:,3)=x(:,3)-min_z;
    x(:,3)=x(:,3)./max(x(:,3));
else
    min_x= min(x(:,1));
    min_y= min(x(:,2));
    min_z= min(x(:,3));
    x2(:,1)=x2(:,1)-min_x;
    x2(:,1)=x2(:,1)./max(x2(:,1));
    x2(:,2)=x2(:,2)-min_y;
    x2(:,2)=x2(:,2)./max(x2(:,2));
    x2(:,3)=x2(:,3)-min_z;
    x2(:,3)=x2(:,3)./max(x2(:,3)); 
    x=x2;
end

end

function [x,n,apt_ind]=reject_design(tDomain,design,varargin)

%varargin ---Dset if consider spatiotemporal rejection.

%x is a design candidate point generated from [0,1]

% note that with mean function the time domain of t must be same as the study region on time


if ~isfield(tDomain,'Dset')
    sptmp = false;
else
    sptmp = true;  
end

 
x=[];


if sptmp
    
    design=design(:,[2 1 3]);  
    designTMP=design;        
     %first check if the generated design is in the study region.
    if isfield(tDomain,'scale') 
         temp=design;
         %for fast sampling, we only scaling the design from [0,1] to the
         %ROI. 
         week =  tDomain.scale{1};
         scaleY=  tDomain.scale{2};       
                  
         torg=round(temp(:,3)*(week(end)+3-(week(1)-3))+(week(1)-3));
         [~,ind_wekd] =min(pdist2(torg,week),[],2);
         ind1=find(ind_wekd==1 |ind_wekd==2);  
         temp(ind1,:)=[];
         design(ind1,:)=[];
         torg(ind1)=[];
         ind_wekd(ind1)=[];
         
         
         temp(:,3)=week(ind_wekd);              
         temp(:,1)=round(temp(:,1).*scaleY(ind_wekd,1)+ scaleY(ind_wekd,2));
         temp(:,2)=round(temp(:,2)*tDomain.scale{3})+1;
         ind2=find(ismember(temp,tDomain.Dset,'rows')); %first rule out the land area.
         
         %start regject method
         if ~isempty(ind2)
             design=design(ind2,:);
             temp=temp(ind2,:);
             torg=torg(ind2);
            
             apt_ind=[];
             for ii1=1:size(temp,1)
                 u=0.5*rand(1);
                 
                 ind_deg= find(ismember(tDomain.Dset,temp(ii1,:),'rows'));
                 if (tDomain.aptRate(ind_deg) >=u)                     
                    x=[x; temp(ii1,:) design(ii1,:)  ind_deg torg(ii1)]; %the last column is the index for retrieve meanfun of the design
                    apt_ind=[apt_ind;ii1]
                 else
                    x=x;
                 end
             end
             
         else
             x=x;
         end          
          
    else 
          
          for ii1=1:size(design,1)
            
            u=rand(1); 
            ind_deg= ismember(tDomain.Dset,design(ii1,:),'rows');
            
            if (tDomain.aptRate(ind_deg) >=u)                
                apt_ind=[apt_ind;ii1];                 
                x=[x;design(ii1,:) ind_deg]; 
            else
                x=x;
            end
        
          end 
     end

    
    
else
    apt_ind=[];
    for ii1=1:size(design,1)

        [~,ind_t]=min(pdist2(design(ii1,3), tDomain.t'));
        

        u=rand(numel(ind_t),1);
        ind_t=ind_t( tDomain.aptRate(ind_t) >=u);

        if ~isempty(ind_t) 
            apt_ind=[apt_ind;ii1]; 
            x=[x;design(ii1,:)];
        else
            x=x;
        end

    end
end
n=size(x,1)  ;
 
  
 
end


function Dset=vain_data_scale_week(Dset, tDomain)
         temp=Dset;         
         week =  tDomain.scale{1};
         scaleY=  tDomain.scale{2};                 


         torg=round(temp(:,3)*(week(end)+3-(week(1)-3))+(week(1)-3));
         [~,ind_wekd] =min(pdist2(torg,week),[],2);
         temp(:,3)=week(ind_wekd);


         temp(:,1)=round(temp(:,1).*scaleY(ind_wekd,1)+ scaleY(ind_wekd,2));
         temp(:,2)=round(temp(:,2)*tDomain.scale{3})+1;        
        for ii1=size(Dset,1)
            ind_deg= find(ismember(tDomain.Dset,temp(ii1,:),'rows'));
        end
         
         Dset=[temp Dset  ind_deg torg];
         
   end





