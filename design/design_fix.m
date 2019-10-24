
function design=design_fix(D_N,varargin)

%Copyright (c) 2019-   Jia Liu 

%weights----is a n*1 vector

%grid_reg---A Degsin set with pre-defined resolution, the default
%resolution for 2D is 20x20, for 3D is 20x20x20.

%dist--methods of distance, see pdist, pdist2, the default dist is euc.

%delta--threthold distance between two points in sampling, the default
%value for min_dist design is 0.06. 

%n--size of design

%nclosepair--number of close pair when calling close pair min_dist design,
%the default value is n/2.

%D_N--number of all the points in the design space; for the deterministic
%and completely random design w/o additional design criterion, D_N points
%to the size of the design, i.e., D_N=n

%dim-- dimensions of the design space

%x-- is the design based minimun distance

%name--the design, the default design is sampling design equvilent to the 
%sampling method.

%a distance that weights each coordinate contribution differently
%only work with Euclidean distance.
 
  
  ip=inputParser;
  ip.FunctionName = 'RANDSAMPDEG';     
  ip.addRequired('D_N', @(x) ~isempty(x) &&isscalar(x) && isfinite(x))
  ip.addParameter('delta',0.06, @(x) isscalar(x) && isfinite(x)) 
  ip.addParameter('n',D_N, @(x) isscalar(x) && isfinite(x))
  ip.addParameter('nclosepair', @(x) isscalar(x) && isfinite(x))
  ip.addParameter('dim',2, @(x) isscalar(x) && isfinite(x)); %dim
  ip.addParameter('grid_reg',[], @(x) isvector(x) && isreal(x));
  ip.addParameter('gridZreg',[], @(x) isvector(x) && isreal(x));
  ip.addParameter('weights',[], @(x) isvector(x) && isreal(x));  
  ip.addParameter('samp_method','lattice',@(x) ischar(x) || iscell(x));  
  ip.addParameter('name','lattice',@(x) ischar(x) || iscell(x));   
  ip.addParameter('dist','euclidean', @(x) ischar(x) || iscell(x));
  ip.addParameter('Pc',5, @(x) isscalar(x) && isfinite(x) ||ismatrix(x)||isvector(x) && isreal(x) );   
  ip.addOptional('design',[], @isstruct);
   
  ip.parse(D_N,varargin{:});
  design=ip.Results.design;
  
  
  % Initialize a design
  if isempty(design)  
      init=true;
  end
  
    % sampling method
  if init || ~ismember('samp_method',ip.UsingDefaults)
      design.samp_method=ip.Results.samp_method;
  end  

   % name of the design
  if  init ||~ismember('name',ip.UsingDefaults)      
      if ~strcmp(design.samp_method,'lattice')&&~ismember(ip.Results.name,{'min_dist' 'rand_lattice' 'close_pair' 'space_fill'})
        design.name=design.samp_method;
      else
        design.name=ip.Results.name;
      end
  end  
    
  % type of distance (minimun distance)
  if  init ||~ismember('dist',ip.UsingDefaults)%&& ...
          %ismember(design.name,{'min_dist' 'space_fill'})
      design.dist=ip.Results.dist;
  end
  %  size of design
  if  init ||~ismember('n',ip.UsingDefaults)
      design.size=ip.Results.n;
  end
 
     % dimension of the design space
  if init || ~ismember('dim',ip.UsingDefaults)
      design.dim=ip.Results.dim;
  end
   % threshold (e.g., minimum distance in min_dist design)
  if init || ~ismember('delta',ip.UsingDefaults)
      if ismember(ip.Results.name,{'min_dist' 'close_pair'})
          design.threshold=ip.Results.delta;
      end
  end  
  
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
  

if  ismember(design.name,{'min_dist' 'space_fill' 'rand_lattice'}) && design.size==D_N
    error('The design size, n, must be smaller than D_N the size of candidate set!')
end
% dimension of the design space
if  ismember('weights',ip.UsingDefaults)
    design.weights=false;
else
    design.weights=true;
end

switch  design.samp_method
    
    case 'Sobol'
        %sobol random sequence in 2 D
        p = sobolset(design.dim,'Skip',1e3,'Leap',1e2); %d---dimensions of design space
        
        Dset = net(p,D_N);
        
    case 'Halton'
        %sobol random sequence in 2 D
        p = haltonset(design.dim,'Skip',1e3,'Leap',1e2); %d---dimensions of design space
        
        Dset = net(p,D_N);
        
    case 'lattice'
        lat_size=nthroot(D_N,design.dim);
        if design.dim == 3
            [X, Y,Z] = meshgrid(grid_reg,grid_reg,gridZreg);
            Dset=[X(:) Y(:) Z(:)];
            
            if length(gridZreg)==lat_size
               lat=num2str(lat_size);
               disp(['The resolution of lattice is ',lat, 'x', lat, 'x', lat, '3D grids.'])
            end
        else
            [X, Y] = meshgrid(grid_reg,grid_reg); %here d is resoluation of the z axis.
            Dset=[X(:) Y(:)];
            if length(grid_reg)==lat_size
                lat=num2str(lat_size);
                disp(['The resolution of lattice is ',lat, 'x', lat,' 2D grids.'])
            end
        end
        
    case 'Random'
        
        Dset=rand(D_N,design.dim);
    case 'Fibonacci lattice'
        if design.dim==3
            Dset = Fibo_design_3D(D_N);
        else
            Dset = Fibo_design(D_N);
        end
end

if strcmp(design.name,'min_dist')
    delta=design.threshold;
    x=datasample(Dset,design.size,'Replace',false);
    xbar = setdiff(Dset,x,'rows');
    
    if  design.weights
        weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
    end
    
    i=1;
    xbar0=xbar;
    nreplace=0;
    while i<= design.size
        x2=x;
        x2(i,:)=[];
        if design.weights
            distS = pdist2(x(i,:),x2,@(Xi,Xj) weuc(Xi,Xj,weights));
        else
            if ismember(design.dist,{'Minkoski' 'mahalanobis' 'seuclidean'})
                distS = pdist2(x(i,:),x2,design.dist,Pc);
            elseif ismember(design.dist,{'squaredeuclidean'...
                    'cityblock' 'chebychev' 'hamming'})
                distS = pdist2(x(i,:),x2,design.dist);
            else
                distS = pdist2(x(i,:),x2);
            end
        end
        
        if min(distS) >= delta
            i=i+1;
            xbar=setdiff(xbar0,x,'rows');
        else
            
            try          
                x(i,:)=datasample(xbar,1,1,'Replace',false);
                
                xbar = setdiff(xbar,x(i,:),'rows');
                nreplace=nreplace+1; 
                 
            catch
                delta=delta*0.95;
                i=1;                
                design.threshold=delta;
            end
        end
    end
elseif strcmp(design.name,'rand_lattice')
    x=datasample(Dset,design.size,'Replace',false);
    
elseif strcmp(design.name,'close_pair')
    delta=design.threshold;
    k=design.nclosepair;
    x= closeParisSamp(design.dim,design.size, k, delta,design.dist);
elseif    strcmp(design.name,'space_fill')
    x=coffee_house_design(Dset,design.size,design.dist);
else
    x=Dset;
end


if ~isempty(x)
    design.x=x;
    if design.dim==2
        if design.size>1
            plot(x(:,1),x(:,2),'ro');
        else
            plot(x(1),x(2),'ro');
        end
        xlabel('x')
        ylabel('y')            
        axis normal
    elseif design.dim==3
        if  design.size>1
            plot3(x(:,1),x(:,2),x(:,3),'ro');
        else
            plot3(x(1),x(2),x(3),'ro');
        end
        xlim([0 1]);        
        xlabel('x')
        ylabel('y') 
        zlabel('time')
         axis equal 
         
    end
    %if ~strcmp(design.samp_method,design.name)
        %title([ design.samp_method ' ' design.name ' design'],'Interpreter','none')
%     if any(strcmp(design.name,{'close_pair','min_dist','space_fill'}))
%         title([ design.name ' design'],'Interpreter','none')
%     else
%         title([ design.samp_method ' design'],'Interpreter','none')
%     end
%     set(gcf, 'Visible', 'off');
    %delete(gca);
end
end