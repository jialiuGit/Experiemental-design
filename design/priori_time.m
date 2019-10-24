

function tDomain = priori_time(n,varargin)

%Copyright (c) 2019-   Jia Liu

%This function aims to include an uneven distribution / layout on time
%domain.

%n---size of the design

%set a convex mean function which has a linear trend of the frequency
%fun, f(t) >0, considering uneven distributed status of the time in the design
%space, m = @(t) a-(t(:)-b).^2;

%n=sum(f(t)) = a2*sum(m) +b2, i.e. a2=max(mt), b2=min(mt), if
%min(min)<0;else b2=0;

ip=inputParser;
ip.FunctionName = 'PRIORI_TIME';  
  
ip.addRequired('n', @(x) ~isempty(x) &&isscalar(x) && isfinite(x));
ip.addParameter('t',linspace(0,1,100), @(x) isvector(x) && isreal(x));
ip.addParameter('a',2, @(x) isscalar(x) && isfinite(x))
ip.addParameter('b',0.5, @(x) isscalar(x) && isfinite(x))
ip.addParameter('c',30, @(x) isscalar(x) && isfinite(x))
ip.addParameter('resource','meanfun',@(x) ischar(x) || iscell(x)); 
ip.addParameter('fileName','/home/jiliu/monotonic_codes/experimental_design/time',@isFunc);
ip.addParameter('p',@(x) isvector(x) && isreal(x));
ip.addParameter('predt',@(x) isvector(x) && isreal(x));
ip.parse(n,varargin{:});

a=ip.Results.a;
b=ip.Results.b; 
c=ip.Results.c; 
t=ip.Results.t;
resource = ip.Results.resource;
fileName = ip.Results.fileName;
p= ip.Results.p;
predt= ip.Results.predt;


switch resource
    
    case 'data'
        
        load(fileName)
        ha=histogram(time,length(unique(time)));
        inProb2=ha.Values./size(time,1)*n; %size of the design


        %check ----sum(inProb2)=n;

        inProb = inProb2./max(inProb2); %rescale to [0 1];

        %rescale time
        tDomain.t=ha.BinEdges./365; %365 whole calender day
        tDomain.aptRate = inProb;
        tDomain.tReg=tDomain.t([1 end]);
        
    case 'meanfun'
 
        m = @(t)a-c*(t(:)-b).^2;
        mt=m(t);
        
        mint=a-c*(1-b)^2;
        
        mt=(mt-mint)/(a-mint);
      
        inProb = mt;
      
        
       
        %figure
        %plot(t,freq);
        %inProb=mt/M;
        %inProb = freq./max(freq); 
        %inProb = freq./sum(freq);
        
        
        tDomain.t=t;
        tDomain.aptRate = inProb;
        tDomain.tReg=tDomain.t([1 end]);
        tDomain.tDelta=max(m(t))-min(m(t));
        
        
    case 'pred_mean'        
        predt=predt-min(predt);%rescale to [0,1]
        tDomain.t= predt./max(predt); 
       
        
        if min(p)<0
            p2=p-min(p);
        end
        
        p2(p<=-1)=0;
        inProb = p2./max(p2);
        tDomain.aptRate = inProb;
        tDomain.tReg=[predt(1) predt(end)];
         
        
        
end
    

    
    

    


