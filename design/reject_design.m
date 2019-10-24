function [designDest,n]=reject_design(tDomain,design, plimt)
%Copyright (c) 2019-   Jia Liu

%x is a design candidate point generated from [0,1]

% note that with mean function the time domain of t must be same as the study region on time


    
if nargin<3
    plimt=1;
end

if ~isfield(tDomain,'Dset')
    sptmp = false;
else
    sptmp = true;  
end

 
x=[];


if sptmp
    
    design=design(:,[2 1 3]); %This is for special case in imaging plotting 
    
     %first check if the generated design is in the study region.
    if isfield(tDomain,'scale') 
         temp=design;
         %for fast sampling, we only scaling the design samplers from [0,1] to the
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
         %only check spatial!
         ind2=find(ismember(temp,tDomain.Dset,'rows')); %first rule out the land area.
         
         %start regject method
         if ~isempty(ind2)
             design=design(ind2,:);
             temp=temp(ind2,:);
             torg=torg(ind2);

             xorg=[temp torg]; 
             for ii1=1:size(temp,1)                      
                     u=plimt*rand(1);
                     ind_deg= find(ismember(tDomain.Dset,temp(ii1,:),'rows'));
                 if (tDomain.aptRate(ind_deg) >=u)
                     
                    x=[x; temp(ii1,:) design(ii1,:)  ind_deg torg(ii1)]; %the last column is the index for retrieve meanfun of the design
                 else
                    x=x;
                 end
             end
            
         else
             x=x;
             xorg=[];
         end
         
         designDest=struct('x',x,'xorg',xorg);        
       

    else 
          for ii1=1:size(design,1)
            
            u=rand(1); 
            ind_deg= ismember(tDomain.Dset,design(ii1,:),'rows');
            if (tDomain.aptRate(ind_deg) >=u)               
                x=[x;design(ii1,:) ind_deg]; %the last column is the index for retrieve meanfun of the design
            else
                x=x;
            end
        
          end 
          designDest=x;
     end

    
    
else
 
 
     
    for ii1=1:size(design,1)

        [~,ind_t]=min(pdist2(design(ii1,3), tDomain.t'));
       
        u=rand(numel(ind_t),1);
        ind_t=ind_t( tDomain.aptRate(ind_t) >=u);

        if ~isempty(ind_t) 

            x=[x;design(ii1,:)];
        else
            x=x;
        end
    
    end
    designDest=x;
end


n=size(x,1) ;
 
    
end