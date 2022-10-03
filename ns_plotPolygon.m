function ns_plotPolygon(vrList,bs,nei,sstn)

col = [0.9,0.9,0.9;0.8,0.8,1;1,0.8,0.9;1,0,0;0,1,0;1,1,0;...
    0.7,0.7,0.7;0,0,1;1,0,1;0,1,1;0.9412,0.4706,0;0.251,0,0.502;];

pb=[1+sstn,1;1,0;1-sstn,-1;...
    sstn,1;0,0;-sstn,-1;-1+sstn,1;-1,0;-1-sstn,-1]*bs;    
vrTmp=cell(size(pb,1),1);
for ii=1:size(pb,1)
    vrTmp{ii}=vrList+repmat(pb(ii,:),size(vrList,1),1);
end

for ii=1:size(pb,1)
    if nei>11
        patch(vrTmp{ii}(:,1),vrTmp{ii}(:,2),'white','FaceAlpha',0.5);
        hold on;
    else
        patch(vrTmp{ii}(:,1),vrTmp{ii}(:,2),col(nei+1,:),'FaceAlpha',0.7);
        hold on;
    end
end
    
end