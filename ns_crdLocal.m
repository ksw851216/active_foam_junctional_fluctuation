function vrF=ns_crdLocal(vrL,bs,sstn)

vrF=vrL;

for j=2:size(vrL,1)
    if vrF(j,2)-vrF(j-1,2)>bs/2
        vrF(j,2)=vrF(j,2)-bs;
        vrF(j,1)=vrF(j,1)-bs*sstn;
    end
    if vrF(j,2)-vrF(j-1,2)<-bs/2
        vrF(j,2)=vrF(j,2)+bs;
        vrF(j,1)=vrF(j,1)+bs*sstn;
    end
    
    if vrF(j,1)-vrF(j-1,1)>bs/2
        vrF(j,1)=vrF(j,1)-bs;
    end
    if vrF(j,1)-vrF(j-1,1)<-bs/2
        vrF(j,1)=vrF(j,1)+bs;
    end
end

end