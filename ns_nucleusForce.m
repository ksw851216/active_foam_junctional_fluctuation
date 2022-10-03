function [eFrc,ncFrc,ncMmt]=ns_nucleusForce(gmp,mcp,fcen,pcen,forn,eMd,ecDis)

eVec=eMd-repmat(fcen,size(eMd,1),1);
eAng=mod(atan2(eVec(:,2),eVec(:,1)),2*pi);
pAng=mod(eAng-forn,2*pi);

eqDis=zeros(size(pAng,1),1);
for ptc=1:size(pAng,1)
    idx=find(gmp.ellAng<pAng(ptc));
    idx=idx(end);
    if idx==size(gmp.ellDis,1)
        dsTmp=[gmp.ellDis(idx);gmp.ellDis(1)];
        angTmp=[gmp.ellAng(idx);gmp.ellAng(1)+2*pi];
    else
        dsTmp=gmp.ellDis(idx:idx+1);
        angTmp=gmp.ellAng(idx:idx+1);
    end
        
    eqDis(ptc)=dsTmp(2)+(dsTmp(1)-dsTmp(2))*(angTmp(2)-pAng(ptc))/(angTmp(2)-angTmp(1));
end


eFrc=zeros(size(eMd,1),2);
ncFrc=zeros(1,2);

for edc=1:size(eMd,1)-1
    if ecDis(edc)<eqDis(edc) || ecDis(edc+1)<eqDis(edc+1)
        
        ptCnt=10;
        ptCrd=[(eMd(edc,1):(eMd(edc+1,1)-eMd(edc,1))/ptCnt:eMd(edc+1,1)).',...
            (eMd(edc,2):(eMd(edc+1,2)-eMd(edc,2))/ptCnt:eMd(edc+1,2)).'];
        eCen=ptCrd(ptCnt/2+1,:);
        
        nVec=[eMd(edc+1,2)-eMd(edc,2),eMd(edc,1)-eMd(edc+1,1)];
        nVec=nVec/norm(nVec);
        if dot(nVec,eCen-pcen)<0
            nVec=-nVec;
        end

        rVec=ptCrd-fcen;        
        rDis=sqrt(sum(rVec.^2,2));

        eLen=norm(ptCrd(end,:)-ptCrd(1,:));
        lenSt=eLen/ptCnt;
        
        eSegVec=ptCrd-repmat(fcen,ptCnt+1,1);
        eSegAng=atan2(eSegVec(:,2),eSegVec(:,1));
        pSegAng=mod(eSegAng-forn,2*pi);
        
        eqSegDis=zeros(size(pSegAng,1),1);
        for ptc=1:size(pSegAng,1)
            idx=find(gmp.ellAng<pSegAng(ptc));
            idx=idx(end);
            if idx==size(gmp.ellDis,1)
                dsTmp=[gmp.ellDis(idx);ellDis(1)];
                angTmp=[gmp.ellAng(idx);ellAng(1)+2*pi];
            else
                dsTmp=gmp.ellDis(idx:idx+1);
                angTmp=gmp.ellAng(idx:idx+1);
            end

            eqSegDis(ptc)=dsTmp(2)+(dsTmp(1)-dsTmp(2))*(angTmp(2)-pSegAng(ptc))/(angTmp(2)-angTmp(1));
        end
        
        lenArr=[lenSt/2;repmat(lenSt,ptCnt-1,1);lenSt/2];
        hvFn=(rDis<eqSegDis);
        
        frc=mcp.nstf.*(eqSegDis/gmp.lsc-rDis/gmp.lsc).*hvFn.*...
            lenArr/gmp.lsc.*nVec;
        lVec=ptCrd-eCen;
        mmt=(lVec(:,1)/gmp.lsc).*frc(:,2)-(lVec(:,2)/gmp.lsc).*frc(:,1);

        frcAll=sum(frc);
        mmtAll=sum(mmt);
        if abs(mmtAll)>0
            pVal=1/2+1/2*mmtAll/det([ptCrd(1,:)-eCen;frcAll]);
        else
            pVal=1/2;
        end
        eFrc(edc,:)=eFrc(edc,:)+pVal*frcAll;
        eFrc(edc+1,:)=eFrc(edc+1,:)+(1-pVal)*frcAll;
        ncFrc=ncFrc-frcAll;
    end
end

rVec=eMd-fcen;
ncMmt=sum(-(rVec(:,1)/gmp.lsc.*eFrc(:,2)-rVec(:,2)/gmp.lsc.*eFrc(:,1)));
        
end