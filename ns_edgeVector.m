function [tv,nv,eln]=ns_edgeVector(eVs,eVf,fCn,bs,sstn)

tv=eVf-eVs;
md=(eVf+eVs)/2;

eln=norm(tv);
tv=tv/eln;

nv=[-tv(2),tv(1)];
fCmp=[fCn;md];
fCmp=ns_crdLocal(fCmp,bs,sstn);
fv=fCmp(2,:)-fCmp(1,:);
if dot(nv,fv)<0
    nv=-nv;
end

end