function vSrt=ns_vrtxIdSort(vLst,vCrd,bs,sstn)    

vCrd=ns_crdLocal(vCrd(vLst,:),bs,sstn); 

cCrd=[mean(vCrd(:,1)),mean(vCrd(:,2))];
angV=atan2(vCrd(:,2)-cCrd(2),vCrd(:,1)-cCrd(1));

[~,idc]=sort(angV);
vSrt=vLst(idc);    

end