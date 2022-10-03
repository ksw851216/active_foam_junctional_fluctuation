function ne=ns_faceNeiCount(fEd,eFa,fId)

neFa=eFa(abs(fEd),:);
neFa=unique(neFa);
neFa=neFa(neFa~=fId);
neFa=neFa(neFa~=0);

ne=size(neFa,1);

end