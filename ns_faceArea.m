function [ar,cen]=ns_faceArea(fEd,eMd,bs,sstn)

fCrd=ns_faceVrtxCrd(fEd,eMd,bs,sstn);

fCrd=[fCrd;fCrd(1,:)];
ar=0;
cen=zeros(1,2);
for vrc=1:size(fCrd,1)-1
    ar = ar + det(fCrd(vrc:vrc+1,:))/2;
    cen = cen + (fCrd(vrc,:)+fCrd(vrc+1,:))*det(fCrd(vrc:vrc+1,:))/6;
end
cen=cen/ar;

end