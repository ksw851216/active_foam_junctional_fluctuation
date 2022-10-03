function ns_plotDomain(edge,face,rg,gmp)

figure;
for fac=1:gmp.nFa
    if face{1}(fac,rg.fi(1))==1
        fEd=face{3}{fac};
        fVr=ns_faceVrtxCrd(fEd,edge{2},gmp.bs,gmp.sstn);
        ns_plotPolygon(fVr,gmp.bs,face{1}(fac,rg.fi(4)),gmp.sstn);
    end
end

axis([0 gmp.bs 0 gmp.bs]);
pbaspect([1 1 1]);

end