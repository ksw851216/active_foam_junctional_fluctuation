function vrtxN=ns_crdDomain(vrtx,rg,gmp)

vrtxN=vrtx;
for ii=1:size(vrtxN,1)
    if vrtxN(ii,rg.vf(3))>gmp.bs
        vrtxN(ii,rg.vf(3))=vrtxN(ii,rg.vf(3))-gmp.bs;
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))-gmp.bs*gmp.sstn;
    elseif vrtxN(ii,rg.vf(3))<0
        vrtxN(ii,rg.vf(3))=vrtxN(ii,rg.vf(3))+gmp.bs;
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))+gmp.bs*gmp.sstn;
    end

    if vrtxN(ii,rg.vi(3))>gmp.bs+gmp.sstn*vrtxN(ii,rg.vf(3))
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))-gmp.bs;
    elseif vrtxN(ii,rg.vi(3))<gmp.sstn*vrtxN(ii,rg.vf(3))
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))+gmp.bs;
    end
end

end