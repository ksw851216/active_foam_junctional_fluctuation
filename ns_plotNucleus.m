function ns_plotNucleus(fCen,orn,gmp)

col = [1,0.8,0.8];

pb=[1+gmp.sstn,1;1,0;1-gmp.sstn,-1;...
    gmp.sstn,1;0,0;-gmp.sstn,-1;-1+gmp.sstn,1;-1,0;-1-gmp.sstn,-1]*gmp.bs;

vrTmp=cell(size(pb,1),1);
for ii=1:size(pb,1)
    vrTmp{ii}=repmat(fCen+pb(ii,:),40,1)+...
        gmp.major*[cos(orn),sin(orn)].*(cos(0:pi/20:2*pi-pi/20)).'+...
        gmp.minor*[-sin(orn),cos(orn)].*(sin(0:pi/20:2*pi-pi/20)).';
end

for ii=1:size(pb,1)
    patch(vrTmp{ii}(:,1),vrTmp{ii}(:,2),col,'FaceAlpha',0.5);
    hold on;
end
    
end