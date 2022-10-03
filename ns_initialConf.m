function [vrtx,edge,face,rg]=ns_initialConf(bs)

[gmp,mcp]=deal(struct);
[gmp.tgAr,gmp.sstn,gmp.bs,gmp.edpc,gmp.asp]=deal(1,0,bs,0.15,0.8);
[mcp.del,mcp.psi,mcp.mu,mcp.taut,mcp.nvf,mcp.nstf,mcp.taun,mcp.taua,mcp.mun,mcp.nrst]=...
    deal(10,1,0,10,0.7,1,0.1,0.1,0,1);
[gmp.nFa,gmp.lsc,gmp.dt]=...
    deal(gmp.bs^2,mcp.psi*sqrt(gmp.tgAr),0.005);
[gmp.major,gmp.minor]=deal(sqrt(mcp.nvf*mcp.psi/pi/gmp.asp),...
    sqrt(mcp.nvf*mcp.psi/pi*gmp.asp));

gmp.ellDis=[gmp.major*(cos(0:pi/30:2*pi-pi/30)).',...
    gmp.minor*(sin(0:pi/30:2*pi-pi/30)).';gmp.major,0];
gmp.ellAng=[mod(atan2(gmp.ellDis(1:end-1,2),gmp.ellDis(1:end-1,1)),2*pi);2*pi];
gmp.ellDis=sqrt(sum(gmp.ellDis.^2,2));


%% Generate an initial configuration.
rng('shuffle');

% randomly generate seed points and compute voronoi Tessellation from that.
sdPt=rand(gmp.nFa,2)*gmp.bs;          
[vrtx,edge,face,rg]=ns_randomVoronoi(gmp,sdPt);
clear sdPt;

gmp.eLnc=ns_edgeLenAll(edge,gmp);
shC=10;
gmp.shEd=0.001*2*sqrt(pi)*shC;

for itc=1:1500
    [vrtx,edge,face]=ns_iteration(vrtx,edge,face,rg,gmp,mcp);

    if mod(itc,5)==0
        eLn=ns_edgeLenAll(edge,gmp);
        t1Id=find(eLn<gmp.shEd);
        if isempty(t1Id)==0
            for jj=1:size(t1Id,1)
                chk=min(face{1}(edge{1}(t1Id(jj),rg.ei(2):rg.ef(2)),rg.fi(5)));
                if chk>3 && gmp.eLnc(t1Id(jj))>eLn(t1Id(jj))
                    [vrtx,edge,face]=...
                        ns_t1Transition(vrtx,edge,face,rg,gmp,mcp,t1Id(jj));
                end
            end
        else
            if shC<90
                shC=shC+1;
                gmp.shEd=0.001*2*sqrt(pi)*shC;
            end
        end        
        gmp.eLnc=eLn;
    end
end

end