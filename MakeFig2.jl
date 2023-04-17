##########################################################################
##########################################################################
##########################################################################
#
# This script generates figure 2 of the paper:
#
# "Spike shape and synaptic-amplitude distribution interact to set ...
# the high-frequency firing-rate response of neuronal populations"
# MJE Richardson (2018) to appear in Physical Review E.
#
# The code is in Julia version 1.0 (https://julialang.org/)
#
##########################################################################
#
# COPYRIGHT: Magnus Richardson (2018).
# This code is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# <http://www.gnu.org/licenses/>
#
##########################################################################
##########################################################################
##########################################################################

using PyPlot
using Random
using SpecialFunctions

function MakeFig2()


    Random.seed!(1)
    K=1000.0

    ######################################################################
    # The neuronal parameters
    ######################################################################
    tau=20
    vT=10
    vre=5
    vth=30

    dT=0.6
    asA=0.2
    asB=0.6
    asC=1.8

    dv=0.0001;

    # finds the synaptic rate for a given output rate for cases A, B and C
    r0aim=5.0/1000

    rA,RsA=Rintorout(r0aim,tau,vT,dT,vre,vth,asA,dv)
    rB,RsB=Rintorout(r0aim,tau,vT,dT,vre,vth,asB,dv)
    rC,RsC=Rintorout(r0aim,tau,vT,dT,vre,vth,asC,dv)

    println("rA=$(round(rA*K;digits=4))Hz and RsA=$(round(RsA;digits=4))kHz")
    println("rB=$(round(rB*K;digits=4))Hz and RsB=$(round(RsB;digits=4))kHz")
    println("rC=$(round(rC*K;digits=4))Hz and RsC=$(round(RsC;digits=4))kHz")

    ######################################################################
    ######################################################################
    # Generate the data for each panel
    ######################################################################
    ######################################################################

    ######################################################################
    # Generate data for figure 2A
    ######################################################################

    # methods fig for the paper
    asE=asC
    RsE=RsC
    tE,RstE,vE,rtE,tEb,rtEsimb,fE,PE,TE=nicemethodsfig(dv,vth,vre,RsE,asE,tau,vT,dT);

    ######################################################################
    # Generate data for figure 2B
    ######################################################################

    # ThIn method for steady-state and modulation

    nfth=33; lfth=0.1/1000; hfth=1000/1000;
    fth=exp.(collect(range(log(lfth),stop=log(hfth),length=nfth)))
    w=2*pi*fth

    rhA=zeros(Complex{Float64},nfth)
    rhB=zeros(Complex{Float64},nfth)
    rhC=zeros(Complex{Float64},nfth)
    for k=1:nfth
        ~,rhA[k],~,~,~=EshotEIFThIn(dv,vth,vre,RsA,asA,tau,vT,dT,fth[k],1)
        ~,rhB[k],~,~,~=EshotEIFThIn(dv,vth,vre,RsB,asB,tau,vT,dT,fth[k],1)
        ~,rhC[k],~,~,~=EshotEIFThIn(dv,vth,vre,RsC,asC,tau,vT,dT,fth[k],1)
    end
    ahA=abs.(rhA);   ahB=abs.(rhB);   ahC=abs.(rhC)
    phA=angle.(rhA); phB=angle.(rhB); phC=angle.(rhC)

    # calculate the frequency domain asymptotics

    # case A as<dT
    r1thasymA=rA*tau*asA./(1im*w*tau*(dT-asA));

    # case B as=dT
    ~,vB,PB,~=EshotEIFThIn0(dv,vth,vre,RsB,asB,tau,vT,dT);
    B1=(PB.*dT.*exp.(-(vT .-vB)/dT)/(tau*rB) .-1);
    I=(vT-vB[1])/dT + sum(B1)*dv/dT;
    c=0.65522*exp(I); # the integral and the 1/exp(1-gamma) term
    r1thasymB=rB*tau*log.(c*1im*w*tau)./(1im*w*tau);

    # case C as>dT
    ~,vC,PC,~=EshotEIFThIn0(dv,vth,vre,RsC,asC,tau,vT,dT);
    Is=sum(dv*PC.*exp.((vC .-vT)/asC))/(rC*tau);
    r1thasymC=rC[1]*tau*Is*gamma(1+dT/asC)./(1im*w*tau).^(dT/asC);

    asymA=abs.(r1thasymA);
    thasymA=-ones(length(w))*pi/2;
    asymB=abs.(r1thasymB); thasymB=-pi/2 .+atan.(pi./(2*log.(c*w*tau)));
    asymC=abs.(r1thasymC); thasymC=-ones(length(w))*pi*dT/(2*asC);

    ######################################################################
    ######################################################################
    # Now plot the figures
    ######################################################################
    ######################################################################

    # set up the figure
    fig2=figure(figsize=(10,3.5));

    # figure parameters

    # useful constants
    K=1000;
    r2d=360/(2*pi)

    lw=0.75
    be=0.125
    ms=3;

    labelfs=8
    panelfs=13
    tickfs=8
    textfs=6

    ######################################################################
    # Plot Figure 2A
    ######################################################################

    pwA=0.25
    phAA=0.125
    leA=0.075

    ax1=PyPlot.axes([leA,0.75,pwA,phAA])
    plot(tE,RstE,"k-",linewidth=lw)
    plot([tE[1],tE[end]],RsC*[1,1],"k:",linewidth=lw)
    ax1[:set_yticks]([0,0.1,0.2])
    axis([TE-2*PE,TE,0,0.2])
    ylabel(L"R_s(t)~"*" (kHz)",fontsize=labelfs)

    ax2=PyPlot.axes([leA,0.45,pwA,0.175])
    plot(tE,vE,linewidth=lw)
    ylabel("Voltage (mV)",fontsize=labelfs)
    ax2[:set_yticks]([0,10,20,30])
    axis([TE-2*PE,TE,-5,30])

    ax3=PyPlot.axes([leA,be,pwA,0.225])
    plot(tE,K*rtE,linewidth=lw,"k")
    plot([tE[1],tE[end]],5*[1,1],"k:",linewidth=lw)
    fill_between(tEb,0,K*rtEsimb,linewidth=lw,color="lightgray")
    plot(tEb,K*rtEsimb,linewidth=lw,color="gray")
    axis([TE-2*PE,TE,0,7.5]);
    ax3[:set_yticks]([0,2.5,5,7.5])
    ylabel("Rate "*L"r(t)"*" (Hz)",fontsize=labelfs)
    xlabel("Time (ms)",fontsize=labelfs)

    for ax in (ax1,ax2,ax3)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig2[:text](0.02,0.925, "A", fontsize=panelfs)

    ######################################################################
    # Plot Figure 2B
    ######################################################################

    pwB=0.225
    phBB=0.325
    leB=0.45

    ax4=PyPlot.axes([leB,0.575,pwB,phBB])
    loglog(K*fth,abs.(rhA/rhA[1]),"b-",linewidth=lw)
    loglog(K*fth,abs.(rhB/rhB[1]),"g-",linewidth=lw)
    loglog(K*fth,abs.(rhC/rhC[1]),"r-",linewidth=lw)
    loglog(K*fth,asymA/ahA[1],"b:",linewidth=lw)
    loglog(K*fth,asymB/ahB[1],"g:",linewidth=lw)
    loglog(K*fth,asymC/ahC[1],"r:",linewidth=lw);
    axis([K*lfth,K*hfth,0.02,1.2]);
    ylabel("Normalised amplitude",fontsize=labelfs)

    ax4[:text](100,0.035,"case (i)",fontsize=textfs)
    ax4[:text](100,0.025,L"a_s<\delta_\mathrm{T}",fontsize=textfs)

    ax4[:text](500,0.12,"case (ii)",fontsize=textfs)
    ax4[:text](800,0.08,L"a_s=\delta_\mathrm{T}",fontsize=textfs)

    ax4[:text](300,0.5,"case (iii)",fontsize=textfs)
    ax4[:text](300,0.35,L"a_s>\delta_\mathrm{T}",fontsize=textfs)

    ax5=PyPlot.axes([leB+0.05,0.68,0.06,0.15])
    loglog(K*fth,K*abs.(rhA),"b-",linewidth=lw)
    loglog(K*fth,K*abs.(rhB),"g-",linewidth=lw)
    loglog(K*fth,K*abs.(rhC),"r-",linewidth=lw)
    loglog(K*fth,K*asymA,"b:",linewidth=lw)
    loglog(K*fth,K*asymB,"g:",linewidth=lw)
    loglog(K*fth,K*asymC,"r:",linewidth=lw)
    ylabel("Amp. (Hz)",fontsize=7)
    ax5[:set_yticks]([1,10,100])
    xlabel("Freq. (Hz)",fontsize=7,labelpad=0)
    ax5[:set_xticks]([1,1000])
    #label_params(axis="x", pad=-1)
    yticks(rotation=90,fontsize=6)
    xticks(fontsize=6)
    axis([1,K*hfth,1.0,100]);
    #xlabel("f (Hz)",fontsize=labelfs)

    ss=findall(K*fth.>100);

    ax6=PyPlot.axes([leB,be,pwB,phBB])
    semilogx(K*fth,phA*r2d,"b-",linewidth=lw)
    semilogx(K*fth,phB*r2d,"g-",linewidth=lw)
    semilogx(K*fth,phC*r2d,"r-",linewidth=lw)
    semilogx(K*fth,thasymA*r2d,"b:",linewidth=lw)
    semilogx(K*fth[ss],thasymB[ss]*r2d,"g:",linewidth=lw)
    semilogx(K*fth,thasymC*r2d,"r:",linewidth=lw)
    axis([K*lfth,K*hfth,-120,0]);
    xlabel("Modulation frequency (Hz)",fontsize=labelfs)
    ylabel("Phase (deg)",fontsize=labelfs)
    ax6[:set_yticks]([0,-30,-60,-90,-120])

    ax6[:text](100,-70,"(i)",fontsize=textfs)
    ax6[:text](200,-52,"(ii)",fontsize=textfs)
    ax6[:text](400,-25,"(iii)",fontsize=textfs)

    for ax in (ax4,ax6)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig2[:text](0.375,0.925, "B", fontsize=panelfs)

    ######################################################################
    # Plot Figure 2C
    ######################################################################

    pwC=0.15
    phCC=0.325
    leC=0.8

    xx=collect(0:0.01:10)
    yy=1*(xx.<1.0)+(xx.>=1.0)./xx
    zz=-90*(xx.<1.0)-(xx.>=1.0)*90 ./xx

    ax7=PyPlot.axes([leC,0.575,pwC,phCC])
    plot(xx,yy,"k-",linewidth=lw)
    plot(asA/dT,1,"bo",markersize=ms)
    plot(asB/dT,1,"go",markersize=ms)
    plot(asC/dT,dT/asC,"ro",markersize=ms)
    fill_between([0,1],[0,0],[1.5,1.5],color="lightgray")
    xlabel(L"a_s / \delta_\mathrm{T}",fontsize=labelfs)
    ylabel("Exponent "*L"\beta",fontsize=labelfs)
    axis([0,4,0,1.5])
    ax7[:text](0.1,0.4, "range of ", fontsize=textfs)
    ax7[:text](0.1,0.25, "diffusion ", fontsize=textfs)
    ax7[:text](0.1,0.1, "approximation", fontsize=textfs)

    ax7[:text](0.3,1.1,"(i)",fontsize=textfs)
    ax7[:text](1.1,1.1,"(ii)",fontsize=textfs)
    ax7[:text](3,0.45,"(iii)",fontsize=textfs)

    ax8=PyPlot.axes([leC,be,pwC,phCC])
    plot(xx,zz,"k-",linewidth=lw)
    plot(asA/dT,-90,"bo",markersize=ms)
    plot(asB/dT,-90,"go",markersize=ms)
    plot(asC/dT,-90*dT/asC,"ro",markersize=ms)
    fill_between([0,1],[0,0],[-120,-120],color="lightgray")
    xlabel(L"a_s / \delta_\mathrm{T}",fontsize=labelfs)
    ylabel("phase (deg)",fontsize=labelfs)
    axis([0,4,-120,0])
    ax8[:set_yticks]([0,-30,-60,-90,-120])

    ax8[:text](0.3,-105,"(i)",fontsize=textfs)
    ax8[:text](1.1,-105,"(ii)",fontsize=textfs)
    ax8[:text](3,-45,"(iii)",fontsize=textfs)

    for ax in (ax7,ax8)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig2[:text](0.725,0.925, "C", fontsize=panelfs)

    ######################################################################
    # Optional save
    ######################################################################

    if true
        savefig("fig2.pdf")
    end

end

#########################################################################
#########################################################################
#########################################################################
# Subfunctions called by the main function
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# finds the steady-state input rates for a particular output rate
#########################################################################
#########################################################################
function Rintorout(r0aim,tau,vT,dT,vre,vth,as,dv)

    # first run
    Rsl=0.01; Rsh=5; Rsn=30;
    RsL=collect(range(Rsl,stop=Rsh,length=Rsn))
    rL=zeros(Rsn)

    for k=1:Rsn
        rL[k],~,~,~=EshotEIFThIn0(dv,vth,vre,RsL[k],as,tau,vT,dT)
    end

    s1=maximum(findall(rL.<r0aim)); R1=RsL[s1]
    s2=minimum(findall(rL.>r0aim)); R2=RsL[s2]

    # second run
    Rsl=R1; Rsh=R2; Rsn=30;
    RsL=collect(range(Rsl,stop=Rsh,length=Rsn))
    rL=zeros(Rsn)

    for k=1:Rsn
        rL[k],~,~,~=EshotEIFThIn0(dv,vth,vre,RsL[k],as,tau,vT,dT)
    end

    s1=maximum(findall(rL.<r0aim)); R1=RsL[s1]; r1=rL[s1]
    s2=minimum(findall(rL.>r0aim)); R2=RsL[s2]; r2=rL[s2]

    Rstar=R1+(r0aim-r1)*(R2-R1)/(r2-r1);
    rstar,~,~,~=EshotEIFThIn0(dv,vth,vre,Rstar,as,tau,vT,dT)

    return rstar,Rstar
end

#########################################################################
#########################################################################
# ThIn for steady-state EIF with one shot-noise process
#########################################################################
#########################################################################
function EshotEIFThIn0(dv,vth,vre,Re,ae,tau,vT,dT)

    # get the fixed points
    vl,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)

    # set up the voltage
    v=collect(vl+dv/2:dv:vth-dv/2)
    F=(dT*exp.((v .-vT)/dT)-v)/tau
    n=length(v)

    # above and below the reset
    krep=minimum(findall(v.>vre))
    krem=krep-1

    # above and below upper fixed point
    nup=minimum(findall(v.>vu))
    num=nup-1

    # above the lower fixed point (NB should be entry 1)
    nlp=minimum(findall(v.>vl))

    p=zeros(n)
    je=zeros(n)

    # useful quantities at the upper fixed point
    dFdvu=(exp((vu .-vT)/dT)-1)/tau
    jeu=1
    pu=(jeu/ae)/(Re+dFdvu)

    #improve the initial conditions around vu - good for stability
    je[num]=jeu-(vu-v[num])*(Re*pu-jeu/ae);
    je[nup]=jeu+(v[nup]-vu)*(Re*pu-jeu/ae);

    #########################################################################
    # integrate vs <-----<< vu ..... vth
    #########################################################################
    for k=num:-1:nlp+1
        p[k]=((v[k]>vre)-je[k])/F[k]
        je[k-1]=je[k]-dv*(Re*p[k]-je[k]/ae)
    end
    p[nlp]=((v[nlp]>vre)-je[nlp])/F[nlp];

    #########################################################################
    # integrate vs ....... vu >>--------> vth
    #########################################################################
    for k=nup:n-1
        p[k]=((v[k]>vre)-je[k])/F[k]
        je[k+1]=je[k]+dv*(Re*p[k]-je[k]/ae)
    end
    p[n]=((v[n]>vre)-je[n])/F[n]

    r=1.0/sum(p*dv)
    P=r*p
    Je=r*je

    return r,v,P,Je

end

#########################################################################
#########################################################################
# Halley's method for lower and upper fixed points
#########################################################################
#########################################################################
function MyfRoots(v1,vT,dT,n)

    v=zeros(n)
    v[1]=v1
    for k=1:n-1
        expk=exp((v[k]-vT)/dT)
        f0=dT*expk-v[k]
        f1=expk-1
        f2=(1/dT)*expk
        v[k+1]=v[k]-2*f0*f1/(2*f1^2-f0*f2)
    end

    return v[end]

end

#########################################################################
#########################################################################
# binning function
#########################################################################
#########################################################################
function binnit(dt,t,tbin,x)

    sbin=Int(ceil(tbin/dt))
    sbin2=Int(sbin/2)
    tb=t[sbin2:sbin:end-sbin2]

    Cx=dt*cumsum(x)
    Cxb=Cx[sbin:sbin:end]
    xb=[Cxb[1];diff(Cxb)]/tbin

    return tb,xb

end

#########################################################################
#########################################################################
# nice methods fig
#########################################################################
#########################################################################
function nicemethodsfig(dv,vth,vre,RsE,asE,tau,vT,dT)

    dtE=0.05
    TE=500
    tE=collect(0:dtE:TE)
    ntE=length(tE)

    fE=10.0/1000
    PE=1.0/fE
    gfE=0.02

    RstE=RsE .+ gfE*sin.(2*pi*fE*tE)

    # get the analytical result
    rE,rhE,~,~,~=EshotEIFThIn(dv,vth,vre,RsE,asE,tau,vT,dT,fE,1)
    ahE=abs(rhE)
    phE=angle(rhE)
    rtE=rE .+ gfE*ahE*sin.(2*pi*fE*tE .+phE)

    # do some sims
    N=2000
    tbinE=2
    nrunsE=10
    ~,vE,rtEsim=EIFshotRtsim(vth,vre,tau,vT,dT,tE,RstE,asE,N)

    for k=2:nrunsE
        ~,vEx,x=EIFshotRtsim(vth,vre,tau,vT,dT,tE,RstE,asE,N)
        rtEsim=rtEsim+x
        vE=[vE vEx]
    end
    tEb,rtEsimb=binnit(dtE,tE,tbinE,rtEsim)
    rtEsimb=rtEsimb/nrunsE

    return tE,RstE,vE,rtE,tEb,rtEsimb,fE,PE,TE

end

#########################################################################
#########################################################################
# ThIn for modulated shot EIF 
#########################################################################
#########################################################################
function EshotEIFThIn(dv1,vth,vre,Re,ae,tau,vT,dT,fkHz,r1thG)

    # f to w
    w=2*pi*fkHz

    # get the fixed points
    vl,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)

    # match grid to fixed points
    su=Int(ceil((vu-vl)/dv1)+1)      # this is the unstable fixed point grid number
    dv=(vu-vl)/(su-1)
    v=collect(vl:dv:vth)
    n=length(v)
    sre=Int(ceil((vre-vl)/dv)+1)     # first point above vre

    # the forcing term and useful functions of it (on half sites)
    vh=v .+dv/2
    fh=(dT*exp.((vh .-vT)/dT)-vh)/tau
    X=dv./fh
    Xw=exp.(-1im*w*X)
    Xe=exp.(-Re*X)
    Xe[su]=0

    vh=v .-dv/2
    fh=(dT*exp.((vh .-vT)/dT)-vh)/tau
    X=dv./fh
    Yw=exp.(1im*w*X)
    Ye=exp.(Re*X)
    Ye[su]=0

    #########################################################################
    # steady-state first
    #########################################################################

    # the value an vu is 1
    je=ones(n)

    # integrate vs ....... vu >>--------> vth
    for k=su:n-1
        je[k+1]=exp(-dv/ae)*(je[k]*Xe[k]+(1-Xe[k]))
    end

    # integrate vs <-----<< vu ..... vth
    for k=su:-1:2
        je[k-1]=exp(dv/ae)*(je[k]*Ye[k]+(k>sre)*(1-Ye[k]));
    end

    # calculation of p using djedv+je/ae=ReP
    djedv=[diff(je);0]/dv
    p=(djedv+je/ae)/Re

    # normalization
    r=ae*Re/sum(je*dv)
    P=r*p
    Je=r*je

    #########################################################################
    # Now the modulation
    #########################################################################

    Aj,Aje=ones(Complex{Float64},n),ones(Complex{Float64},n)
    Ej,Eje=zeros(Complex{Float64},n),zeros(Complex{Float64},n)

    #########################################################################
    # integrate from vlb to vs
    # vs ...... vu ------> vth
    #
    # for the A and E trial solns at vu we have
    # Aj=1 Aje=1
    # Ej=0 Eje=0
    #########################################################################

    for k=su:n-1

        Aj[k+1]=Aj[k]*Xw[k]+Aje[k]*(1-Xw[k])
        Aje[k+1]=exp(-dv/ae)*(Aje[k]*Xe[k] + Aj[k]*(1-Xe[k]))

        Ej[k+1]=Ej[k]*Xw[k]+Eje[k]*(1-Xw[k])
        Eje[k+1]=exp(-dv/ae)*(Eje[k]*Xe[k] + Ej[k]*(1-Xe[k])) + dv*P[k]

    end

    # at threshold require rj=a*Aj=1 and ej=aa*Aj+Ej=0
    a=1/Aj[n];          rj=a*Aj;     rje=a*Aje;
    aa=-Ej[n]/Aj[n];    ej=aa*Aj+Ej; eje=aa*Aje+Eje;

    #########################################################################
    # now integrate from vu to vs
    # vs <----- vu ..... vth
    #########################################################################

    r1thGA=r1thG*0.999
    r1thGB=r1thG*1.111

    jA,jB=zeros(Complex{Float64},n),zeros(Complex{Float64},n)
    jeA,jeB=zeros(Complex{Float64},n),zeros(Complex{Float64},n)

    jA[su]=r1thGA*rj[su]+ej[su]
    jB[su]=r1thGB*rj[su]+ej[su]

    jeA[su]=r1thGA*rje[su]+eje[su]
    jeB[su]=r1thGB*rje[su]+eje[su]

    for k=su:-1:2

        jA[k-1]=jA[k]*Yw[k] + jeA[k]*(1-Yw[k]) - r1thGA*(k==sre)
        jeA[k-1]=exp(dv/ae)*jeA[k]*Ye[k] + (jA[k])*(1-Ye[k]) - dv*P[k]

        jB[k-1]=jB[k]*Yw[k] + jeB[k]*(1-Yw[k]) - r1thGB*(k==sre)
        jeB[k-1]=exp(dv/ae)*jeB[k]*Ye[k] + (jB[k])*(1-Ye[k]) - dv*P[k]

        rj[k-1]=rj[k]*Yw[k]+rje[k]*(1-Yw[k])-(k==sre)
        rje[k-1]=exp(dv/ae)*rje[k]*Ye[k]+(rj[k])*(1-Ye[k])

        ej[k-1]=ej[k]*Yw[k]+eje[k]*(1-Yw[k])
        eje[k-1]=exp(dv/ae)*eje[k]*Ye[k]+ej[k]*(1-Ye[k])-dv*P[k]

    end

    re=(r1thGA*jB[1]-r1thGB*jA[1])/(jB[1]-jA[1])

    j=zeros(Complex{Float64},n)
    je=zeros(Complex{Float64},n)

    # check this below
    j[1:su]=(jA[1:su]*jB[1]-jB[1:su]*jA[1])/(jB[1]-jA[1])
    je[1:su]=(jeA[1:su]*jB[1]-jeB[1:su]*jA[1])/(jB[1]-jA[1])
    j[su+1:end]=rj[su+1:end]*re+ej[su+1:end]
    je[su+1:end]=rje[su+1:end]*re+eje[su+1:end]

    return r,re,v,j,je
end

#########################################################################
#########################################################################
# response to a time course t, R(t)
# uses a second-order RKG
#########################################################################
#########################################################################
function EIFshotRtsim(vth,vre,tau,vT,dT,t,Rt,ae,N)

    dt=t[2]-t[1]
    ns=length(t)

    Xt=Rt*dt
    S=zeros(ns) # spike time series
    v=zeros(ns) # voltage holder

    for n=1:N

        # create the noise vector
        se=ae*(rand(ns).<Xt).*randexp(ns)

        # initialise
        v=zeros(ns)
        v[1]=1.0

        for k=1:ns-1

            v1=v[k]
            F1=(dT*exp((v1-vT)/dT)-v1)/tau

            v2=v[k] + dt*F1/2
            F2=(dT*exp((v2-vT)/dT)-v2)/tau

            v[k+1] = v[k] + dt*F2 + se[k]

            if v[k+1]>vth
                v[k]=50;        # decorative spike
                S[k+1]=S[k+1]+1
                v[k+1]=vre
            end

        end
    end

    r=S/(dt*N)

    return   t,v,r

end
