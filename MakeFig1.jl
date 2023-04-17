##########################################################################
##########################################################################
##########################################################################
#
# This script generates figure 1 of the paper:
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
using Distributions

function MakeFig1()

    ######################################################################
    # The neuronal parameters
    ######################################################################

    tau=20
    vT,dT=10,0.6
    vre,vth=5,30

    asA=0.2
    asB=0.6
    asC=1.8

    Random.seed!(2);

    #####################################################################
    ######################################################################
    # Generate the data for each panel
    ######################################################################
    ######################################################################

    ######################################################################
    # Generate data for figure 1A
    ######################################################################

    vl,vu=fRoots.((0.0,vT+dT),vT,dT,10)

    # voltage for the trapping potential
    dv=0.05
    v=collect(-12:dv:vT+6*dT)
    n=length(v)

    f=(dT*exp.((v .-vT)/dT) -v)/tau
    phi=-cumsum(f*dv)
    sl=Int(ceil(n*(vl-v[1])/(v[end]-v[1])))
    su=Int(ceil(n*(vu-v[1])/(v[end]-v[1])))

    phi=phi .-phi[sl]
    phil=0
    phiu=phi[su];

    ######################################################################
    # Generate data for figure 1B (and 1C)
    ######################################################################

    # parameters for simulations
    T=5000
    dt=0.1;

    # numerical parameters for Fig C
    dv=0.01;

    # This is mainly for panel C but RsA etc are needed for Fig B too
    Rsl=0.01; Rsh=3; Rsn=100;
    RsL=collect(range(Rsl,stop=Rsh,length=Rsn))

    rAL=zeros(Rsn)
    rBL=zeros(Rsn)
    rCL=zeros(Rsn)

    for k=1:Rsn
        rAL[k],~,~,~=EshotEIFThIn0(dv,vth,vre,RsL[k],asA,tau,vT,dT)
        rBL[k],~,~,~=EshotEIFThIn0(dv,vth,vre,RsL[k],asB,tau,vT,dT)
        rCL[k],~,~,~=EshotEIFThIn0(dv,vth,vre,RsL[k],asC,tau,vT,dT)
    end

    # the firing rate for the example
    r0aim=5.0/1000.0;

    sA1=maximum(findall(rAL.<r0aim)); sA2=minimum(findall(rAL.>r0aim));
    R1=RsL[sA1]; R2=RsL[sA2]; r1=rAL[sA1]; r2=rAL[sA2];
    RsA=R1+(r0aim-r1)*(R2-R1)/(r2-r1);
    rA=r1+(RsA-R1)*(r2-r1)/(R2-R1);

    sB1=maximum(findall(rBL.<r0aim)); sB2=minimum(findall(rBL.>r0aim));
    R1=RsL[sB1]; R2=RsL[sB2]; r1=rBL[sB1]; r2=rBL[sB2];
    RsB=R1+(r0aim-r1)*(R2-R1)/(r2-r1);
    rB=r1+(RsB-R1)*(r2-r1)/(R2-R1);

    sC1=maximum(findall(rCL.<r0aim)); sC2=minimum(findall(rCL.>r0aim));
    R1=RsL[sC1]; R2=RsL[sC2]; r1=rCL[sC1]; r2=rCL[sC2];
    RsC=R1+(r0aim-r1)*(R2-R1)/(r2-r1);
    rC=r1+(RsC-R1)*(r2-r1)/(R2-R1);

    t,vA,rAsim=EIFsim(vth,vre,tau,vT,dT,RsA,asA,T,dt)
    t,vB,rBsim=EIFsim(vth,vre,tau,vT,dT,RsB,asB,T,dt)
    t,vC,rCsim=EIFsim(vth,vre,tau,vT,dT,RsC,asC,T,dt);

    ######################################################################
    # Generate data for Figure 1D data
    ######################################################################

    rA,vD,PA,JsA=EshotEIFThIn0(dv,vth,vre,RsA,asA,tau,vT,dT)
    rB,vD,PB,JsB=EshotEIFThIn0(dv,vth,vre,RsB,asB,tau,vT,dT)
    rC,vD,PC,JsC=EshotEIFThIn0(dv,vth,vre,RsC,asC,tau,vT,dT);

    # now the integrand for Js

    # these are the full values
    IA=PA.*exp.(vD/asA)
    IB=PB.*exp.(vD/asB)
    IC=PC.*exp.(vD/asC);

    ######################################################################
    ######################################################################
    # Now plot the figures
    ######################################################################
    ######################################################################

    # set up the figure
    fig1=figure(figsize=(10,3.5))

    # figure parameters
    clf()

    K=1000.0 # useful for plots

    lw=0.75
    be=0.125
    ms=3;
    labelfs=8
    panelfs=13
    tickfs=8
    textfs=6

    ######################################################################
    # Plot Figure 1A
    ######################################################################

    pwA=0.2
    phA=0.125
    leA=0.1

    ax1=PyPlot.axes([leA,0.75,pwA,phA])
    plot(v,f,"k-",linewidth=lw)
    plot([v[1],20],[0,0],"k:",linewidth=lw)
    plot(vl,0,"ko",markersize=ms)
    plot(vu,0,"ko",markersize=ms)
    axis([0,20,minimum(f)*1.1,-minimum(f)*1.1])
    ylabel("Force f(v)",fontsize=labelfs)
    #title("v_T=$vT dT=$dT mV");

    ax2=PyPlot.axes([leA,0.55,pwA,phA])
    plot(v,phi,"k",linewidth=lw)
    plot(vl,phil,"ko",markersize=ms)
    plot(vu,phiu,"ko",markersize=ms)
    axis([0,20,phil-1,phiu+1]);
    ylabel("Potential well",fontsize=labelfs)


    ax1[:text](12.5,0.1,L"v_\mathrm{u}",fontsize=textfs)
    ax1[:text](0.5,0.1,L"v_\mathrm{s}",fontsize=textfs)

    for ax in (ax1,ax2)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig1[:text](0.02,0.925, "A", fontsize=panelfs)

    ######################################################################
    # Plot Figure 1B
    ######################################################################

    pwB=0.2
    phB=0.2
    leB=0.4

    ax3=PyPlot.axes([leB,0.7,pwB,phB])
    plot(t,vA,"b",linewidth=lw)
    ylabel("Voltage (mV)",fontsize=labelfs)
    axis([0,1000,0,20]);
    #title("rA=$(round(rAsim*K),2)Hz"])

    ax4=PyPlot.axes([leB,0.4125,pwB,phB])
    plot(t,vB,"g",linewidth=lw)
    ylabel("Voltage (mV)",fontsize=labelfs)
    axis([0,1000,0,20]);
    #title("rB=$(round(rBsim*K),2)Hz"])

    ax5=PyPlot.axes([leB,be,pwB,phB])
    plot(t,vC,"r",linewidth=lw)
    ylabel("Voltage (mV)",fontsize=labelfs)
    xlabel("Time (ms)",fontsize=labelfs)
    axis([0,1000,0,20]);
    #title("rC=$(round(rCsim*K),2)Hz"])

    ax3[:text](100,16,"case (i)",fontsize=textfs)
    ax3[:text](100,13,L"a_s<\delta_\mathrm{T}",fontsize=textfs)

    ax4[:text](300,16,"case (ii)",fontsize=textfs)
    ax4[:text](300,13,L"a_s=\delta_\mathrm{T}",fontsize=textfs)

    ax5[:text](200,16,"case (iii)",fontsize=textfs)
    ax5[:text](200,13,L"a_s>\delta_\mathrm{T}",fontsize=textfs)


    for ax in (ax3,ax4,ax5)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig1[:text](0.325,0.925, "B", fontsize=panelfs)

    ######################################################################
    # Plot Figure 1C
    ######################################################################

    pwC=0.27
    phC=0.7
    leC=0.69

    ax6=PyPlot.axes([leC,be,pwC,phC])
    plot(RsL,K*rAL,"b-",linewidth=lw)
    plot(RsL,K*rBL,"g-",linewidth=lw)
    plot(RsL,K*rCL,"r-",linewidth=lw)
    plot(RsA,K*rA,"b*",markersize=ms)
    plot(RsB,K*rB,"g*",markersize=ms)
    plot(RsC,K*rC,"r*",markersize=ms)
    plot([Rsl,Rsh],K*r0aim*[1,1],"k:",linewidth=lw)
    axis([0,maximum(RsL),0,20]);
    xlabel("Synaptic rate R (kHz)",fontsize=labelfs);
    ylabel("Firing rate r (Hz)",fontsize=labelfs);

    ax6[:text](2.2,8.5,L"R_s=2.1\mathrm{kHz}",fontsize=textfs,rotation=65)
    ax6[:text](0.65,9,L"R_s=0.59\mathrm{kHz}",fontsize=textfs,rotation=70)
    ax6[:text](0.2,9,L"R_s=0.14\mathrm{kHz}",fontsize=textfs,rotation=75)

    ax6[:text](2.5,21.5,"case (i)",fontsize=textfs)
    ax6[:text](0.8,21.5,"case (ii)",fontsize=textfs)
    ax6[:text](0.2,21.5,"case (iii)",fontsize=textfs)

    ax6[:text](2.5,20.5,L"a_s<\delta_\mathrm{T}",fontsize=textfs)
    ax6[:text](0.8,20.5,L"a_s=\delta_\mathrm{T}",fontsize=textfs)
    ax6[:text](0.2,20.5,L"a_s>\delta_\mathrm{T}",fontsize=textfs)

    ax6[:spines]["top"][:set_visible](false)
    ax6[:spines]["right"][:set_visible](false)
    ax6[:tick_params](axis="both",labelsize=tickfs)

    fig1[:text](0.625,0.925, "C", fontsize=panelfs)

    ######################################################################
    # Plot Figure 1C
    ######################################################################

    pwD=pwA
    phD=phA
    leA=0.1

    ax7=PyPlot.axes([leA,0.325,pwD,phD])
    plot(vD,PA,"b-",linewidth=lw)
    plot(vD,PB,"g-",linewidth=lw)
    plot(vD,PC,"r-",linewidth=lw)
    ax7[:set_yticks]([0,0.1,0.2,0.3])
    axis([0,20,0,0.31])
    ylabel("Density P(v)",fontsize=labelfs)

    ax7[:text](11,0.05,L"P\propto \exp(-(v-v_\mathrm{T})/\delta_\mathrm{T})",fontsize=textfs)

    s15=minimum(findall(vD.>15))

    ax8=PyPlot.axes([leA,be,pwD,phD])
    plot(vD,IA/IA[s15],"b",linewidth=lw)
    plot(vD,IB/maximum(IB),"g",linewidth=lw)
    plot(vD,IC/maximum(IC),"r",linewidth=lw)
    axis([0,20,0,1.1])
    ylabel(L"{\sim}P(v)e^{v/a_s}",fontsize=labelfs)
    xlabel("Voltage (mV)",fontsize=labelfs)

    ax8[:text](15,0.5,"case (i)",fontsize=textfs)
    ax8[:text](6.5,0.3,"case (ii)",fontsize=textfs)
    ax8[:text](1.5,0.8,"case (iii)",fontsize=textfs)

    for ax in (ax7,ax8)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig1[:text](0.02,0.45, "D", fontsize=panelfs)

    ######################################################################
    # Optional save
    ######################################################################

    if false
        savefig("fig1.pdf")
    end

end

#########################################################################
#########################################################################
#########################################################################
# Subfunctions called by main function
#########################################################################
#########################################################################
#########################################################################

#########################################################################
#########################################################################
# Halley's method for lower and upper fixed points
#########################################################################
#########################################################################
function fRoots(v1,vT,dT,n)

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
# ThIn for steady-state EIF with one shot-noise process
#########################################################################
#########################################################################
function EshotEIFThIn0(dv,vth,vre,Re,ae,tau,vT,dT)

    # get the fixed points
    vl,vu=fRoots.((0.0,vT+dT),vT,dT,10)

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
    dFdvu=(exp((vu-vT)/dT)-1)/tau
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
# ThIn for steady-state EIF with one shot-noise process
#########################################################################
#########################################################################
function EIFsim(vth,vre,tau,vT,dT,Rs,as,T,dt)

    # useful constants
    dtbytau=dt/tau
    Rsdt=Rs*dt

    # prepare the arrays and initialise
    n=Int(ceil(T/dt))
    t=dt*collect(1:n)
    v=zeros(n)
    v[1]=vre
    c=0.0 # spike counter

    # generate the noise vector
    sn=as*rand(Bernoulli(Rsdt),n).*randexp(n,1)

    for k=1:n-1
        v[k+1]=v[k]+dt*(dT*exp((v[k]-vT)/dT)-v[k])/tau + sn[k]

        if v[k+1]>vth
            v[k]=vth
            v[k+1]=vre
            c=c+1.0
        end
    end

    r=c/T;

    return t,v,r

end
