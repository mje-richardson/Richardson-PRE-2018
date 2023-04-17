##########################################################################
##########################################################################
##########################################################################
#
# This script generates figure 3 of the paper:
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

function MakeFig3()

    # useful constants
    K=1000;
    r2d=360/(2*pi)

    ######################################################################
    # The neuronal parameters
    ######################################################################
    tau=20
    vT=10
    vre=5
    vth=15

    # factors of 3
    dT=0.6
    asA=0.2
    asB=0.6
    asC=1.8

    ######################################################################
    ######################################################################
    # Get the data
    ######################################################################
    ######################################################################

    # some choices for the level of discretization

    #dt=0.01;  nz=500; dv=0.001 # instantaneous
    dt=0.01;  nz=1000; dv=0.001 # about 20secs
    #dt=0.002; nz=1000; dv=0.0001 # about 2mins
    #dt=0.002; nz=2000; dv=0.0001 # about 3mins
    #dt=0.002; nz=5000; dv=0.0001 # about 8mins

    # finds the synaptic rate for a given output rate for cases A, B and C
    r0aim=5.0/1000

    rA,RsA=Rintorout(r0aim,tau,vT,dT,vre,vth,asA,dv)
    rB,RsB=Rintorout(r0aim,tau,vT,dT,vre,vth,asB,dv)
    rC,RsC=Rintorout(r0aim,tau,vT,dT,vre,vth,asC,dv)

    println("rA=$(round(rA*K;digits=4))Hz and RsA=$(round(RsA;digits=4))kHz")
    println("rB=$(round(rB*K;digits=4))Hz and RsB=$(round(RsB;digits=4))kHz")
    println("rC=$(round(rC*K;digits=4))Hz and RsC=$(round(RsC;digits=4))kHz")

    # ThIn method for steady-state and modulation

    nfth=33; lfth=0.1/1000; hfth=1000/1000;
    fth=exp.(collect(range(log(lfth),stop=log(hfth),length=nfth)))  # exp spacing
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
    phA=angle.(rhA); phB=angle.(rhB); phC=angle.(rhC);

    # match the amplitudes so steady response the same
    gfA=0.025*ahA[1]/ahA[1]
    gfB=0.025*ahA[1]/ahB[1]
    gfC=0.025*ahA[1]/ahC[1]

    # get the stimulus: double sampled for the RK4
    t2,RsAt2=thestimulus(dt/2,RsA,gfA)
    t2,RsBt2=thestimulus(dt/2,RsB,gfB)
    t2,RsCt2=thestimulus(dt/2,RsC,gfC);

    t,rAt,RsAt=numint(t2,RsAt2,asA,tau,dT,vT,vre,vth,nz)
    t,rBt,RsBt=numint(t2,RsBt2,asB,tau,dT,vT,vre,vth,nz)
    t,rCt,RsCt=numint(t2,RsCt2,asC,tau,dT,vT,vre,vth,nz)
    T=t[end];

    ######################################################################
    ######################################################################
    # Make the figure
    ######################################################################
    ######################################################################

    # set up the figure
    fig1=figure(figsize=(10,4));

    # figure parameters
    clf()

    lw=0.75
    be=0.125
    ms=3;

    labelfs=8
    panelfs=13
    tickfs=8
    textfs=6

    lw=0.75
    be=0.125

    ######################################################################
    # Plot fig 3a
    ######################################################################

    leA=0.075
    pw=0.85
    phA=0.2

    ax1=PyPlot.axes([leA,0.725,pw,phA])
    plot(t,RsAt,"b-",linewidth=lw)
    plot(t,RsBt,"g-",linewidth=lw)
    plot(t,RsCt,"r-",linewidth=lw)
    axis([0,550,0,2.5])
    ylabel("Stimulus R(t) (kHz)",fontsize=labelfs)

    ax1[:text](20,1.7,"case (i) "   *L"R_s=2.11"*"kHz",fontsize=textfs)
    ax1[:text](20,0.8,"case (ii) "  *L"R_s=0.59"*"kHz",fontsize=textfs)
    ax1[:text](20,0.25,"case (iii) "*L"R_s=0.14"*"kHz",fontsize=textfs)

    ax1[:text](141,1.5,"Sharpening Gaussian impulses",fontsize=textfs)
    ax1[:text](295,1.5,"square pulse ",fontsize=textfs)
    ax1[:text](420,1.5,"100Hz high gamma",fontsize=textfs)
    ax1[:text](500,1.5,"200Hz ripple",fontsize=textfs)

    ax1[:spines]["top"][:set_visible](false)
    ax1[:spines]["right"][:set_visible](false)
    ax1[:tick_params](axis="both",labelsize=tickfs)

    fig1[:text](0.01,0.925, "A", fontsize=panelfs)

    ######################################################################
    # Plot fig 3b
    ######################################################################

    phB=0.2

    ax2=PyPlot.axes([leA,0.45,pw,phB])
    plot(t,K*rAt,"b-",linewidth=lw)
    plot(t,K*rBt,"g-",linewidth=lw)
    plot(t,K*rCt,"r-",linewidth=lw)
    axis([0,550,0,6.5])
    plot([0,T],K*rA*[1,1],"b:",linewidth=lw)
    plot([0,T],K*rB*[1,1],"g:",linewidth=lw)
    plot([0,T],K*rC*[1,1],"r:",linewidth=lw);
    xlabel("Time (ms)",fontsize=labelfs);
    ylabel("Response r(t) (Hz)",fontsize=labelfs);

    plot([18,51.15],2.5*[1,1],"k:",linewidth=lw)
    ax2[:text](55,2.25,"33ms",fontsize=textfs)

    ax2[:text](48,1.0,"(i)",fontsize=textfs)
    ax2[:text](30,1.0,"(ii)",fontsize=textfs)
    ax2[:text](12,1.0,"(iii)",fontsize=textfs)

    ax2[:spines]["top"][:set_visible](false)
    ax2[:spines]["right"][:set_visible](false)
    ax2[:tick_params](axis="both",labelsize=tickfs)

    fig1[:text](0.01,0.65, "B", fontsize=panelfs)

    ######################################################################
    # Plot fig 3c
    ######################################################################

    pwC=0.175
    Cgap=0.225
    phC=0.2

    lC1=leA
    lC2=lC1+Cgap
    lC3=lC2+Cgap
    lC4=lC3+Cgap

    ax3=PyPlot.axes([lC1,0.1,pwC,phC])
    plot(t,K*rAt,"b-",linewidth=lw)
    plot(t,K*rBt,"g-",linewidth=lw)
    plot(t,K*rCt,"r-",linewidth=lw)
    plot(t .-40,K*rAt,"b-",linewidth=lw)
    plot(t .-40,K*rBt,"g-",linewidth=lw)
    plot(t .-40,K*rCt,"r-",linewidth=lw)
    plot(t .-80,K*rAt,"b-",linewidth=lw)
    plot(t .-80,K*rBt,"g-",linewidth=lw)
    plot(t .-80,K*rCt,"r-",linewidth=lw)
    axis([135,180,4.8,7]);
    xlabel("Time (ms)",fontsize=labelfs);
    ylabel("r(t) (Hz)",fontsize=labelfs);
    ax3[:text](141,6.9,"Gaussian impulse responses",fontsize=textfs)

    ax4=PyPlot.axes([lC2,0.1,pwC,phC])
    plot(t,K*rAt,"b-",linewidth=lw)
    plot(t,K*rBt,"g-",linewidth=lw)
    plot(t,K*rCt,"r-",linewidth=lw);
    axis([285,310,4.8,6]);
    xlabel("Time (ms)",fontsize=labelfs);
    ylabel("r(t) (Hz)",fontsize=labelfs);
    ax4[:text](294,5.2,"square-pulse responses",fontsize=textfs)

    ax5=PyPlot.axes([lC3,0.1,pwC,phC])
    plot(t,K*rAt,"b-",linewidth=lw)
    plot(t,K*rBt,"g-",linewidth=lw)
    plot(t,K*rCt,"r-",linewidth=lw);
    axis([415,475,3.5,6.5]);
    xlabel("Time (ms)",fontsize=labelfs);
    ylabel("r(t) (Hz)",fontsize=labelfs);
    ax5[:text](420,6.5,"100Hz high gamma responses",fontsize=textfs)

    ax6=PyPlot.axes([lC4,0.1,pwC,phC])
    plot(t,K*rAt,"b-",linewidth=lw)
    plot(t,K*rBt,"g-",linewidth=lw)
    plot(t,K*rCt,"r-",linewidth=lw);
    axis([490,550,3.5,6.5]);
    xlabel("Time (ms)",fontsize=labelfs);
    ylabel("r(t) (Hz)",fontsize=labelfs);
    ax6[:text](500,6.5,"200Hz ripple responses",fontsize=textfs)

    for ax in (ax3,ax4,ax5,ax6)
        ax[:spines]["top"][:set_visible](false)
        ax[:spines]["right"][:set_visible](false)
        ax[:tick_params](axis="both",labelsize=tickfs)
    end

    fig1[:text](0.01,0.3, "C", fontsize=panelfs)

    ######################################################################
    # Optional save
    ######################################################################

    if false
        savefig("fig3.pdf")
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
    RsL=range(Rsl,stop=Rsh,length=Rsn)
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

#########################################################################
#########################################################################
# Solves for an arbitrary waveform using an RK4 algorithm
#########################################################################
#########################################################################
function numint(t2,Rst2,as,tau,dT,vT,vre,vth,nz)

    z,gz,vz,dvdz=ztransformB(tau,dT,vT,vth,nz)
    dz=z[2]-z[1]

    t=t2[1:2:end]
    Rst=Rst2[1:2:end]
    RstH=Rst2[2:2:end] # half step needed for the RK4

    dt=t[2]-t[1]
    nt=length(t)
    dz=z[2]-z[1]
    r=zeros(nt)

    # initial conditions and initialisation
    muv=1
    sigv=0.1
    Pz=exp.(-(vz .-muv).^2/(2*sigv^2)).*dvdz
    Pz=Pz./sum(dz*Pz)

    # smooth the insertion at vre with a Gaussian (in voltage space)
    # this has little effect on result but improves convergence
    muvre=vre
    sigvre=0.1
    hzre=exp.(-(vz .-muvre).^2/(2*sigvre^2)).*dvdz
    hzre=hzre/sum(dz*hzre)

    F1=exp.((vT .-vz)/as)
    G1=1.0 ./F1

    # loop over time
    for k=1:nt-1

        Jz,Kz,dJdz=getarray(z,F1,G1,gz,dz,Pz,Rst[k])
        k1=-dt*dJdz + dt*hzre*Jz[end]

        Jz,Kz,dJdz=getarray(z,F1,G1,gz,dz,Pz+k1/2,RstH[k])
        k2=-dt*dJdz + dt*hzre*Jz[end]

        Jz,Kz,dJdz=getarray(z,F1,G1,gz,dz,Pz+k2/2,RstH[k])
        k3=-dt*dJdz + dt*hzre*Jz[end]

        Jz,Kz,dJdz=getarray(z,F1,G1,gz,dz,Pz+k3,Rst[k+1])
        k4=-dt*dJdz + dt*hzre*Jz[end]

        Pz=Pz+k1/6+k2/3+k3/3+k4/6
        r[k]=Jz[end]

    end

    r[end]=r[end-1]

    return t,r,Rst

end

#########################################################################
#########################################################################
# get arrays
#########################################################################
#########################################################################
function getarray(z,F1,G1,gz,dz,Pz,Ret)

    Kz=dz*Ret*F1.*(cumsum(Pz.*G1)-Pz[1]*G1/2-Pz.*G1/2)
    Jz=Kz + gz.*Pz
    dJdz=[(Jz[2]-Jz[1]);(Jz[3:end]-Jz[1:end-2])/2;Jz[end]-Jz[end-1]]/dz;

    return Jz,Kz,dJdz

end

#########################################################################
#########################################################################
# the stimulus which needs to be multiplied by a gain factor
#########################################################################
#########################################################################
function thestimulus(dt,Rs,gf)

    t=collect(0:dt:570)
    nt=length(t)

    x=zeros(nt)

    # number of extra spikes for the Gaussian
    ng=10

    # Gaussian pulse 1
    tgauss=140; FWHM=4; siggau=FWHM/sqrt(8*log(2))
    x=x + ng*exp.(-(t .-tgauss).^2/(2*siggau^2))/sqrt(2*pi*siggau^2)

    # Gaussian pulse 2
    tgauss=190; FWHM=2; siggau=FWHM/sqrt(8*log(2))
    x=x + ng*exp.(-(t .-tgauss).^2/(2*siggau^2))/sqrt(2*pi*siggau^2)

    # Gaussian pulse 3
    tgauss=240; FWHM=1; siggau=FWHM/sqrt(8*log(2))
    x=x + ng*exp.(-(t .-tgauss).^2/(2*siggau^2))/sqrt(2*pi*siggau^2)

    # step
    x=x + 2*(t.>290).*(t.<350)

    # chirp 1 - fast gamma
    tch1=445; sch1=10; fch1=100/1000;  ampch1=5;
    x=x+ampch1*cos.(2*pi*fch1*t).*exp.(-(t .-tch1).^2/(2*sch1^2))

    # chirp 2 - ripples
    tch2=520; sch2=10; fch2=200/1000; ampch2=7;
    x=x+ampch2*cos.(2*pi*fch2*t).*exp.(-(t .-tch2).^2/(2*sch2^2))

    # Multiply by prefactor and add steady-state
    Rst=Rs .+ gf*x
    Rst[findall(t.<5)] .=0

    return t,Rst

end

#########################################################################
#########################################################################
# z transform to get v for a constant array over z
#########################################################################
#########################################################################
function ztransformB(tau,dT,vT,vth,nz)

    # get the bounds for v
    vl,vu=MyfRoots.((0.0,vT+dT),vT,dT,10)
    vh=vth

    # add ve to v and sample at dv for the spline
    v=collect(range(vl,stop=vh,length=nz))
    psi=exp.((v .-vT)/dT)
    zv=(psi .+v/vT)./(psi .+v/vT .+1)    # z on the v grid
    z=collect(range(zv[1],stop=zv[end],length=nz))
    dz=z[2]-z[1]

    # do one step using Halley method for finding roots
    vz=vT*ones(nz)
    for k=1:50
        psi=exp.((vz .-vT)/dT)
        H=psi .+vz/vT .-z./(1 .-z)
        dHdv=psi/dT .+1/vT
        d2Hdv2=psi/dT^2
        vz=vz-(2*H.*dHdv)./(2*dHdv.^2-H.*d2Hdv2)
    end

    if isnan(vz[end])
        vz[end]=vth
    end

    if sum(isnan.(vz))>0
        println("ERROR! there are some values of z that are NaN")
    end

    # now the other quantities we need on the z grid
    psi=exp.((vz .-vT)/dT)
    dzdv=(psi/dT .+1/vT)./(psi .+vz/vT .+1).^2
    dvdz=1.0 ./dzdv
    gz=dzdv.*(dT*psi -vz)/tau

    return z,gz,vz,dvdz

end
