HEADER 'Compton Scattering';
!!  Very high energy Compton scattering including nonlinear scattering.
ALLOCATE MP=1000000;
SET   photon=1, electron=2, positron=3,
   mm=1e-3, micron=1e-6, nm=1e-9, mu0=4*Pi*1e-7, psec=1e-12*Cvel,
   nsec=1e-9*Cvel;
! define variables for electron beam
SET ee=40D6,  gamma=ee/Emass,  an=5D5,  sigz=0.01*psec, 
   betax=0.8*mm, betay=0.8*mm, emitx=2D-8, emity=2D-8,
   sigx=Sqrt(emitx*betax), sigy=Sqrt(emity*betay),sige = 0.25,
   ntcut=3.0;
! define variables for laser
SET  laserwl=0.8*micron, lambar=laserwl/(2*Pi), omegal=Hbarc/lambar, 
   rlx=100*micron, rly=100*micron, sigt=1*psec,   !! Rayleigh length and pulse length
   pulseE=0.8,        !!  pulse energy in Joule
   powerd=pulseE/[sigt*laserwl*Sqrt(rlx*rly)/(4*Cvel)],
   !aNaught = 0.1,
   !powerd=(10000)*(aNaught^2)*1.37*10^18/(0.8)^2,
   xisq=powerd*mu0*Cvel*(lambar/Emass)^2,   xi=Sqrt(xisq),
   eta=omegal*ee/Emass^2, lambda=4*eta,
   angle=0.0 ;
SET MsgLevel=1;
BEAM  RIGHT, KIND=electron, NP=500000, AN=an, E0=ee,
   TXYS=(0,0,0,0),  GCUTT=ntcut,
   BETA=(betax,betay), EMIT=(emitx,emity), SIGT=sigz, SPIN=(0,0,-1), SIGE=sige;
LASER LEFT, WAVEL=laserwl, POWERD=powerd,
      TXYS=(0,0,0,0),
      E3=(0,-Sin(angle),-Cos(angle)), E1=(1,0,0), 
      RAYLEIGH=(rlx,rly), SIGT=sigt, GCUTT=ntcut, STOKES=(0,1,0) ;
LASERQED  COMPTON, LINEARPOL, NPH=0, XIMAX=1.1*xi, LAMBDAMAX=1.1*lambda ,
       ENHANCE=1, PMAX=0.5 ;
SET MsgLevel=0;  FLAG OFF ECHO;
SET Smesh=sigt/3;
SET emax=1.001*ee, wmax=emax;
SET  it=0;
PRINT CPUTIME;
PUSH  Time=(-ntcut*(sigt+sigz),ntcut*(sigt+sigz),1000);
      IF Mod(it,50)=0;
        PRINT it, FORMAT=(F6.0,'-th time step'); PRINT STAT, SHORT;
      ENDIF;
     SET it=it+1;
ENDPUSH;
PRINT CPUTIME;
!  Pull all particles to the back to the focal point
DRIFT S=0.7;
WRITE BEAM, KIND=photon, FILE='pt4mmJet.dat';
PRINT STAT;
!PLOT  HIST, RIGHT, KIND=electron, H=En/1D9, HSCALE=(0,emax/1e9,50),
!
!
!        TITLE='Right-Going Electron Energy Spectrum;',
!        HTITLE='E0e1 (GeV); X X      ;';
!PLOT  HIST, KIND=photon, H=En/1D9, HSCALE=(0,emax/1e9,50),
!        TITLE='All Photon Energy Spectrum;',
!        HTITLE='E0G1 (GeV); XGX      ;'  ;
!PLOT  SCAT, KIND=photon, RIGHT,
!        H=Y, V=X,
!        HSCALE=(-1*micron,1*micron), VSCALE=(-1*micron,1*micron),
!        TITLE='Photon Energy vs. Helicity;',
!        HTITLE='E0G1 (GeV); XGX      ;',
!        VTITLE='X021; X X;';
STOP;
END;
