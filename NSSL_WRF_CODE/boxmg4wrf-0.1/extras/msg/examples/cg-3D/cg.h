      integer MaxNumProc, maxndimx, maxndimy, maxndimz
      integer MaxStencilSize, MaxHalo
      integer k0, kw, knw, ksw, ks
      integer kp,kpw,kps,kb,kpsw,kpnw,kbw,kbnw,kbn,kbne
      integer kbe,kbse,kbs,kbsw
      integer nk,MSGptrn1,MSGptrn2
      

      parameter (MSGptrn1=1,MSGptrn2=2)
      parameter (k0=1,kw=2,ks=3,ksw=4,knw=5)

      parameter (kp=1,kpw=2,kps=3,kb=4,kpsw=5,kpnw=6,kbw=7,kbnw=8,
     &           kbn=9,kbne=10,kbe=11,kbse=12,kbs=13,kbsw=14)
      parameter (MaxNumProc=16)
      parameter (maxndimx=138,maxndimy=138,maxndimz=138)
      parameter (MaxStencilSize=14)
      parameter (MaxHalo=2)

