! ===============
  module Type_Def
! ===============
!
! Cellc
! -----
     type Cellc
        SEQUENCE
        integer        :: nl(8)       ! その点の全体番号
     end type Cellc

  end module Type_Def
!
!
! ===========
  module com0
! ===========
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER(NF=20,NF2=19,Nunit=21)
      PARAMETER(EPS_MAX=1.D100,EPS_MIN=-1.D100)

      PARAMETER(NPMAX=27)
!------CONSTANT
!
        real(8) PAI
!
        integer NON,NOE,NT,NGP,nlmax,non3,iterst
! 
        real(8) YOUNG,ANU,gratio
!--------ITE
!
        integer ITE_PSOR,ITE_PSOR_DISP
!
        real(8) OMG_PSOR
!--------MATRIX
!
        integer NX,NY,NZ
        real(8) DNDXAI(8,NPMAX),DNDETA(8,NPMAX),DNDZTA(8,NPMAX)
!
        real(8) SF(8,NPMAX),GW(NPMAX),SGP(6,8,3)
!--------FILE
!
!
!--------BOUNDARY
!
        integer LP(6,4)
        integer nnlmax
        real(8),allocatable:: px(:),py(:),pz(:)
         integer NB(6,3)

  end module com0
! ===========
  module com1
! ===========
      use Type_def
!        real(4),allocatable:: CNOD(:,:),DNOD(:,:)
        real(4),allocatable:: CNODX(:),CNODY(:),CNODZ(:),CNODX1(:),CNODY1(:),CNODZ1(:)
        real(8),allocatable:: A(:),DENS(:),pxyz(:)
!        real(8),allocatable:: A(:),AOLD(:),resA(:)
        type(Cellc) ,allocatable:: NEC(:)
     integer,allocatable:: NODFAN(:,:,:),NECFAN(:,:,:),NECNCO(:,:)
        real(8) NEN(8)
        real(8),allocatable:: DIS(:,:),EP(:,:),&
&                    SIG(:,:),FF(:),&
&                    S_SIG(:),earea(:)
!
        real(8),allocatable:: AD(:),AU(:,:)
!
        integer,allocatable:: IU(:,:),NZU(:)
!
        real(8),allocatable:: STR(:,:)
        integer NLD
        integer,allocatable:: NEL(:),NPL(:),NFL(:)
!-------FIX
       integer,allocatable:: NFIX(:),ND(:)
       integer NTFIX
!.......Material
          integer noem
         integer,allocatable:: NECM(:)

  end module com1
!////////////////////////////////////////////////////////////
!::::::::::::::::
 module com2

          integer NTSEL,NTSNEG,NTPABS
         integer,allocatable::NSEL(:),MSNEG(:),MEREF(:),NETOS(:),METOS(:,:)
        real(8),allocatable::SDDENS(:),DENSDEV(:)
        real(8) DENSMIN

!      XD=10.D0
!      YD=10.D0
!      ZD=10.D0

!      NX=int(xd)+1
!      NY=int(yd)+1
!      NZ=int(zd)+1

!      allocate( NODFAN(nx,ny,nz),NECFAN(nx,ny,nz))
!      allocate( CNOD(nx*ny*nz,3),NEC(nx*ny*nz,8))

 end module com2
!:::::::::::::::

! ===========
  module com3
! ===========
!qqqq
        real(8)a
!        real(8)pai
        real(8)ori(3),tra(3),q(4,3),v(3,2),surf(3),pnl(3),pnr(3),p3d(3),veclr(3),uveclr(3),cvecl(3),cvecr(3)
!	psl():perpendicular slope length
	real(8)uel(3),uer(3)
        real(8)rub,rhb,rdb12,rdb34,rcb1,rcb2,rad,radt,ul,ul2,usq,ucub,cc1,cc2,ans,ans2,ata,ata2,ccz,cx,cy,cz,ptrox,ptroy,ptroz,pp1,pp2,pp3,pp2o
!	rub=sul/sbl,rhb=hts/sbl,rdb12=d12/sbl,rdb34=d34/sbl,rcb1=cc1/sbl,rcb2=cc2/sbl
!	rub=1+rdb12(medial: + )+rdb34(lateral: + )
	integer nls
        integer,allocatable::melem(:,:,:),nhelem(:,:,:),mshader2(:,:,:),mround(:,:),nshade(:,:,:)
        integer,allocatable::mshader(:,:,:),numnls(:)
        real(8),allocatable::aelem(:,:,:),delem(:,:,:),helem(:,:,:),shadec(:,:,:),shader(:,:,:),shade(:,:,:)
        real(8),allocatable::pelem(:,:),den(:),densnls(:)
!                                                         density of near loded surface
  end module com3
!============

        program main
        use com0
        use com1
        use com2
        use com3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	real(8)trpzl(4,3),trpzu(4,3),psl(4),trpzl2(4,3),trpzu2(4,3),psl2(4),vec(3),p1(5),p2(5)
!
	write(*,*)' mkg.f90 '
       DENSMIN=1.D-30
        pai = DACOS(-1.D0)
        pai2=pai/2
        anso=120.d0!original neck-shaft angle
        ans=180.d0-anso
        ans=pai*ans/180.d0
        ans2=120.d0
        ans2=180.d0-ans2
        ans2=pai*ans2/180.d0
	cosans=dcos(ans)

        rad=25.d0
	radt=rad*0.7d0
	radm=13.0d0/25.d0*rad!radius of femoral head
!        sa=30.d0!loaded surface angle
	sah=60.d0!loaded globe(femoral head) opening angle.(half)
	als=180.d0-sah!angle of loaded surface
	sa2=20.d0
	sa=180.0*dasin(radm*dsin(sah*pai/180.d0)/rad)/pai!when sah=60.0degree it must be 25.65degree(0.447radian)
	write(*,*)'sa= ',sa
!	if(ans0-sa.lt.90.d0)stop 'sa is too big '
!********	Trochanter major variable	***************
!q2
	rhb=12.0d0/9.d0!Ratio of ht by sides:hts/(0.5*(sbl+sul))
	rub=10.0d0/8.0d0!Ratio of sul by sbl:sul/sbl
        rdb34=1.d0/8.0d0!Ratio of length between left end of sul and sbl by sbl:d12/sbl
!	rdb34=0.0d0
!	rub+rdb12-rdb34=1
!	rdb34=rub+rdb12-1.d0!Ratio of length between right end of sul and sbl by sbl:d34/sbl
        rdb12=rdb34-rub+1.d0
!        rdb12=rdb34+1.0d0-rub
	rcb1=4.0d0/8.d0!Ratio of length between lower end of sbl and side_center(dcs) by sbl:cc1/sbl 
	rcb2=5.d0/10.d0!Ratio of length between left_side and hight_center(cc2) by sbl:cc2/sbl
	httr=0.2*rad!hight of trochanter major
	httr2=0.1*rad!hight2 of trochanter major
	ra12=8.0d0/8.0d0!ratio of area2/area1
	ra21=1.0d0/ra12 !ratio of area1/area2
!	trpzl(4,3),trpzu(4,3):lower trapezoid,upper trapezoid:1=UpperLeft,2=LowerLeft,3=LowerRight,4=UpperRight
	as12=90.d0!:slope angle of left side
	as34=90.d0!:slope angle of right side
	ab41=90.d0!:slope angle of upper base
	ab23=90.d0!:slope angle of lower base
	p1(1)=0.d0
	p1(2)=0.d0
	p1(3)=0.d0
	p1(4)=0.d0
	p1(5)=0.d0

	p2(1)=0.d0
	p2(2)=0.d0
	p2(3)=0.d0
	p2(4)=0.d0
	p2(5)=1.0d0


        csa=90.d0-sa!complementary angle of surface angle
	csa2=90.d0-sa2
        sa=sa*pai/180.d0
        sah=sah*pai/180.d0
        als=als*pai/180.d0
        csa=csa*pai/180.d0
        sa2=sa2*pai/180.d0
        csa2=csa2*pai/180.d0
	as12=as12*pai/180.d0
	as34=as34*pai/180.d0
	ab41=ab41*pai/180.d0
	ab23=ab23*pai/180.d0

        ata=0.d0
	ata2=180.d0
        ata=pai*ata/180.d0 
	ata2=pai*ata2/180.d0!great trochanteric ante_torsion angle corresponding to neck ata
	pp1=1.0d0!force1/area
	pp2=0.05d0!force2/area
	pp3=1.0d0
	pp2o=pp2

        nnx=int(rad)*10
        ddeg=0.5*pai/nnx

        tall=rad*1.0d0
!	tall=0.d0
      YOUNG=100.D0
      ANU=0.3D0

      gratio=10.d0/rad!in standerdizing for paraview , radius becomes 10.0


	room=rad/20.d0
	if(room.gt.2.d0 .or. room.le.1.d0)room=2.d0
!	room=3.0d0
	iroom=int(room)

	hc=rad*dcos(sa)+radm!hight of head center
	toph=hc*dcos(ans)
	toph=hc+radm!hight of toph
	if(rad.gt.toph)toph=rad
	cy=hc*dsin(ans)+radm
	cy=cy+room
	if(rad.gt.cy)cy=rad
        cy=cy+room
	cx=hc*dsin(ans)*dsin(ata)+radm
        if(rad.gt.cx)cx=rad
	cx=cx+room
        cz=tall
	cx=real(int(cx))
	cy=real(int(cy))

        write(*,*)'cx1,cy1,cz1= '
        write(*,*)cx,cy,cz

        irad1=int((cx)*2)
        irad2=int((cy)*2)
        itall=int(tall+toph+room)+1
        nx=irad1
        ny=irad2
        nz=itall

	xd=real(nx)
	yd=real(ny)
	zd=real(nz)

	cx=real(nx/2)
	cy=real(ny/2)
	cz=tall

        write(*,*)'nx,ny,nz= ',nx,ny,nz

        write(*,*)'cx,cy,cz= '
        write(*,*)cx,cy,cz


      allocate(shadec(nx,ny,nz),shader(nx,ny,nz),shade(nx,ny,nz),aelem(nx,ny,nz),mround(nx,ny))
      allocate(mshader(nx,ny,nz),nshade(nx,ny,nz))
      allocate(NODFAN((nx+1),(ny+1),(nz+1)),NECFAN(nx,ny,nz),NECNCO(nx*ny*nz,3))
      allocate( CNODX((nx+1)*(ny+1)*(nz+1)),CNODY((nx+1)*(ny+1)*(nz+1)))
      allocate( CNODZ((nx+1)*(ny+1)*(nz+1)),NEC(nx*ny*nz),DENS(nx*ny*nz),pelem(nx*ny*nz,3),earea(nx*ny*nz))
      allocate( CNODX1((nx+1)*(ny+1)*(nz+1)),CNODY1((nx+1)*(ny+1)*(nz+1)))
      allocate( CNODZ1((nx+1)*(ny+1)*(nz+1)))
      allocate(den((NX+1)*(NY+1)*(NZ+1)),px(nx*ny*nz),py(nx*ny*nz),pz(nx*ny*nz),NEL(nx*ny*nz),pxyz(nx*ny*nz))

        pelem=0.d0
        nshade=0

      NUM=0
      DO k=1,NZ+1
      DO j=1,NY+1
      DO i=1,NX+1
      NUM=NUM+1
      NODFAN(i,j,k)=NUM
      ENDDO
      ENDDO
      ENDDO
      NON=NUM
          non3=non*3
!        allocate (a(non3),AOLD(non3),resA(non3),NZU(non3),FF(non3))
!        allocate (a(non3),NZU(non3),FF(non3))
!        allocate (DIS(non,3),cnod(non,3),DNOD(non,3))
        allocate (DIS(non,3))
      STEPX=XD/NX
      STEPY=YD/NY
      STEPZ=ZD/NZ
      DO k=1,NZ+1
      DO j=1,NY+1
      DO i=1,NX+1
      CNODX(NODFAN(i,j,k))=STEPX*(i-1)
      CNODY(NODFAN(i,j,k))=STEPY*(j-1)
      CNODZ(NODFAN(i,j,k))=STEPZ*(k-1)
      CNODX1(NODFAN(i,j,k))=CNODX(NODFAN(i,j,k))*gratio
      CNODY1(NODFAN(i,j,k))=CNODY(NODFAN(i,j,k))*gratio
      CNODZ1(NODFAN(i,j,k))=CNODZ(NODFAN(i,j,k))*gratio
      ENDDO
      ENDDO
      ENDDO

      NUM=0
      DO k=1,NZ
      DO j=1,NY
      DO i=1,NX
      NUM=NUM+1
      NECFAN(i,j,k)=NUM
       NECNCO(NUM,1)=i
       NECNCO(NUM,2)=j
       NECNCO(NUM,3)=k
!
!−−−−計算結果出力
!
!      CALL OUTPUT
!
      NEC(NUM)%NL(1)=NODFAN(i,j,k)
      NEC(NUM)%NL(2)=NODFAN(i+1,j,k)
      NEC(NUM)%NL(3)=NODFAN(i+1,j+1,k)
      NEC(NUM)%NL(4)=NODFAN(i,j+1,k)
      NEC(NUM)%NL(5)=NODFAN(i,j,k+1)
      NEC(NUM)%NL(6)=NODFAN(i+1,j,k+1)
      NEC(NUM)%NL(7)=NODFAN(i+1,j+1,k+1)
      NEC(NUM)%NL(8)=NODFAN(i,j+1,k+1)
      ENDDO
      ENDDO
      ENDDO
       NOE=NUM
        write(*,*)'NOE,NON= ',noe,non

        do k=1,nz
        do j=1,ny
        do i=1,nx
        dens(necfan(i,j,k))=densmin
        enddo
        enddo
        enddo

        shade=0.d0
        shadec=0.d0
        shader=0.d0
        mshade=0
        mshader=0
	aelem=0
        allocate (melem(nx,ny,nz),delem(nx,ny,nz),helem(nx,ny,nz),mshader2(nx,ny,nz))
        mshader2=0
! ::::: quarterglobe ::::
	vol=0.d0
	ndev=10
	ul=1.d0/ndev
	ul2=0.5d0*ul
	usq=ul*ul
	ucub=ul*ul*ul

	nrad=int((rad-ul2)/ul)
	nnx=nrad

	Do i=0,nnx
	xn=ul2+ul*i
	if(xn.gt.radt)exit
	yn=rad*rad-xn*xn
	if(yn.ge.0.d0)then
	yn=sqrt(yn)
	else
	stop 'yn negative!'
	endif
	nny=int((yn-ul2)/ul)
	xx=xn
	xr=cx+xx
	xr1=cx-xx
	nxr=int(xr)
	nxr1=int(xr1)
	if(xr.eq.real(nxr))nxr=nxr-1
	nxr=nxr+1
	if(nxr.eq.0)nxr=1
	if(xr1.eq.real(nxr1))nxr1=nxr1-1
	nxr1=nxr1+1
	if(nxr.eq.0)nxr1=1
	  Do j=0,nny
	  yy=ul2+ul*j
	  yr=cy+yy
	  nyr=int(yr)
	  if(yr.eq.real(nyr))nyr=nyr-1
	  nyr=nyr+1
	  if(nyr.eq.0)nyr=1
	  yr1=cy-yy
	  nyr1=int(yr1)
	  if(yr1.eq.real(nyr1))nyr1=nyr1-1
	  nyr1=nyr1+1
	  if(nyr1.eq.0)nyr=1
	  zz=rad*rad-xx*xx-yy*yy
	  if(zz.ge.0.d0)then
	  zz=sqrt(zz)
	  else
	write(*,*)'zz negative in i,j= ',i,j
	  stop 'zz negative!'
	  endif
	  zr=zz+cz
	  nzr=int(zr)
	  if(zr.eq.real(nzr))nzr=nzr-1
	  nzr=nzr+1
	  if(nzr.eq.0)nzr=1
	  nzz=int(zz)
	  zzn=real(nzz)
	  if(zz.gt.zzn)then
	  ht=zz-zzn
!	  dens(necfan(nxr,nyr,nzr))=dens(necfan(nxr,nyr,nzr))+ht*usq
!	  dens(necfan(nxr1,nyr,nzr))=dens(necfan(nxr1,nyr,nzr))+ht*usq
	  dens(necfan(nxr,nyr1,nzr))=dens(necfan(nxr,nyr1,nzr))+ht*usq
	  dens(necfan(nxr1,nyr1,nzr))=dens(necfan(nxr1,nyr1,nzr))+ht*usq
	  vol=vol+ht*usq
	  if(nzr.ge.1)then
	    do k=1,nzr-1
!	    dens(necfan(nxr,nyr,k))=dens(necfan(nxr,nyr,k))+usq!(usq*1)
!	    dens(necfan(nxr1,nyr,k))=dens(necfan(nxr1,nyr,k))+usq!(usq*1)
	    dens(necfan(nxr,nyr1,k))=dens(necfan(nxr,nyr1,k))+usq!(usq*1)
	    dens(necfan(nxr1,nyr1,k))=dens(necfan(nxr1,nyr1,k))+usq!(usq*1)
	    vol=vol+usq
	    enddo! k=1,nzz
	  endif
	  else
	    do k=1,nzr
!	    dens(necfan(nxr,nyr,k))=dens(necfan(nxr,nyr,k))+usq!(usq*1)
!	    dens(necfan(nxr1,nyr,k))=dens(necfan(nxr1,nyr,k))+usq!(usq*1)
	    dens(necfan(nxr,nyr1,k))=dens(necfan(nxr,nyr1,k))+usq!(usq*1)
	    dens(necfan(nxr1,nyr1,k))=dens(necfan(nxr1,nyr1,k))+usq!(usq*1)
	    vol=vol+usq
	    enddo! k=1,nzz
	  endif
	  Enddo! j=0,ny
	Enddo!i=0,nx


	nradt=int((radt-ul2)/ul)
	nnxt=nradt

	Do i=0,nnxt
	xn=ul2+ul*i
	yn=radt*radt-xn*xn
	if(yn.ge.0.d0)then
	yn=sqrt(yn)
	else
	stop 'yn negative!'
	endif
	nny=int((yn-ul2)/ul)
	xx=xn
	xr=cx+xx
	xr1=cx-xx
	nxr=int(xr)
	nxr1=int(xr1)
	if(xr.eq.real(nxr))nxr=nxr-1
	nxr=nxr+1
	if(nxr.eq.0)nxr=1
	if(xr1.eq.real(nxr1))nxr1=nxr1-1
	nxr1=nxr1+1
	if(nxr.eq.0)nxr1=1
	  Do j=0,nny
	  yy=ul2+ul*j
	  yr=cy+yy
	  nyr=int(yr)
	  if(yr.eq.real(nyr))nyr=nyr-1
	  nyr=nyr+1
	  if(nyr.eq.0)nyr=1
	  yr1=cy-yy
	  nyr1=int(yr1)
	  if(yr1.eq.real(nyr1))nyr1=nyr1-1
	  nyr1=nyr1+1
	  if(nyr1.eq.0)nyr=1
	  zz=radt*radt-xx*xx-yy*yy
	  if(zz.ge.0.d0)then
	  zz=sqrt(zz)
	  else
	write(*,*)'zz negative in i,j= ',i,j
	  stop 'zz negative!'
	  endif
	  zr=zz+cz
	  nzr=int(zr)
	  if(zr.eq.real(nzr))nzr=nzr-1
	  nzr=nzr+1
	  if(nzr.eq.0)nzr=1
	  nzz=int(zz)
	  zzn=real(nzz)
	  if(zz.gt.zzn)then
	  ht=zz-zzn
	  dens(necfan(nxr,nyr,nzr))=dens(necfan(nxr,nyr,nzr))+ht*usq
	  dens(necfan(nxr1,nyr,nzr))=dens(necfan(nxr1,nyr,nzr))+ht*usq
!	  dens(necfan(nxr,nyr1,nzr))=dens(necfan(nxr,nyr1,nzr))+ht*usq
!	  dens(necfan(nxr1,nyr1,nzr))=dens(necfan(nxr1,nyr1,nzr))+ht*usq
	  vol=vol+ht*usq
	  if(nzr.ge.1)then
	    do k=1,nzr-1
	    dens(necfan(nxr,nyr,k))=dens(necfan(nxr,nyr,k))+usq!(usq*1)
	    dens(necfan(nxr1,nyr,k))=dens(necfan(nxr1,nyr,k))+usq!(usq*1)
!	    dens(necfan(nxr,nyr1,k))=dens(necfan(nxr,nyr1,k))+usq!(usq*1)
!	    dens(necfan(nxr1,nyr1,k))=dens(necfan(nxr1,nyr1,k))+usq!(usq*1)
	    vol=vol+usq
	    enddo! k=1,nzz
	  endif
	  else
	    do k=1,nzr
	    dens(necfan(nxr,nyr,k))=dens(necfan(nxr,nyr,k))+usq!(usq*1)
	    dens(necfan(nxr1,nyr,k))=dens(necfan(nxr1,nyr,k))+usq!(usq*1)
!	    dens(necfan(nxr,nyr1,k))=dens(necfan(nxr,nyr1,k))+usq!(usq*1)
!	    dens(necfan(nxr1,nyr1,k))=dens(necfan(nxr1,nyr1,k))+usq!(usq*1)
	    vol=vol+usq
	    enddo! k=1,nzz
	  endif
	  Enddo! j=0,ny
	Enddo!i=0,nnx2
	  vol=vol*2
!	  write(*,*)'vol,4*pai*rad*rad*rad/3= ',vol,4*pai*rad*rad*rad/3

!************** mkhead *********************
!:::::::::::::::: lower bowl ::::::::::::::::::::

	rads=radm*dsin(sah)
	nxs=int((rads-ul2)/ul)
	nxs1=nxs+1
	nxm=int((radm-ul2)/ul)
	chz=rad*dcos(sa)+radm*dcos(sah)!z coordinate of femoral head center
	nchz=(chz-ul2)/ul!until under chz
	nchz1=nchz+1
	Do ix=0,nxs
	xn=ul2+ul*ix
	yns=rads*rads-xn*xn
	ynm=radm*radm-xn*xn
	  if(yns.ge.0.d0)then
	  yns=sqrt(yns)
	  else
	  stop 'yns negative'
	  endif!if(yns.ge.0.d0)
	  if(ynm.ge.0.d0)then
	  ynm=sqrt(ynm)
	  else
	  stop 'ynm negative'
	  endif!if(ynm.ge.0.d0)
	nys=int((yns-ul2)/ul)
	nys1=nys+1
	nym=int((ynm-ul2)/ul)
	    do iy=0,nys
	      yn=ul2+iy*ul
	      znr=rad*rad-xn*xn-yn*yn
	      if(znr.ge.0.d0)then
	      znr=sqrt(znr)
	      else
	  write(*,*) 'znr negative in lower bowl1 ix,iy=',ix,iy
	      stop
	      endif!if(znr.ge.0.d0)
	      nzr=int((znr+ul2)/ul)!from over znr
		do iz=nzr,nchz
		  zn=ul2+iz*ul
		  ori(1)=xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		enddo!iz=nzr,nchz
	    enddo! iy=0,nys
	nys1=nys+1
	    do iy=nys1,nym
	      yn=ul2+iy*ul
	      znr=radm*radm-xn*xn-yn*yn
	      if(znr.ge.0.d0)then
	      znr=sqrt(znr)
	      else
	  write(*,*) 'znr negative in lower bowl2 ix,iy=',ix,iy
	      stop
	      endif!if(znr.ge.0.d0)
	      znr=chz-znr
	      nzr=int((znr+ul2)/ul)!from over znr
		do iz=nzr,nchz
		  zn=ul2+iz*ul
		  ori(1)=xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		enddo!iz=nzr,nchz
	    enddo! iy=nys1,nym
	Enddo! ix=0,nxs
	Do ix=nxs1,nxm
	xn=ul2+ul*ix
	ynm=radm*radm-xn*xn
	  if(ynm.ge.0.d0)then
	  ynm=sqrt(ynm)
	  else
	  stop 'ynm negative'
	  endif!if(ynm.ge.0.d0)
	nym=int((ynm-ul2)/ul)
	    do iy=0,nym
	      yn=ul2+iy*ul
	      znr=radm*radm-xn*xn-yn*yn
	      if(znr.ge.0.d0)then
	      znr=sqrt(znr)
	      else
	  write(*,*) 'znr negative in lower bowl3 ix,iy=',ix,iy
	      stop
	      endif!if(znr.ge.0.d0)
	      znr=chz-znr
	      nzr=int((znr+ul2)/ul)!from over znr
		do iz=nzr,nchz
		  zn=ul2+iz*ul
		  ori(1)=xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	mshader(iix,iiy,iiz)=10
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
		  ori(1)=xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		enddo!iz=nzr,nchz
	    enddo! iy=nys1,nym
	Enddo! ix=nxs1,nxm

!:::::::::::::::: upper hat ::::::::::::::
	Do ix=0,nxm
	xn=ul2+ul*ix
	ynm=radm*radm-xn*xn
	  if(ynm.ge.0.d0)then
	  ynm=sqrt(ynm)
	  else
	  stop 'ynm negative'
	  endif!if(ynm.ge.0.d0)
	nym=int((ynm-ul2)/ul)
	    do iy=0,nym
	      yn=ul2+iy*ul
	      znr=radm*radm-xn*xn-yn*yn
	      if(znr.ge.0.d0)then
	      znr=sqrt(znr)
	      else
	  write(*,*) 'znr negative in lower bowl3 ix,iy=',ix,iy
	      stop
	      endif!if(znr.ge.0.d0)
	      znr=chz+znr
	      nzr=int((znr-ul2)/ul)!from over znr
		do iz=nchz1,nzr
		  zn=ul2+iz*ul
		  ori(1)=xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		  ori(1)=-xn
		  ori(2)=-yn
		  ori(3)=zn
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub
	mshader(iix,iiy,iiz)=10
		enddo! iz=nchz1,nzr
	    enddo! iy=nys1,nym
	Enddo! ix=0,nxm

!************** mkround ********************
!                  area of loaded surface of the femoral head
!        PAI = DACOS(-1.D0)
!        rad=10.d0
! latitude:equater to pole;north to south , longitude:west to east
	ori(1)=0.d0
	ori(2)=0.d0
	ori(3)=chz
        call trans(ori,tra,ans,ata)
	chxtra=tra(1)!x coordinate of translated femoral head center
	chytra=tra(2)!y coordinate of translated femoral head center
	chztra=tra(3)!z coordinate of translated femoral head center

        ndlat=int(als*radm/ul)
        ndlon=int(2.d0*pai*radm/ul)
!        sa=90.d0
!        if(sa.ge.90)then
        nlatn=ndlat
!        endif
        dlatn=als/ndlat
        dlon=2*pai/ndlon
        garea=0.d0
        do nn=1,ndlon
        dlonnn=dlon*nn
        coslon=dcos(dlonnn)
        sinlon=dsin(dlonnn)
        rad2=radm*dsin(als)
!
        q(1,1)=rad2!rad*cos(0.)
        q(1,2)=0.d0!rad*sin(0.)
        q(1,3)=radm*dcos(als)
!
        q(4,1)=rad2*dcos(dlon)
        q(4,2)=rad2*dsin(dlon)
        q(4,3)=q(1,3)
!
        rad2=rad*dcos(csa)
!        nlatn1=nlatn-1
        do i=1,nlatn
        parea=0.d0
        rad1=rad2
        rad2=radm*dsin(als-dlatn*i)
                do ic=1,3
                q(2,ic)=q(1,ic)
                q(3,ic)=q(4,ic)
                enddo!ic=1,3
!
                q(1,1)=rad2
                q(1,2)=0.d0
                q(1,3)=radm*dcos(als-dlatn*i)
                q(4,1)=rad2*dcos(dlon)
                q(4,2)=rad2*dsin(dlon)
                q(4,3)=q(1,3)
          do ic=1,3
           v(ic,1)=q(2,ic)-q(1,ic)
           v(ic,2)=q(3,ic)-q(1,ic)
          enddo
          v11=0.d0
          do ic=1,3
           v11=v11+v(ic,1)*v(ic,1)
          enddo
          v22=0.d0
          do ic=1,3
           v22=v22+v(ic,2)*v(ic,2)
          enddo
          v12=0.d0
          do ic=1,3
           v12=v12+v(ic,1)*v(ic,2)
          enddo
          parea1=v11*v22-v12*v12
         if(parea1 .ge. 0.d0)then
        parea1=0.5d0*sqrt(parea1)
          garea=garea+parea1
         else
        write(*,*)'negative triangle, 1 of nlatn= ',1,'of ',i
          stop'parea missing;negative!'
         endif
        parea2=0.d0
          do ic=1,3
           v(ic,1)=q(4,ic)-q(3,ic)
           v(ic,2)=q(1,ic)-q(3,ic)
          enddo
          v11=0.d0
          do ic=1,3
           v11=v11+v(ic,1)*v(ic,1)
          enddo
          v22=0.d0
          do ic=1,3
           v22=v22+v(ic,2)*v(ic,2)
          enddo
          v12=0.d0
          do ic=1,3
           v12=v12+v(ic,1)*v(ic,2)
          enddo
          parea2=v11*v22-v12*v12
         if(parea2 .ge. 0.d0)then
        parea2=0.5d0*sqrt(parea2)
          garea=garea+parea2
         else
        write(*,*)'negative triangle, 2 of nlatn= ',2,'of ',i
          stop'parea missing;negative!'
         endif
        parea=parea1+parea2
!        write(*,*)'parea,ntse= ',parea,ntse
!        ase(NTSE)=parea
        do ic=1,3
        ori(ic)=0.5d0*(q(1,ic)+q(3,ic))
        enddo
        radq=sqrt(ori(1)*ori(1)+ori(2)*ori(2))
        ori(1)=radq*coslon
        ori(2)=radq*sinlon
	ori(3)=ori(3)+chz
!        write(*,*)'calling trans'
        call trans(ori,tra,ans,ata)
!        write(*,*)'return from trans'
        ix=int(tra(1)+cx)+1
        iy=int(tra(2)+cy)+1
        iz=int(tra(3)+cz)+1
        if(ix.gt.nx)ix=nx
        if(iy.gt.ny)iy=ny
        if(iz.gt.nz)iz=nz
        if(ix.lt.1)ix=1
        if(iy.lt.1)iy=1
        if(iz.lt.1)iz=1
	tra(1)=chxtra-tra(1)
	tra(2)=chytra-tra(2)
	tra(3)=chztra-tra(3)

        abstra=sqrt(tra(1)*tra(1)+tra(2)*tra(2)+tra(3)*tra(3))
!	abstra=sqrt((tra(1)-chxtra)*(tra(1)-chxtra)+(tra(2)-chytra)*(tra(2)-chytra)+(tra(3)-chztra)*(tra(3)-chztra))
        shader(ix,iy,iz)=shader(ix,iy,iz)+parea
        mshader(ix,iy,iz)=1
        nshade(ix,iy,iz)=nshade(ix,iy,iz)+1
        pxyz(NECFAN(ix,iy,iz))=pp1
        do ic=1,3
        pelem(NECFAN(ix,iy,iz),ic)=pelem(NECFAN(ix,iy,iz),ic)+tra(ic)/abstra
        enddo! ic=1,3
        enddo! i=1,nlatn

        enddo!nn=1,ndlon
!

!************** mk trochanter maj. **************
!	rhb=1.d0!Ratio of hts by sides:hts/(0.5*(sbl+sul))
!	rub=1.d0!Ratio of sul by sbl:sul/sbl
!	rdb12=0.d0!Ratio of length between lower end of sul and sbl by sbl:d12/sbl
!	rcb1=0.5d0!Ratio of length between lower end of sbl and side_center(cc1) by sbl:dcs/sbl 
!	rcb2=0.5d0!Ratio of length between left_side and hight_center(cc2) by sbl:cc2/sbl
!	area=0.5*(sbl+sul)*hts=0.5*sbl*sbl*(1+rub)*rhb
!	trpzl(4,3):lower trapezoid,trpzu(4,3):upper trapezoid,1=UpperLeft,2=LowerLeft,3=LowerRight,4=UpperRight
!	as12:slope angle of left side
!	as34:slope angle of right side
!	ab41:slope angle of upper base
!	ab23:slope angle of lower base
	sarea=pai*rads*rads
!	sarea=sarea*0.8d0
	sbl=sqrt(sarea*4.d0/((1.d0+rub)*(1.d0+rub)*rhb))
	write(*,*)'sbl= ',sbl
	sul=sbl*rub
	hts=0.5d0*(sbl+sul)*rhb
	d12=sbl*rdb12
	d34=sbl*rdb34
	cc1=sbl*rcb1!trochanter load center;x coordinate
	cc2=hts*rcb2!trochanter load center;y coordinate

	trpzl2(1,1)=d12
	trpzl2(2,1)=0.d0
	trpzl2(3,1)=sbl
	trpzl2(4,1)=d12+sul
	trpzl2(1,2)=hts
	trpzl2(2,2)=0.d0
	trpzl2(3,2)=0.d0
	trpzl2(4,2)=hts
	ccz=0.d0
	trpzl2(1,3)=ccz
	trpzl2(2,3)=ccz
	trpzl2(3,3)=ccz
	trpzl2(4,3)=ccz

	sleng=0.5d0*(trpzu(3,1)-trpzu(2,1)+trpzu(4,1)-trpzu(1,1))*(1.d0-ra12)

	trpzu(1,2)=trpzu(1,2)
	trpzu(2,2)=trpzu(2,2)
	trpzu(3,2)=trpzu(3,2)
	trpzu(4,2)=trpzu(4,2)

	trpzu(1,3)=trpzl2(1,3)
	trpzu(2,3)=trpzu(1,3)
	trpzu(3,3)=trpzu(1,3)
	trpzu(4,3)=trpzu(1,3)

!	trpzl(1,1)=d12
!	trpzl(2,1)=0.d0
!	trpzl(3,1)=sbl
!	trpzl(4,1)=d12+sul
!	trpzl(1,2)=hts
!	trpzl(2,2)=0.d0
!	trpzl(3,2)=0.d0
!	trpzl(4,2)=hts
!	ccz=0.d0
!	trpzl(1,3)=ccz
!	trpzl(2,3)=ccz
!	trpzl(3,3)=ccz
!	trpzl(4,3)=ccz

	sl12=sqrt(d12*d12+hts*hts)
	if(sl12.le.0.d0)stop'sl12 is negative!'
	csla12=hts/sl12
	tsla12=d12/hts
	if(abs(dcos(as12)).gt.0.d0)then
	 psl(1)=httr/dtan(as12)!psl12
	 th12=psl(1)/csla12
	else
	 psl(1)=0.d0
	 th12=0.d0
	endif

	sl34=sqrt(d34*d34+hts*hts)
	if(sl34.le.0.d0)stop'sl12 is negative!'
	csla34=hts/sl34
	tsla34=d34/hts
	if(abs(dcos(as34)).gt.0.d0)then
	 psl(3)=httr/dtan(as34)!psl34
	 th34=psl(3)/csla34
	else
	 psl(3)=0.d0
	 th34=0.d0
	endif

	if(abs(dcos(ab23)).gt.0.)then
	 psl(2)=httr/dtan(ab23)!psl23
	else
	 psl(2)=0.d0
	endif

	if(abs(dcos(ab41)).gt.0.d0)then
	 psl(4)=httr/dtan(ab41)!psl41
	else
	 psl(4)=0.d0
	endif
!q1

!	trpzu(1,1)=trpzl(1,1)+th12+psl(4)*tsla12
	trpzl(1,1)=trpzu(1,1)-th12-psl(4)*tsla12
!	trpzu(2,1)=trpzl(2,1)+th12-psl(2)*tsla12
	trpzl(2,1)=trpzu(2,1)-th12+psl(2)*tsla12
!	trpzu(3,1)=trpzl(3,1)-th34+psl(2)*tsla34
	trpzl(3,1)=trpzu(3,1)+th34-psl(2)*tsla34
!	trpzu(4,1)=trpzl(4,1)-th34-psl(4)*tsla34
	trpzl(4,1)=trpzu(4,1)+th34+psl(4)*tsla34

!	trpzu(1,2)=trpzl(1,2)-psl(4)
!	trpzu(2,2)=trpzl(2,2)+psl(2)
!	trpzu(3,2)=trpzl(3,2)+psl(2)
!	trpzu(4,2)=trpzl(4,2)-psl(4)
	trpzl(1,2)=trpzu(1,2)+psl(4)
	trpzl(2,2)=trpzu(2,2)-psl(2)
	trpzl(3,2)=trpzu(3,2)-psl(2)
	trpzl(4,2)=trpzu(4,2)+psl(4)
!	ccz=rad
!	trpzu(1,3)=httr
!	trpzu(2,3)=httr
!	trpzu(3,3)=httr
!	trpzu(4,3)=httr
	trpzl(1,3)=0.d0
	trpzl(2,3)=0.d0
	trpzl(3,3)=0.d0
	trpzl(4,3)=0.d0

!	sleng=0.5d0*(trpzu(3,1)-trpzu(2,1)+trpzu(4,1)-trpzu(1,1))*(1.d0-ra12)
!	trpzl2(1,1)=trpzu(1,1)+sleng
!	trpzl2(2,1)=trpzu(2,1)+sleng
!	trpzl2(3,1)=trpzu(3,1)
!	trpzl2(4,1)=trpzu(4,1)

!	trpzl2(1,2)=trpzu(1,2)
!	trpzl2(2,2)=trpzu(2,2)
!	trpzl2(3,2)=trpzu(3,2)
!	trpzl2(4,2)=trpzu(4,2)

!	trpzl2(1,3)=trpzu(1,3)
!	trpzl2(2,3)=trpzl2(1,3)
!	trpzl2(3,3)=trpzl2(1,3)
!	trpzl2(4,3)=trpzl2(1,3)

	do nslide=1,4
	  do i=1,3
	vec(i)=trpzu(nslide,i)-trpzl(nslide,i)
	  enddo!i=1,3
	  do i=1,3
	trpzu2(nslide,i)=trpzl2(nslide,i)+vec(i)*httr2/vec(3)
	  enddo!i=1,3
	enddo!nslide=1,4

	do i=1,4
	psl2(i)=psl(i)*httr2/httr
	enddo! i=1,4

	do j=1,4
	do i=1,3
!	write(*,*)'i,j= ',i,j
!	write(*,*)'trpzl(j,i),trpzu(j,i)= ',trpzl(j,i),trpzu(j,i)
	enddo!i=1,3
!	write(*,*)' '
	enddo!j=1,4

	do j=1,4
	do i=1,3
!	write(*,*)'i,j= ',i,j
!	write(*,*)'trpzl2(j,i),trpzu2(j,i)= ',trpzl2(j,i),trpzu2(j,i)
	enddo!i=1,3
!	write(*,*)' '
	enddo!j=1,4

!!!!	call mktr(trpzl,trpzu,psl,p1)
	pp2=pp3
!!!!	call mktr(trpzl2,trpzu2,psl2,p2)
!	goto 1000
!1000	continue


!        write(*,*)'garea,garea*2= ',garea,garea*2
!        write(*,*)'globe(4*pai*rad**2) = ',4.d0*pai*rad*rad

        gosa=4.d0*pai*rad*rad-garea*2
!        write(*,*)'gosa= ',gosa
!        ase(NTSE)=parea

!1000 continue

!***************mkround **********************
        do k=1,nz
        do j=1,ny
        do i=1,nx
        if(mshader(i,j,k).ne.0)then
        dens(necfan(i,j,k))=dens(necfan(i,j,k))+shade(i,j,k)
        endif
        if(mshader(i,j,k).eq.1)then
        do ic=1,3
        pelem(NECFAN(i,j,k),ic)=pelem(NECFAN(i,j,k),ic)/nshade(i,j,k)
        enddo! ic=1,3
        endif
        enddo
        enddo
        enddo

        nls=0
        do k=1,nz
        do j=1,ny
        do i=1,nx
	if(dens(necfan(i,j,k)).gt.1.0d0)dens(necfan(i,j,k))=1.0d0
        if(mshader(i,j,k).eq.10) then
        nls=nls+1
!        dens(necfan(i,j,k))=10.d0
        endif
!qqqq
        if(mshader2(i,j,k).eq.1) then
!        write(*,*)'mshader= ',mshader(i,j,k)
!        write(*,*)'nshade= ',nshade(i,j,k)
!        dens(necfan(i,j,k))=15.d0
        endif
!        if(mshader(i,j,k).eq.1) then
!        dens(necfan(i,j,k))=10.d0
!        endif
!        if(dens(necfan(i,j,k)).gt.0.99d0) then
!        dens(necfan(i,j,k))=1.d0
!        endif
        if(mshader(i,j,k).eq.1) then
        NLD=NLD+1
        NEL(NLD)=necfan(i,j,k)
        earea(NLD)=shader(i,j,k)
        px(NLD)=pxyz(NECFAN(i,j,k))*pelem(NECFAN(i,j,k),1)
        py(NLD)=pxyz(NECFAN(i,j,k))*pelem(NECFAN(i,j,k),2)
        pz(NLD)=pxyz(NECFAN(i,j,k))*pelem(NECFAN(i,j,k),3)
!        write(*,*)'NLD,p(NLD)*earea(in main))= ',nld,sqrt(px(NLD)**2+py(NLD)**2+pz(NLD)**2)*earea(nld)
        endif

        enddo! i=1,nx
        enddo! j=1,ny
        enddo! k=1,nz
!	nls:near loaded surface
	allocate(densnls(nls),numnls(nls))
	nnls=0
        do k=1,nz
        do j=1,ny
        do i=1,nx
        if(mshader(i,j,k).eq.10) then
	nnls=nnls+1
        numnls(nnls)=necfan(i,j,k)
	densnls(nnls)=dens(necfan(i,j,k))
        endif
        enddo! i=1,nx
        enddo! j=1,ny
        enddo! k=1,nz
!*******************************************************
        call output


        iterst=1
      OPEN(NF,FILE='DAT3D.TXT',STATUS='UNKNOWN')
      WRITE(NF,*)'DATA in '
      WRITE(NF,*)gratio
      WRITE(NF,*)iterst
      WRITE(NF,*)NX
      WRITE(NF,*)NY
      WRITE(NF,*)NZ
      WRITE(NF,*)DENSMIN
!      write(*,*)' NOE = ',NOE

      NDX=(NX-NCX)/2
      WRITE(NF,*)NOE
      DO I=1,NOE
!        DENS(I)=0.001
         WRITE(NF,*)DENS(I)
      ENDDO
      WRITE(NF,*)(NX+1)*(NY+1)
      NUM=0
      DO j=1,NY+1
      DO i=1,NX+1
      NUM=NUM+1
      IF(i.eq.int(cx+1).and.j.eq.int(cy+1))THEN
!      IF(NUM.eq.1)THEN
         WRITE(NF,'(2I10)')NODFAN(i,j,1),4
        write(*,*)'ix,iy(in Fix 4 point:x,y,z_fixed) = ',i,j
!        write(*,*)'       cx,cy=  ',cx,cy
      ELSE
         WRITE(NF,'(2I10)')NODFAN(i,j,1),3
      ENDIF
      ENDDO
      ENDDO

      WRITE(NF,*)NLD
      WRITE(*,*)NLD
        area=0.d0
      DO I=1,NLD
         WRITE(NF,*)NEL(I),earea(i),px(i),py(i),pz(i)
!        pp=sqrt(px(i)*px(i)+py(i)*py(i)+pz(i)*pz(i))
!        write(*,*)'pp=',pp
!        write(*,*)'earea(i)*pp= ',earea(i)*pp
        area=area+earea(i)
      ENDDO
        write(*,*)'area*2= ',area*2
!        write(*,*)'globe(4*pai*rad**2) = ',4.d0*pai*rad*rad
!      WRITE(NF,*)px,py,pz

!
	write(NF,'(I10)')nls
!	write(*,*)'nls= ',nls
	
      DO i=1,nls
	write(NF,*)numnls(i),densnls(i)
!	write(*,*)'numnls(i),densnls(i)= ',numnls(i),densnls(i)
      ENDDO

      WRITE(NF,*)YOUNG,ANU

      CLOSE(NF)
        write(*,*)'Normal Ending!'
        stop
        end program main
!
!    ==================================
     subroutine trans(ori,tra,ans,ata)
!    ==================================
        real(8)::ori(3),tra(3),ans,ata
        real(8) tmp(3)
        real(8) radxy,radyz

        sin2=sin(ans)
        cos2=cos(ans)
        radyz=ori(2)*ori(2)+ori(3)*ori(3)
!        write(*,*)'radyz= ',radyz
        radyz=sqrt(radyz)
        if(radyz .eq. 0.0d0)then
!        write(*,*)'radyz=0 in i= ',i
        tmp(1)=ori(1)
        tmp(2)=ori(2)
        tmp(3)=ori(3)
        else
        cos1=ori(2)/radyz
        sin1=ori(3)/radyz
        cos12=cos1*cos2-sin1*sin2
        sin12=sin1*cos2+cos1*sin2
        tmp(1)=ori(1)
        tmp(2)=radyz*cos12
        tmp(3)=radyz*sin12
        endif

        sin4=sin(ata)
        cos4=cos(ata)
        radxy=tmp(1)*tmp(1)+tmp(2)*tmp(2)
        radxy=sqrt(radxy)
        if(radxy .eq. 0.0d0)then
        tra(1)=tmp(1)
        tra(2)=tmp(2)
        tra(3)=tmp(3)
        else
        cos3=tmp(1)/radxy
        sin3=tmp(2)/radxy
        sin34=sin3*cos4+cos3*sin4
        cos34=cos3*cos4-sin3*sin4
        tra(1)=radxy*cos34
        tra(2)=radxy*sin34
        tra(3)=tmp(3)
        endif
     end subroutine trans

!-----------------------------------------------------------
      SUBROUTINE OUTPUT
!-----------------------------------------------------------
        use com0
        use com1
        use com2
        use com3

      NX1=NX+1
      NY1=NY+1
      NZ1=NZ+1
      DO k=1,NZ1
      kk=k
      if(kk.eq.NZ1)kk=NZ
      DO j=1,NY1
      jj=j
      if(jj.eq.NY1)jj=NY
      DO i=1,NX1
      ii=i
      if(ii.eq.NX1)ii=NX
      den(NODFAN(i,j,k))=DENS(NECFAN(ii,jj,kk))
      ENDDO
      ENDDO
      ENDDO

!     Base Plate
        do j=1,NY1
        do i=1,NX1
        den(NODFAN(i,j,1))=2.0d0
!        den(NODFAN(i,j,1))=0.5d0
        enddo! i=1,NX1
        enddo! j=1,NY1
!OUTPUT
!
        NDATA=8*NOE+NOE
!
!     call outdatb(NON,NOE,NDATA, CNODX, CNODY, CNODZ,NEC,S_SIG,dens)

!
	write(*,*)'gratio= ',gratio
    write(*,*)'now in output'
    OPEN(15,FILE='Initial.vtk')
    WRITE(15,'(''# vtk DataFile Version 2.0'')')
    WRITE(15,'(''Data'')')
    WRITE(15,'(''ASCII'')')
!    WRITE(15,'(''           '')')
    WRITE(15,'(''DATASET UNSTRUCTURED_GRID'')')
    WRITE(15,'(''POINTS '',i10,'' float'')') NON
    DO i=1,NON
       WRITE(15,'(3e15.7)') CNODX(I)*gratio,CNODY(I)*gratio,CNODZ(I)*gratio
    ENDDO
!
        NDATA=8*NOE+NOE
    WRITE(15,'(''CELLS '',2i10)') NOE, NDATA
    DO i=1,NOE
!      ne = ele(i)%non
      WRITE(15,'(9i10)') 8,NEC(I)%NL(1)-1,NEC(I)%NL(2)-1,NEC(I)%NL(3)-1,NEC(I)%NL(4)-1,&
      NEC(I)%NL(5)-1,NEC(I)%NL(6)-1,NEC(I)%NL(7)-1,NEC(I)%NL(8)-1
    ENDDO
!
    WRITE(15,'(''CELL_TYPES '',i10)') NOE
    DO i=1,NOE
       WRITE(15,*) 12
    ENDDO
!
!    WRITE(15,'(''POINT_DATA '',i10)') NON
    WRITE(15,'(''POINT_DATA '',i10)') NON
!
! Velocity
! density
!----------
!    WRITE(15,'(''VECTORS vectors float'')')
    WRITE(15,'(''SCALARS DENSITY float'')')
    WRITE(15,'(''LOOKUP_TABLE default'')')
!    DO i=1,NON
!       WRITE(15,'(3e15.7)') Sx(i), Sy(i), Sz(i)
!    ENDDO
    DO i=1,NON
       WRITE(15,'(e15.7)') den(i)
    ENDDO
!
!
    CLOSE(15)
!
!
      RETURN
      END SUBROUTINE OUTPUT

                function vnear(a)result(v)!seisuuka
                        real(8),intent(in)::a
                        real(8) v,aa
                        aa=gehalf(a)
                        if(abs(aa-a) .lt. 0.0001)then
                        v=aa
                        else
                        v=a
                        endif

                endfunction vnear

                function gehalf(a)result(v)!sisyagonyuu
                        real(8),intent(in)::a
                        real(8) v,aa
                if(a .gt. 0.0d0)then
                        aa=dble(int(a))
                        if(a-aa .ge. 0.5)aa=aa+1
                else
                        aa=dble(int(a))-1
                        if(a-aa .gt. 0.5)aa=aa+1
                endif
        v=aa
                endfunction gehalf
!*********************************Start	line process
	subroutine lprocv
        use com0
        use com1
        use com2
        use com3
	real(8)ucub2
!	line process
	  do i=1,3
	veclr(i)=pnr(i)-pnl(i)
	  enddo! i=1,3
	abslr=sqrt(veclr(1)*veclr(1)+veclr(2)*veclr(2)+veclr(3)*veclr(3))
	if(abslr.eq.0.)stop'abslr.eq.0. in lprocv'
!	write(*,*)'abslr= ',abslr
	dnlr=abslr/ul
	ndlr=int(dnlr)
	reslr=dnlr-real(ndlr)
	reslr2=0.5d0*reslr
	reslr22=0.5d0*reslr2
	  do i=1,3
	uveclr(i)=veclr(i)/dnlr
	  enddo! i=1,3

!	first point
	rp3d=reslr22
	  do i=1,3
	p3d(i)=pnl(i)+rp3d*uveclr(i)
	  enddo! i=1,3
	ucub2=ul*ul*ul*reslr2
!	write(*,*)'Calling pprocv first point'
!**********************************call	point process
	call pprocv(ucub2)
!	mid points
	ucub2=ul*ul*ul
	if(ndlr.ge.1)then
	  do ipoin=1,ndlr
	rp3d=reslr2+0.5d0+real(ipoin-1)
	    do i=1,3
	p3d(i)=pnl(i)+rp3d*uveclr(i)
	    enddo! i=1,3

!	write(*,*)'Calling pprocv in mid point ipoin= ',ipoin
!	***** Call Point process *****
!	write(*,*)'ipoin= ',ipoin
	call pprocv(ucub2)


	  enddo! ipoin=1,ndlr
	endif!ndlr.ge.1
!	last point
	rp3d=reslr2+real(ndlr)+reslr22
	  do i=1,3
	p3d(i)=pnl(i)+rp3d*uveclr(i)
	  enddo! i=1,3
	ucub2=ul*ul*ul*reslr2
!	write(*,*)'Calling pprocv in last point'
!	***** Call Point process *****
	call pprocv(ucub2)


!***********************************End	line process
	return
	end subroutine lprocv


        subroutine lproca(hta,nmshader)
        use com0
        use com1
        use com2
        use com3

	real(8)hta,uarea
	integer nmshader
!        write(*,*)'hta,ul(in lproca)= ',hta,ul
!       line process
!	line process
	  do i=1,3
	veclr(i)=pnr(i)-pnl(i)
!	write(*,*)'i,veclr,pnr,pnl= ',i,veclr(i),pnr(i),pnl(i)
	  enddo! i=1,3
	abslr=sqrt(veclr(1)*veclr(1)+veclr(2)*veclr(2)+veclr(3)*veclr(3))
	if(abslr.eq.0.)stop'abslr.eq.0. in lproca'

	dnlr=abslr/ul
	ndlr=int(dnlr)
	reslr=dnlr-real(ndlr)
	reslr2=0.5d0*reslr
	reslr22=0.5d0*reslr2
	  do i=1,3
	uveclr(i)=veclr(i)/dnlr
	  enddo! i=1,3

!	first point
	rp3d=reslr22
	  do i=1,3
	p3d(i)=pnl(i)+rp3d*uveclr(i)
	  enddo! i=1,3
	uarea=hta*reslr2*ul
!**********************************call	point process
!        write(*,*)'uarea= ',uarea
	call pproca(uarea,nmshader)

!	mid points
	uarea=hta*ul
	if(ndlr.ge.1)then
	  do ipoin=1,ndlr
	rp3d=reslr2+0.5+real(ipoin-1)
	    do i=1,3
	p3d(i)=pnl(i)+rp3d*uveclr(i)
	    enddo! i=1,3

!	***** Call Point process *****
!        write(*,*)'uarea= ',uarea
	call pproca(uarea,nmshader)


	  enddo! ipoin=1,ndlr
	endif!ndlr.ge.1
!	last point
	rp3d=reslr2+real(ndlr)+reslr22
	  do i=1,3
	p3d(i)=pnl(i)+rp3d*uveclr(i)
	  enddo! i=1,3

	uarea=hta*reslr2*ul
!	***** Call Point process *****
!        write(*,*)'uarea= ',uarea
	call pproca(uarea,nmshader)

        return
        end subroutine lproca



!**********************************Start 	point process
	subroutine pprocv(ucub2)
        use com0
        use com1
        use com2
        use com3
	real(8)ucub2
!	point process for volume calculation
	ccz=radt+p3d(3)
	nccz=(ccz-ul2)/ul
!	write(*,*)'ucub= ',ucub
!	nhts=int((hts-ul2)/ul)
!	dsh=sbl-d12
!	Do ih=0,nhts
!	  hi=ul2+ul*ih
	  yn=p3d(2)-cc2!	rcb2=0.5d0!Ratio of length between left_side and hight_center(cc2) by sbl:cc2/sbl
!	  do is=nsil,nsih
!	    si=ul2+ul*is
	    xn=p3d(1)-cc1!	rcb1=0.5d0!Ratio of length between lower end of sbl and side_center(cc1) by sbl:dcls/sbl 
	    znr=radt*radt-xn*xn-yn*yn
	    if(znr.ge.0.d0)then
	      znr=sqrt(znr)
	    else
	    znr=0.d0
!q
!	write(*,*)'cc1,cc2,rad,xn,yn,znr= ',cc1,cc2,rad,xn,yn,znr
!	write(*,*)'p3d= ',p3d(1),p3d(2),p3d(3)
!	      write(*,*)'znr negative in pprocv '
!	      stop
	    endif!if(znr.ge.0.do)then
	    nzr=(znr+ul2)/ul
!	write(*,*)'nzr,nccz,nccz-nzr= ',nzr,nccz,nccz-nzr
	      do iz=nzr,nccz
		  zn=ul2+iz*ul
		  ori(1)=xn
		  ori(2)=yn
		  ori(3)=zn
        call trans(ori,tra,ans2,ata2)
!        write(*,*)'return from trans'
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
	dens(necfan(iix,iiy,iiz))=dens(necfan(iix,iiy,iiz))+ucub2
!qqqq
!        mshader(iix,iiy,iiz)=10
	if(zn.gt.radt)mshader(iix,iiy,iiz)=10
	      enddo! iz=nzr,nccz
!	  enddo! is=nsil,nsih
!	Enddo! ih=0,nhts
	return
	end subroutine pprocv


	subroutine pproca(usq2,nmshader)
!	point process for area calculation
        use com0
        use com1
        use com2
        use com3

	real(8)usq2
	integer nmshader

!	  do i=1,3
!	ori(i)=p3d(i)
!	  enddo! i=1,3
	ori(1)=p3d(1)-cc1
	ori(2)=p3d(2)-cc2
	ori(3)=p3d(3)+radt
        call trans(ori,tra,ans2,ata2)
        iix=int(tra(1)+cx)+1
        iiy=int(tra(2)+cy)+1
        iiz=int(tra(3)+cz)+1
        if(iix.gt.nx)iix=nx
        if(iiy.gt.ny)iiy=ny
        if(iiz.gt.nz)iiz=nz
        if(iix.lt.1)iix=1
        if(iiy.lt.1)iiy=1
        if(iiz.lt.1)iiz=1
        shader(iix,iiy,iiz)=shader(iix,iiy,iiz)+usq2
        mshader(iix,iiy,iiz)=nmshader
        mshader2(iix,iiy,iiz)=1
        nshade(iix,iiy,iiz)=nshade(iix,iiy,iiz)+1
        pxyz(NECFAN(iix,iiy,iiz))=pp2

        pelem(NECFAN(iix,iiy,iiz),1)=pelem(NECFAN(iix,iiy,iiz),1)+ptrox
        pelem(NECFAN(iix,iiy,iiz),2)=pelem(NECFAN(iix,iiy,iiz),2)+ptroy
        pelem(NECFAN(iix,iiy,iiz),3)=pelem(NECFAN(iix,iiy,iiz),3)+ptroz


	return
	end subroutine pproca

