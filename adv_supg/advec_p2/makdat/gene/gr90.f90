!===============================================================
program mesh 
!===============================================================
      implicit double precision (a-h , o-z)

     integer,parameter ::   md1 = 400000,   md2 = 730000, md3 = 20000,   &
                   md4 = 20000,   md5 = 20000

   double precision, allocatable ::     xyz1(:,:),         &
                   xyz2(:,:),     lmax(:),     knod(:)
   double precision    ::           rn(20),        xyza(3),       xyzb(3)
   integer, allocatable :: nodc1(:,:),nnod(:,:),nodc2(:,:)

      character*50         inpfile,mesfile

  allocate(  xyz1(3,md1), xyz2(3,md1), lmax(md1), knod(md5),   &
                nodc1(20,md3),nnod(md4,md1),nodc2(8,md2)  ) 


      open (9,file='gr.file',status='unknown')
      read (9,'(a)') inpfile
      read (9,'(a)') mesfile

      open(50, file = inpfile)
!cube72.mesファイルの作成
      open(60, file = mesfile)

      call inmesh  (5,       isw,      icoo,     icon1,    icon, &
                   nmax1,    iemax1,   xyz1,     nodc1,    eps,  &
                   md1,      md3)

      nod2   = 0
      ie2    = 0
      idum   = 0
      nmax2  = 0
      iemax2 = 0
      idmax  = 0

      do 1000  loop = 1, iemax1

      call step01  (5,      isw,     ie1,      idivx,    idivy,  &
                  idivz,  itype,   nodx,     nody,     nodz,     &
                  dsx,      dsy,      dsz)

      call step02  (isw,    nodx,    nody,    nodz,    nod2,  &
                  rn,     idivx,   idivy,   idivz,   icoo,    &
                  icon1,  xyz1,    xyz2,    ie1,     nodc1,   &
                  dsx,    dsy,     dsz ,    md1,     md3)

      call step03  (isw,    idivx,   idivy,   idivz,   idum,  &
                   ie2,    nmax2,   idmax,   nodx,    itype,  &
                  nodc2,  iemax2,  md2,     xyz2,    md1)

      nmax2  = nod2
      iemax2 = ie2
      idmax  = idum

 1000 continue
      if(nmax2 .gt. md1)  then
       write(0,*)   'nmax2 .gt. md1'
       write(0,*)   'nmax2 =', nmax2
       write(0,*)   '  md1 =', md1
      stop
      end if
      if(iemax2 .gt. md2)  then
       write(0,*)   'iemax2 .gt. md2'
       write(0,*)   'iemax2 =', iemax2
      write(0,*)   '   md2 =', md2
      stop
      end if

      call step04  (nmax2,   nnod,     xyza,     xyzb,     xyz2,  & 
                 icoo,    eps,      kmax,     knod,     lmax,     & 
                 md1,      md4,     md5,      21)
 write(*,*)"1"

      call step05  (nmax2,   lmax,    nnod,     icoo,     md4,  &
                  iemax2,  icon,    nodc2,    md1,      md2)
 write(*,*)"2"

      call step06  (nmax2,   kmax,    knod,    icoo,     xyz1,  &
                  xyz2,    iemax2,  icon,    nodc2,   nmax1,    &
                  md1,      md2,      md5)
 write(*,*)"3"

      call renumb  (xyz1,     nodc2,    nmax1,    iemax2,   md1,  &
                  md2,      isw,      icoo,     icon)
 write(*,*)"4"

      call oumesh  (20,      isw,     nmax1,    iemax2,   xyz1,  &
                  nodc2,   icoo,    icon,     md1,      md2)
 write(*,*)"5"

      call echeck  (nmax1,   iemax2,  xyz1,     nodc2,    icoo,   &   
              icon,     isw,     md1,    md2)
      stop
      end

!--------------------------------------------------------------

      subroutine inmesh (ir,   isw,    icoo,   icon1,   icon,  &
                       nmax, iemax,  xyz,    nodc,    eps,     &
                       md1,      md3)

      implicit double precision (a-h , o-z)

      dimension     xyz(3,md1),    nodc(20,md3)

      write(0,*)"in inmesh"

      read(50,*)   eps
      read(50,*)   isw 

      write(0,*)"eps,isw OK",isw

      if(isw .eq. 1)  then
       icoo  = 2
       icon1 = 8
       icon  = 3
  !これが適用される
      else if(isw .eq. 2)  then
       icoo  = 2
       icon1 = 8
       icon  = 4

      else
       icoo  = 3
       icon1 = 20
       icon  = 8
      end if

       write(0,*)"icoo OK",icoo
       write(0,*)"icon OK",icon

      read(50,*)   nmax,  iemax
      if(nmax .gt. md1)  then
       write(0,*)   'nmax .gt. md1'
       write(0,*)   'nmax =', nmax
       write(0,*)   ' md1 =', md1
      stop
      end if
      if(iemax .gt. md3)  then
       write(0,*)   'iemax .gt. md2'
       write(0,*)   'iemax =', iemax
       write(0,*)   '  md3 =', md3
      stop
      end if

       write(0,*)"nmax OK",nmax
       write(0,*)"iemax OK",iemax
       write(0,*)"ir OK",ir

      read(50,*) (nod,  (xyz(io,nod),  io = 1, icoo),  & 
                 i = 1, nmax)
       write(0,*)"nod OK",nod
      read(50,*) (ie,  (nodc(in,ie),  in = 1, icon1), &
                   i = 1, iemax)

  100 format(i5,2f10.0)
  200 format(9i5)
       write(0,*)   'eps OK'
!     read(50,*)   eps
!     read(50,*)   eps

      return
      end

!--------------------------------------------------------------

      subroutine step01 (ir,     isw,    ie1,    idivx,  idivy,  & 
                      idivz,  itype,  nodx,   nody,   nodz,      & 
                      dsx,    dsy,    dsz)      

      implicit double precision (a-h , o-z)

      if(isw .eq. 1)  then
       read(50,*)   ie1,  idivx,  idivy,  itype
       idivz = 1
       nodz = 1             

      else if(isw .eq. 2)  then
       read(50,*)   ie1,  idivx,  idivy,  idummy
       idivz = 1
       nodz = 1
  
      else  
       read(50,*)   ie1,  idivx,  idivy,  idivz
       nodz = idivz + 1
       divz = float(idivz)
       dsz = 2.0 / divz
      end if
 
      nodx = idivx + 1
      nody = idivy + 1
      divx = float(idivx)
      divy = float(idivy)
      dsx = 2.0 / divx
      dsy = 2.0 / divy

      return
      end

!--------------------------------------------------------------

      subroutine step02 (isw,    nodx,   nody,   nodz,   nod2, & 
                       rn,     idivx,  idivy,  idivz,  icoo,   & 
                       icon1,  xyz1,   xyz2,   ie1,    nodc1,  & 
                       dsx,    dsy,    dsz,    md1,    md3)  

      implicit double precision (a-h , o-z)

      dimension     xyz1(3,md1),   xyz2(3,md1),   rn(20), & 
                  nodc1(20,md3)

      sz = - 1.0
      do 100  iz = 1, nodz
      sy = - 1.0
      do 200  iy = 1, nody
      sx = - 1.0
      do 300  ix = 1, nodx

      nod2 = nod2 + 1

      if(isw .eq. 1  .or.  isw .eq. 2)  then
      rn(1) = 0.25 * (1.0 - sx) * (1.0 - sy) * (- 1.0 - sx - sy)
      rn(3) = 0.25 * (1.0 + sx) * (1.0 - sy) * (- 1.0 + sx - sy)
      rn(5) = 0.25 * (1.0 + sx) * (1.0 + sy) * (- 1.0 + sx + sy)
      rn(7) = 0.25 * (1.0 - sx) * (1.0 + sy) * (- 1.0 - sx + sy)
      rn(2) = 0.5 * (1.0 - sx * sx) * (1.0 - sy)
      rn(4) = 0.5 * (1.0 - sy * sy) * (1.0 + sx)
      rn(6) = 0.5 * (1.0 - sx * sx) * (1.0 + sy)
      rn(8) = 0.5 * (1.0 - sy * sy) * (1.0 - sx)

      else
      end if

      do 10  io = 1, icoo
      xyz2(io,nod2) = 0.0
   10 continue

      do 11  in = 1, icon1
      do 11  io = 1, icoo
      inie1 = nodc1(in,ie1)
      xyz2(io,nod2) = xyz2(io,nod2) + xyz1(io,inie1) * rn(in)
   11 continue

      sx = sx + dsx
      if(ix .eq. idivx)  sx = 1.0
  300 continue
      sy = sy + dsy
      if(iy .eq. idivy)  sy = 1.0
  200 continue
      sz = sz + dsz
      if(iz .eq. idivz)  sy = 1.0
  100 continue

      return
      end
 
!--------------------------------------------------------------

      subroutine step03 (isw,  idivx,  idivy,  idivz,  idum,   &    
                   ie2,  nmax2,  idmax,  nodx,   itype,        &    
                   nodc2,   iemax2,  md2,   xyz2,   md1)

      implicit double precision (a-h , o-z)

      dimension     nodc2(8,md2),  xyz2(3,md1)
      do 100  iez = 1, idivz
      do 100  iey = 1, idivy
      inui = iey - iey / 2 * 2 
      do 100  iex = 1, idivx

      if(isw .eq. 1)  then
       idum = idum + 1
       ie2 = ie2 + 1
       nc1 = nmax2 + idum - idmax + iey - 1
       nc2 = nc1 + 1
       nc3 = nc2 + nodx 
       nc4 = nc3 - 1   
       if(itype .eq. 1)  then
        nodc2(1,ie2) = nc1
        nodc2(2,ie2) = nc3
        nodc2(3,ie2) = nc4
        ie2 = ie2 + 1
        nodc2(1,ie2) = nc3
        nodc2(2,ie2) = nc1
        nodc2(3,ie2) = nc2
       else if(itype .eq. 2)  then
        nodc2(1,ie2) = nc2
        nodc2(2,ie2) = nc4
        nodc2(3,ie2) = nc1
        ie2 = ie2 + 1
        nodc2(1,ie2) = nc4
        nodc2(2,ie2) = nc2
        nodc2(3,ie2) = nc3
       else if(itype .eq. 3)  then
        idum2 = idum - idmax - (iey - 1) * idivx
        idivx2 = idivx / 2
        if(idum2 .le. idivx2)  then
         nodc2(1,ie2) = nc1
         nodc2(2,ie2) = nc3
         nodc2(3,ie2) = nc4
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc3
         nodc2(2,ie2) = nc1
         nodc2(3,ie2) = nc2
        else
         nodc2(1,ie2) = nc2
         nodc2(2,ie2) = nc4
         nodc2(3,ie2) = nc1
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc4
         nodc2(2,ie2) = nc2
         nodc2(3,ie2) = nc3
        end if
       else if(itype .eq. 4)  then
        idum2 = idum - idmax - (iey - 1) * idivx
        idivx2 = idivx / 2
        if(idum2 .le. idivx2)  then
         nodc2(1,ie2) = nc2
         nodc2(2,ie2) = nc4
         nodc2(3,ie2) = nc1
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc4
         nodc2(2,ie2) = nc2
         nodc2(3,ie2) = nc3
        else
         nodc2(1,ie2) = nc1
         nodc2(2,ie2) = nc3
         nodc2(3,ie2) = nc4
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc3
         nodc2(2,ie2) = nc1
         nodc2(3,ie2) = nc2
        end if
       else if(itype .eq. 5)  then
        x1 = xyz2(1,nc1)
        x2 = xyz2(1,nc2)
        x3 = xyz2(1,nc3)
        x4 = xyz2(1,nc4)
        y1 = xyz2(2,nc1)
        y2 = xyz2(2,nc2)
        y3 = xyz2(2,nc3)
        y4 = xyz2(2,nc4)
        xy13 = (x1 - x3)**2 + (y1 - y3)**2
        xy24 = (x2 - x4)**2 + (y2 - y4)**2
        xy = xy24 - xy13
        if(xy .gt. 1e-5)  then
         nodc2(1,ie2) = nc1
         nodc2(2,ie2) = nc3
         nodc2(3,ie2) = nc4
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc3
         nodc2(2,ie2) = nc1
         nodc2(3,ie2) = nc2
        else
         nodc2(1,ie2) = nc2
         nodc2(2,ie2) = nc4
         nodc2(3,ie2) = nc1
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc4
         nodc2(2,ie2) = nc2
         nodc2(3,ie2) = nc3
        end if
       else if(itype .eq. 6)  then
        if(inui .eq. 1)  then
         nodc2(1,ie2) = nc1
         nodc2(2,ie2) = nc3
         nodc2(3,ie2) = nc4
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc3
         nodc2(2,ie2) = nc1
         nodc2(3,ie2) = nc2
         inui = 0
        else
!write(100,*) ie2
         nodc2(1,ie2) = nc2
         nodc2(2,ie2) = nc4
         nodc2(3,ie2) = nc1
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc4
         nodc2(2,ie2) = nc2
         nodc2(3,ie2) = nc3
         inui = 1
        end if
       else if(itype .eq. 7)  then
        x1 = xyz2(1,nc1)
        x2 = xyz2(1,nc2)
        x3 = xyz2(1,nc3)
        x4 = xyz2(1,nc4)
        y1 = xyz2(2,nc1)
        y2 = xyz2(2,nc2)
        y3 = xyz2(2,nc3)
        y4 = xyz2(2,nc4)
        xxx1 = (x1 - x3)**2 + (y1 - y3)**2
        xxx2 = (x2 - x4)**2 + (y2 - y4)**2
        xyxy = xxx2 - xxx1
        if(xyxy.gt.1e-05)  then
         nodc2(1,ie2) = nc1
         nodc2(2,ie2) = nc2
         nodc2(3,ie2) = nc3
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc3
         nodc2(2,ie2) = nc4
         nodc2(3,ie2) = nc1
        else
         nodc2(1,ie2) = nc1
         nodc2(2,ie2) = nc2
         nodc2(3,ie2) = nc4
         ie2 = ie2 + 1
         nodc2(1,ie2) = nc3
         nodc2(2,ie2) = nc4
         nodc2(3,ie2) = nc2
        end if
       end if

      else if(isw .eq. 2)  then
       ie2 = ie2 + 1
       nc1 = nmax2 + ie2 - iemax2 + iey - 1
       nc2 = nc1 + 1
       nc3 = nc2 + nodx 
       nc4 = nc3 - 1   
       nodc2(1,ie2) = nc1
       nodc2(2,ie2) = nc2
       nodc2(3,ie2) = nc3
       nodc2(4,ie2) = nc4

      else
      end if

  100 continue

      return
      end

!--------------------------------------------------------------

      subroutine step04 (nmax2,  nnod,   xyza,   xyzb,    xyz2,   &     
                 icoo,   eps,    kmax,   knod,   lmax,            &     
                 md1,    md4,    md5,    iw)

      implicit double precision (a-h , o-z)

      dimension     xyz2(3,md1),   nnod(md4,md1),  &  
                lmax(md1),     knod(md5),          &  
                xyza(3),       xyzb(3)

      km = 0
      do 100  nod1 = 1, nmax2-1

      lm = 0
      do 200  l = 1, md4 
      nnod(l,nod1) = 0
  200 continue

      do 300  nod2 = nod1+1, nmax2

      do 400  io = 1, icoo
      xyza(io) = xyz2(io,nod1)
      xyzb(io) = xyz2(io,nod2)
  400 continue

      do 500  io = 1, icoo
      xyzab = abs(xyza(io) - xyzb(io))
      if(xyzab .gt. eps)  goto 300  
  500 continue

      lm = lm + 1
      km = km + 1
      nnod(lm,nod1) = nod2
      knod(km)      = nod2

  300 continue

      lmax(nod1) = lm
      if(lm .gt. md4)  then     
       write(0,*)   'LMAX .gt. MD4'
       write(0,*)   'lmax =', lm
       write(0,*)   ' md4 =', md4
      stop     
      end if

      if(lm .eq. 0)  goto 100
!      write(0,*)   nod1,  (nnod(l,nod1),  l = 1, lm)

  100 continue
      kmax = km
      if(kmax .gt. md5)   then
       write(0,*)   'KMAX .gt. MD5'
       write(0,*)   'kmax =',kmax
       write(0,*)   ' md5 =',md5
      stop
      end if

      return
      end

!--------------------------------------------------------------

      subroutine step05 (nmax2,  lmax,   nnod,   icoo,   md4,  & 
                      iemax2, icon,   nodc2,  md1,    md2)
 
      implicit double precision (a-h , o-z)

      dimension     nodc2(8,md2),  lmax(md1),    nnod(md4,md1)

      do 100  nod1 = 1, nmax2
      lm = lmax(nod1)

      do 200  l = 1, lm
      nod2 = nnod(l,nod1)

      do 300  ie = 1, iemax2
      do 400  in = 1, icon
      if(nod2 .ne. nodc2(in,ie))  goto 400 
      nodc2(in,ie) = nod1
  400 continue
  300 continue

  200 continue
  100 continue

      return
      end

!--------------------------------------------------------------

      subroutine step06 (nmax2,  kmax,   knod,   icoo,   xyz1,     &   
                   xyz2,   iemax2, icon,   nodc2,  nmax1,          &   
                   md1,    md2,    md5)

      implicit double precision (a-h , o-z)

      dimension     knod(md5),     xyz1(3,md1),   xyz2(3,md2),  & 
                 nodc2(8,md2) 

      nod1 = 0
      do 100  nod2 = 1, nmax2

      do 200  km = 1, kmax
      if(nod2 .eq. knod(km))  goto 400 
  200 continue

      nod1 = nod1 + 1
      do 300  io = 1, icoo
      xyz1(io,nod1) = xyz2(io,nod2)
  300 continue
    
  400 continue
      do 500  ie = 1, iemax2
      do 600  in = 1, icon
      if(nod2 .ne. nodc2(in,ie))  goto 600
      nodc2(in,ie) = nod1
      goto 500 
  600 continue
  500 continue

  100 continue
      nmax1 = nod1
 
      return
      end

!--------------------------------------------------------------
!cube72.mesへの出力文
      subroutine oumesh (iw,     isw,    nmax,   iemax,  xyz,  &   
                    nodc,   icoo,   icon,   md1,    md2)
 
      implicit double precision (a-h , o-z)

      dimension     xyz(3,md1),    nodc(8,md2)

!こいつが使われている
      if(isw .eq. 1)  then
      write(60,603)   nmax,  iemax
      write(60,602)   (nod,  (xyz(io,nod),  io = 1, icoo),  & 
                    nod = 1, nmax)
      write(60,603)   (ie,  (nodc(in,ie),  in = 1, icon),   & 
                   ie = 1, iemax)
!     write(60,700)   (i,i= 1,iemax) 

      else if(isw .eq. 2)  then
      write(60,603)   nmax,  iemax
      write(60,601)   (nod,  (xyz(io,nod),  io = 1, icoo),   & 
                    nod = 1, nmax)
      write(60,613)   (ie,(nodc(in,ie),  in = 1, icon),  & 
                    ie = 1, iemax)

      else
      write(60,604)   nmax,  iemax
      write(60,602)   (nod,  (xyz(io,nod),  io = 1, icoo),   &
                     nod = 1, nmax)
      write(60,603)   ((nodc(in,ie),  in = 1, icon), &
                      ie = 1, iemax)
      end if

  601 format(i9,2f15.6)
  602 format(i9,2f15.5)
  603 format(4i9)
  613 format(5i9)
  604 format(2i9)
  710 format(4i9)
  700 format(2(i9,'        2.1e+04','       0.3'))

      return
      end

!--------------------------------------------------------------

      subroutine echeck (nmax,   iemax,  xyz,    nodc,   icoo, & 
                       icon,   isw,    md1,    md2)

      implicit double precision (a-h , o-z)

      dimension     xyz(3,md1),    nodc(8,md2)

      if(isw .eq. 1)  then
      do 100  ie = 1, iemax
      nod1 = nodc(1,ie)
      nod2 = nodc(2,ie)
      nod3 = nodc(3,ie)
      x1 = xyz(1,nod1)
      x2 = xyz(1,nod2)
      x3 = xyz(1,nod3)
      y1 = xyz(2,nod1)
      y2 = xyz(2,nod2)
      y3 = xyz(2,nod3)
      area = (y1 - y3) * (x1 - x2) - (x1 - x3) * (y1 - y2)
      if(area .le. 0.0)  then
       write(0,*)   'area =',area  
       write(0,*)   '  ie =',ie
       write(0,*)   'nod1 =',nod1,x1,y1
       write(0,*)   'nod2 =',nod2,x2,y2
       write(0,*)   'nod3 =',nod3,x3,y3
      stop
      end if
  100 continue
      else if(isw .eq. 2)  then
      do 200  ie = 1, iemax
      nod1 = nodc(1,ie)
      nod2 = nodc(2,ie)
      nod3 = nodc(3,ie)
      nod4 = nodc(4,ie)
      x1 = xyz(1,nod1)
      x2 = xyz(1,nod2)
      x3 = xyz(1,nod3)
      x4 = xyz(1,nod4)
      y1 = xyz(2,nod1)
      y2 = xyz(2,nod2)
      y3 = xyz(2,nod3)
      y4 = xyz(2,nod4)
      xa1 = x1
      xa2 = x2
      xa3 = x4
      ya1 = y1
      ya2 = y2
      ya3 = y4
      xb1 = x3
      xb2 = x4
      xb3 = x2
      yb1 = y3
      yb2 = y4
      yb3 = y2
      areaa=(ya1 - ya3)*(xa1 - xa2)-(xa1 - xa3)*(ya1 - ya2)
      areab=(yb1 - yb3)*(xb1 - xb2)-(xb1 - xb3)*(yb1 - yb2)
      if(areaa .le. 0.0  .or.  areab .le. 0.0)  then
       write(0,*)   'area .le. 0.0'
       write(0,*)   'ie =',ie
      stop
      end if
  200 continue
      else
      end if

      write(0,*)  'element check  OK '

      return
      end

!--------------------------------------------------------------

      subroutine renumb (xyz,    nodc,   nmax,   iemax,   md1,  & 
                      md2,    isw,    icoo,   icon)

      implicit double precision (a-h , o-z)

      dimension     xyz(3,md1),    nodc(8,md2)

      if(isw .eq. 3)  then

       do 10  nod1 =1, nmax-1

        nodold = nod1
        do 20  nod2 = nod1+1, nmax
        xyz3 = abs(xyz(3,nodold) - xyz(3,nod2))
        if(xyz3 .lt. 1e-5)  then
         xyz2 = abs(xyz(2,nodold) - xyz(2,nod2))
         if(xyz2 .lt. 1e-5)  then
          if(xyz(1,nod2) .lt. xyz(1,nodold))  then
           nodold=nod2
          end if
         else if(xyz(2,nod2) .lt. xyz(2,nodold))  then
          nodold=nod2
         end if
        else if(xyz(3,nod2).lt.xyz(3,nodold)) then
         nodold=nod2
        end if
   20   continue

        if(nodold.eq.nod1) go to 10

        do 25  io = 1, icoo
         xyzold        = xyz(io,nod1)
         xyz(io,nod1)   = xyz(io,nodold)
         xyz(io,nodold) = xyzold
   25   continue 

        do 30  ie = 1, iemax
         do 40  in = 1, icon
          if(nodc(in,ie) .eq. nod1)  then
           nodc(in,ie) = nodold
          else if(nodc(in,ie) .eq. nodold)  then
           nodc(in,ie) = nod1
          end if
   40    continue
   30   continue

   10  continue

      else

       do 100  nod1 =1, nmax-1

        nodold = nod1
        do 200  nod2 = nod1+1, nmax
        xyz1 = abs(xyz(1,nodold) - xyz(1,nod2))
        if(xyz1 .lt. 1e-5)  then
         if(xyz(2,nod2) .gt. xyz(2,nodold))  then
          nodold=nod2
         end if
        else if(xyz(1,nod2) .lt. xyz(1,nodold))  then
         nodold=nod2
        end if
  200   continue

        if(nodold.eq.nod1) go to 100

        do 250  io = 1, icoo
         xyzold        = xyz(io,nod1)
         xyz(io,nod1)   = xyz(io,nodold)
         xyz(io,nodold) = xyzold
  250   continue 

        do 300  ie = 1, iemax
         do 400  in = 1, icon
          if(nodc(in,ie) .eq. nod1)  then
           nodc(in,ie) = nodold
          else if(nodc(in,ie) .eq. nodold)  then
           nodc(in,ie) = nod1
          end if
  400    continue
  300   continue

  100  continue
  
      end if

      return
      end
