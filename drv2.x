      subroutine drv2

      include "real8.h"
      character jti*1,kfi*1,dummy*132,fname*132
      include "basicd.h"
      include "const.h"
      include "files.h"
      include 'cdkwr.h'

c  hobi  height of burst for evaluation (ft)
c  cep   circular error probable for weapon (km)
c  r95   95% of damage radius lies within this area (km)
c  gname file name containing ground range, altitude and yield data
c  ivnb  starting vn #
c  ivne  ending vn #
c  ivns  step in vn #
c  jti   target type
c  kfb   k-factor begin value (typically 0)
c        the end value will be min(ivn-1,9)
c  kfs   k-factor step (typicall 1)
c  az    azimuth iin degrees from dgz to target

      namelist /plst/ ivnb,ivne,ivns,jti,kfb,kfs,yld,hobi,r95,cep,az

      namelist /covlst/ sig_cr,sig_dr,sig_hob,
     *                  rho_cr_dr,rho_cr_hob,rho_dr_hob

      lin  =  1
      lout =  6
      ldbg = 62
      ldbg = -1

      call acon

      ivnb =  1 
      ivne = 30
      ivns =  1

      jti  = 'p'

      kfb  = 0
      kfs  = 1

      yld  = 1000.0d0
      hobi = 0.0d0

      r95  = 5.0d0
      cep  = 0.0d0

      az   = 0.0d0

      iflg = 2

      sig_cr    = 0.0d0
      sig_dr    = 0.0d0
      sig_hob   = 0.0d0
      rho_cr_dr = 0.0d0
      rho_cr_hob= 0.0d0
      rho_dr_hob= 0.0d0

      fname = 'drv2.nml'

      do iarg=1,iargc()
         jarg = iarg + 1

         call getarg(iarg,dummy)

         dummy = dummy(1:lnblnk(dummy))

         if (dummy.eq."-i") then
            call getarg(jarg,fname)
         endif
      enddo

      open (unit=lin,status='old',file=fname)
      read(lin,nml=plst)
      read(lin,nml=covlst,err=997,end=997)
 997  continue
      close (unit=lin)

c  convert sig_cr, sig_dr from km to feet
      sig_cr = sig_cr * ckm2ft
      sig_dr = sig_dr * ckm2ft

      open (unit=lout,status='unknown',file='lout.out')

      if ((ldbg.gt.0).and.(ldbg.ne.lout)) then
         open (unit=ldbg,status='unknown',file='dbg.out')
      endif

      open (unit=10,status='unknown',file='drv2.out')

      do ivn=ivnb,ivne,ivns
         do kf=kfb,min(ivn-1,9),kfs
            kfi = char(kf + 48)

            d   = 0.0d0
            wr  = 0.0d0
            pod = 0.0d0

            if (sig_cr .gt. 0.0d0 .or. sig_dr .gt. 0.0d0) then
               call pdcov(ivn,jti,kfi,yld,hobi,r95,sig_cr,sig_dr,
     *                    sig_hob,rho_cr_dr,rho_cr_hob,rho_dr_hob,
     *                    d,wr,pod,iflg,az)
            else
               call pdcalc(ivn,jti,kfi,yld,hobi,r95,cep,
     *                     d,wr,pod,iflg,az)
            endif

            write(10,20)ivn,kf,avn,d,wr,pod
         enddo
         write(10,*)
      enddo
      close (unit=10)

      iflg = 6

      open (unit=11,status='unknown',file='iflg6.out')

      do ivn=ivnb,ivne,ivns
         do kf=kfb,min(ivn-1,9),kfs
             kfi = char(kf + 48)

            if (kf.eq.kfb)icnt = 0
            podi = 0.01
            do while (podi.le.1.0d0)

               pod = podi
               if (sig_cr .gt. 0.0d0 .or. sig_dr .gt. 0.0d0) then
                  call pdcov(ivn,jti,kfi,yld,hobi,r95,sig_cr,sig_dr,
     *                       sig_hob,rho_cr_dr,rho_cr_hob,rho_dr_hob,
     *                       d,wr,pod,iflg,az)
               else
                  call pdcalc(ivn,jti,kfi,yld,hobi,r95,cep,
     *                        d,wr,pod,iflg,az)
               endif

               write(11,30)ivn,kf,avn,d,wr,pod,icnt
               icnt = icnt + 1

               podi = podi + 0.01
            enddo
            write(11,*)''
         enddo
         write(11,*)''
      enddo
      close (unit=11)

      close (unit=lout)
      if ((ldbg.gt.0).and.(ldbg.ne.lout))close(unit = ldbg)

      return

 20   format(2(1x,i3),4(1x,f15.5))
 30   format(2(1x,i3),4(1x,f15.5),1x,i5)

      end
