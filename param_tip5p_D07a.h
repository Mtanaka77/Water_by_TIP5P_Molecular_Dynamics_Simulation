!  param_tip5p_D07a.h
!
      integer(C_INT) np0,nq0,npq0,npq5,nq1,            &
                     ip0,mesh,mintpol,mx,my,mz,        &
                     mx1,my1,mz1,kstart,brillouin
      logical        if_xyz1,if_xyz2,if_obsv
!
! /home2, or /lv01 for photon
!     character     praefixs*29,praefixc*24,praefixe*24, &
!                   praefixi*24,suffix2*2,suffix1*2,suffix0*1
      character*4   praefixs*30,praefixc*30,praefixe*30, &
                    praefixi*30,suffix2*2,suffix1*2,suffix0*1
!
      real(C_DOUBLE) epsilon
      real(C_DOUBLE) t_init,t_wipe_sta,t_wipe_end
      parameter  (epsilon=88.d0/88.d0)  ! 273 K, eps=88. 
      parameter  (t_init=1000.d0,t_wipe_sta=1700.d0,  &
                  t_wipe_end=4700.d0)
!
      parameter  (kstart=0,suffix2='0a', & ! 120a is created
                           suffix1='0a', & ! TIP701__0 
                           suffix0='1')    ! 
!     parameter  (kstart=2,suffix2='0b', & ! 120b, kstart=2
!                          suffix1='0a', & ! TIP701__1
!                          suffix0='1')    !        +++ 
!     parameter  (kstart=2,suffix2='0c', & ! 120d, kstart=2
!                          suffix1='0b', & ! TIP701__1
!                          suffix0='1')    ! 
! /home2, /lv01                          +++++++ short 
!     parameter (praefixs='/home/tanakam/MPI_wat5/TIP507', & ! LXwat3
!                praefixi='/data/sht/tanakam/tip507',      &
!                praefixc='/data/sht/tanakam/tip507',      &
!                praefixe='/data/sht/tanakam/tip507')
      parameter (praefixs='/home/mtanaka/MPI_wat5A/TIP707', & 
                 praefixi='/home/mtanaka/MPI_wat5A/tip707', &
                 praefixc='/home/mtanaka/MPI_wat5A/tip707', &
                 praefixe='/home/mtanaka/MPI_wat5A/tip707')
      parameter (if_xyz1=.false., if_xyz2=.true.,          &
                 if_obsv=.false.)  !! if .false,, save to a file 
!
! moldyn  L > 32 Ang
!  Total number is nq0; for mh3: nq0=6210,np0=216
!                              +++++ null
!     parameter  (nq0=6210,np0=216) ! 6210 atoms by 5-body CS1
      parameter  (nq0=8640,np0=0)   ! 8640 atoms by 5-body CS1
!
      parameter  (npq5=nq0+np0)
      parameter  (npq0=nq0/5+np0)
      parameter  (nq1=nq0/5)
!
      parameter  (ip0=3,mesh=32,mintpol=4*50048)
      parameter  (brillouin=1)
      parameter  (mx=mesh,my=mesh,mz=mesh,mx1=mx+1,my1=my+1,mz1=mz+1)
!
!  forces.h
!
      real(C_DOUBLE) AA1,AA2,AA3,AA4,AA5,PPP,sqrtpi
      parameter ( AA1 = 0.254829592d0, AA2 = -0.284496736d0 )
      parameter ( AA3 = 1.421413741d0, AA4 = -1.453152027d0 )
      parameter ( AA5 = 1.061405429d0, PPP =  0.3275911d0   )
      parameter ( sqrtpi = sqrt(3.14159265358979323846d0)   ) 

