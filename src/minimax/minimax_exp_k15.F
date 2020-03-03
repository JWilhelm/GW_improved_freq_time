!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2020  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Routines to calculate the minimax coefficients in order to
!>        approximate 1/x as a sum over exponential functions
!>        1/x ~ SUM_{i}^{K} w_i EXP(-a_i * x) for x belonging to [1:Rc].
!>        This module contains coefficients for minimax approximations with 1 <= k <= 15.
!> \par History
!>      11.2012 created [Mauro Del Ben]
!>      02.2015 changed module name minimax to minimax_k15 to replace this module
!>              with new minimax module, added function to get minimax error
!>              [Patrick Seewald]
! **************************************************************************************************
MODULE minimax_exp_k15
   USE kinds,                           ONLY: dp
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'minimax_exp_k15'

   PUBLIC :: get_minimax_coeff_k15, check_range_k15, get_minimax_numerical_error, get_minimax_coeff_k35

CONTAINS

! **************************************************************************************************
!> \brief Check that the range for the minimax approximation is not too small
!>        for the chosen number of integration point, or that the number of
!>        integration point is not exceeding the maximum allowed by the actual
!>        implementation.
!> \param k ...
!> \param Rc ...
!> \param ierr ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE check_range_k15(k, Rc, ierr)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      INTEGER                                            :: ierr

      REAL(KIND=dp)                                      :: R_min

      ierr = 0
      R_min = 0.0_dp

      SELECT CASE (k)
      CASE (1)
         R_min = 2.0000000000000000E+00_dp
      CASE (2)
         R_min = 2.0000000000000000E+00_dp
      CASE (3)
         R_min = 2.0000000000000000E+00_dp
      CASE (4)
         R_min = 3.0000000000000000E+00_dp
      CASE (5)
         R_min = 4.0000000000000000E+00_dp
      CASE (6)
         R_min = 5.0000000000000000E+00_dp
      CASE (7)
         R_min = 7.0000000000000000E+00_dp
      CASE (8)
         R_min = 1.0000000000000000E+01_dp
      CASE (9)
         R_min = 1.0000000000000000E+01_dp
      CASE (10)
         R_min = 1.0000000000000000E+01_dp
      CASE (11)
         R_min = 1.0000000000000000E+01_dp
      CASE (12)
         R_min = 1.0000000000000000E+01_dp
      CASE (13)
         R_min = 1.0000000000000000E+01_dp
      CASE (14)
         R_min = 1.0000000000000000E+01_dp
      CASE (15)
         R_min = 2.0000000000000000E+01_dp
      CASE DEFAULT
         ierr = -1
      END SELECT

      IF (Rc < R_min .AND. ierr == 0) ierr = 1

   END SUBROUTINE check_range_k15

! **************************************************************************************************
!> \brief Sample numerical error and return its maximum.
!> \param Rc ...
!> \param aw ...
!> \return ...
!> \author Patrick Seewald
! **************************************************************************************************
   FUNCTION get_minimax_numerical_error(Rc, aw) RESULT(numerr)
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(:), INTENT(IN)            :: aw
      REAL(KIND=dp)                                      :: numerr

      INTEGER, PARAMETER                                 :: npoints = 100

      INTEGER                                            :: j, k, np
      REAL(KIND=dp)                                      :: mm, x

      k = SIZE(aw)/2
      numerr = 0.0_dp
      DO np = 1, npoints
         x = 1.0_dp + (Rc - 1.0_dp)*(npoints - np)/(npoints - 1)
         mm = 0.0_dp
         DO j = 1, k
            mm = mm + aw(j + k)*EXP(-aw(j)*x)
         ENDDO
         numerr = MAX(numerr, ABS(mm - 1.0_dp/x))
      ENDDO

   END FUNCTION get_minimax_numerical_error

! **************************************************************************************************
!> \brief Get minimax coefficients: k15 implementation (coefficients up to k=15 terms).
!>        All a_i and w_i have been fitted with a 12th order polynomial as
!>        a function of Rc for each value of K.
!> \param k ...
!> \param Rc ...
!> \param aw ...
!> \param mm_error ...
!> \author Mauro Del Ben
! **************************************************************************************************
   SUBROUTINE get_minimax_coeff_k15(k, Rc, aw, mm_error)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: Rc
      REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
      REAL(KIND=dp), INTENT(OUT), OPTIONAL               :: mm_error

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_minimax_coeff_k15', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icoef, kkk
      REAL(KIND=dp)                                      :: L_b, R_max, R_min, R_transf, U_b
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef

      CALL timeset(routineN, handle)
      R_min = 0.0_dp
      R_max = 0.0_dp

      SELECT CASE (k)
      CASE (1)
         R_min = 2.0000000000000000E+00_dp
         R_max = 8.6669999999999998E+00_dp
      CASE (2)
         R_min = 2.0000000000000000E+00_dp
         R_max = 4.1539999999999999E+01_dp
      CASE (3)
         R_min = 2.0000000000000000E+00_dp
         R_max = 1.4680000000000001E+02_dp
      CASE (4)
         R_min = 3.0000000000000000E+00_dp
         R_max = 4.3610000000000002E+02_dp
      CASE (5)
         R_min = 4.0000000000000000E+00_dp
         R_max = 1.1540000000000000E+03_dp
      CASE (6)
         R_min = 5.0000000000000000E+00_dp
         R_max = 2.8070000000000000E+03_dp
      CASE (7)
         R_min = 7.0000000000000000E+00_dp
         R_max = 6.3730000000000000E+03_dp
      CASE (8)
         R_min = 1.0000000000000000E+01_dp
         R_max = 1.3749000000000000E+04_dp
      CASE (9)
         R_min = 1.0000000000000000E+01_dp
         R_max = 2.8387000000000000E+04_dp
      CASE (10)
         R_min = 1.0000000000000000E+01_dp
         R_max = 5.6502000000000000E+04_dp
      CASE (11)
         R_min = 1.0000000000000000E+01_dp
         R_max = 1.0890000000000000E+05_dp
      CASE (12)
         R_min = 1.0000000000000000E+01_dp
         R_max = 2.0420000000000000E+05_dp
      CASE (13)
         R_min = 1.0000000000000000E+01_dp
         R_max = 3.7370000000000000E+05_dp
      CASE (14)
         R_min = 1.0000000000000000E+01_dp
         R_max = 6.6910000000000000E+05_dp
      CASE (15)
         R_min = 2.0000000000000000E+01_dp
         R_max = 1.1750000000000000E+06_dp
      CASE DEFAULT
         CPABORT("")
      END SELECT

      fit_coef = 0.0
      aw = 0.0
      IF (Rc < R_min) THEN
         CPABORT("")
      ELSE
         SELECT CASE (k)
         CASE (1)
            IF (Rc >= R_max) THEN
               aw(1) = 4.4649260629947779E-01_dp
               aw(2) = 1.4290997869805788E+00_dp
            ELSE
               CALL get_fit_coef_k1(k, fit_coef, L_b, U_b)
            END IF
         CASE (2)
            IF (Rc >= R_max) THEN
               aw(1) = 9.3106172568836743E-02_dp
               aw(2) = 1.0301539335400713E+00_dp
               aw(3) = 2.9765554323276328E-01_dp
               aw(4) = 1.9917443679019897E+00_dp
            ELSE
               CALL get_fit_coef_k2(k, fit_coef, L_b, U_b)
            END IF
         CASE (3)
            IF (Rc >= R_max) THEN
               aw(1) = 2.6349272848143715E-02_dp
               aw(2) = 2.8996855484815137E-01_dp
               aw(3) = 1.5533661852788561E+00_dp
               aw(4) = 8.4233696465358193E-02_dp
               aw(5) = 5.5403643315333395E-01_dp
               aw(6) = 2.3557688429224690E+00_dp
            ELSE
               CALL get_fit_coef_k3(k, fit_coef, L_b, U_b)
            END IF
         CASE (4)
            IF (Rc >= R_max) THEN
               aw(1) = 8.8685367946663256E-03_dp
               aw(2) = 9.7565835027097439E-02_dp
               aw(3) = 5.1691899901700877E-01_dp
               aw(4) = 2.0308311260118175E+00_dp
               aw(5) = 2.8350982414551848E-02_dp
               aw(6) = 1.8628933861057068E-01_dp
               aw(7) = 7.6510507307920383E-01_dp
               aw(8) = 2.6286812528458148E+00_dp
            ELSE
               CALL get_fit_coef_k4(k, fit_coef, L_b, U_b)
            END IF
         CASE (5)
            IF (Rc >= R_max) THEN
               aw(1) = 3.3529112764591387E-03_dp
               aw(2) = 3.6885523808630270E-02_dp
               aw(3) = 1.9524288957651931E-01_dp
               aw(4) = 7.5508129953934144E-01_dp
               aw(5) = 2.4730352459564653E+00_dp
               aw(6) = 1.0718601062497496E-02_dp
               aw(7) = 7.0423935703868623E-02_dp
               aw(8) = 2.8840120533084923E-01_dp
               aw(9) = 9.4333101049811663E-01_dp
               aw(10) = 2.8485313369244247E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k5_sh(k, fit_coef, L_b, U_b)
               ELSE
                  CALL get_fit_coef_k5_lo(k, fit_coef, L_b, U_b)
               END IF
            END IF
         CASE (6)
            IF (Rc >= R_max) THEN
               aw(1) = 1.3802019452414500E-03_dp
               aw(2) = 1.5183614381572768E-02_dp
               aw(3) = 8.0361511855253240E-02_dp
               aw(4) = 3.1025560863277402E-01_dp
               aw(5) = 9.9661148626887119E-01_dp
               aw(6) = 2.8871396539466740E+00_dp
               aw(7) = 4.4122353526579831E-03_dp
               aw(8) = 2.8989220321940696E-02_dp
               aw(9) = 1.1867791139299246E-01_dp
               aw(10) = 3.8613951254853507E-01_dp
               aw(11) = 1.0973265611814316E+00_dp
               aw(12) = 3.0334106111925840E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k6_sh(k, fit_coef, L_b, U_b)
               ELSE
                  CALL get_fit_coef_k6_lo(k, fit_coef, L_b, U_b)
               END IF
            END IF
         CASE (7)
            IF (Rc >= R_max) THEN
               aw(1) = 6.0680046054254941E-04_dp
               aw(2) = 6.6754145503086018E-03_dp
               aw(3) = 3.5330085311658307E-02_dp
               aw(4) = 1.3636790644450009E-01_dp
               aw(5) = 4.3692597185507065E-01_dp
               aw(6) = 1.2378541669818244E+00_dp
               aw(7) = 3.2781922981409610E+00_dp
               aw(8) = 1.9398222440315960E-03_dp
               aw(9) = 1.2744981074293062E-02_dp
               aw(10) = 5.2173800858709714E-02_dp
               aw(11) = 1.6963212116290538E-01_dp
               aw(12) = 4.7830405076036075E-01_dp
               aw(13) = 1.2328765453712285E+00_dp
               aw(14) = 3.1933984156189350E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k7_sh(k, fit_coef, L_b, U_b)
               ELSE
                  CALL get_fit_coef_k7_lo(k, fit_coef, L_b, U_b)
               END IF
            END IF
         CASE (8)
            IF (Rc >= R_max) THEN
               aw(1) = 2.8127089949768415E-04_dp
               aw(2) = 3.0942621869177655E-03_dp
               aw(3) = 1.6376552972548756E-02_dp
               aw(4) = 6.3208135461484877E-02_dp
               aw(5) = 2.0243640179262407E-01_dp
               aw(6) = 5.7158548817956678E-01_dp
               aw(7) = 1.4769970948293096E+00_dp
               aw(8) = 3.6498865034832644E+00_dp
               aw(9) = 8.9916798455210248E-04_dp
               aw(10) = 5.9076940528835931E-03_dp
               aw(11) = 2.4183992214449269E-02_dp
               aw(12) = 7.8619729999942961E-02_dp
               aw(13) = 2.2140056400660871E-01_dp
               aw(14) = 5.6481870666490019E-01_dp
               aw(15) = 1.3539633386811780E+00_dp
               aw(16) = 3.3347117690388681E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k8_sh(k, fit_coef, L_b, U_b)
               ELSE
                  CALL get_fit_coef_k8_lo(k, fit_coef, L_b, U_b)
               END IF
            END IF
         CASE (9)
            IF (Rc >= R_max) THEN
               aw(1) = 1.3620787587729222E-04_dp
               aw(2) = 1.4984233172990983E-03_dp
               aw(3) = 7.9304843136686811E-03_dp
               aw(4) = 3.0608852526075443E-02_dp
               aw(5) = 9.8023317512591052E-02_dp
               aw(6) = 2.7660071140675380E-01_dp
               aw(7) = 7.1177472290822164E-01_dp
               aw(8) = 1.7131253759530025E+00_dp
               aw(9) = 4.0050154343014466E+00_dp
               aw(10) = 4.3542990561554510E-04_dp
               aw(11) = 2.8608520316584030E-03_dp
               aw(12) = 1.1711291382515341E-02_dp
               aw(13) = 3.8071387509506595E-02_dp
               aw(14) = 1.0718758733492406E-01_dp
               aw(15) = 2.7293433084542407E-01_dp
               aw(16) = 6.4600134991152369E-01_dp
               aw(17) = 1.4634218815677373E+00_dp
               aw(18) = 3.4614670015318758E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k9_sh(k, fit_coef, L_b, U_b)
               ELSE
                  CALL get_fit_coef_k9_lo(k, fit_coef, L_b, U_b)
               END IF
            END IF
         CASE (10)
            IF (Rc >= R_max) THEN
               aw(1) = 6.8443901388305083E-05_dp
               aw(2) = 7.5295159558786751E-04_dp
               aw(3) = 3.9850355998364926E-03_dp
               aw(4) = 1.5380799076269465E-02_dp
               aw(5) = 4.9255465262658012E-02_dp
               aw(6) = 1.3897073828118570E-01_dp
               aw(7) = 3.5731192893101876E-01_dp
               aw(8) = 8.5579284635745334E-01_dp
               aw(9) = 1.9457922045681717E+00_dp
               aw(10) = 4.3457510815955152E+00_dp
               aw(11) = 2.1880174937623183E-04_dp
               aw(12) = 1.4375664500890244E-03_dp
               aw(13) = 5.8848743519528662E-03_dp
               aw(14) = 1.9130624082222252E-02_dp
               aw(15) = 5.3858487442246458E-02_dp
               aw(16) = 1.3708846821179432E-01_dp
               aw(17) = 3.2364445020328991E-01_dp
               aw(18) = 7.2229153747049757E-01_dp
               aw(19) = 1.5633339391232095E+00_dp
               aw(20) = 3.5765361929512172E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k10_sh(k, fit_coef, L_b, U_b)
               ELSE
                  CALL get_fit_coef_k10_lo(k, fit_coef, L_b, U_b)
               END IF
            END IF
         CASE (11)
            IF (Rc >= R_max) THEN
               aw(1) = 3.5503484659246813E-05_dp
               aw(2) = 3.9057395758855312E-04_dp
               aw(3) = 2.0671330057150483E-03_dp
               aw(4) = 7.9783846601041101E-03_dp
               aw(5) = 2.5549886450088238E-02_dp
               aw(6) = 7.2085122754821276E-02_dp
               aw(7) = 1.8530557881445084E-01_dp
               aw(8) = 4.4334723019604549E-01_dp
               aw(9) = 1.0024332764988158E+00_dp
               aw(10) = 2.1748065663060481E+00_dp
               aw(11) = 4.6738220578828598E+00_dp
               aw(12) = 1.1349768781191745E-04_dp
               aw(13) = 7.4570001674067511E-04_dp
               aw(14) = 3.0526245415056301E-03_dp
               aw(15) = 9.9235004487624248E-03_dp
               aw(16) = 2.7937353788559310E-02_dp
               aw(17) = 7.1104122225773142E-02_dp
               aw(18) = 1.6777050240847788E-01_dp
               aw(19) = 3.7320560837397676E-01_dp
               aw(20) = 7.9414366264691305E-01_dp
               aw(21) = 1.6552701124330567E+00_dp
               aw(22) = 3.6820059233978353E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k11_sh(k, fit_coef, L_b, U_b)
               ELSE
                  IF (Rc <= 1000.0_dp) THEN
                     CALL get_fit_coef_k11_mi(k, fit_coef, L_b, U_b)
                  ELSE
                     CALL get_fit_coef_k11_lo(k, fit_coef, L_b, U_b)
                  END IF
               END IF
            END IF
         CASE (12)
            IF (Rc >= R_max) THEN
               aw(1) = 1.8934012296256358E-05_dp
               aw(2) = 2.0829314615342741E-04_dp
               aw(3) = 1.1024023520154727E-03_dp
               aw(4) = 4.2548734477736863E-03_dp
               aw(5) = 1.3625745873543612E-02_dp
               aw(6) = 3.8442710765665453E-02_dp
               aw(7) = 9.8818287645666872E-02_dp
               aw(8) = 2.3636401719642108E-01_dp
               aw(9) = 5.3373454579640744E-01_dp
               aw(10) = 1.1508219943995057E+00_dp
               aw(11) = 2.4001213099764143E+00_dp
               aw(12) = 4.9906311509941723E+00_dp
               aw(13) = 6.0528329999281435E-05_dp
               aw(14) = 3.9768193515047327E-04_dp
               aw(15) = 1.6279651361197933E-03_dp
               aw(16) = 5.2922029903677377E-03_dp
               aw(17) = 1.4898953114611908E-02_dp
               aw(18) = 3.7918950606510460E-02_dp
               aw(19) = 8.9457745237352296E-02_dp
               aw(20) = 1.9884427477508254E-01_dp
               aw(21) = 4.2144796502526766E-01_dp
               aw(22) = 8.6198574283723861E-01_dp
               aw(23) = 1.7404432216810901E+00_dp
               aw(24) = 3.7794416490913081E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k12_sh(k, fit_coef, L_b, U_b)
               ELSE
                  IF (Rc <= 1000.0_dp) THEN
                     CALL get_fit_coef_k12_mi(k, fit_coef, L_b, U_b)
                  ELSE
                     CALL get_fit_coef_k12_lo(k, fit_coef, L_b, U_b)
                  END IF
               END IF
            END IF
         CASE (13)
            IF (Rc >= R_max) THEN
               aw(1) = 1.0347391661194427E-05_dp
               aw(2) = 1.1383169144875790E-04_dp
               aw(3) = 6.0246014985531054E-04_dp
               aw(4) = 2.3252777266436543E-03_dp
               aw(5) = 7.4464345169693259E-03_dp
               aw(6) = 2.1008804344928592E-02_dp
               aw(7) = 5.4003219260973512E-02_dp
               aw(8) = 1.2916211533151045E-01_dp
               aw(9) = 2.9156437505279720E-01_dp
               aw(10) = 6.2769461815368310E-01_dp
               aw(11) = 1.3003150067055842E+00_dp
               aw(12) = 2.6217712688276107E+00_dp
               aw(13) = 5.2973357440769391E+00_dp
               aw(14) = 3.3078584812160216E-05_dp
               aw(15) = 2.1733219160086539E-04_dp
               aw(16) = 8.8967889768721382E-04_dp
               aw(17) = 2.8921756053065669E-03_dp
               aw(18) = 8.1422349627770872E-03_dp
               aw(19) = 2.0722487658696522E-02_dp
               aw(20) = 4.8886434214697162E-02_dp
               aw(21) = 1.0864184956111259E-01_dp
               aw(22) = 2.3003323004923831E-01_dp
               aw(23) = 4.6829507330399833E-01_dp
               aw(24) = 9.2620553596664335E-01_dp
               aw(25) = 1.8198087510110939E+00_dp
               aw(26) = 3.8700490617962435E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k13_sh(k, fit_coef, L_b, U_b)
               ELSE
                  IF (Rc <= 1000.0_dp) THEN
                     CALL get_fit_coef_k13_mi(k, fit_coef, L_b, U_b)
                  ELSE
                     CALL get_fit_coef_k13_lo(k, fit_coef, L_b, U_b)
                  END IF
               END IF
            END IF
         CASE (14)
            IF (Rc >= R_max) THEN
               aw(1) = 5.7793103608868073E-06_dp
               aw(2) = 6.3578213256604128E-05_dp
               aw(3) = 3.3649100160858427E-04_dp
               aw(4) = 1.2987332474693207E-03_dp
               aw(5) = 4.1590436292552775E-03_dp
               aw(6) = 1.1734004029327220E-02_dp
               aw(7) = 3.0162212092092591E-02_dp
               aw(8) = 7.2139138815286746E-02_dp
               aw(9) = 1.6282875663814733E-01_dp
               aw(10) = 3.5040070791304617E-01_dp
               aw(11) = 7.2459730019137303E-01_dp
               aw(12) = 1.4504311395006662E+00_dp
               aw(13) = 2.8398376292459662E+00_dp
               aw(14) = 5.5949042429760905E+00_dp
               aw(15) = 1.8475323460824520E-05_dp
               aw(16) = 1.2138616450563895E-04_dp
               aw(17) = 4.9691078029321441E-04_dp
               aw(18) = 1.6153616951963146E-03_dp
               aw(19) = 4.5476671955009150E-03_dp
               aw(20) = 1.1574075510211794E-02_dp
               aw(21) = 2.7304150630042847E-02_dp
               aw(22) = 6.0675597068478228E-02_dp
               aw(23) = 1.2843669732401369E-01_dp
               aw(24) = 2.6114026945752189E-01_dp
               aw(25) = 5.1372708252492216E-01_dp
               aw(26) = 9.8714811599006669E-01_dp
               aw(27) = 1.8941325485709735E+00_dp
               aw(28) = 3.9547772395130827E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k14_sh(k, fit_coef, L_b, U_b)
               ELSE
                  IF (Rc <= 1000.0_dp) THEN
                     CALL get_fit_coef_k14_mi(k, fit_coef, L_b, U_b)
                  ELSE
                     CALL get_fit_coef_k14_lo(k, fit_coef, L_b, U_b)
                  END IF
               END IF
            END IF
         CASE (15)
            IF (Rc >= R_max) THEN
               aw(1) = 3.2916678993378631E-06_dp
               aw(2) = 3.6211649942309474E-05_dp
               aw(3) = 1.9165204138720421E-04_dp
               aw(4) = 7.3970738152952427E-04_dp
               aw(5) = 2.3688276459558682E-03_dp
               aw(6) = 6.6832262130187142E-03_dp
               aw(7) = 1.7179191903350474E-02_dp
               aw(8) = 4.1087365085774624E-02_dp
               aw(9) = 9.2737907068114855E-02_dp
               aw(10) = 1.9954458193741717E-01_dp
               aw(11) = 4.1243407495664974E-01_dp
               aw(12) = 8.2392876030355866E-01_dp
               aw(13) = 1.6008068838318577E+00_dp
               aw(14) = 3.0544270604919670E+00_dp
               aw(15) = 5.8841569972355829E+00_dp
               aw(16) = 1.0522817667026174E-05_dp
               aw(17) = 6.9136785573218804E-05_dp
               aw(18) = 2.8302083846862936E-04_dp
               aw(19) = 9.2004648989894703E-04_dp
               aw(20) = 2.5901722444491987E-03_dp
               aw(21) = 6.5921351784279962E-03_dp
               aw(22) = 1.5551321239944434E-02_dp
               aw(23) = 3.4557794033658169E-02_dp
               aw(24) = 7.3145459688388037E-02_dp
               aw(25) = 1.4866796426944323E-01_dp
               aw(26) = 2.9202494628235121E-01_dp
               aw(27) = 5.5775841334330534E-01_dp
               aw(28) = 1.0451182163084207E+00_dp
               aw(29) = 1.9640378140134063E+00_dp
               aw(30) = 4.0343871534325908E+00_dp
            ELSE
               IF (Rc <= 100.0_dp) THEN
                  CALL get_fit_coef_k15_sh(k, fit_coef, L_b, U_b)
               ELSE
                  IF (Rc <= 1000.0_dp) THEN
                     CALL get_fit_coef_k15_mi(k, fit_coef, L_b, U_b)
                  ELSE
                     CALL get_fit_coef_k15_lo(k, fit_coef, L_b, U_b)
                  END IF
               END IF
            END IF
         CASE DEFAULT
            CPABORT("")
         END SELECT
      END IF

      IF (Rc < R_max .AND. Rc >= R_min) THEN
         R_transf = SQRT(L_b/Rc)
         R_transf = SQRT(R_transf)*(SQRT(U_b/(U_b - L_b)) - SQRT(Rc/(U_b - L_b)))
         DO kkk = 1, 2*k
            aw(kkk) = fit_coef(0, kkk)
            DO icoef = 1, 12
               aw(kkk) = aw(kkk) + fit_coef(icoef, kkk)*R_transf**icoef
            END DO
         END DO
      END IF

      IF (PRESENT(mm_error)) mm_error = get_minimax_numerical_error(Rc, aw)

      CALL timestop(handle)

   END SUBROUTINE get_minimax_coeff_k15

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k1(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 8.6669999999999998E+00_dp
      fit_coef(0, 1) = 4.4649260681804581E-01_dp
      fit_coef(1, 1) = 3.0675836870877398E-06_dp
      fit_coef(2, 1) = 1.2472664651690426E+00_dp
      fit_coef(3, 1) = -1.4325378312702204E+00_dp
      fit_coef(4, 1) = -2.8852101654931545E-01_dp
      fit_coef(5, 1) = 1.3224636333703371E+01_dp
      fit_coef(6, 1) = -5.6553448268455575E+01_dp
      fit_coef(7, 1) = 1.5090215216505641E+02_dp
      fit_coef(8, 1) = -2.8212708107138189E+02_dp
      fit_coef(9, 1) = 3.6970660684742296E+02_dp
      fit_coef(10, 1) = -3.2374043561240165E+02_dp
      fit_coef(11, 1) = 1.6996735688896140E+02_dp
      fit_coef(12, 1) = -4.0395921268987330E+01_dp
      fit_coef(0, 2) = 1.4290997884622449E+00_dp
      fit_coef(1, 2) = 5.6745594674266576E-06_dp
      fit_coef(2, 2) = 2.4380748907132181E+00_dp
      fit_coef(3, 2) = -2.8016418715477878E+00_dp
      fit_coef(4, 2) = 3.7584116167213760E-01_dp
      fit_coef(5, 2) = 2.3327242171726457E+01_dp
      fit_coef(6, 2) = -1.0520185473804629E+02_dp
      fit_coef(7, 2) = 2.8832898650596110E+02_dp
      fit_coef(8, 2) = -5.4924137182200832E+02_dp
      fit_coef(9, 2) = 7.3027624158467188E+02_dp
      fit_coef(10, 2) = -6.4709409025644425E+02_dp
      fit_coef(11, 2) = 3.4309801277522899E+02_dp
      fit_coef(12, 2) = -8.2224836383626666E+01_dp

   END SUBROUTINE get_fit_coef_k1

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k2(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 4.1539999999999999E+01_dp
      fit_coef(0, 1) = 9.3106187850842653E-02_dp
      fit_coef(1, 1) = -5.0629856861437123E-05_dp
      fit_coef(2, 1) = 1.0015082423014863E+00_dp
      fit_coef(3, 1) = -3.5061593037011978E+00_dp
      fit_coef(4, 1) = 1.4483017017771591E+01_dp
      fit_coef(5, 1) = -4.5087114171411322E+01_dp
      fit_coef(6, 1) = 1.0677856101159048E+02_dp
      fit_coef(7, 1) = -1.8217614209987144E+02_dp
      fit_coef(8, 1) = 2.1105105002638726E+02_dp
      fit_coef(9, 1) = -1.5567397052444477E+02_dp
      fit_coef(10, 1) = 6.5143729915539296E+01_dp
      fit_coef(11, 1) = -1.1144201390642060E+01_dp
      fit_coef(12, 1) = -4.0413697217822347E-01_dp
      fit_coef(0, 2) = 1.0301539533889652E+00_dp
      fit_coef(1, 2) = -1.8120259734086499E-04_dp
      fit_coef(2, 2) = 3.9565096665651787E+00_dp
      fit_coef(3, 2) = -1.3814853129007110E+01_dp
      fit_coef(4, 2) = 5.7316987025660488E+01_dp
      fit_coef(5, 2) = -1.7656913969609036E+02_dp
      fit_coef(6, 2) = 4.1570794568371764E+02_dp
      fit_coef(7, 2) = -7.0182518512857257E+02_dp
      fit_coef(8, 2) = 8.0388072043543650E+02_dp
      fit_coef(9, 2) = -5.9083346503424480E+02_dp
      fit_coef(10, 2) = 2.5295462667462738E+02_dp
      fit_coef(11, 2) = -4.9543705304205133E+01_dp
      fit_coef(12, 2) = 1.0078618878729646E+00_dp
      fit_coef(0, 3) = 2.9765555120179177E-01_dp
      fit_coef(1, 3) = -1.0353208854227029E-04_dp
      fit_coef(2, 3) = 2.3095283302989413E+00_dp
      fit_coef(3, 3) = -8.0576570729831154E+00_dp
      fit_coef(4, 3) = 3.3739368092954251E+01_dp
      fit_coef(5, 3) = -1.0475419759213936E+02_dp
      fit_coef(6, 3) = 2.4672759964823950E+02_dp
      fit_coef(7, 3) = -4.1482385671796834E+02_dp
      fit_coef(8, 3) = 4.6585822070233695E+02_dp
      fit_coef(9, 3) = -3.2366345592266214E+02_dp
      fit_coef(10, 3) = 1.1849771979400440E+02_dp
      fit_coef(11, 3) = -1.0906004349060501E+01_dp
      fit_coef(12, 3) = -3.7598063728027054E+00_dp
      fit_coef(0, 4) = 1.9917444037245748E+00_dp
      fit_coef(1, 4) = -1.5897415044355091E-04_dp
      fit_coef(2, 4) = 3.2506199533242417E+00_dp
      fit_coef(3, 4) = -1.1373065869782103E+01_dp
      fit_coef(4, 4) = 4.6558059922949951E+01_dp
      fit_coef(5, 4) = -1.4224165879889640E+02_dp
      fit_coef(6, 4) = 3.3658512396870015E+02_dp
      fit_coef(7, 4) = -5.7666672132123392E+02_dp
      fit_coef(8, 4) = 6.8845868913664435E+02_dp
      fit_coef(9, 4) = -5.5227838411317748E+02_dp
      fit_coef(10, 4) = 2.8109979885115376E+02_dp
      fit_coef(11, 4) = -8.0922550511268881E+01_dp
      fit_coef(12, 4) = 9.8057271611259083E+00_dp

   END SUBROUTINE get_fit_coef_k2

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k3(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 1.4680000000000001E+02_dp
      fit_coef(0, 1) = 2.6348727105229584E-02_dp
      fit_coef(1, 1) = 2.0153944068439938E-04_dp
      fit_coef(2, 1) = 6.1530657742665373E-01_dp
      fit_coef(3, 1) = -3.1917423654295867E+00_dp
      fit_coef(4, 1) = 1.8476493141601150E+01_dp
      fit_coef(5, 1) = -7.6512288788660001E+01_dp
      fit_coef(6, 1) = 2.3572982749223701E+02_dp
      fit_coef(7, 1) = -5.2565813705739686E+02_dp
      fit_coef(8, 1) = 8.2467549297012602E+02_dp
      fit_coef(9, 1) = -8.8144565993107483E+02_dp
      fit_coef(10, 1) = 6.0887047274283691E+02_dp
      fit_coef(11, 1) = -2.4457518923875350E+02_dp
      fit_coef(12, 1) = 4.3343959760454148E+01_dp
      fit_coef(0, 2) = 2.8996604313967794E-01_dp
      fit_coef(1, 2) = 9.1966653451787525E-04_dp
      fit_coef(2, 2) = 2.9603964889982954E+00_dp
      fit_coef(3, 2) = -1.5436155715607415E+01_dp
      fit_coef(4, 2) = 8.9046449784370381E+01_dp
      fit_coef(5, 2) = -3.6819211532781554E+02_dp
      fit_coef(6, 2) = 1.1368031324897308E+03_dp
      fit_coef(7, 2) = -2.5382696326038131E+03_dp
      fit_coef(8, 2) = 3.9825571880586440E+03_dp
      fit_coef(9, 2) = -4.2549854059515601E+03_dp
      fit_coef(10, 2) = 2.9378687444731368E+03_dp
      fit_coef(11, 2) = -1.1796924244872741E+03_dp
      fit_coef(12, 2) = 2.0901178723785654E+02_dp
      fit_coef(0, 3) = 1.5533604484128483E+00_dp
      fit_coef(1, 3) = 2.0797203873888733E-03_dp
      fit_coef(2, 3) = 6.8823846111885398E+00_dp
      fit_coef(3, 3) = -3.6081265871868659E+01_dp
      fit_coef(4, 3) = 2.0382721580134501E+02_dp
      fit_coef(5, 3) = -8.3588368834778782E+02_dp
      fit_coef(6, 3) = 2.5589042797513875E+03_dp
      fit_coef(7, 3) = -5.6774399573846031E+03_dp
      fit_coef(8, 3) = 8.8818039188117691E+03_dp
      fit_coef(9, 3) = -9.4791778618478966E+03_dp
      fit_coef(10, 3) = 6.5412200642397502E+03_dp
      fit_coef(11, 3) = -2.6248240660426704E+03_dp
      fit_coef(12, 3) = 4.6461771734428112E+02_dp
      fit_coef(0, 4) = 8.4232388261126792E-02_dp
      fit_coef(1, 4) = 4.8315351483134558E-04_dp
      fit_coef(2, 4) = 1.4901688878913604E+00_dp
      fit_coef(3, 4) = -7.7278009557035778E+00_dp
      fit_coef(4, 4) = 4.5200180529539338E+01_dp
      fit_coef(5, 4) = -1.8811770415936255E+02_dp
      fit_coef(6, 4) = 5.8400962109149202E+02_dp
      fit_coef(7, 4) = -1.3122325202133370E+03_dp
      fit_coef(8, 4) = 2.0704762349353055E+03_dp
      fit_coef(9, 4) = -2.2210949812031968E+03_dp
      fit_coef(10, 4) = 1.5374220056673389E+03_dp
      fit_coef(11, 4) = -6.1817894452963594E+02_dp
      fit_coef(12, 4) = 1.0958974382960224E+02_dp
      fit_coef(0, 5) = 5.5403377016655386E-01_dp
      fit_coef(1, 5) = 9.6266405874694016E-04_dp
      fit_coef(2, 5) = 3.2463696099948760E+00_dp
      fit_coef(3, 5) = -1.7046533273203099E+01_dp
      fit_coef(4, 5) = 9.6159163168231984E+01_dp
      fit_coef(5, 5) = -3.9387879060912519E+02_dp
      fit_coef(6, 5) = 1.2057094045527570E+03_dp
      fit_coef(7, 5) = -2.6721670291204064E+03_dp
      fit_coef(8, 5) = 4.1758049164403919E+03_dp
      fit_coef(9, 5) = -4.4574085233127789E+03_dp
      fit_coef(10, 5) = 3.0806227002555916E+03_dp
      fit_coef(11, 5) = -1.2391922190720968E+03_dp
      fit_coef(12, 5) = 2.1996432808098118E+02_dp
      fit_coef(0, 6) = 2.3557652649415379E+00_dp
      fit_coef(1, 6) = 1.2866271731043475E-03_dp
      fit_coef(2, 6) = 4.2915899775667716E+00_dp
      fit_coef(3, 6) = -2.2590368917521662E+01_dp
      fit_coef(4, 6) = 1.2467969331840358E+02_dp
      fit_coef(5, 6) = -5.0586944026595319E+02_dp
      fit_coef(6, 6) = 1.5294941344789261E+03_dp
      fit_coef(7, 6) = -3.3597570045760758E+03_dp
      fit_coef(8, 6) = 5.2196087711777982E+03_dp
      fit_coef(9, 6) = -5.5365578831549756E+03_dp
      fit_coef(10, 6) = 3.7972016920254159E+03_dp
      fit_coef(11, 6) = -1.5146892655309941E+03_dp
      fit_coef(12, 6) = 2.6667582714296555E+02_dp

   END SUBROUTINE get_fit_coef_k3

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k4(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 4.3610000000000002E+02_dp
      fit_coef(0, 1) = 8.8672649711285020E-03_dp
      fit_coef(1, 1) = 5.7946595950817889E-04_dp
      fit_coef(2, 1) = 3.5675945534374653E-01_dp
      fit_coef(3, 1) = -2.0716213353930191E+00_dp
      fit_coef(4, 1) = 1.2657470729413074E+01_dp
      fit_coef(5, 1) = -4.9885844724999998E+01_dp
      fit_coef(6, 1) = 1.4055190870387446E+02_dp
      fit_coef(7, 1) = -2.8346908048131127E+02_dp
      fit_coef(8, 1) = 4.0424697790393645E+02_dp
      fit_coef(9, 1) = -3.9722796612734163E+02_dp
      fit_coef(10, 1) = 2.5527240184308903E+02_dp
      fit_coef(11, 1) = -9.6356867177556367E+01_dp
      fit_coef(12, 1) = 1.6170840941701190E+01_dp
      fit_coef(0, 2) = 9.7559203877924675E-02_dp
      fit_coef(1, 2) = 3.0161149229635965E-03_dp
      fit_coef(2, 2) = 1.8835046917379701E+00_dp
      fit_coef(3, 2) = -1.1044051251373379E+01_dp
      fit_coef(4, 2) = 6.6851851055246527E+01_dp
      fit_coef(5, 2) = -2.6245846800700389E+02_dp
      fit_coef(6, 2) = 7.3917115340985299E+02_dp
      fit_coef(7, 2) = -1.4908588623752182E+03_dp
      fit_coef(8, 2) = 2.1239922172484485E+03_dp
      fit_coef(9, 2) = -2.0828446618499129E+03_dp
      fit_coef(10, 2) = 1.3350612779873063E+03_dp
      fit_coef(11, 2) = -5.0262753942439002E+02_dp
      fit_coef(12, 2) = 8.4153648005965323E+01_dp
      fit_coef(0, 3) = 5.1690155310834918E-01_dp
      fit_coef(1, 3) = 7.9435328867281975E-03_dp
      fit_coef(2, 3) = 5.1805612337915026E+00_dp
      fit_coef(3, 3) = -3.1004834727213211E+01_dp
      fit_coef(4, 3) = 1.8433410572917765E+02_dp
      fit_coef(5, 3) = -7.2586403084370136E+02_dp
      fit_coef(6, 3) = 2.0342579068012742E+03_dp
      fit_coef(7, 3) = -4.0678145171858919E+03_dp
      fit_coef(8, 3) = 5.7545086314708260E+03_dp
      fit_coef(9, 3) = -5.6192254926256765E+03_dp
      fit_coef(10, 3) = 3.5949228753416078E+03_dp
      fit_coef(11, 3) = -1.3527284950245175E+03_dp
      fit_coef(12, 3) = 2.2651488475207483E+02_dp
      fit_coef(0, 4) = 2.0307980329330455E+00_dp
      fit_coef(1, 4) = 1.5103721465188732E-02_dp
      fit_coef(2, 4) = 1.0170027292275309E+01_dp
      fit_coef(3, 4) = -6.1629536377936731E+01_dp
      fit_coef(4, 4) = 3.6339135279187411E+02_dp
      fit_coef(5, 4) = -1.4370905260427323E+03_dp
      fit_coef(6, 4) = 4.0240409210934931E+03_dp
      fit_coef(7, 4) = -8.0281015706711596E+03_dp
      fit_coef(8, 4) = 1.1327126412565874E+04_dp
      fit_coef(9, 4) = -1.1018179056966479E+04_dp
      fit_coef(10, 4) = 7.0109171513047604E+03_dp
      fit_coef(11, 4) = -2.6210011647009460E+03_dp
      fit_coef(12, 4) = 4.3584541357166989E+02_dp
      fit_coef(0, 5) = 2.8347760489671650E-02_dp
      fit_coef(1, 5) = 1.4649670575243566E-03_dp
      fit_coef(2, 5) = 8.8299650888912973E-01_dp
      fit_coef(3, 5) = -5.0888392853266016E+00_dp
      fit_coef(4, 5) = 3.1075839685772152E+01_dp
      fit_coef(5, 5) = -1.2128602025520448E+02_dp
      fit_coef(6, 5) = 3.3937414729640119E+02_dp
      fit_coef(7, 5) = -6.8302048253288706E+02_dp
      fit_coef(8, 5) = 9.7534981752611350E+02_dp
      fit_coef(9, 5) = -9.6151742177938490E+02_dp
      fit_coef(10, 5) = 6.2037394636212662E+02_dp
      fit_coef(11, 5) = -2.3512801932664019E+02_dp
      fit_coef(12, 5) = 3.9611471687460110E+01_dp
      fit_coef(0, 6) = 1.8628153335173334E-01_dp
      fit_coef(1, 6) = 3.5520276122272556E-03_dp
      fit_coef(2, 6) = 2.3086588899079099E+00_dp
      fit_coef(3, 6) = -1.3800728356647832E+01_dp
      fit_coef(4, 6) = 8.2177612215645453E+01_dp
      fit_coef(5, 6) = -3.2302705970393248E+02_dp
      fit_coef(6, 6) = 9.0586751432164715E+02_dp
      fit_coef(7, 6) = -1.8072280789056376E+03_dp
      fit_coef(8, 6) = 2.5395651527460632E+03_dp
      fit_coef(9, 6) = -2.4571249439954981E+03_dp
      fit_coef(10, 6) = 1.5570737965484213E+03_dp
      fit_coef(11, 6) = -5.8102885219535062E+02_dp
      fit_coef(12, 6) = 9.6652055418531859E+01_dp
      fit_coef(0, 7) = 7.6509143586699357E-01_dp
      fit_coef(1, 7) = 6.2314075869334632E-03_dp
      fit_coef(2, 7) = 4.2963285111812013E+00_dp
      fit_coef(3, 7) = -2.6274891972880781E+01_dp
      fit_coef(4, 7) = 1.5420391696751719E+02_dp
      fit_coef(5, 7) = -6.1288005128817349E+02_dp
      fit_coef(6, 7) = 1.7216536923227670E+03_dp
      fit_coef(7, 7) = -3.4485704001993868E+03_dp
      fit_coef(8, 7) = 4.8977145964963656E+03_dp
      fit_coef(9, 7) = -4.8034792371771946E+03_dp
      fit_coef(10, 7) = 3.0816756873370769E+03_dp
      fit_coef(11, 7) = -1.1601064628670051E+03_dp
      fit_coef(12, 7) = 1.9390631611709998E+02_dp
      fit_coef(0, 8) = 2.6286642321590428E+00_dp
      fit_coef(1, 8) = 7.7803605643329997E-03_dp
      fit_coef(2, 8) = 5.3941542397285547E+00_dp
      fit_coef(3, 8) = -3.3087210185862091E+01_dp
      fit_coef(4, 8) = 1.9273310816176871E+02_dp
      fit_coef(5, 8) = -7.6265175017497995E+02_dp
      fit_coef(6, 8) = 2.1209378243539750E+03_dp
      fit_coef(7, 8) = -4.1903565462442602E+03_dp
      fit_coef(8, 8) = 5.8461185745331686E+03_dp
      fit_coef(9, 8) = -5.6134540471156033E+03_dp
      fit_coef(10, 8) = 3.5238957022458790E+03_dp
      fit_coef(11, 8) = -1.3009219235427045E+03_dp
      fit_coef(12, 8) = 2.1404715392962453E+02_dp

   END SUBROUTINE get_fit_coef_k4

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k5_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 1.2552725467386876E-02_dp
      fit_coef(1, 1) = 4.6649769538903360E-02_dp
      fit_coef(2, 1) = 9.2959477678781308E-02_dp
      fit_coef(3, 1) = 1.0149142608425701E-01_dp
      fit_coef(4, 1) = 2.8480683386527764E-02_dp
      fit_coef(5, 1) = 1.1113054160128813E-01_dp
      fit_coef(6, 1) = -6.0921078972010467E-01_dp
      fit_coef(7, 1) = 1.2596196983676409E+00_dp
      fit_coef(8, 1) = -1.8451880523904862E+00_dp
      fit_coef(9, 1) = 1.7957653403032408E+00_dp
      fit_coef(10, 1) = -1.0556136831977423E+00_dp
      fit_coef(11, 1) = 3.3910276380778348E-01_dp
      fit_coef(12, 1) = -4.6181724196745676E-02_dp
      fit_coef(0, 2) = 8.6549382218897089E-02_dp
      fit_coef(1, 2) = 2.4941345083475788E-01_dp
      fit_coef(2, 2) = 4.9245083501834641E-01_dp
      fit_coef(3, 2) = 5.3644207042808922E-01_dp
      fit_coef(4, 2) = 1.2428540956728855E-01_dp
      fit_coef(5, 2) = 6.7557445154066409E-01_dp
      fit_coef(6, 2) = -3.7485327413027623E+00_dp
      fit_coef(7, 2) = 8.1470783598203376E+00_dp
      fit_coef(8, 2) = -1.2089572776244864E+01_dp
      fit_coef(9, 2) = 1.1847479598929507E+01_dp
      fit_coef(10, 2) = -7.1090620100857240E+00_dp
      fit_coef(11, 2) = 2.3701867666711376E+00_dp
      fit_coef(12, 2) = -3.3987884749260155E-01_dp
      fit_coef(0, 3) = 3.3100121276962868E-01_dp
      fit_coef(1, 3) = 6.4559516531419092E-01_dp
      fit_coef(2, 3) = 1.1695370338909203E+00_dp
      fit_coef(3, 3) = 1.1909933210765451E+00_dp
      fit_coef(4, 3) = 3.2031142163339787E-01_dp
      fit_coef(5, 3) = 1.9886975419904467E+00_dp
      fit_coef(6, 3) = -9.7799471531969644E+00_dp
      fit_coef(7, 3) = 2.2252882892705870E+01_dp
      fit_coef(8, 3) = -3.7167325466273596E+01_dp
      fit_coef(9, 3) = 4.2278212835435340E+01_dp
      fit_coef(10, 3) = -2.9956267208984151E+01_dp
      fit_coef(11, 3) = 1.1968794532343871E+01_dp
      fit_coef(12, 3) = -2.0781674372557291E+00_dp
      fit_coef(0, 4) = 1.0393831770652862E+00_dp
      fit_coef(1, 4) = 1.2828501087122686E+00_dp
      fit_coef(2, 4) = 2.0825379734819873E+00_dp
      fit_coef(3, 4) = 1.7998297820064115E+00_dp
      fit_coef(4, 4) = 2.5080914149094569E-01_dp
      fit_coef(5, 4) = 4.4983607972730759E+00_dp
      fit_coef(6, 4) = -1.8966342225496536E+01_dp
      fit_coef(7, 4) = 4.5287344234453357E+01_dp
      fit_coef(8, 4) = -7.1920958793917023E+01_dp
      fit_coef(9, 4) = 6.9428467733194353E+01_dp
      fit_coef(10, 4) = -3.7831237207121376E+01_dp
      fit_coef(11, 4) = 1.0235977584153179E+01_dp
      fit_coef(12, 4) = -9.4268597782017161E-01_dp
      fit_coef(0, 5) = 2.9584698863596302E+00_dp
      fit_coef(1, 5) = 2.1184323245195742E+00_dp
      fit_coef(2, 5) = 3.1709055258739989E+00_dp
      fit_coef(3, 5) = 2.2881184795555805E+00_dp
      fit_coef(4, 5) = 1.3681290961089937E-01_dp
      fit_coef(5, 5) = 5.7626338647647000E+00_dp
      fit_coef(6, 5) = -2.2759359631501187E+01_dp
      fit_coef(7, 5) = 5.5231152856101538E+01_dp
      fit_coef(8, 5) = -8.4310336806818569E+01_dp
      fit_coef(9, 5) = 7.9823159301719699E+01_dp
      fit_coef(10, 5) = -4.7210036737799030E+01_dp
      fit_coef(11, 5) = 1.6850346492810033E+01_dp
      fit_coef(12, 5) = -2.8940063886358156E+00_dp
      fit_coef(0, 6) = 3.4210297969183165E-02_dp
      fit_coef(1, 6) = 1.1998863470153487E-01_dp
      fit_coef(2, 6) = 2.4069816123621954E-01_dp
      fit_coef(3, 6) = 2.6041301008819967E-01_dp
      fit_coef(4, 6) = 6.3759797012621147E-02_dp
      fit_coef(5, 6) = 2.9954791472086834E-01_dp
      fit_coef(6, 6) = -1.6121273118458419E+00_dp
      fit_coef(7, 6) = 3.4105969041090662E+00_dp
      fit_coef(8, 6) = -5.1176148344713379E+00_dp
      fit_coef(9, 6) = 5.1110976139035680E+00_dp
      fit_coef(10, 6) = -3.1089038444203969E+00_dp
      fit_coef(11, 6) = 1.0468712585814062E+00_dp
      fit_coef(12, 6) = -1.5183313871493764E-01_dp
      fit_coef(0, 7) = 1.3026228141842919E-01_dp
      fit_coef(1, 7) = 2.9032232823923010E-01_dp
      fit_coef(2, 7) = 5.4877493441677716E-01_dp
      fit_coef(3, 7) = 5.9286988770199878E-01_dp
      fit_coef(4, 7) = 1.6832665275575420E-01_dp
      fit_coef(5, 7) = 7.6364041194666443E-01_dp
      fit_coef(6, 7) = -4.1778527336132214E+00_dp
      fit_coef(7, 7) = 8.2782930216328303E+00_dp
      fit_coef(8, 7) = -1.0896573954399400E+01_dp
      fit_coef(9, 7) = 9.5485112538544357E+00_dp
      fit_coef(10, 7) = -5.1102322033174294E+00_dp
      fit_coef(11, 7) = 1.4792267331631863E+00_dp
      fit_coef(12, 7) = -1.7439726860494367E-01_dp
      fit_coef(0, 8) = 4.0415838257065068E-01_dp
      fit_coef(1, 8) = 5.1255268146080202E-01_dp
      fit_coef(2, 8) = 7.9992895750622572E-01_dp
      fit_coef(3, 8) = 6.6970622776170818E-01_dp
      fit_coef(4, 8) = 3.1464540387258828E-02_dp
      fit_coef(5, 8) = 2.7391878500879105E+00_dp
      fit_coef(6, 8) = -1.1671256087331150E+01_dp
      fit_coef(7, 8) = 3.2335565030282325E+01_dp
      fit_coef(8, 8) = -6.2446813951907615E+01_dp
      fit_coef(9, 8) = 7.6210993792312763E+01_dp
      fit_coef(10, 8) = -5.5813236665948722E+01_dp
      fit_coef(11, 8) = 2.2568205101339789E+01_dp
      fit_coef(12, 8) = -3.8998743959104960E+00_dp
      fit_coef(0, 9) = 1.1232181639604371E+00_dp
      fit_coef(1, 9) = 7.5524302843049163E-01_dp
      fit_coef(2, 9) = 1.0162656347519898E+00_dp
      fit_coef(3, 9) = 5.2091134380973092E-01_dp
      fit_coef(4, 9) = 6.2773571552543109E-02_dp
      fit_coef(5, 9) = 4.7009056670980359E-01_dp
      fit_coef(6, 9) = 9.6931752144443684E-01_dp
      fit_coef(7, 9) = -8.1767165915400302E+00_dp
      fit_coef(8, 9) = 3.1589453385128493E+01_dp
      fit_coef(9, 9) = -6.0876284446976854E+01_dp
      fit_coef(10, 9) = 6.0604503383962623E+01_dp
      fit_coef(11, 9) = -3.0336494827398063E+01_dp
      fit_coef(12, 9) = 6.0926762476688889E+00_dp
      fit_coef(0, 10) = 3.0638169445090688E+00_dp
      fit_coef(1, 10) = 8.9165891212009063E-01_dp
      fit_coef(2, 10) = 1.1487286532265266E+00_dp
      fit_coef(3, 10) = 5.1064617564126968E-01_dp
      fit_coef(4, 10) = -4.7446553971131300E-01_dp
      fit_coef(5, 10) = 3.6370240067130397E+00_dp
      fit_coef(6, 10) = -1.5479102145878249E+01_dp
      fit_coef(7, 10) = 4.6883317751508919E+01_dp
      fit_coef(8, 10) = -8.9069835651344420E+01_dp
      fit_coef(9, 10) = 1.1482810068630945E+02_dp
      fit_coef(10, 10) = -9.5886031308673950E+01_dp
      fit_coef(11, 10) = 4.5677182175500711E+01_dp
      fit_coef(12, 10) = -9.3187455185522623E+00_dp

   END SUBROUTINE get_fit_coef_k5_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k5_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.1540000000000000E+03_dp
      fit_coef(0, 1) = 3.3529086066115803E-03_dp
      fit_coef(1, 1) = -8.2268999028435328E-06_dp
      fit_coef(2, 1) = 3.3611781342021224E-02_dp
      fit_coef(3, 1) = -1.3832592990876919E-01_dp
      fit_coef(4, 1) = 6.8079696115130983E-01_dp
      fit_coef(5, 1) = -2.6681612815812645E+00_dp
      fit_coef(6, 1) = 8.2290959869689022E+00_dp
      fit_coef(7, 1) = -1.9074468376984107E+01_dp
      fit_coef(8, 1) = 3.2127429353482611E+01_dp
      fit_coef(9, 1) = -3.7859356041773175E+01_dp
      fit_coef(10, 1) = 2.9446783510204408E+01_dp
      fit_coef(11, 1) = -1.3538459324077703E+01_dp
      fit_coef(12, 1) = 2.7808973445415268E+00_dp
      fit_coef(0, 2) = 3.6885509656511273E-02_dp
      fit_coef(1, 2) = -4.6135404124615802E-05_dp
      fit_coef(2, 2) = 1.8694797387922529E-01_dp
      fit_coef(3, 2) = -7.7044086984549065E-01_dp
      fit_coef(4, 2) = 3.7458845491550794E+00_dp
      fit_coef(5, 2) = -1.4597914775662142E+01_dp
      fit_coef(6, 2) = 4.4865271280504238E+01_dp
      fit_coef(7, 2) = -1.0374981069891059E+02_dp
      fit_coef(8, 2) = 1.7445816708209344E+02_dp
      fit_coef(9, 2) = -2.0533721999642813E+02_dp
      fit_coef(10, 2) = 1.5956721649659403E+02_dp
      fit_coef(11, 2) = -7.3312269788822348E+01_dp
      fit_coef(12, 2) = 1.5050733022098333E+01_dp
      fit_coef(0, 3) = 1.9524285077618225E-01_dp
      fit_coef(1, 3) = -1.3983711740977594E-04_dp
      fit_coef(2, 3) = 5.5888718949926475E-01_dp
      fit_coef(3, 3) = -2.3088502072129957E+00_dp
      fit_coef(4, 3) = 1.0999754812308348E+01_dp
      fit_coef(5, 3) = -4.2542674867782623E+01_dp
      fit_coef(6, 3) = 1.2991866971851243E+02_dp
      fit_coef(7, 3) = -2.9910040841843386E+02_dp
      fit_coef(8, 3) = 5.0131493127323301E+02_dp
      fit_coef(9, 3) = -5.8857681580426390E+02_dp
      fit_coef(10, 3) = 4.5649949573815150E+02_dp
      fit_coef(11, 3) = -2.0942533223478736E+02_dp
      fit_coef(12, 3) = 4.2945558512158080E+01_dp
      fit_coef(0, 4) = 7.5508121759060887E-01_dp
      fit_coef(1, 4) = -3.1605312969997662E-04_dp
      fit_coef(2, 4) = 1.2522786699714503E+00_dp
      fit_coef(3, 4) = -5.1813437090592105E+00_dp
      fit_coef(4, 4) = 2.4364506762575651E+01_dp
      fit_coef(5, 4) = -9.3776364610573921E+01_dp
      fit_coef(6, 4) = 2.8515029307910078E+02_dp
      fit_coef(7, 4) = -6.5452025447927747E+02_dp
      fit_coef(8, 4) = 1.0946041909024659E+03_dp
      fit_coef(9, 4) = -1.2829510036669149E+03_dp
      fit_coef(10, 4) = 9.9375334913507118E+02_dp
      fit_coef(11, 4) = -4.5544368131405361E+02_dp
      fit_coef(12, 4) = 9.3324529570079108E+01_dp
      fit_coef(0, 5) = 2.4730351053385511E+00_dp
      fit_coef(1, 5) = -5.6369602817021478E-04_dp
      fit_coef(2, 5) = 2.2229968075692117E+00_dp
      fit_coef(3, 5) = -9.2054846165557898E+00_dp
      fit_coef(4, 5) = 4.2976267592818459E+01_dp
      fit_coef(5, 5) = -1.6496009719447150E+02_dp
      fit_coef(6, 5) = 5.0038928803322085E+02_dp
      fit_coef(7, 5) = -1.1466294337635131E+03_dp
      fit_coef(8, 5) = 1.9151887594969048E+03_dp
      fit_coef(9, 5) = -2.2425568445298004E+03_dp
      fit_coef(10, 5) = 1.7357499822754255E+03_dp
      fit_coef(11, 5) = -7.9505099488495580E+02_dp
      fit_coef(12, 5) = 1.6284234597227066E+02_dp
      fit_coef(0, 6) = 1.0718594277376004E-02_dp
      fit_coef(1, 6) = -2.0665714271050207E-05_dp
      fit_coef(2, 6) = 8.4582172811458273E-02_dp
      fit_coef(3, 6) = -3.4797663076359847E-01_dp
      fit_coef(4, 6) = 1.7165661645244641E+00_dp
      fit_coef(5, 6) = -6.7291693231375893E+00_dp
      fit_coef(6, 6) = 2.0768418975747345E+01_dp
      fit_coef(7, 6) = -4.8168030617751583E+01_dp
      fit_coef(8, 6) = 8.1164964031411273E+01_dp
      fit_coef(9, 6) = -9.5679411237784350E+01_dp
      fit_coef(10, 6) = 7.4441710956106434E+01_dp
      fit_coef(11, 6) = -3.4234363834212587E+01_dp
      fit_coef(12, 6) = 7.0335229296269830E+00_dp
      fit_coef(0, 7) = 7.0423918791211709E-02_dp
      fit_coef(1, 7) = -6.0238742257850590E-05_dp
      fit_coef(2, 7) = 2.4113679202510896E-01_dp
      fit_coef(3, 7) = -9.9588383189867091E-01_dp
      fit_coef(4, 7) = 4.7554782250087468E+00_dp
      fit_coef(5, 7) = -1.8401716891304922E+01_dp
      fit_coef(6, 7) = 5.6243855385486498E+01_dp
      fit_coef(7, 7) = -1.2955708976498997E+02_dp
      fit_coef(8, 7) = 2.1724259904024584E+02_dp
      fit_coef(9, 7) = -2.5514786204737450E+02_dp
      fit_coef(10, 7) = 1.9794444495810336E+02_dp
      fit_coef(11, 7) = -9.0825962860161653E+01_dp
      fit_coef(12, 7) = 1.8627265430579317E+01_dp
      fit_coef(0, 8) = 2.8840117189215214E-01_dp
      fit_coef(1, 8) = -1.3225944743723339E-04_dp
      fit_coef(2, 8) = 5.2242393022997591E-01_dp
      fit_coef(3, 8) = -2.1627434438440303E+00_dp
      fit_coef(4, 8) = 1.0121981765560314E+01_dp
      fit_coef(5, 8) = -3.8888859154975577E+01_dp
      fit_coef(6, 8) = 1.1806445276665656E+02_dp
      fit_coef(7, 8) = -2.7069811465697614E+02_dp
      fit_coef(8, 8) = 4.5233627491539119E+02_dp
      fit_coef(9, 8) = -5.2982955497378055E+02_dp
      fit_coef(10, 8) = 4.1019408642078287E+02_dp
      fit_coef(11, 8) = -1.8792358285758019E+02_dp
      fit_coef(12, 8) = 3.8496219600202558E+01_dp
      fit_coef(0, 9) = 9.4333095808344458E-01_dp
      fit_coef(1, 9) = -2.1878257569565257E-04_dp
      fit_coef(2, 9) = 8.5869507294830094E-01_dp
      fit_coef(3, 9) = -3.5589224290290429E+00_dp
      fit_coef(4, 9) = 1.6493170506220842E+01_dp
      fit_coef(5, 9) = -6.3130491079340260E+01_dp
      fit_coef(6, 9) = 1.9101930146767009E+02_dp
      fit_coef(7, 9) = -4.3694612969128224E+02_dp
      fit_coef(8, 9) = 7.2886349059237500E+02_dp
      fit_coef(9, 9) = -8.5258297326024899E+02_dp
      fit_coef(10, 9) = 6.5938598174532456E+02_dp
      fit_coef(11, 9) = -3.0184732852622824E+02_dp
      fit_coef(12, 9) = 6.1796237778410379E+01_dp
      fit_coef(0, 10) = 2.8485312741055728E+00_dp
      fit_coef(1, 10) = -2.6597242771707791E-04_dp
      fit_coef(2, 10) = 1.0421996494576409E+00_dp
      fit_coef(3, 10) = -4.3207432293774151E+00_dp
      fit_coef(4, 10) = 1.9972564302073831E+01_dp
      fit_coef(5, 10) = -7.6371765475451213E+01_dp
      fit_coef(6, 10) = 2.3088260828483425E+02_dp
      fit_coef(7, 10) = -5.2780572532914982E+02_dp
      fit_coef(8, 10) = 8.8002099877005821E+02_dp
      fit_coef(9, 10) = -1.0290334029188657E+03_dp
      fit_coef(10, 10) = 7.9563294851641638E+02_dp
      fit_coef(11, 10) = -3.6413961414797035E+02_dp
      fit_coef(12, 10) = 7.4536963875157937E+01_dp

   END SUBROUTINE get_fit_coef_k5_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k6_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 1.0105361899470566E-02_dp
      fit_coef(1, 1) = 3.9366983227627925E-02_dp
      fit_coef(2, 1) = 7.8053377439735977E-02_dp
      fit_coef(3, 1) = 8.6837580919162383E-02_dp
      fit_coef(4, 1) = 4.6771807433750277E-02_dp
      fit_coef(5, 1) = -6.6734672326170835E-02_dp
      fit_coef(6, 1) = 9.9926171459107288E-02_dp
      fit_coef(7, 1) = -4.1644827063229539E-01_dp
      fit_coef(8, 1) = 7.6706989454888019E-01_dp
      fit_coef(9, 1) = -8.7152052671237190E-01_dp
      fit_coef(10, 1) = 6.5955681044267345E-01_dp
      fit_coef(11, 1) = -2.9363049092199256E-01_dp
      fit_coef(12, 1) = 5.6414782302946012E-02_dp
      fit_coef(0, 2) = 6.3047382400045412E-02_dp
      fit_coef(1, 2) = 2.1292104950028926E-01_dp
      fit_coef(2, 2) = 4.1580501626077238E-01_dp
      fit_coef(3, 2) = 4.5381507206227650E-01_dp
      fit_coef(4, 2) = 2.3261948283202027E-01_dp
      fit_coef(5, 2) = -3.7333612007818145E-01_dp
      fit_coef(6, 2) = 5.8891901066912178E-01_dp
      fit_coef(7, 2) = -2.2435149407329531E+00_dp
      fit_coef(8, 2) = 4.0761036364131602E+00_dp
      fit_coef(9, 2) = -4.5995236732489717E+00_dp
      fit_coef(10, 2) = 3.4552630197713277E+00_dp
      fit_coef(11, 2) = -1.5313736452266173E+00_dp
      fit_coef(12, 2) = 2.9408638427612482E-01_dp
      fit_coef(0, 3) = 2.1176071064496055E-01_dp
      fit_coef(1, 3) = 5.4982332295197556E-01_dp
      fit_coef(2, 3) = 1.0068197532238754E+00_dp
      fit_coef(3, 3) = 1.0684972989083634E+00_dp
      fit_coef(4, 3) = 6.0915448276900608E-01_dp
      fit_coef(5, 3) = -1.0238203491794160E+00_dp
      fit_coef(6, 3) = 2.1584333978631443E+00_dp
      fit_coef(7, 3) = -8.0611200001081951E+00_dp
      fit_coef(8, 3) = 1.5009961267213102E+01_dp
      fit_coef(9, 3) = -1.6828339663434313E+01_dp
      fit_coef(10, 3) = 1.2012563864685774E+01_dp
      fit_coef(11, 3) = -4.9849632599680698E+00_dp
      fit_coef(12, 3) = 9.0294945051842401E-01_dp
      fit_coef(0, 4) = 5.9366013171339316E-01_dp
      fit_coef(1, 4) = 1.1012914782020269E+00_dp
      fit_coef(2, 4) = 1.8080642501136919E+00_dp
      fit_coef(3, 4) = 1.7378820099827390E+00_dp
      fit_coef(4, 4) = 9.4465841905941816E-01_dp
      fit_coef(5, 4) = -7.9254597628432621E-01_dp
      fit_coef(6, 4) = 8.0775877014080266E-01_dp
      fit_coef(7, 4) = -2.8523919108461318E+00_dp
      fit_coef(8, 4) = -5.3307961129947901E-01_dp
      fit_coef(9, 4) = 6.7453388823676228E+00_dp
      fit_coef(10, 4) = -6.9248289720281742E+00_dp
      fit_coef(11, 4) = 2.9021292958698686E+00_dp
      fit_coef(12, 4) = -4.5434573537973694E-01_dp
      fit_coef(0, 5) = 1.5164074934712080E+00_dp
      fit_coef(1, 5) = 1.8998922859528067E+00_dp
      fit_coef(2, 5) = 2.7894259870776525E+00_dp
      fit_coef(3, 5) = 2.2836065839718458E+00_dp
      fit_coef(4, 5) = 1.2416212609604089E+00_dp
      fit_coef(5, 5) = -1.9153393783150017E+00_dp
      fit_coef(6, 5) = 7.9420474748342826E+00_dp
      fit_coef(7, 5) = -2.4492358033298018E+01_dp
      fit_coef(8, 5) = 4.7450493814300060E+01_dp
      fit_coef(9, 5) = -6.6407651969522661E+01_dp
      fit_coef(10, 5) = 6.0384585557333999E+01_dp
      fit_coef(11, 5) = -3.0372584465966241E+01_dp
      fit_coef(12, 5) = 6.3523080267784335E+00_dp
      fit_coef(0, 6) = 3.7104671974427847E+00_dp
      fit_coef(1, 6) = 2.8916704696341959E+00_dp
      fit_coef(2, 6) = 3.9056283047263785E+00_dp
      fit_coef(3, 6) = 2.7310607551007684E+00_dp
      fit_coef(4, 6) = 1.4374042741229791E+00_dp
      fit_coef(5, 6) = -3.4614723670541747E+00_dp
      fit_coef(6, 6) = 1.5897273148597373E+01_dp
      fit_coef(7, 6) = -4.2970295002223452E+01_dp
      fit_coef(8, 6) = 7.7894359559436722E+01_dp
      fit_coef(9, 6) = -9.5319740659866412E+01_dp
      fit_coef(10, 6) = 7.1950730534128169E+01_dp
      fit_coef(11, 6) = -2.9556628041991047E+01_dp
      fit_coef(12, 6) = 5.0101543449317205E+00_dp
      fit_coef(0, 7) = 2.6903751790128751E-02_dp
      fit_coef(1, 7) = 1.0173916967403172E-01_dp
      fit_coef(2, 7) = 2.0130716226577558E-01_dp
      fit_coef(3, 7) = 2.2172005830916267E-01_dp
      fit_coef(4, 7) = 1.1797631841567084E-01_dp
      fit_coef(5, 7) = -1.7271060154154347E-01_dp
      fit_coef(6, 7) = 2.7144317536179635E-01_dp
      fit_coef(7, 7) = -1.1089252106255005E+00_dp
      fit_coef(8, 7) = 2.0389731343247885E+00_dp
      fit_coef(9, 7) = -2.3121180806820467E+00_dp
      fit_coef(10, 7) = 1.7413776939963503E+00_dp
      fit_coef(11, 7) = -7.7162742625008818E-01_dp
      fit_coef(12, 7) = 1.4782756356123158E-01_dp
      fit_coef(0, 8) = 8.6888406935520149E-02_dp
      fit_coef(1, 8) = 2.4888995444204018E-01_dp
      fit_coef(2, 8) = 4.7199012197419227E-01_dp
      fit_coef(3, 8) = 5.1105786073145709E-01_dp
      fit_coef(4, 8) = 2.5740296918257327E-01_dp
      fit_coef(5, 8) = -4.3538515482693113E-01_dp
      fit_coef(6, 8) = 6.0029937052834403E-01_dp
      fit_coef(7, 8) = -2.3033251784048994E+00_dp
      fit_coef(8, 8) = 4.3994596232924978E+00_dp
      fit_coef(9, 8) = -5.1923409302629455E+00_dp
      fit_coef(10, 8) = 4.0222988922179050E+00_dp
      fit_coef(11, 8) = -1.8176576142253549E+00_dp
      fit_coef(12, 8) = 3.5423906095949409E-01_dp
      fit_coef(0, 9) = 2.3221985873092740E-01_dp
      fit_coef(1, 9) = 4.3438739423022149E-01_dp
      fit_coef(2, 9) = 7.0152946783034409E-01_dp
      fit_coef(3, 9) = 6.7672001661556891E-01_dp
      fit_coef(4, 9) = 4.8622153708717586E-01_dp
      fit_coef(5, 9) = -8.1041365351394823E-01_dp
      fit_coef(6, 9) = 2.6701912318720957E+00_dp
      fit_coef(7, 9) = -8.7205715741371570E+00_dp
      fit_coef(8, 9) = 1.3607358717959904E+01_dp
      fit_coef(9, 9) = -1.1724874566877528E+01_dp
      fit_coef(10, 9) = 5.7575285125779478E+00_dp
      fit_coef(11, 9) = -1.4357902050465099E+00_dp
      fit_coef(12, 9) = 1.1555583770987149E-01_dp
      fit_coef(0, 10) = 5.7977445413011319E-01_dp
      fit_coef(1, 10) = 6.7547720201746209E-01_dp
      fit_coef(2, 10) = 8.9708210777254471E-01_dp
      fit_coef(3, 10) = 6.3516703724773627E-01_dp
      fit_coef(4, 10) = 1.3236506065579925E-01_dp
      fit_coef(5, 10) = 1.1242277490562433E+00_dp
      fit_coef(6, 10) = -4.2391903965270812E+00_dp
      fit_coef(7, 10) = 1.3658960134998143E+01_dp
      fit_coef(8, 10) = -2.7055678231066494E+01_dp
      fit_coef(9, 10) = 2.7613753748428412E+01_dp
      fit_coef(10, 10) = -1.3209392906437861E+01_dp
      fit_coef(11, 10) = 1.7306395341708665E+00_dp
      fit_coef(12, 10) = 4.1257627070756253E-01_dp
      fit_coef(0, 11) = 1.3735367200768180E+00_dp
      fit_coef(1, 11) = 9.1204205070843380E-01_dp
      fit_coef(2, 11) = 1.0566274747133615E+00_dp
      fit_coef(3, 11) = 4.5662075463077345E-01_dp
      fit_coef(4, 11) = 4.9184167220410918E-01_dp
      fit_coef(5, 11) = -3.3759011651764603E+00_dp
      fit_coef(6, 11) = 1.7420644860890455E+01_dp
      fit_coef(7, 11) = -5.0539046535940329E+01_dp
      fit_coef(8, 11) = 1.0164085731266914E+02_dp
      fit_coef(9, 11) = -1.3287580093754113E+02_dp
      fit_coef(10, 11) = 1.0416384226984417E+02_dp
      fit_coef(11, 11) = -4.4207454491849361E+01_dp
      fit_coef(12, 11) = 7.8051805666697041E+00_dp
      fit_coef(0, 12) = 3.3571561549545903E+00_dp
      fit_coef(1, 12) = 1.0554507466925072E+00_dp
      fit_coef(2, 12) = 1.1790184249005025E+00_dp
      fit_coef(3, 12) = 4.8131299480494466E-01_dp
      fit_coef(4, 12) = -1.4373971957160242E-01_dp
      fit_coef(5, 12) = 7.0579667305022653E-01_dp
      fit_coef(6, 12) = -2.8433888321215206E+00_dp
      fit_coef(7, 12) = 1.4837347864587512E+01_dp
      fit_coef(8, 12) = -3.6747200470183884E+01_dp
      fit_coef(9, 12) = 5.8955148084383744E+01_dp
      fit_coef(10, 12) = -5.7353458415070598E+01_dp
      fit_coef(11, 12) = 2.9868904883232268E+01_dp
      fit_coef(12, 12) = -6.3780881018464042E+00_dp

   END SUBROUTINE get_fit_coef_k6_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k6_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 2.8070000000000000E+03_dp
      fit_coef(0, 1) = 1.3801918030773162E-03_dp
      fit_coef(1, 1) = -1.2020528483486374E-05_dp
      fit_coef(2, 1) = 2.3317704806718915E-02_dp
      fit_coef(3, 1) = -1.1872432794031121E-01_dp
      fit_coef(4, 1) = 6.7174538480429080E-01_dp
      fit_coef(5, 1) = -2.7971289614009751E+00_dp
      fit_coef(6, 1) = 8.6652048337367535E+00_dp
      fit_coef(7, 1) = -1.9404181258270555E+01_dp
      fit_coef(8, 1) = 3.0766165134677419E+01_dp
      fit_coef(9, 1) = -3.3559777770469893E+01_dp
      fit_coef(10, 1) = 2.3905145262192065E+01_dp
      fit_coef(11, 1) = -9.9977817398860651E+00_dp
      fit_coef(12, 1) = 1.8600935931530318E+00_dp
      fit_coef(0, 2) = 1.5183557920944469E-02_dp
      fit_coef(1, 2) = -7.0300568913718165E-05_dp
      fit_coef(2, 2) = 1.3444945392494670E-01_dp
      fit_coef(3, 2) = -6.8712467440053460E-01_dp
      fit_coef(4, 2) = 3.8439092624947708E+00_dp
      fit_coef(5, 2) = -1.5943863755337579E+01_dp
      fit_coef(6, 2) = 4.9281124351364589E+01_dp
      fit_coef(7, 2) = -1.1020727704358747E+02_dp
      fit_coef(8, 2) = 1.7458545156698841E+02_dp
      fit_coef(9, 2) = -1.9032304924079617E+02_dp
      fit_coef(10, 2) = 1.3551010987771855E+02_dp
      fit_coef(11, 2) = -5.6654664500478205E+01_dp
      fit_coef(12, 2) = 1.0537734527387647E+01_dp
      fit_coef(0, 3) = 8.0361343459697018E-02_dp
      fit_coef(1, 3) = -2.2644879609720777E-04_dp
      fit_coef(2, 3) = 4.2400787838782916E-01_dp
      fit_coef(3, 3) = -2.1793513661368733E+00_dp
      fit_coef(4, 3) = 1.1986721121277252E+01_dp
      fit_coef(5, 3) = -4.9550061013597045E+01_dp
      fit_coef(6, 3) = 1.5269170262558245E+02_dp
      fit_coef(7, 3) = -3.4093496564056039E+02_dp
      fit_coef(8, 3) = 5.3967525217022001E+02_dp
      fit_coef(9, 3) = -5.8799761545490060E+02_dp
      fit_coef(10, 3) = 4.1842708147977544E+02_dp
      fit_coef(11, 3) = -1.7483665022530386E+02_dp
      fit_coef(12, 3) = 3.2499971977161401E+01_dp
      fit_coef(0, 4) = 3.1025521547243406E-01_dp
      fit_coef(1, 4) = -5.5670514020829420E-04_dp
      fit_coef(2, 4) = 1.0283519940275256E+00_dp
      fit_coef(3, 4) = -5.3051948737708079E+00_dp
      fit_coef(4, 4) = 2.8857705580976056E+01_dp
      fit_coef(5, 4) = -1.1904005396455449E+02_dp
      fit_coef(6, 4) = 3.6607040216746668E+02_dp
      fit_coef(7, 4) = -8.1652968181369386E+02_dp
      fit_coef(8, 4) = 1.2917726681162258E+03_dp
      fit_coef(9, 4) = -1.4068028264311984E+03_dp
      fit_coef(10, 4) = 1.0006527523899351E+03_dp
      fit_coef(11, 4) = -4.1792217327242378E+02_dp
      fit_coef(12, 4) = 7.7650266966660254E+01_dp
      fit_coef(0, 5) = 9.9661072333548284E-01_dp
      fit_coef(1, 5) = -1.1171889827697565E-03_dp
      fit_coef(2, 5) = 2.0461119945840913E+00_dp
      fit_coef(3, 5) = -1.0580601944377532E+01_dp
      fit_coef(4, 5) = 5.7146266925973983E+01_dp
      fit_coef(5, 5) = -2.3540296082342908E+02_dp
      fit_coef(6, 5) = 7.2293490323806668E+02_dp
      fit_coef(7, 5) = -1.6114328355080609E+03_dp
      fit_coef(8, 5) = 2.5483891682119961E+03_dp
      fit_coef(9, 5) = -2.7745316619650062E+03_dp
      fit_coef(10, 5) = 1.9729645919543370E+03_dp
      fit_coef(11, 5) = -8.2377007268050909E+02_dp
      fit_coef(12, 5) = 1.5301212059298075E+02_dp
      fit_coef(0, 6) = 2.8871384041721235E+00_dp
      fit_coef(1, 6) = -1.8663713535930727E-03_dp
      fit_coef(2, 6) = 3.4015318790079183E+00_dp
      fit_coef(3, 6) = -1.7613395948721490E+01_dp
      fit_coef(4, 6) = 9.4741249061261925E+01_dp
      fit_coef(5, 6) = -3.8994245150367010E+02_dp
      fit_coef(6, 6) = 1.1965993777329293E+03_dp
      fit_coef(7, 6) = -2.6661783300231391E+03_dp
      fit_coef(8, 6) = 4.2155136266861591E+03_dp
      fit_coef(9, 6) = -4.5888706525865791E+03_dp
      fit_coef(10, 6) = 3.2626368424723773E+03_dp
      fit_coef(11, 6) = -1.3620271028068767E+03_dp
      fit_coef(12, 6) = 2.5294961462574193E+02_dp
      fit_coef(0, 7) = 4.4122094885822010E-03_dp
      fit_coef(1, 7) = -3.0506026172807132E-05_dp
      fit_coef(2, 7) = 5.9262798671013028E-02_dp
      fit_coef(3, 7) = -3.0161167176282039E-01_dp
      fit_coef(4, 7) = 1.7086225790243108E+00_dp
      fit_coef(5, 7) = -7.1111130868609358E+00_dp
      fit_coef(6, 7) = 2.2029502145209662E+01_dp
      fit_coef(7, 7) = -4.9330933530700591E+01_dp
      fit_coef(8, 7) = 7.8203553981541233E+01_dp
      fit_coef(9, 7) = -8.5283995941844253E+01_dp
      fit_coef(10, 7) = 6.0734937581336254E+01_dp
      fit_coef(11, 7) = -2.5396331019128979E+01_dp
      fit_coef(12, 7) = 4.7244104190757472E+00_dp
      fit_coef(0, 8) = 2.8989148486487894E-02_dp
      fit_coef(1, 8) = -9.5531400887653847E-05_dp
      fit_coef(2, 8) = 1.7940835754435794E-01_dp
      fit_coef(3, 8) = -9.2139806209629660E-01_dp
      fit_coef(4, 8) = 5.0793149159385784E+00_dp
      fit_coef(5, 8) = -2.0996406588101028E+01_dp
      fit_coef(6, 8) = 6.4728578095110564E+01_dp
      fit_coef(7, 8) = -1.4454498390418735E+02_dp
      fit_coef(8, 8) = 2.2882213635238523E+02_dp
      fit_coef(9, 8) = -2.4934889854822558E+02_dp
      fit_coef(10, 8) = 1.7747811343984634E+02_dp
      fit_coef(11, 8) = -7.4175788381531561E+01_dp
      fit_coef(12, 8) = 1.3791699026484896E+01_dp
      fit_coef(0, 9) = 1.1867775101904349E-01_dp
      fit_coef(1, 9) = -2.3002545083394702E-04_dp
      fit_coef(2, 9) = 4.2351088978816531E-01_dp
      fit_coef(3, 9) = -2.1868332609774375E+00_dp
      fit_coef(4, 9) = 1.1863193969553315E+01_dp
      fit_coef(5, 9) = -4.8911381713192242E+01_dp
      fit_coef(6, 9) = 1.5033689264253886E+02_dp
      fit_coef(7, 9) = -3.3524655332010315E+02_dp
      fit_coef(8, 9) = 5.3030221224968204E+02_dp
      fit_coef(9, 9) = -5.7746463929962727E+02_dp
      fit_coef(10, 9) = 4.1070292367082004E+02_dp
      fit_coef(11, 9) = -1.7150997861389590E+02_dp
      fit_coef(12, 9) = 3.1862954979470953E+01_dp
      fit_coef(0, 10) = 3.8613921710382521E-01_dp
      fit_coef(1, 10) = -4.4141730884221628E-04_dp
      fit_coef(2, 10) = 8.0440006468030922E-01_dp
      fit_coef(3, 10) = -4.1653901689708244E+00_dp
      fit_coef(4, 10) = 2.2403190787145370E+01_dp
      fit_coef(5, 10) = -9.2211974413485876E+01_dp
      fit_coef(6, 10) = 2.8295925585741463E+02_dp
      fit_coef(7, 10) = -6.3046689253660929E+02_dp
      fit_coef(8, 10) = 9.9682217293895803E+02_dp
      fit_coef(9, 10) = -1.0850822438050727E+03_dp
      fit_coef(10, 10) = 7.7146124564531726E+02_dp
      fit_coef(11, 10) = -3.2204769431435153E+02_dp
      fit_coef(12, 10) = 5.9807896790490901E+01_dp
      fit_coef(0, 11) = 1.0973261216294632E+00_dp
      fit_coef(1, 11) = -6.7305790071552931E-04_dp
      fit_coef(2, 11) = 1.2191563472414737E+00_dp
      fit_coef(3, 11) = -6.3236622978399160E+00_dp
      fit_coef(4, 11) = 3.3838641149158832E+01_dp
      fit_coef(5, 11) = -1.3912998330249843E+02_dp
      fit_coef(6, 11) = 4.2651665499720167E+02_dp
      fit_coef(7, 11) = -9.4985408069387211E+02_dp
      fit_coef(8, 11) = 1.5014096221052871E+03_dp
      fit_coef(9, 11) = -1.6340457638756484E+03_dp
      fit_coef(10, 11) = 1.1615527639524207E+03_dp
      fit_coef(11, 11) = -4.8480157740581882E+02_dp
      fit_coef(12, 11) = 9.0015704436915854E+01_dp
      fit_coef(0, 12) = 3.0334100910855679E+00_dp
      fit_coef(1, 12) = -8.0128677103265039E-04_dp
      fit_coef(2, 12) = 1.4492945609869972E+00_dp
      fit_coef(3, 12) = -7.5204430315974555E+00_dp
      fit_coef(4, 12) = 4.0192517443549171E+01_dp
      fit_coef(5, 12) = -1.6520813224035805E+02_dp
      fit_coef(6, 12) = 5.0634186437136140E+02_dp
      fit_coef(7, 12) = -1.1274845628371224E+03_dp
      fit_coef(8, 12) = 1.7820754347486288E+03_dp
      fit_coef(9, 12) = -1.9394275221726205E+03_dp
      fit_coef(10, 12) = 1.3785794074607632E+03_dp
      fit_coef(11, 12) = -5.7535942912689143E+02_dp
      fit_coef(12, 12) = 1.0682548193028630E+02_dp

   END SUBROUTINE get_fit_coef_k6_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k7_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+00_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 8.4890821794056659E-03_dp
      fit_coef(1, 1) = 3.3909965321659449E-02_dp
      fit_coef(2, 1) = 6.7634541203882806E-02_dp
      fit_coef(3, 1) = 7.6353562644566234E-02_dp
      fit_coef(4, 1) = 3.5948563267733900E-02_dp
      fit_coef(5, 1) = -2.5771176000865662E-02_dp
      fit_coef(6, 1) = -4.0945070980132330E-02_dp
      fit_coef(7, 1) = -4.3343260753829285E-02_dp
      fit_coef(8, 1) = 1.3815432901925578E-01_dp
      fit_coef(9, 1) = -1.7538024437982933E-01_dp
      fit_coef(10, 1) = 1.6543062150575619E-01_dp
      fit_coef(11, 1) = -9.1169573964789280E-02_dp
      fit_coef(12, 1) = 2.0264349482170308E-02_dp
      fit_coef(0, 2) = 5.0116642286956788E-02_dp
      fit_coef(1, 2) = 1.8365009326344528E-01_dp
      fit_coef(2, 2) = 3.5986371694373215E-01_dp
      fit_coef(3, 2) = 3.9798008193331352E-01_dp
      fit_coef(4, 2) = 1.7835093713215566E-01_dp
      fit_coef(5, 2) = -1.3327325046872260E-01_dp
      fit_coef(6, 2) = -2.2419539361558347E-01_dp
      fit_coef(7, 2) = -1.5459559592831851E-01_dp
      fit_coef(8, 2) = 5.7498920152238775E-01_dp
      fit_coef(9, 2) = -7.2077559903241084E-01_dp
      fit_coef(10, 2) = 6.9718992787938672E-01_dp
      fit_coef(11, 2) = -3.9787065520994253E-01_dp
      fit_coef(12, 2) = 9.0654023942071268E-02_dp
      fit_coef(0, 3) = 1.5366998360970857E-01_dp
      fit_coef(1, 3) = 4.7204846669452077E-01_dp
      fit_coef(2, 3) = 8.8325935893055196E-01_dp
      fit_coef(3, 3) = 9.5894564770013224E-01_dp
      fit_coef(4, 3) = 4.2827132793975348E-01_dp
      fit_coef(5, 3) = -3.2332267433019096E-01_dp
      fit_coef(6, 3) = -6.1277314101509384E-01_dp
      fit_coef(7, 3) = -3.1904891130696167E-01_dp
      fit_coef(8, 3) = 1.7347388703763822E+00_dp
      fit_coef(9, 3) = -2.4531109689297721E+00_dp
      fit_coef(10, 3) = 2.2903712151885607E+00_dp
      fit_coef(11, 3) = -1.2225979103964606E+00_dp
      fit_coef(12, 3) = 2.6689160566637776E-01_dp
      fit_coef(0, 4) = 3.9160355914825828E-01_dp
      fit_coef(1, 4) = 9.3781039645429920E-01_dp
      fit_coef(2, 4) = 1.6046991663294892E+00_dp
      fit_coef(3, 4) = 1.6443527376353326E+00_dp
      fit_coef(4, 4) = 8.0944769705145370E-01_dp
      fit_coef(5, 4) = -5.1143473303498077E-01_dp
      fit_coef(6, 4) = -5.2679642212467515E-01_dp
      fit_coef(7, 4) = -1.8364447956089240E+00_dp
      fit_coef(8, 4) = 2.5683156855002216E+00_dp
      fit_coef(9, 4) = 1.5501948487951303E-01_dp
      fit_coef(10, 4) = -1.8725579182282932E+00_dp
      fit_coef(11, 4) = 1.2221919802984111E+00_dp
      fit_coef(12, 4) = -2.7544768102671480E-01_dp
      fit_coef(0, 5) = 9.1783363747242741E-01_dp
      fit_coef(1, 5) = 1.6309461105431255E+00_dp
      fit_coef(2, 5) = 2.4937510867743007E+00_dp
      fit_coef(3, 5) = 2.2928527582942770E+00_dp
      fit_coef(4, 5) = 9.1128547426324147E-01_dp
      fit_coef(5, 5) = 9.6182376136104375E-01_dp
      fit_coef(6, 5) = -6.3804461030632824E+00_dp
      fit_coef(7, 5) = 1.5896262996437345E+01_dp
      fit_coef(8, 5) = -3.2697959752155285E+01_dp
      fit_coef(9, 5) = 4.0204694276308771E+01_dp
      fit_coef(10, 5) = -2.6993883785520257E+01_dp
      fit_coef(11, 5) = 9.1873159371192017E+00_dp
      fit_coef(12, 5) = -1.2193412409631701E+00_dp
      fit_coef(0, 6) = 2.0447798807240689E+00_dp
      fit_coef(1, 6) = 2.5647815805469976E+00_dp
      fit_coef(2, 6) = 3.5222500643076748E+00_dp
      fit_coef(3, 6) = 2.7811523063687598E+00_dp
      fit_coef(4, 6) = 1.2018538663100924E+00_dp
      fit_coef(5, 6) = -5.4703655436328180E-01_dp
      fit_coef(6, 6) = 3.2225185493071979E+00_dp
      fit_coef(7, 6) = -1.4146455603976277E+01_dp
      fit_coef(8, 6) = 3.4172972270754862E+01_dp
      fit_coef(9, 6) = -5.7518986665041560E+01_dp
      fit_coef(10, 6) = 5.7935519876217747E+01_dp
      fit_coef(11, 6) = -3.0472542712278077E+01_dp
      fit_coef(12, 6) = 6.4684369280250085E+00_dp
      fit_coef(0, 7) = 4.4889127849133947E+00_dp
      fit_coef(1, 7) = 3.6844329124185262E+00_dp
      fit_coef(2, 7) = 4.6629325988561545E+00_dp
      fit_coef(3, 7) = 3.2053175028746703E+00_dp
      fit_coef(4, 7) = 1.1894024823027469E+00_dp
      fit_coef(5, 7) = -5.5664546189271846E-01_dp
      fit_coef(6, 7) = 4.4797254977157817E+00_dp
      fit_coef(7, 7) = -1.3162135892918478E+01_dp
      fit_coef(8, 7) = 2.6920639733012177E+01_dp
      fit_coef(9, 7) = -3.7885356014248515E+01_dp
      fit_coef(10, 7) = 3.0404958289755267E+01_dp
      fit_coef(11, 7) = -1.2048234826847173E+01_dp
      fit_coef(12, 7) = 1.7536434043357469E+00_dp
      fit_coef(0, 8) = 2.2333642667863175E-02_dp
      fit_coef(1, 8) = 8.7650336191702150E-02_dp
      fit_coef(2, 8) = 1.7407545680352357E-01_dp
      fit_coef(3, 8) = 1.9523723539373550E-01_dp
      fit_coef(4, 8) = 9.1078382735736693E-02_dp
      fit_coef(5, 8) = -6.5238131379128483E-02_dp
      fit_coef(6, 8) = -1.0790743432696713E-01_dp
      fit_coef(7, 8) = -1.0055162956331706E-01_dp
      fit_coef(8, 8) = 3.3608939631307738E-01_dp
      fit_coef(9, 8) = -4.2744783871398162E-01_dp
      fit_coef(10, 8) = 4.0569757577173882E-01_dp
      fit_coef(11, 8) = -2.2517859275370616E-01_dp
      fit_coef(12, 8) = 5.0308943083685061E-02_dp
      fit_coef(0, 9) = 6.5162011076589435E-02_dp
      fit_coef(1, 9) = 2.1474236673046285E-01_dp
      fit_coef(2, 9) = 4.1104112490299882E-01_dp
      fit_coef(3, 9) = 4.4722258512292351E-01_dp
      fit_coef(4, 9) = 1.8925685414469240E-01_dp
      fit_coef(5, 9) = -1.5238925014408933E-01_dp
      fit_coef(6, 9) = -2.8181109625687978E-01_dp
      fit_coef(7, 9) = 4.2073721801476338E-03_dp
      fit_coef(8, 9) = 2.9605389145708244E-01_dp
      fit_coef(9, 9) = -3.8389360909607761E-01_dp
      fit_coef(10, 9) = 4.5264861145927554E-01_dp
      fit_coef(11, 9) = -3.0203559346875924E-01_dp
      fit_coef(12, 9) = 7.4860217823464928E-02_dp
      fit_coef(0, 10) = 1.5346948637055580E-01_dp
      fit_coef(1, 10) = 3.6854433133610354E-01_dp
      fit_coef(2, 10) = 6.2892668627471893E-01_dp
      fit_coef(3, 10) = 6.5477471935299725E-01_dp
      fit_coef(4, 10) = 3.3294850052423303E-01_dp
      fit_coef(5, 10) = -2.7649376642801776E-01_dp
      fit_coef(6, 10) = -8.0513348586680189E-03_dp
      fit_coef(7, 10) = -2.0751899699120990E+00_dp
      fit_coef(8, 10) = 5.2487682729558243E+00_dp
      fit_coef(9, 10) = -6.5715600622772241E+00_dp
      fit_coef(10, 10) = 4.9805286462168388E+00_dp
      fit_coef(11, 10) = -2.1493517756117457E+00_dp
      fit_coef(12, 10) = 4.0028581387747814E-01_dp
      fit_coef(0, 11) = 3.4692296309076204E-01_dp
      fit_coef(1, 11) = 5.7228042968532444E-01_dp
      fit_coef(2, 11) = 8.0923016588837160E-01_dp
      fit_coef(3, 11) = 6.8825242004374665E-01_dp
      fit_coef(4, 11) = 3.0900931504334961E-01_dp
      fit_coef(5, 11) = 5.2224154947033075E-01_dp
      fit_coef(6, 11) = -2.6356703498487173E+00_dp
      fit_coef(7, 11) = 8.5966604379962970E+00_dp
      fit_coef(8, 11) = -2.2765680947013255E+01_dp
      fit_coef(9, 11) = 3.4684627103439190E+01_dp
      fit_coef(10, 11) = -2.9504573991969686E+01_dp
      fit_coef(11, 11) = 1.3352471045720314E+01_dp
      fit_coef(12, 11) = -2.5338463087754426E+00_dp
      fit_coef(0, 12) = 7.5471792316070130E-01_dp
      fit_coef(1, 12) = 8.1698667304635597E-01_dp
      fit_coef(2, 12) = 9.6453392023721507E-01_dp
      fit_coef(3, 12) = 5.8219771460039760E-01_dp
      fit_coef(4, 12) = 5.6834288665372115E-02_dp
      fit_coef(5, 12) = 9.5294461790630214E-01_dp
      fit_coef(6, 12) = -1.9967581073724616E+00_dp
      fit_coef(7, 12) = 3.9290913995308978E+00_dp
      fit_coef(8, 12) = 7.8240450388508198E-01_dp
      fit_coef(9, 12) = -1.8452411993447775E+01_dp
      fit_coef(10, 12) = 2.9544447334169273E+01_dp
      fit_coef(11, 12) = -1.9185547926257918E+01_dp
      fit_coef(12, 12) = 4.6420265475532112E+00_dp
      fit_coef(0, 13) = 1.6044221851675884E+00_dp
      fit_coef(1, 13) = 1.0402801401561750E+00_dp
      fit_coef(2, 13) = 1.0858841570371729E+00_dp
      fit_coef(3, 13) = 4.1610384890729457E-01_dp
      fit_coef(4, 13) = 3.6067842390926969E-01_dp
      fit_coef(5, 13) = -2.5175945757607381E+00_dp
      fit_coef(6, 13) = 1.3559516170001507E+01_dp
      fit_coef(7, 13) = -3.8412585801701184E+01_dp
      fit_coef(8, 13) = 7.5725440762246222E+01_dp
      fit_coef(9, 13) = -9.4629438016079902E+01_dp
      fit_coef(10, 13) = 6.8485370488181275E+01_dp
      fit_coef(11, 13) = -2.5761663952186389E+01_dp
      fit_coef(12, 13) = 3.8061288493554284E+00_dp
      fit_coef(0, 14) = 3.6234684753265820E+00_dp
      fit_coef(1, 14) = 1.1911013659903149E+00_dp
      fit_coef(2, 14) = 1.2089364306597108E+00_dp
      fit_coef(3, 14) = 4.6442125027798664E-01_dp
      fit_coef(4, 14) = -3.1777576152208559E-01_dp
      fit_coef(5, 14) = 2.0241094124978116E+00_dp
      fit_coef(6, 14) = -8.2031751444577026E+00_dp
      fit_coef(7, 14) = 2.9168093356344343E+01_dp
      fit_coef(8, 14) = -6.1746832221017200E+01_dp
      fit_coef(9, 14) = 8.7168512210629885E+01_dp
      fit_coef(10, 14) = -7.6686519988709378E+01_dp
      fit_coef(11, 14) = 3.7037473909730942E+01_dp
      fit_coef(12, 14) = -7.4610849374495922E+00_dp

   END SUBROUTINE get_fit_coef_k7_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k7_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 6.3730000000000000E+03_dp
      fit_coef(0, 1) = 6.0677583366177139E-04_dp
      fit_coef(1, 1) = 1.2203511941973786E-05_dp
      fit_coef(2, 1) = 1.5542282929846361E-02_dp
      fit_coef(3, 1) = -8.9625369583599124E-02_dp
      fit_coef(4, 1) = 5.4937244690257003E-01_dp
      fit_coef(5, 1) = -2.3465783422261945E+00_dp
      fit_coef(6, 1) = 7.2618002882762989E+00_dp
      fit_coef(7, 1) = -1.5964210654587358E+01_dp
      fit_coef(8, 1) = 2.4547214366628502E+01_dp
      fit_coef(9, 1) = -2.5732576011672748E+01_dp
      fit_coef(10, 1) = 1.7495704575487210E+01_dp
      fit_coef(11, 1) = -6.9489399028728807E+00_dp
      fit_coef(12, 1) = 1.2232329673889215E+00_dp
      fit_coef(0, 2) = 6.6752714913863925E-03_dp
      fit_coef(1, 2) = 7.0822129335160658E-05_dp
      fit_coef(2, 2) = 9.1986377885775919E-02_dp
      fit_coef(3, 2) = -5.3405428642774699E-01_dp
      fit_coef(4, 2) = 3.2354397768485041E+00_dp
      fit_coef(5, 2) = -1.3775260985114405E+01_dp
      fit_coef(6, 2) = 4.2538770413228939E+01_dp
      fit_coef(7, 2) = -9.3398515752623808E+01_dp
      fit_coef(8, 2) = 1.4350316188715925E+02_dp
      fit_coef(9, 2) = -1.5036166589456491E+02_dp
      fit_coef(10, 2) = 1.0220160925073999E+02_dp
      fit_coef(11, 2) = -4.0584863952257315E+01_dp
      fit_coef(12, 2) = 7.1433481247893695E+00_dp
      fit_coef(0, 3) = 3.5329627891532003E-02_dp
      fit_coef(1, 3) = 2.2595479303685819E-04_dp
      fit_coef(2, 3) = 3.0068670299306055E-01_dp
      fit_coef(3, 3) = -1.7613684261571489E+00_dp
      fit_coef(4, 3) = 1.0488351739421789E+01_dp
      fit_coef(5, 3) = -4.4542175766068205E+01_dp
      fit_coef(6, 3) = 1.3710066819725051E+02_dp
      fit_coef(7, 3) = -3.0052458124580761E+02_dp
      fit_coef(8, 3) = 4.6158839612526646E+02_dp
      fit_coef(9, 3) = -4.8381807226730933E+02_dp
      fit_coef(10, 3) = 3.2906228868490689E+02_dp
      fit_coef(11, 3) = -1.3076538707988345E+02_dp
      fit_coef(12, 3) = 2.3031862239615609E+01_dp
      fit_coef(0, 4) = 1.3636675797272446E-01_dp
      fit_coef(1, 4) = 5.6647704292489109E-04_dp
      fit_coef(2, 4) = 7.6618173910360343E-01_dp
      fit_coef(3, 4) = -4.5142657914579676E+00_dp
      fit_coef(4, 4) = 2.6580955043085275E+01_dp
      fit_coef(5, 4) = -1.1272540497018385E+02_dp
      fit_coef(6, 4) = 3.4623696291949358E+02_dp
      fit_coef(7, 4) = -7.5826687162428038E+02_dp
      fit_coef(8, 4) = 1.1645610278715124E+03_dp
      fit_coef(9, 4) = -1.2210276266592973E+03_dp
      fit_coef(10, 4) = 8.3085597800660651E+02_dp
      fit_coef(11, 4) = -3.3033932676724800E+02_dp
      fit_coef(12, 4) = 5.8211140143063389E+01_dp
      fit_coef(0, 5) = 4.3692354129416011E-01_dp
      fit_coef(1, 5) = 1.1976571152026670E-03_dp
      fit_coef(2, 5) = 1.6377604996306871E+00_dp
      fit_coef(3, 5) = -9.6866759364374868E+00_dp
      fit_coef(4, 5) = 5.6612700356628949E+01_dp
      fit_coef(5, 5) = -2.3984102923713382E+02_dp
      fit_coef(6, 5) = 7.3561273922776581E+02_dp
      fit_coef(7, 5) = -1.6099560282310861E+03_dp
      fit_coef(8, 5) = 2.4723564756113196E+03_dp
      fit_coef(9, 5) = -2.5927151309452665E+03_dp
      fit_coef(10, 5) = 1.7647692924297401E+03_dp
      fit_coef(11, 5) = -7.0188939174830352E+02_dp
      fit_coef(12, 5) = 1.2372418269084275E+02_dp
      fit_coef(0, 6) = 1.2378497399176813E+00_dp
      fit_coef(1, 6) = 2.1799404763566546E-03_dp
      fit_coef(2, 6) = 3.0033304053549266E+00_dp
      fit_coef(3, 6) = -1.7809339259291960E+01_dp
      fit_coef(4, 6) = 1.0356146048512531E+02_dp
      fit_coef(5, 6) = -4.3841186307805572E+02_dp
      fit_coef(6, 6) = 1.3433098010233427E+03_dp
      fit_coef(7, 6) = -2.9385308911410843E+03_dp
      fit_coef(8, 6) = 4.5121527260780776E+03_dp
      fit_coef(9, 6) = -4.7323142457902477E+03_dp
      fit_coef(10, 6) = 3.2217506991135392E+03_dp
      fit_coef(11, 6) = -1.2816488835906916E+03_dp
      fit_coef(12, 6) = 2.2596889581599186E+02_dp
      fit_coef(0, 7) = 3.2781853114633690E+00_dp
      fit_coef(1, 7) = 3.4389001435898561E-03_dp
      fit_coef(2, 7) = 4.7592429785757018E+00_dp
      fit_coef(3, 7) = -2.8265375133025159E+01_dp
      fit_coef(4, 7) = 1.6386545527527471E+02_dp
      fit_coef(5, 7) = -6.9336876574524956E+02_dp
      fit_coef(6, 7) = 2.1232162962264451E+03_dp
      fit_coef(7, 7) = -4.6431587196422133E+03_dp
      fit_coef(8, 7) = 7.1291035126709503E+03_dp
      fit_coef(9, 7) = -7.4773722512143786E+03_dp
      fit_coef(10, 7) = 5.0911609048640648E+03_dp
      fit_coef(11, 7) = -2.0255881694388165E+03_dp
      fit_coef(12, 7) = 3.5717987126526839E+02_dp
      fit_coef(0, 8) = 1.9397592535656612E-03_dp
      fit_coef(1, 8) = 3.1224136945126820E-05_dp
      fit_coef(2, 8) = 3.9794210955105026E-02_dp
      fit_coef(3, 8) = -2.2945121723050055E-01_dp
      fit_coef(4, 8) = 1.4081093618036822E+00_dp
      fit_coef(5, 8) = -6.0128883949026397E+00_dp
      fit_coef(6, 8) = 1.8615856974219142E+01_dp
      fit_coef(7, 8) = -4.0950738259707443E+01_dp
      fit_coef(8, 8) = 6.2995147189900791E+01_dp
      fit_coef(9, 8) = -6.6049282238872578E+01_dp
      fit_coef(10, 8) = 4.4906773129267386E+01_dp
      fit_coef(11, 8) = -1.7833720418073618E+01_dp
      fit_coef(12, 8) = 3.1386820276994531E+00_dp
      fit_coef(0, 9) = 1.2744789430102283E-02_dp
      fit_coef(1, 9) = 9.4709522213266069E-05_dp
      fit_coef(2, 9) = 1.2561186892205739E-01_dp
      fit_coef(3, 9) = -7.3485105104893944E-01_dp
      fit_coef(4, 9) = 4.3871454899236513E+00_dp
      fit_coef(5, 9) = -1.8627306443837497E+01_dp
      fit_coef(6, 9) = 5.7355825078700910E+01_dp
      fit_coef(7, 9) = -1.2570281346740651E+02_dp
      fit_coef(8, 9) = 1.9298761023815177E+02_dp
      fit_coef(9, 9) = -2.0218882506002689E+02_dp
      fit_coef(10, 9) = 1.3746187852198693E+02_dp
      fit_coef(11, 9) = -5.4608420832822929E+01_dp
      fit_coef(12, 9) = 9.6158489721068801E+00_dp
      fit_coef(0, 10) = 5.2173333294315491E-02_dp
      fit_coef(1, 10) = 2.3055569105486745E-04_dp
      fit_coef(2, 10) = 3.1278409064824336E-01_dp
      fit_coef(3, 10) = -1.8448833978165269E+00_dp
      fit_coef(4, 10) = 1.0840143395546759E+01_dp
      fit_coef(5, 10) = -4.5958901600709979E+01_dp
      fit_coef(6, 10) = 1.4110750881102726E+02_dp
      fit_coef(7, 10) = -3.0897514713921822E+02_dp
      fit_coef(8, 10) = 4.7453270093475732E+02_dp
      fit_coef(9, 10) = -4.9757949866356881E+02_dp
      fit_coef(10, 10) = 3.3861207241883631E+02_dp
      fit_coef(11, 10) = -1.3464036769091453E+02_dp
      fit_coef(12, 10) = 2.3727726896949317E+01_dp
      fit_coef(0, 11) = 1.6963116889689617E-01_dp
      fit_coef(1, 11) = 4.6901741492581534E-04_dp
      fit_coef(2, 11) = 6.4447000621439654E-01_dp
      fit_coef(3, 11) = -3.8181538304305924E+00_dp
      fit_coef(4, 11) = 2.2242335828744238E+01_dp
      fit_coef(5, 11) = -9.4193777329169180E+01_dp
      fit_coef(6, 11) = 2.8872166699409729E+02_dp
      fit_coef(7, 11) = -6.3173717103578508E+02_dp
      fit_coef(8, 11) = 9.7012917536585030E+02_dp
      fit_coef(9, 11) = -1.0174610499996350E+03_dp
      fit_coef(10, 11) = 6.9265148866175264E+02_dp
      fit_coef(11, 11) = -2.7552649004097742E+02_dp
      fit_coef(12, 11) = 4.8575054818192108E+01_dp
      fit_coef(0, 12) = 4.7830241431219889E-01_dp
      fit_coef(1, 12) = 8.0536633187941779E-04_dp
      fit_coef(2, 12) = 1.1161391199045863E+00_dp
      fit_coef(3, 12) = -6.6319638481849212E+00_dp
      fit_coef(4, 12) = 3.8412707149315935E+01_dp
      fit_coef(5, 12) = -1.6252589178932661E+02_dp
      fit_coef(6, 12) = 4.9759933805308168E+02_dp
      fit_coef(7, 12) = -1.0881301553940079E+03_dp
      fit_coef(8, 12) = 1.6707481111324018E+03_dp
      fit_coef(9, 12) = -1.7524442770987675E+03_dp
      fit_coef(10, 12) = 1.1932577680169172E+03_dp
      fit_coef(11, 12) = -4.7477746229497609E+02_dp
      fit_coef(12, 12) = 8.3723024185087226E+01_dp
      fit_coef(0, 13) = 1.2328742147617662E+00_dp
      fit_coef(1, 13) = 1.1464277359986964E-03_dp
      fit_coef(2, 13) = 1.5973172261127475E+00_dp
      fit_coef(3, 13) = -9.5083265310456007E+00_dp
      fit_coef(4, 13) = 5.4875950711264657E+01_dp
      fit_coef(5, 13) = -2.3203772329691336E+02_dp
      fit_coef(6, 13) = 7.0989929741104845E+02_dp
      fit_coef(7, 13) = -1.5517453595579404E+03_dp
      fit_coef(8, 13) = 2.3823114822950174E+03_dp
      fit_coef(9, 13) = -2.4989178165580306E+03_dp
      fit_coef(10, 13) = 1.7017491667373313E+03_dp
      fit_coef(11, 13) = -6.7719904341330835E+02_dp
      fit_coef(12, 13) = 1.1943636331258227E+02_dp
      fit_coef(0, 14) = 3.1933956854310037E+00_dp
      fit_coef(1, 14) = 1.3428436858287465E-03_dp
      fit_coef(2, 14) = 1.8732506344005615E+00_dp
      fit_coef(3, 14) = -1.1155445162652324E+01_dp
      fit_coef(4, 14) = 6.4329861295470778E+01_dp
      fit_coef(5, 14) = -2.7197046053760539E+02_dp
      fit_coef(6, 14) = 8.3192845835379330E+02_dp
      fit_coef(7, 14) = -1.8182993373416557E+03_dp
      fit_coef(8, 14) = 2.7914396608814409E+03_dp
      fit_coef(9, 14) = -2.9280874922129306E+03_dp
      fit_coef(10, 14) = 1.9940610145207027E+03_dp
      fit_coef(11, 14) = -7.9354809788071543E+02_dp
      fit_coef(12, 14) = 1.3996131946484607E+02_dp

   END SUBROUTINE get_fit_coef_k7_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k8_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 7.3303678789083827E-03_dp
      fit_coef(1, 1) = 1.9068147607663823E-02_dp
      fit_coef(2, 1) = 2.4523188619260956E-02_dp
      fit_coef(3, 1) = 1.7843299989374534E-02_dp
      fit_coef(4, 1) = 5.4214701717497675E-03_dp
      fit_coef(5, 1) = -2.4694578878648259E-03_dp
      fit_coef(6, 1) = -3.0284493938072380E-03_dp
      fit_coef(7, 1) = -1.6609340715173272E-05_dp
      fit_coef(8, 1) = 5.8996552393493821E-04_dp
      fit_coef(9, 1) = -5.5212702213677864E-05_dp
      fit_coef(10, 1) = 3.4862398750072344E-04_dp
      fit_coef(11, 1) = -3.9218355438884901E-04_dp
      fit_coef(12, 1) = 1.2569394142797241E-04_dp
      fit_coef(0, 2) = 4.1899841856663761E-02_dp
      fit_coef(1, 2) = 1.0302234723218924E-01_dp
      fit_coef(2, 2) = 1.3023821601879626E-01_dp
      fit_coef(3, 2) = 9.3099102625574220E-02_dp
      fit_coef(4, 2) = 2.7404221544145291E-02_dp
      fit_coef(5, 2) = -1.3085499573245986E-02_dp
      fit_coef(6, 2) = -1.5582408546316897E-02_dp
      fit_coef(7, 2) = 5.2970248514313688E-04_dp
      fit_coef(8, 2) = 2.2609393067795900E-03_dp
      fit_coef(9, 2) = 4.6279627025403949E-04_dp
      fit_coef(10, 2) = 1.4146159928236672E-03_dp
      fit_coef(11, 2) = -1.9800386255810288E-03_dp
      fit_coef(12, 2) = 6.6847497768620658E-04_dp
      fit_coef(0, 3) = 1.2090283026330555E-01_dp
      fit_coef(1, 3) = 2.6358046342322966E-01_dp
      fit_coef(2, 3) = 3.2173608677475352E-01_dp
      fit_coef(3, 3) = 2.2555617000622660E-01_dp
      fit_coef(4, 3) = 6.4242254728770876E-02_dp
      fit_coef(5, 3) = -3.3531313583337982E-02_dp
      fit_coef(6, 3) = -3.7471068968577388E-02_dp
      fit_coef(7, 3) = 3.6816061058378741E-03_dp
      fit_coef(8, 3) = 6.5687067201372204E-03_dp
      fit_coef(9, 3) = -2.6384863361432163E-03_dp
      fit_coef(10, 3) = 6.4385556221874555E-03_dp
      fit_coef(11, 3) = -6.0511747469316332E-03_dp
      fit_coef(12, 3) = 1.8525016625111832E-03_dp
      fit_coef(0, 4) = 2.8592429103323302E-01_dp
      fit_coef(1, 4) = 5.1841367625987889E-01_dp
      fit_coef(2, 4) = 5.9115684453736916E-01_dp
      fit_coef(3, 4) = 4.0105593630503783E-01_dp
      fit_coef(4, 4) = 1.1948149198856918E-01_dp
      fit_coef(5, 4) = -5.2658401069846560E-02_dp
      fit_coef(6, 4) = -7.3898713202086491E-02_dp
      fit_coef(7, 4) = 9.1566644314021242E-03_dp
      fit_coef(8, 4) = -1.1867605759483740E-02_dp
      fit_coef(9, 4) = 4.8318451940700932E-02_dp
      fit_coef(10, 4) = -3.5186576456035019E-02_dp
      fit_coef(11, 4) = 7.7443536084135598E-03_dp
      fit_coef(12, 4) = 6.1488354568711538E-04_dp
      fit_coef(0, 5) = 6.2290281260653402E-01_dp
      fit_coef(1, 5) = 8.9580902133545759E-01_dp
      fit_coef(2, 5) = 9.2609808946198835E-01_dp
      fit_coef(3, 5) = 5.8224600441757135E-01_dp
      fit_coef(4, 5) = 1.7839482993919215E-01_dp
      fit_coef(5, 5) = -4.6849100928790438E-02_dp
      fit_coef(6, 5) = -8.6070770575879779E-02_dp
      fit_coef(7, 5) = -1.8640277231677663E-02_dp
      fit_coef(8, 5) = 1.1471460387476767E-02_dp
      fit_coef(9, 5) = -4.1058254270079884E-02_dp
      fit_coef(10, 5) = 1.1306530929746011E-01_dp
      fit_coef(11, 5) = -8.8348931981328011E-02_dp
      fit_coef(12, 5) = 2.3298598302007475E-02_dp
      fit_coef(0, 6) = 1.2950331297092936E+00_dp
      fit_coef(1, 6) = 1.4208263619906618E+00_dp
      fit_coef(2, 6) = 1.3168706585134538E+00_dp
      fit_coef(3, 6) = 7.3735161688744788E-01_dp
      fit_coef(4, 6) = 2.1308021650008105E-01_dp
      fit_coef(5, 6) = -2.2140636924009825E-02_dp
      fit_coef(6, 6) = -6.5523715696762824E-02_dp
      fit_coef(7, 6) = -5.5227429757704954E-03_dp
      fit_coef(8, 6) = -3.6319163794681919E-02_dp
      fit_coef(9, 6) = -2.1057002605029421E-02_dp
      fit_coef(10, 6) = 6.7929722094043590E-02_dp
      fit_coef(11, 6) = -2.9821898442364107E-02_dp
      fit_coef(12, 6) = 2.5847668701003404E-03_dp
      fit_coef(0, 7) = 2.6150915219181412E+00_dp
      fit_coef(1, 7) = 2.0925013488168447E+00_dp
      fit_coef(2, 7) = 1.7517994308499847E+00_dp
      fit_coef(3, 7) = 8.5858648940090110E-01_dp
      fit_coef(4, 7) = 2.1877243019911102E-01_dp
      fit_coef(5, 7) = 9.2380438422983800E-03_dp
      fit_coef(6, 7) = -3.5233332041510933E-02_dp
      fit_coef(7, 7) = 5.9579661965596109E-02_dp
      fit_coef(8, 7) = -1.5530877010597718E-01_dp
      fit_coef(9, 7) = 1.1970417799859782E-01_dp
      fit_coef(10, 7) = -7.3859936805166476E-02_dp
      fit_coef(11, 7) = 3.9401452538859211E-02_dp
      fit_coef(12, 7) = -8.0701941615812720E-03_dp
      fit_coef(0, 8) = 5.2886811819505075E+00_dp
      fit_coef(1, 8) = 2.8791408561019143E+00_dp
      fit_coef(2, 8) = 2.2282496597143453E+00_dp
      fit_coef(3, 8) = 9.6463274698100221E-01_dp
      fit_coef(4, 8) = 2.0804339162869817E-01_dp
      fit_coef(5, 8) = 3.8463667107493361E-02_dp
      fit_coef(6, 8) = 6.2738138193085286E-03_dp
      fit_coef(7, 8) = 1.0100983178342585E-01_dp
      fit_coef(8, 8) = -1.8879341637967198E-01_dp
      fit_coef(9, 8) = 1.2577637243618989E-01_dp
      fit_coef(10, 8) = -4.4240946155634182E-02_dp
      fit_coef(11, 8) = -2.7654276880568853E-03_dp
      fit_coef(12, 8) = 7.1269217913952432E-03_dp
      fit_coef(0, 9) = 1.9153697185214646E-02_dp
      fit_coef(1, 9) = 4.9243996218940246E-02_dp
      fit_coef(2, 9) = 6.3063039313140626E-02_dp
      fit_coef(3, 9) = 4.5675586741478570E-02_dp
      fit_coef(4, 9) = 1.3789973749172879E-02_dp
      fit_coef(5, 9) = -6.3486543694890422E-03_dp
      fit_coef(6, 9) = -7.7438232796799697E-03_dp
      fit_coef(7, 9) = 4.2857574556254048E-05_dp
      fit_coef(8, 9) = 1.4170649293170209E-03_dp
      fit_coef(9, 9) = -6.5915991195257174E-05_dp
      fit_coef(10, 9) = 8.5657276913500048E-04_dp
      fit_coef(11, 9) = -1.0012786918386677E-03_dp
      fit_coef(12, 9) = 3.2416699208711279E-04_dp
      fit_coef(0, 10) = 5.2487499010059821E-02_dp
      fit_coef(1, 10) = 1.2020367140285117E-01_dp
      fit_coef(2, 10) = 1.4885896496539552E-01_dp
      fit_coef(3, 10) = 1.0446916703848536E-01_dp
      fit_coef(4, 10) = 2.9452247957851872E-02_dp
      fit_coef(5, 10) = -1.5096946389090963E-02_dp
      fit_coef(6, 10) = -1.6978769798053302E-02_dp
      fit_coef(7, 10) = 1.5621840217608987E-03_dp
      fit_coef(8, 10) = 1.1498062866988877E-03_dp
      fit_coef(9, 10) = 1.5932508408093571E-03_dp
      fit_coef(10, 10) = 1.2044209864296811E-03_dp
      fit_coef(11, 10) = -2.2799861036947064E-03_dp
      fit_coef(12, 10) = 8.0794913771834831E-04_dp
      fit_coef(0, 11) = 1.1214831732759308E-01_dp
      fit_coef(1, 11) = 2.0374016088705749E-01_dp
      fit_coef(2, 11) = 2.3269755141858203E-01_dp
      fit_coef(3, 11) = 1.5862165696203001E-01_dp
      fit_coef(4, 11) = 4.5536019398208713E-02_dp
      fit_coef(5, 11) = -2.4771937060547151E-02_dp
      fit_coef(6, 11) = -2.8674783757891342E-02_dp
      fit_coef(7, 11) = 3.5158594299659887E-03_dp
      fit_coef(8, 11) = 9.3535039711395342E-03_dp
      fit_coef(9, 11) = -3.0604392897587429E-03_dp
      fit_coef(10, 11) = -6.0143657962914525E-04_dp
      fit_coef(11, 11) = 4.2950330778982355E-04_dp
      fit_coef(12, 11) = 1.9255873491051108E-05_dp
      fit_coef(0, 12) = 2.3171781907279218E-01_dp
      fit_coef(1, 12) = 3.1091079731599292E-01_dp
      fit_coef(2, 12) = 3.0389179989885584E-01_dp
      fit_coef(3, 12) = 1.8532205840012622E-01_dp
      fit_coef(4, 12) = 6.2306830325449335E-02_dp
      fit_coef(5, 12) = -7.6106859087436904E-03_dp
      fit_coef(6, 12) = -3.6599453891127606E-02_dp
      fit_coef(7, 12) = 1.5509386794819603E-03_dp
      fit_coef(8, 12) = -4.0853132735310875E-02_dp
      fit_coef(9, 12) = 6.6595641401483052E-02_dp
      fit_coef(10, 12) = -1.8834300269643112E-02_dp
      fit_coef(11, 12) = -1.5460640921107466E-02_dp
      fit_coef(12, 12) = 8.1050164047865251E-03_dp
      fit_coef(0, 13) = 4.6839273718001340E-01_dp
      fit_coef(1, 13) = 4.4848685344087941E-01_dp
      fit_coef(2, 13) = 3.6456593426754630E-01_dp
      fit_coef(3, 13) = 1.7154603893334677E-01_dp
      fit_coef(4, 13) = 5.0029581620847113E-02_dp
      fit_coef(5, 13) = 1.7043310105799386E-02_dp
      fit_coef(6, 13) = 1.3624711461006474E-02_dp
      fit_coef(7, 13) = -4.0467589755335792E-02_dp
      fit_coef(8, 13) = 6.0917156869027100E-02_dp
      fit_coef(9, 13) = -1.7391850945484083E-01_dp
      fit_coef(10, 13) = 1.9540773247888954E-01_dp
      fit_coef(11, 13) = -8.8667454786059030E-02_dp
      fit_coef(12, 13) = 1.3511223411471382E-02_dp
      fit_coef(0, 14) = 9.2519886811346685E-01_dp
      fit_coef(1, 14) = 6.0160298038546722E-01_dp
      fit_coef(2, 14) = 4.1486918310874005E-01_dp
      fit_coef(3, 14) = 1.3735307260002044E-01_dp
      fit_coef(4, 14) = 1.8815781506297060E-02_dp
      fit_coef(5, 14) = 3.0299118585638541E-02_dp
      fit_coef(6, 14) = 2.2993724328518487E-02_dp
      fit_coef(7, 14) = 6.5423665107552431E-02_dp
      fit_coef(8, 14) = -1.4716829076991533E-01_dp
      fit_coef(9, 14) = 1.9132722353109541E-01_dp
      fit_coef(10, 14) = -2.2566102318270043E-01_dp
      fit_coef(11, 14) = 1.4512891962640706E-01_dp
      fit_coef(12, 14) = -3.4869303190922442E-02_dp
      fit_coef(0, 15) = 1.8180761389140441E+00_dp
      fit_coef(1, 15) = 7.3545411157338625E-01_dp
      fit_coef(2, 15) = 4.5382578616281255E-01_dp
      fit_coef(3, 15) = 1.0815904182157396E-01_dp
      fit_coef(4, 15) = -5.0564441213593594E-03_dp
      fit_coef(5, 15) = 3.1009296969077782E-02_dp
      fit_coef(6, 15) = 3.8543619111406462E-02_dp
      fit_coef(7, 15) = 4.9577607054581169E-02_dp
      fit_coef(8, 15) = -6.0901972233113555E-02_dp
      fit_coef(9, 15) = 4.4164548268671569E-02_dp
      fit_coef(10, 15) = -1.6420936751689871E-02_dp
      fit_coef(11, 15) = -1.9128788572305898E-02_dp
      fit_coef(12, 15) = 1.2610902336520378E-02_dp
      fit_coef(0, 16) = 3.8678537708718710E+00_dp
      fit_coef(1, 16) = 8.3755983422995561E-01_dp
      fit_coef(2, 16) = 5.0744880083562161E-01_dp
      fit_coef(3, 16) = 1.1252857834008548E-01_dp
      fit_coef(4, 16) = -1.3783956985565306E-02_dp
      fit_coef(5, 16) = 2.6967047539614656E-02_dp
      fit_coef(6, 16) = 4.1533124322160005E-02_dp
      fit_coef(7, 16) = 4.2761766824422134E-02_dp
      fit_coef(8, 16) = -3.1110880765649707E-02_dp
      fit_coef(9, 16) = 1.0013876986386597E-02_dp
      fit_coef(10, 16) = 2.7350239856939598E-02_dp
      fit_coef(11, 16) = -3.3885317157732221E-02_dp
      fit_coef(12, 16) = 1.0044298911104206E-02_dp

   END SUBROUTINE get_fit_coef_k8_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k8_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.3749000000000000E+04_dp
      fit_coef(0, 1) = 2.8123802981039521E-04_dp
      fit_coef(1, 1) = 1.5926645040140255E-05_dp
      fit_coef(2, 1) = 1.0329030933196252E-02_dp
      fit_coef(3, 1) = -6.4007908632844474E-02_dp
      fit_coef(4, 1) = 4.1161308446326278E-01_dp
      fit_coef(5, 1) = -1.7624808185774974E+00_dp
      fit_coef(6, 1) = 5.3986045573131776E+00_dp
      fit_coef(7, 1) = -1.1672753089064145E+01_dp
      fit_coef(8, 1) = 1.7591571465892009E+01_dp
      fit_coef(9, 1) = -1.8028596170855813E+01_dp
      fit_coef(10, 1) = 1.1956816057777168E+01_dp
      fit_coef(11, 1) = -4.6228922203591667E+00_dp
      fit_coef(12, 1) = 7.9068650631410775E-01_dp
      fit_coef(0, 2) = 3.0940667380651658E-03_dp
      fit_coef(1, 2) = 9.4603372213690319E-05_dp
      fit_coef(2, 2) = 6.2371593714423038E-02_dp
      fit_coef(3, 2) = -3.9061079256549203E-01_dp
      fit_coef(4, 2) = 2.4800148089165819E+00_dp
      fit_coef(5, 2) = -1.0590097222199683E+01_dp
      fit_coef(6, 2) = 3.2358288023437467E+01_dp
      fit_coef(7, 2) = -6.9843326031928697E+01_dp
      fit_coef(8, 2) = 1.0512898764379077E+02_dp
      fit_coef(9, 2) = -1.0764631696630022E+02_dp
      fit_coef(10, 2) = 7.1348387712368961E+01_dp
      fit_coef(11, 2) = -2.7573670968312566E+01_dp
      fit_coef(12, 2) = 4.7147055876112587E+00_dp
      fit_coef(0, 3) = 1.6375908158103483E-02_dp
      fit_coef(1, 3) = 3.1160349661731832E-04_dp
      fit_coef(2, 3) = 2.0927456992109450E-01_dp
      fit_coef(3, 3) = -1.3273785824948792E+00_dp
      fit_coef(4, 3) = 8.2688701711323898E+00_dp
      fit_coef(5, 3) = -3.5220567469598883E+01_dp
      fit_coef(6, 3) = 1.0701331897299248E+02_dp
      fit_coef(7, 3) = -2.2989760166900771E+02_dp
      fit_coef(8, 3) = 3.4493893446724235E+02_dp
      fit_coef(9, 3) = -3.5248615129301214E+02_dp
      fit_coef(10, 3) = 2.3334466251650863E+02_dp
      fit_coef(11, 3) = -9.0116552857752083E+01_dp
      fit_coef(12, 3) = 1.5403005349166568E+01_dp
      fit_coef(0, 4) = 6.3206454085883088E-02_dp
      fit_coef(1, 4) = 8.1154192163640840E-04_dp
      fit_coef(2, 4) = 5.5144861056801797E-01_dp
      fit_coef(3, 4) = -3.5260842601451339E+00_dp
      fit_coef(4, 4) = 2.1693237689818268E+01_dp
      fit_coef(5, 4) = -9.2258467452491558E+01_dp
      fit_coef(6, 4) = 2.7922864280490649E+02_dp
      fit_coef(7, 4) = -5.9806470556712168E+02_dp
      fit_coef(8, 4) = 8.9553395470219289E+02_dp
      fit_coef(9, 4) = -9.1394179906769205E+02_dp
      fit_coef(10, 4) = 6.0454016322906523E+02_dp
      fit_coef(11, 4) = -2.3336202357129835E+02_dp
      fit_coef(12, 4) = 3.9877696019475444E+01_dp
      fit_coef(0, 5) = 2.0243266526935697E-01_dp
      fit_coef(1, 5) = 1.8020515157672840E-03_dp
      fit_coef(2, 5) = 1.2344900509541614E+00_dp
      fit_coef(3, 5) = -7.9367169391085666E+00_dp
      fit_coef(4, 5) = 4.8422935474594851E+01_dp
      fit_coef(5, 5) = -2.0571215987444512E+02_dp
      fit_coef(6, 5) = 6.2097865514974285E+02_dp
      fit_coef(7, 5) = -1.3272491727851836E+03_dp
      fit_coef(8, 5) = 1.9845109548223063E+03_dp
      fit_coef(9, 5) = -2.0233745953849727E+03_dp
      fit_coef(10, 5) = 1.3376052915228147E+03_dp
      fit_coef(11, 5) = -5.1616417908155108E+02_dp
      fit_coef(12, 5) = 8.8189346087826550E+01_dp
      fit_coef(0, 6) = 5.7157821118623775E-01_dp
      fit_coef(1, 6) = 3.5076317937191783E-03_dp
      fit_coef(2, 6) = 2.4171981047624476E+00_dp
      fit_coef(3, 6) = -1.5600811633970423E+01_dp
      fit_coef(4, 6) = 9.4627104824817778E+01_dp
      fit_coef(5, 6) = -4.0167179845150463E+02_dp
      fit_coef(6, 6) = 1.2102930856782566E+03_dp
      fit_coef(7, 6) = -2.5828405060244945E+03_dp
      fit_coef(8, 6) = 3.8576577717746441E+03_dp
      fit_coef(9, 6) = -3.9303584111835289E+03_dp
      fit_coef(10, 6) = 2.5970876357091988E+03_dp
      fit_coef(11, 6) = -1.0019128608117558E+03_dp
      fit_coef(12, 6) = 1.7115745096843858E+02_dp
      fit_coef(0, 7) = 1.4769845791046401E+00_dp
      fit_coef(1, 7) = 6.0305390650928353E-03_dp
      fit_coef(2, 7) = 4.1735972990944301E+00_dp
      fit_coef(3, 7) = -2.7010363741735947E+01_dp
      fit_coef(4, 7) = 1.6316505318502820E+02_dp
      fit_coef(5, 7) = -6.9218832414970177E+02_dp
      fit_coef(6, 7) = 2.0829974158076584E+03_dp
      fit_coef(7, 7) = -4.4403073409985855E+03_dp
      fit_coef(8, 7) = 6.6265908711188122E+03_dp
      fit_coef(9, 7) = -6.7478402025357373E+03_dp
      fit_coef(10, 7) = 4.4572806839765526E+03_dp
      fit_coef(11, 7) = -1.7191849878725106E+03_dp
      fit_coef(12, 7) = 2.9365492248421521E+02_dp
      fit_coef(0, 8) = 3.6498674343702047E+00_dp
      fit_coef(1, 8) = 9.1861328914861536E-03_dp
      fit_coef(2, 8) = 6.3746908376074938E+00_dp
      fit_coef(3, 8) = -4.1325382231508740E+01_dp
      fit_coef(4, 8) = 2.4901070747707607E+02_dp
      fit_coef(5, 8) = -1.0559668335879119E+03_dp
      fit_coef(6, 8) = 3.1752073241924800E+03_dp
      fit_coef(7, 8) = -6.7638081975788209E+03_dp
      fit_coef(8, 8) = 1.0088934681089724E+04_dp
      fit_coef(9, 8) = -1.0269971775238855E+04_dp
      fit_coef(10, 8) = 6.7823044591615362E+03_dp
      fit_coef(11, 8) = -2.6155916262455330E+03_dp
      fit_coef(12, 8) = 4.4673488306700750E+02_dp
      fit_coef(0, 9) = 8.9908324754226225E-04_dp
      fit_coef(1, 9) = 4.1044297500840920E-05_dp
      fit_coef(2, 9) = 2.6592960271921522E-02_dp
      fit_coef(3, 9) = -1.6475343310038174E-01_dp
      fit_coef(4, 9) = 1.0601282924777382E+00_dp
      fit_coef(5, 9) = -4.5352614218846243E+00_dp
      fit_coef(6, 9) = 1.3894117699928159E+01_dp
      fit_coef(7, 9) = -3.0070697693720529E+01_dp
      fit_coef(8, 9) = 4.5371751202142313E+01_dp
      fit_coef(9, 9) = -4.6549283112947627E+01_dp
      fit_coef(10, 9) = 3.0899823469725881E+01_dp
      fit_coef(11, 9) = -1.1955338603808340E+01_dp
      fit_coef(12, 9) = 2.0459285416163961E+00_dp
      fit_coef(0, 10) = 5.9074268743020538E-03_dp
      fit_coef(1, 10) = 1.2917982248413781E-04_dp
      fit_coef(2, 10) = 8.6643119919001477E-02_dp
      fit_coef(3, 10) = -5.4867858026705008E-01_dp
      fit_coef(4, 10) = 3.4301753057411557E+00_dp
      fit_coef(5, 10) = -1.4614938748597057E+01_dp
      fit_coef(6, 10) = 4.4473633864726359E+01_dp
      fit_coef(7, 10) = -9.5626734151356388E+01_dp
      fit_coef(8, 10) = 1.4351698350787132E+02_dp
      fit_coef(9, 10) = -1.4665263177842959E+02_dp
      fit_coef(10, 10) = 9.7069655720786628E+01_dp
      fit_coef(11, 10) = -3.7481250052554522E+01_dp
      fit_coef(12, 10) = 6.4052741628428116E+00_dp
      fit_coef(0, 11) = 2.4183310977353760E-02_dp
      fit_coef(1, 11) = 3.2874969051224556E-04_dp
      fit_coef(2, 11) = 2.2379265645203872E-01_dp
      fit_coef(3, 11) = -1.4327470925490047E+00_dp
      fit_coef(4, 11) = 8.7976488483593407E+00_dp
      fit_coef(5, 11) = -3.7406648187178646E+01_dp
      fit_coef(6, 11) = 1.1314765967557075E+02_dp
      fit_coef(7, 11) = -2.4223242835221623E+02_dp
      fit_coef(8, 11) = 3.6262537126551155E+02_dp
      fit_coef(9, 11) = -3.7003484902257611E+02_dp
      fit_coef(10, 11) = 2.4474926980111050E+02_dp
      fit_coef(11, 11) = -9.4473072063769564E+01_dp
      fit_coef(12, 11) = 1.6143361281904095E+01_dp
      fit_coef(0, 12) = 7.8618264891410383E-02_dp
      fit_coef(1, 12) = 7.0637070558397471E-04_dp
      fit_coef(2, 12) = 4.8525758636061883E-01_dp
      fit_coef(3, 12) = -3.1258030882321890E+00_dp
      fit_coef(4, 12) = 1.9013088969838954E+01_dp
      fit_coef(5, 12) = -8.0742838543906203E+01_dp
      fit_coef(6, 12) = 2.4349690163060086E+02_dp
      fit_coef(7, 12) = -5.2005299814464468E+02_dp
      fit_coef(8, 12) = 7.7720399577073124E+02_dp
      fit_coef(9, 12) = -7.9217171146497253E+02_dp
      fit_coef(10, 12) = 5.2358376975962551E+02_dp
      fit_coef(11, 12) = -2.0202246440009665E+02_dp
      fit_coef(12, 12) = 3.4515024598323109E+01_dp
      fit_coef(0, 13) = 2.2139783616846259E-01_dp
      fit_coef(1, 13) = 1.3144046612244267E-03_dp
      fit_coef(2, 13) = 9.0893314563834471E-01_dp
      fit_coef(3, 13) = -5.8797685805649715E+00_dp
      fit_coef(4, 13) = 3.5538337045891197E+01_dp
      fit_coef(5, 13) = -1.5078315089899843E+02_dp
      fit_coef(6, 13) = 4.5381740391672264E+02_dp
      fit_coef(7, 13) = -9.6759884283658630E+02_dp
      fit_coef(8, 13) = 1.4442703065970643E+03_dp
      fit_coef(9, 13) = -1.4708851169077066E+03_dp
      fit_coef(10, 13) = 9.7168186776942412E+02_dp
      fit_coef(11, 13) = -3.7480558358376271E+02_dp
      fit_coef(12, 13) = 6.4023927170066514E+01_dp
      fit_coef(0, 14) = 5.6481431024685758E-01_dp
      fit_coef(1, 14) = 2.1175592913405812E-03_dp
      fit_coef(2, 14) = 1.4715067019515555E+00_dp
      fit_coef(3, 14) = -9.5481410751299673E+00_dp
      fit_coef(4, 14) = 5.7450639914728498E+01_dp
      fit_coef(5, 14) = -2.4358417666644485E+02_dp
      fit_coef(6, 14) = 7.3209348863097784E+02_dp
      fit_coef(7, 14) = -1.5589336783178264E+03_dp
      fit_coef(8, 14) = 2.3247369231211455E+03_dp
      fit_coef(9, 14) = -2.3660752817146231E+03_dp
      fit_coef(10, 14) = 1.5624107971537085E+03_dp
      fit_coef(11, 14) = -6.0251181283667290E+02_dp
      fit_coef(12, 14) = 1.0290483462401868E+02_dp
      fit_coef(0, 15) = 1.3539573260602098E+00_dp
      fit_coef(1, 15) = 2.8953150568411467E-03_dp
      fit_coef(2, 15) = 2.0183462102822127E+00_dp
      fit_coef(3, 15) = -1.3121809745507548E+01_dp
      fit_coef(4, 15) = 7.8730592243301871E+01_dp
      fit_coef(5, 15) = -3.3365664063421633E+02_dp
      fit_coef(6, 15) = 1.0019300382950731E+03_dp
      fit_coef(7, 15) = -2.1317727140515958E+03_dp
      fit_coef(8, 15) = 3.1770143228209445E+03_dp
      fit_coef(9, 15) = -3.2321383186557791E+03_dp
      fit_coef(10, 15) = 2.1337098849426625E+03_dp
      fit_coef(11, 15) = -8.2267510237563181E+02_dp
      fit_coef(12, 15) = 1.4049188918686588E+02_dp
      fit_coef(0, 16) = 3.3347047868343913E+00_dp
      fit_coef(1, 16) = 3.3620771239965489E-03_dp
      fit_coef(2, 16) = 2.3453111222189276E+00_dp
      fit_coef(3, 16) = -1.5253659417099389E+01_dp
      fit_coef(4, 16) = 9.1468857073866545E+01_dp
      fit_coef(5, 16) = -3.8760249737899488E+02_dp
      fit_coef(6, 16) = 1.1637180704439147E+03_dp
      fit_coef(7, 16) = -2.4755718153726775E+03_dp
      fit_coef(8, 16) = 3.6888887814157529E+03_dp
      fit_coef(9, 16) = -3.7525428895000127E+03_dp
      fit_coef(10, 16) = 2.4771015308842470E+03_dp
      fit_coef(11, 16) = -9.5503423745582688E+02_dp
      fit_coef(12, 16) = 1.6309118629666474E+02_dp

   END SUBROUTINE get_fit_coef_k8_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k9_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 6.4552940952308244E-03_dp
      fit_coef(1, 1) = 1.6981103140462521E-02_dp
      fit_coef(2, 1) = 2.1954315188636421E-02_dp
      fit_coef(3, 1) = 1.6041342676653195E-02_dp
      fit_coef(4, 1) = 4.8873501014382193E-03_dp
      fit_coef(5, 1) = -2.2607549832820450E-03_dp
      fit_coef(6, 1) = -2.5816889517990111E-03_dp
      fit_coef(7, 1) = -3.3177472147924311E-04_dp
      fit_coef(8, 1) = 6.7633586217191212E-04_dp
      fit_coef(9, 1) = 5.0154491931341720E-04_dp
      fit_coef(10, 1) = -7.5477037926869911E-04_dp
      fit_coef(11, 1) = 4.3750564946997106E-04_dp
      fit_coef(12, 1) = -1.0884574145659750E-04_dp
      fit_coef(0, 2) = 3.6158523047046677E-02_dp
      fit_coef(1, 2) = 9.1448302892550112E-02_dp
      fit_coef(2, 2) = 1.1643393031728949E-01_dp
      fit_coef(3, 2) = 8.3846218207491621E-02_dp
      fit_coef(4, 2) = 2.4929838698097138E-02_dp
      fit_coef(5, 2) = -1.2000438048214228E-02_dp
      fit_coef(6, 2) = -1.3233969051869923E-02_dp
      fit_coef(7, 2) = -1.6781088808427976E-03_dp
      fit_coef(8, 2) = 3.7175799646608039E-03_dp
      fit_coef(9, 2) = 2.3200967656010554E-03_dp
      fit_coef(10, 2) = -3.7296821680115493E-03_dp
      fit_coef(11, 2) = 2.1900296497794772E-03_dp
      fit_coef(12, 2) = -5.4966296541064880E-04_dp
      fit_coef(0, 3) = 1.0021175251987058E-01_dp
      fit_coef(1, 3) = 2.3290658157250768E-01_dp
      fit_coef(2, 3) = 2.8826132533071419E-01_dp
      fit_coef(3, 3) = 2.0349781078404069E-01_dp
      fit_coef(4, 3) = 5.8186389475606060E-02_dp
      fit_coef(5, 3) = -3.0264996373773254E-02_dp
      fit_coef(6, 3) = -3.0958058136578551E-02_dp
      fit_coef(7, 3) = -3.3385124438745062E-03_dp
      fit_coef(8, 3) = 9.8293593154800095E-03_dp
      fit_coef(9, 3) = 3.3092123173594926E-03_dp
      fit_coef(10, 3) = -6.9541522350748128E-03_dp
      fit_coef(11, 3) = 4.3124239191782493E-03_dp
      fit_coef(12, 3) = -1.1315599631624136E-03_dp
      fit_coef(0, 4) = 2.2411697245587953E-01_dp
      fit_coef(1, 4) = 4.5426701810076259E-01_dp
      fit_coef(2, 4) = 5.3412584957488496E-01_dp
      fit_coef(3, 4) = 3.6826576404725364E-01_dp
      fit_coef(4, 4) = 1.0602701865363849E-01_dp
      fit_coef(5, 4) = -5.4577502966425870E-02_dp
      fit_coef(6, 4) = -5.8813664835258211E-02_dp
      fit_coef(7, 4) = -6.1835169673331971E-03_dp
      fit_coef(8, 4) = 1.9342657781255358E-02_dp
      fit_coef(9, 4) = 1.5715487208569516E-02_dp
      fit_coef(10, 4) = -2.9256796795427198E-02_dp
      fit_coef(11, 4) = 1.8078340581926571E-02_dp
      fit_coef(12, 4) = -4.4135712366918118E-03_dp
      fit_coef(0, 5) = 4.5965700239101925E-01_dp
      fit_coef(1, 5) = 7.7718614768161642E-01_dp
      fit_coef(2, 5) = 8.4369237846467748E-01_dp
      fit_coef(3, 5) = 5.5382304325752407E-01_dp
      fit_coef(4, 5) = 1.6629027091834961E-01_dp
      fit_coef(5, 5) = -6.3214899288962637E-02_dp
      fit_coef(6, 5) = -8.7377144006966334E-02_dp
      fit_coef(7, 5) = -1.9259962748352861E-02_dp
      fit_coef(8, 5) = 5.7845664647801637E-03_dp
      fit_coef(9, 5) = 4.9430272241004450E-02_dp
      fit_coef(10, 5) = -3.0204735992589018E-02_dp
      fit_coef(11, 5) = 1.1560750816998883E-03_dp
      fit_coef(12, 5) = 2.2620193192966114E-03_dp
      fit_coef(0, 6) = 9.0150523193514287E-01_dp
      fit_coef(1, 6) = 1.2281870858285731E+00_dp
      fit_coef(2, 6) = 1.2075268204406777E+00_dp
      fit_coef(3, 6) = 7.2753542801068738E-01_dp
      fit_coef(4, 6) = 2.1703768564838269E-01_dp
      fit_coef(5, 6) = -4.9393451753158686E-02_dp
      fit_coef(6, 6) = -7.2910241488968047E-02_dp
      fit_coef(7, 6) = -7.2551142234446420E-02_dp
      fit_coef(8, 6) = 8.0982810239361422E-02_dp
      fit_coef(9, 6) = -1.2121381989765857E-01_dp
      fit_coef(10, 6) = 1.5600379928593311E-01_dp
      fit_coef(11, 6) = -8.4966639309547370E-02_dp
      fit_coef(12, 6) = 1.6060557175561714E-02_dp
      fit_coef(0, 7) = 1.7180040458559993E+00_dp
      fit_coef(1, 7) = 1.8243250896998631E+00_dp
      fit_coef(2, 7) = 1.6163577010542527E+00_dp
      fit_coef(3, 7) = 8.7029777481360915E-01_dp
      fit_coef(4, 7) = 2.4164335157245592E-01_dp
      fit_coef(5, 7) = -2.1288880285781901E-02_dp
      fit_coef(6, 7) = -4.1346282161013680E-02_dp
      fit_coef(7, 7) = -5.2892337289649802E-02_dp
      fit_coef(8, 7) = 1.2566897761471612E-02_dp
      fit_coef(9, 7) = -3.1238115025798451E-02_dp
      fit_coef(10, 7) = 8.8197509236568605E-03_dp
      fit_coef(11, 7) = 3.4268320342607482E-02_dp
      fit_coef(12, 7) = -1.7500888643553550E-02_dp
      fit_coef(0, 8) = 3.2202910229369826E+00_dp
      fit_coef(1, 8) = 2.5584863603841099E+00_dp
      fit_coef(2, 8) = 2.0603420885505477E+00_dp
      fit_coef(3, 8) = 9.8186936771714084E-01_dp
      fit_coef(4, 8) = 2.4292352348458562E-01_dp
      fit_coef(5, 8) = 1.0777326484119814E-02_dp
      fit_coef(6, 8) = 3.4596476543145072E-03_dp
      fit_coef(7, 8) = -2.5878994552358304E-02_dp
      fit_coef(8, 8) = -2.8942176969688937E-02_dp
      fit_coef(9, 8) = 1.7377801484123011E-02_dp
      fit_coef(10, 8) = -5.3614378940648612E-02_dp
      fit_coef(11, 8) = 5.9949568881798150E-02_dp
      fit_coef(12, 8) = -1.7982770689698030E-02_dp
      fit_coef(0, 9) = 6.1060543796063369E+00_dp
      fit_coef(1, 9) = 3.4046538736762950E+00_dp
      fit_coef(2, 9) = 2.5438587661473857E+00_dp
      fit_coef(3, 9) = 1.0825640947785635E+00_dp
      fit_coef(4, 9) = 2.3188547903431525E-01_dp
      fit_coef(5, 9) = 4.0838177926133055E-02_dp
      fit_coef(6, 9) = 5.8580319059672786E-02_dp
      fit_coef(7, 9) = -2.5983271132225795E-02_dp
      fit_coef(8, 9) = 1.2322432042298232E-02_dp
      fit_coef(9, 9) = -5.5917635514737024E-02_dp
      fit_coef(10, 9) = 2.7316598068796355E-02_dp
      fit_coef(11, 9) = 1.8442499628175768E-03_dp
      fit_coef(12, 9) = -1.7781863984523483E-03_dp
      fit_coef(0, 10) = 1.6794417808730849E-02_dp
      fit_coef(1, 10) = 4.3812238879122306E-02_dp
      fit_coef(2, 10) = 5.6433388526163281E-02_dp
      fit_coef(3, 10) = 4.1089213098926332E-02_dp
      fit_coef(4, 10) = 1.2450115917928132E-02_dp
      fit_coef(5, 10) = -5.8158386932097572E-03_dp
      fit_coef(6, 10) = -6.5845058729652479E-03_dp
      fit_coef(7, 10) = -8.4117835741441777E-04_dp
      fit_coef(8, 10) = 1.7550704199386239E-03_dp
      fit_coef(9, 10) = 1.2453478991127647E-03_dp
      fit_coef(10, 10) = -1.9039470025569718E-03_dp
      fit_coef(11, 10) = 1.1072558306069514E-03_dp
      fit_coef(12, 10) = -2.7613956295636030E-04_dp
      fit_coef(0, 11) = 4.4208436478371291E-02_dp
      fit_coef(1, 11) = 1.0636068480303659E-01_dp
      fit_coef(2, 11) = 1.3297982738437528E-01_dp
      fit_coef(3, 11) = 9.4166685284241089E-02_dp
      fit_coef(4, 11) = 2.7104946950991506E-02_dp
      fit_coef(5, 11) = -1.3725955515035670E-02_dp
      fit_coef(6, 11) = -1.4441472385392985E-02_dp
      fit_coef(7, 11) = -1.8001873375527177E-03_dp
      fit_coef(8, 11) = 4.3120736020466057E-03_dp
      fit_coef(9, 11) = 2.3403116157563499E-03_dp
      fit_coef(10, 11) = -3.9978958575902293E-03_dp
      fit_coef(11, 11) = 2.3671743796162751E-03_dp
      fit_coef(12, 11) = -5.9737357121890616E-04_dp
      fit_coef(0, 12) = 8.7964729281669957E-02_dp
      fit_coef(1, 12) = 1.7861597760308465E-01_dp
      fit_coef(2, 12) = 2.1020751169557236E-01_dp
      fit_coef(3, 12) = 1.4439439121052172E-01_dp
      fit_coef(4, 12) = 3.9618119968394774E-02_dp
      fit_coef(5, 12) = -2.3084280501379426E-02_dp
      fit_coef(6, 12) = -2.1408162571924363E-02_dp
      fit_coef(7, 12) = -1.6462884108872455E-03_dp
      fit_coef(8, 12) = 1.0736589853496247E-02_dp
      fit_coef(9, 12) = -4.2541437088059438E-03_dp
      fit_coef(10, 12) = -9.3740060359603919E-05_dp
      fit_coef(11, 12) = 1.2300444443144750E-03_dp
      fit_coef(12, 12) = -5.1803647309730399E-04_dp
      fit_coef(0, 13) = 1.6822075290435659E-01_dp
      fit_coef(1, 13) = 2.6777049222376398E-01_dp
      fit_coef(2, 13) = 2.7956089737268164E-01_dp
      fit_coef(3, 13) = 1.8086459490735671E-01_dp
      fit_coef(4, 13) = 5.6336421826288820E-02_dp
      fit_coef(5, 13) = -2.0705833559352838E-02_dp
      fit_coef(6, 13) = -3.4790488764663681E-02_dp
      fit_coef(7, 13) = -2.0331621831735812E-03_dp
      fit_coef(8, 13) = -1.1735723784168657E-02_dp
      fit_coef(9, 13) = 6.1673712195088486E-02_dp
      fit_coef(10, 13) = -6.7206330987665763E-02_dp
      fit_coef(11, 13) = 3.3093275404480781E-02_dp
      fit_coef(12, 13) = -6.5493688196066237E-03_dp
      fit_coef(0, 14) = 3.1816462977746751E-01_dp
      fit_coef(1, 14) = 3.8275682615881462E-01_dp
      fit_coef(2, 14) = 3.3806848941519918E-01_dp
      fit_coef(3, 14) = 1.8434916370789287E-01_dp
      fit_coef(4, 14) = 5.9808120500728569E-02_dp
      fit_coef(5, 14) = 4.0170881027690826E-03_dp
      fit_coef(6, 14) = -9.8506182208028918E-03_dp
      fit_coef(7, 14) = -4.1552347574310020E-02_dp
      fit_coef(8, 14) = 4.0058352212536501E-02_dp
      fit_coef(9, 14) = -1.0194992831605759E-01_dp
      fit_coef(10, 14) = 1.6261296510646620E-01_dp
      fit_coef(11, 14) = -1.0779408129071036E-01_dp
      fit_coef(12, 14) = 2.6181546202421097E-02_dp
      fit_coef(0, 15) = 5.9260717069203162E-01_dp
      fit_coef(1, 15) = 5.2236714064676937E-01_dp
      fit_coef(2, 15) = 3.8810235742485544E-01_dp
      fit_coef(3, 15) = 1.5987501622242634E-01_dp
      fit_coef(4, 15) = 3.8676464470426337E-02_dp
      fit_coef(5, 15) = 2.2397530474472072E-02_dp
      fit_coef(6, 15) = 2.9241369156544045E-02_dp
      fit_coef(7, 15) = -2.8237942949318490E-02_dp
      fit_coef(8, 15) = 2.5790590509430979E-02_dp
      fit_coef(9, 15) = -7.6958891581302161E-02_dp
      fit_coef(10, 15) = 3.2763833296732159E-02_dp
      fit_coef(11, 15) = 3.0093663681386516E-02_dp
      fit_coef(12, 15) = -1.8130008571038608E-02_dp
      fit_coef(0, 16) = 1.0894766158803928E+00_dp
      fit_coef(1, 16) = 6.6875957514002238E-01_dp
      fit_coef(2, 16) = 4.2769216255296599E-01_dp
      fit_coef(3, 16) = 1.2558051545388288E-01_dp
      fit_coef(4, 16) = 1.1302801231878518E-02_dp
      fit_coef(5, 16) = 3.1745607467806851E-02_dp
      fit_coef(6, 16) = 3.5806902619106150E-02_dp
      fit_coef(7, 16) = 4.4306027813559613E-02_dp
      fit_coef(8, 16) = -1.0289320673076535E-01_dp
      fit_coef(9, 16) = 1.4936934424649737E-01_dp
      fit_coef(10, 16) = -1.9115632729633894E-01_dp
      fit_coef(11, 16) = 1.1798136645873131E-01_dp
      fit_coef(12, 16) = -2.5669389351849693E-02_dp
      fit_coef(0, 17) = 2.0166575300282217E+00_dp
      fit_coef(1, 17) = 7.9409079746788347E-01_dp
      fit_coef(2, 17) = 4.6030756808902995E-01_dp
      fit_coef(3, 17) = 1.0086447786948728E-01_dp
      fit_coef(4, 17) = -6.7793884282790666E-03_dp
      fit_coef(5, 17) = 3.1537489875332711E-02_dp
      fit_coef(6, 17) = 5.2729459288171666E-02_dp
      fit_coef(7, 17) = 6.4791698755050928E-03_dp
      fit_coef(8, 17) = 2.1264307513768375E-02_dp
      fit_coef(9, 17) = -5.0214257290051294E-02_dp
      fit_coef(10, 17) = 5.5541231917168123E-02_dp
      fit_coef(11, 17) = -4.9625274322895069E-02_dp
      fit_coef(12, 17) = 1.7390208326760141E-02_dp
      fit_coef(0, 18) = 4.0940777457341033E+00_dp
      fit_coef(1, 18) = 9.0230847484016108E-01_dp
      fit_coef(2, 18) = 5.1722017559623490E-01_dp
      fit_coef(3, 18) = 1.0860069828810702E-01_dp
      fit_coef(4, 18) = -1.3407150393913490E-02_dp
      fit_coef(5, 18) = 2.8202598978792404E-02_dp
      fit_coef(6, 18) = 5.5571528490456060E-02_dp
      fit_coef(7, 18) = 1.1439449932901843E-03_dp
      fit_coef(8, 18) = 4.1437472839587285E-02_dp
      fit_coef(9, 18) = -6.3654618071805502E-02_dp
      fit_coef(10, 18) = 6.9641516666139436E-02_dp
      fit_coef(11, 18) = -4.3730228367233441E-02_dp
      fit_coef(12, 18) = 9.7082624411000788E-03_dp

   END SUBROUTINE get_fit_coef_k9_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k9_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 2.8387000000000000E+04_dp
      fit_coef(0, 1) = 1.3617241890792169E-04_dp
      fit_coef(1, 1) = 1.7628977768292345E-05_dp
      fit_coef(2, 1) = 6.8456013917834527E-03_dp
      fit_coef(3, 1) = -4.3895761081829160E-02_dp
      fit_coef(4, 1) = 2.9175828400673420E-01_dp
      fit_coef(5, 1) = -1.2306713809247385E+00_dp
      fit_coef(6, 1) = 3.6801538863592786E+00_dp
      fit_coef(7, 1) = -7.7402871595381262E+00_dp
      fit_coef(8, 1) = 1.1345092910279265E+01_dp
      fit_coef(9, 1) = -1.1320431537786655E+01_dp
      fit_coef(10, 1) = 7.3192031213325555E+00_dp
      fit_coef(11, 1) = -2.7615646170479575E+00_dp
      fit_coef(12, 1) = 4.6126151831636680E-01_dp
      fit_coef(0, 2) = 1.4982089284441154E-03_dp
      fit_coef(1, 2) = 1.0653572616093308E-04_dp
      fit_coef(2, 2) = 4.2011508058427197E-02_dp
      fit_coef(3, 2) = -2.7370993918916198E-01_dp
      fit_coef(4, 2) = 1.7939751803412998E+00_dp
      fit_coef(5, 2) = -7.5589803739294439E+00_dp
      fit_coef(6, 2) = 2.2557089057857102E+01_dp
      fit_coef(7, 2) = -4.7359259599571885E+01_dp
      fit_coef(8, 2) = 6.9309690166002625E+01_dp
      fit_coef(9, 2) = -6.9067026631983438E+01_dp
      fit_coef(10, 2) = 4.4604006663335653E+01_dp
      fit_coef(11, 2) = -1.6812960687505747E+01_dp
      fit_coef(12, 2) = 2.8059777999526059E+00_dp
      fit_coef(0, 3) = 7.9297670024187777E-03_dp
      fit_coef(1, 3) = 3.5670578026461845E-04_dp
      fit_coef(2, 3) = 1.4395693117221828E-01_dp
      fit_coef(3, 3) = -9.5715933171411527E-01_dp
      fit_coef(4, 3) = 6.1575946587200887E+00_dp
      fit_coef(5, 3) = -2.6009500954650985E+01_dp
      fit_coef(6, 3) = 7.7302143694085373E+01_dp
      fit_coef(7, 3) = -1.6146363554826667E+02_dp
      fit_coef(8, 3) = 2.3520423885679529E+02_dp
      fit_coef(9, 3) = -2.3348174715185502E+02_dp
      fit_coef(10, 3) = 1.5031746608470456E+02_dp
      fit_coef(11, 3) = -5.6519157967148246E+01_dp
      fit_coef(12, 3) = 9.4136465031151371E+00_dp
      fit_coef(0, 4) = 3.0606944390022209E-02_dp
      fit_coef(1, 4) = 9.4936819479453946E-04_dp
      fit_coef(2, 4) = 3.8897471198892203E-01_dp
      fit_coef(3, 4) = -2.6197251988016608E+00_dp
      fit_coef(4, 4) = 1.6650781475997896E+01_dp
      fit_coef(5, 4) = -7.0453305665193994E+01_dp
      fit_coef(6, 4) = 2.0877294584192904E+02_dp
      fit_coef(7, 4) = -4.3466952414824391E+02_dp
      fit_coef(8, 4) = 6.3140792713339783E+02_dp
      fit_coef(9, 4) = -6.2527262377050522E+02_dp
      fit_coef(10, 4) = 4.0172259183111106E+02_dp
      fit_coef(11, 4) = -1.5078038623522440E+02_dp
      fit_coef(12, 4) = 2.5075879780026515E+01_dp
      fit_coef(0, 5) = 9.8018953522175933E-02_dp
      fit_coef(1, 5) = 2.1720322644464903E-03_dp
      fit_coef(2, 5) = 8.9930042879493677E-01_dp
      fit_coef(3, 5) = -6.1094240846085883E+00_dp
      fit_coef(4, 5) = 3.8522442672428134E+01_dp
      fit_coef(5, 5) = -1.6318345003883445E+02_dp
      fit_coef(6, 5) = 4.8266158600585260E+02_dp
      fit_coef(7, 5) = -1.0027451259770581E+03_dp
      fit_coef(8, 5) = 1.4536819184010467E+03_dp
      fit_coef(9, 5) = -1.4370066077396380E+03_dp
      fit_coef(10, 5) = 9.2184153436145550E+02_dp
      fit_coef(11, 5) = -3.4555746263788461E+02_dp
      fit_coef(12, 5) = 5.7407542232062980E+01_dp
      fit_coef(0, 6) = 2.7659185448391987E-01_dp
      fit_coef(1, 6) = 4.4094406641574943E-03_dp
      fit_coef(2, 6) = 1.8398642403823151E+00_dp
      fit_coef(3, 6) = -1.2577480695844354E+01_dp
      fit_coef(4, 6) = 7.8862889072705556E+01_dp
      fit_coef(5, 6) = -3.3434171658322168E+02_dp
      fit_coef(6, 6) = 9.8773755594540955E+02_dp
      fit_coef(7, 6) = -2.0489450065104556E+03_dp
      fit_coef(8, 6) = 2.9660415791608848E+03_dp
      fit_coef(9, 6) = -2.9282468772949642E+03_dp
      fit_coef(10, 6) = 1.8764113820048669E+03_dp
      fit_coef(11, 6) = -7.0273375035027368E+02_dp
      fit_coef(12, 6) = 1.1665511108284656E+02_dp
      fit_coef(0, 7) = 7.1175857336521264E-01_dp
      fit_coef(1, 7) = 8.0417831942619766E-03_dp
      fit_coef(2, 7) = 3.3752364445702137E+00_dp
      fit_coef(3, 7) = -2.3180805628186413E+01_dp
      fit_coef(4, 7) = 1.4475705309197897E+02_dp
      fit_coef(5, 7) = -6.1407495752819284E+02_dp
      fit_coef(6, 7) = 1.8127278994439107E+03_dp
      fit_coef(7, 7) = -3.7561578028747299E+03_dp
      fit_coef(8, 7) = 5.4315325157499246E+03_dp
      fit_coef(9, 7) = -5.3571860391096889E+03_dp
      fit_coef(10, 7) = 3.4300556791845847E+03_dp
      fit_coef(11, 7) = -1.2837010342914398E+03_dp
      fit_coef(12, 7) = 2.1297330497738872E+02_dp
      fit_coef(0, 8) = 1.7130988759895882E+00_dp
      fit_coef(1, 8) = 1.3197929578399978E-02_dp
      fit_coef(2, 8) = 5.5632864218169500E+00_dp
      fit_coef(3, 8) = -3.8336822165302195E+01_dp
      fit_coef(4, 8) = 2.3871040314392795E+02_dp
      fit_coef(5, 8) = -1.0130802214079603E+03_dp
      fit_coef(6, 8) = 2.9890590755362450E+03_dp
      fit_coef(7, 8) = -6.1888221737685317E+03_dp
      fit_coef(8, 8) = 8.9422907047984954E+03_dp
      fit_coef(9, 8) = -8.8137663823846488E+03_dp
      fit_coef(10, 8) = 5.6398489234791505E+03_dp
      fit_coef(11, 8) = -2.1096575748430814E+03_dp
      fit_coef(12, 8) = 3.4985649609391197E+02_dp
      fit_coef(0, 9) = 4.0049762380285738E+00_dp
      fit_coef(1, 9) = 1.9523136227240855E-02_dp
      fit_coef(2, 9) = 8.2524358089376726E+00_dp
      fit_coef(3, 9) = -5.6989983733349561E+01_dp
      fit_coef(4, 9) = 3.5421093521195178E+02_dp
      fit_coef(5, 9) = -1.5036827242000363E+03_dp
      fit_coef(6, 9) = 4.4352294689617192E+03_dp
      fit_coef(7, 9) = -9.1786236202633227E+03_dp
      fit_coef(8, 9) = 1.3255721756847372E+04_dp
      fit_coef(9, 9) = -1.3059405760710591E+04_dp
      fit_coef(10, 9) = 8.3534120377636300E+03_dp
      fit_coef(11, 9) = -3.1236948984785458E+03_dp
      fit_coef(12, 9) = 5.1787990413580712E+02_dp
      fit_coef(0, 10) = 4.3533793103490108E-04_dp
      fit_coef(1, 10) = 4.5700924266841888E-05_dp
      fit_coef(2, 10) = 1.7699021404652904E-02_dp
      fit_coef(3, 10) = -1.1339163743096678E-01_dp
      fit_coef(4, 10) = 7.5317562249438585E-01_dp
      fit_coef(5, 10) = -3.1670388958213715E+00_dp
      fit_coef(6, 10) = 9.4488857516948261E+00_dp
      fit_coef(7, 10) = -1.9855290060320243E+01_dp
      fit_coef(8, 10) = 2.9102458261293140E+01_dp
      fit_coef(9, 10) = -2.9053092441811895E+01_dp
      fit_coef(10, 10) = 1.8796923562697987E+01_dp
      fit_coef(11, 10) = -7.0973082087264441E+00_dp
      fit_coef(12, 10) = 1.1862930929920559E+00_dp
      fit_coef(0, 11) = 2.8605561196639808E-03_dp
      fit_coef(1, 11) = 1.4711961159333449E-04_dp
      fit_coef(2, 11) = 5.9181533881525895E-02_dp
      fit_coef(3, 11) = -3.9234994207703527E-01_dp
      fit_coef(4, 11) = 2.5337109677677274E+00_dp
      fit_coef(5, 11) = -1.0698664293959414E+01_dp
      fit_coef(6, 11) = 3.1866930157362180E+01_dp
      fit_coef(7, 11) = -6.6687974631903131E+01_dp
      fit_coef(8, 11) = 9.7258475471561397E+01_dp
      fit_coef(9, 11) = -9.6610121349263864E+01_dp
      fit_coef(10, 11) = 6.2222012749124062E+01_dp
      fit_coef(11, 11) = -2.3400982843737399E+01_dp
      fit_coef(12, 11) = 3.8982563637438763E+00_dp
      fit_coef(0, 12) = 1.1710521653280541E-02_dp
      fit_coef(1, 12) = 3.8298132150440530E-04_dp
      fit_coef(2, 12) = 1.5720490880478075E-01_dp
      fit_coef(3, 12) = -1.0604827594360497E+00_dp
      fit_coef(4, 12) = 6.7293547068829866E+00_dp
      fit_coef(5, 12) = -2.8478199973116869E+01_dp
      fit_coef(6, 12) = 8.4349514492679333E+01_dp
      fit_coef(7, 12) = -1.7552622863849189E+02_dp
      fit_coef(8, 12) = 2.5489178159805036E+02_dp
      fit_coef(9, 12) = -2.5238005448515145E+02_dp
      fit_coef(10, 12) = 1.6214051051889587E+02_dp
      fit_coef(11, 12) = -6.0856023343944223E+01_dp
      fit_coef(12, 12) = 1.0120675732982292E+01_dp
      fit_coef(0, 13) = 3.8069679293470479E-02_dp
      fit_coef(1, 13) = 8.5030453244577039E-04_dp
      fit_coef(2, 13) = 3.5319383125255671E-01_dp
      fit_coef(3, 13) = -2.4058326640152856E+00_dp
      fit_coef(4, 13) = 1.5130832210336226E+01_dp
      fit_coef(5, 13) = -6.4118678617210037E+01_dp
      fit_coef(6, 13) = 1.8951627589384418E+02_dp
      fit_coef(7, 13) = -3.9344807722056498E+02_dp
      fit_coef(8, 13) = 5.7003275418870396E+02_dp
      fit_coef(9, 13) = -5.6318587947825517E+02_dp
      fit_coef(10, 13) = 3.6110912251064480E+02_dp
      fit_coef(11, 13) = -1.3530674224009485E+02_dp
      fit_coef(12, 13) = 2.2470490007659631E+01_dp
      fit_coef(0, 14) = 1.0718425215372540E-01_dp
      fit_coef(1, 14) = 1.6606498612547411E-03_dp
      fit_coef(2, 14) = 6.9558041141151850E-01_dp
      fit_coef(3, 14) = -4.7698051853919008E+00_dp
      fit_coef(4, 14) = 2.9821463877897131E+01_dp
      fit_coef(5, 14) = -1.2648112443458319E+02_dp
      fit_coef(6, 14) = 3.7339814874668650E+02_dp
      fit_coef(7, 14) = -7.7395557459000918E+02_dp
      fit_coef(8, 14) = 1.1195459460253383E+03_dp
      fit_coef(9, 14) = -1.1045564613814279E+03_dp
      fit_coef(10, 14) = 7.0739903710729982E+02_dp
      fit_coef(11, 14) = -2.6480250102979971E+02_dp
      fit_coef(12, 14) = 4.3940376198319392E+01_dp
      fit_coef(0, 15) = 2.7292854447821768E-01_dp
      fit_coef(1, 15) = 2.8818089651899471E-03_dp
      fit_coef(2, 15) = 1.2148125826442298E+00_dp
      fit_coef(3, 15) = -8.3719688971021089E+00_dp
      fit_coef(4, 15) = 5.2119215804340925E+01_dp
      fit_coef(5, 15) = -2.2119562023064080E+02_dp
      fit_coef(6, 15) = 6.5252980054868908E+02_dp
      fit_coef(7, 15) = -1.3509643415136143E+03_dp
      fit_coef(8, 15) = 1.9519497444270371E+03_dp
      fit_coef(9, 15) = -1.9238403235944149E+03_dp
      fit_coef(10, 15) = 1.2310185427991041E+03_dp
      fit_coef(11, 15) = -4.6046906916991361E+02_dp
      fit_coef(12, 15) = 7.6360956106162789E+01_dp
      fit_coef(0, 16) = 6.4599249663823077E-01_dp
      fit_coef(1, 16) = 4.4100315779845725E-03_dp
      fit_coef(2, 16) = 1.8680668872908257E+00_dp
      fit_coef(3, 16) = -1.2921923578119841E+01_dp
      fit_coef(4, 16) = 8.0194113618089176E+01_dp
      fit_coef(5, 16) = -3.4051136254025738E+02_dp
      fit_coef(6, 16) = 1.0040346645389629E+03_dp
      fit_coef(7, 16) = -2.0769701859400898E+03_dp
      fit_coef(8, 16) = 2.9983623465977071E+03_dp
      fit_coef(9, 16) = -2.9529215447188330E+03_dp
      fit_coef(10, 16) = 1.8882605454973300E+03_dp
      fit_coef(11, 16) = -7.0592135980760236E+02_dp
      fit_coef(12, 16) = 1.1701028544478041E+02_dp
      fit_coef(0, 17) = 1.4634101645424438E+00_dp
      fit_coef(1, 17) = 5.8372190193024661E-03_dp
      fit_coef(2, 17) = 2.4803142380725367E+00_dp
      fit_coef(3, 17) = -1.7197461327455997E+01_dp
      fit_coef(4, 17) = 1.0652223495438764E+02_dp
      fit_coef(5, 17) = -4.5244169342907389E+02_dp
      fit_coef(6, 17) = 1.3337313060706883E+03_dp
      fit_coef(7, 17) = -2.7575829556743033E+03_dp
      fit_coef(8, 17) = 3.9787771328072104E+03_dp
      fit_coef(9, 17) = -3.9165784027159407E+03_dp
      fit_coef(10, 17) = 2.5034305875991531E+03_dp
      fit_coef(11, 17) = -9.3557062065768866E+02_dp
      fit_coef(12, 17) = 1.5502979180126792E+02_dp
      fit_coef(0, 18) = 3.4614534817764171E+00_dp
      fit_coef(1, 18) = 6.7354426924417570E-03_dp
      fit_coef(2, 18) = 2.8636831995229377E+00_dp
      fit_coef(3, 18) = -1.9864414118527225E+01_dp
      fit_coef(4, 18) = 1.2299824240049257E+02_dp
      fit_coef(5, 18) = -5.2245197249800015E+02_dp
      fit_coef(6, 18) = 1.5400592407663903E+03_dp
      fit_coef(7, 18) = -3.1838942159350613E+03_dp
      fit_coef(8, 18) = 4.5934235919896164E+03_dp
      fit_coef(9, 18) = -4.5212045254281238E+03_dp
      fit_coef(10, 18) = 2.8896730980666093E+03_dp
      fit_coef(11, 18) = -1.0798434840197956E+03_dp
      fit_coef(12, 18) = 1.7892665893271254E+02_dp

   END SUBROUTINE get_fit_coef_k9_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k10_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 5.7695951023970724E-03_dp
      fit_coef(1, 1) = 1.5304265823073267E-02_dp
      fit_coef(2, 1) = 1.9872484454133494E-02_dp
      fit_coef(3, 1) = 1.4567559207979195E-02_dp
      fit_coef(4, 1) = 4.4496400713523688E-03_dp
      fit_coef(5, 1) = -2.0505040635942838E-03_dp
      fit_coef(6, 1) = -2.4017738586966987E-03_dp
      fit_coef(7, 1) = -2.1712350568040649E-04_dp
      fit_coef(8, 1) = 5.6984531809305442E-04_dp
      fit_coef(9, 1) = 3.6792198394363707E-04_dp
      fit_coef(10, 1) = -4.8212959819146314E-04_dp
      fit_coef(11, 1) = 2.3577845786711992E-04_dp
      fit_coef(12, 1) = -5.1828587653862264E-05_dp
      fit_coef(0, 2) = 3.1885037358175702E-02_dp
      fit_coef(1, 2) = 8.2163350099799248E-02_dp
      fit_coef(2, 2) = 1.0528892661444159E-01_dp
      fit_coef(3, 2) = 7.6258768872869365E-02_dp
      fit_coef(4, 2) = 2.2829853937072226E-02_dp
      fit_coef(5, 2) = -1.0880783303710033E-02_dp
      fit_coef(6, 2) = -1.2408166352033598E-02_dp
      fit_coef(7, 2) = -9.9928746104626666E-04_dp
      fit_coef(8, 2) = 2.9811486982521352E-03_dp
      fit_coef(9, 2) = 1.8293482458015533E-03_dp
      fit_coef(10, 2) = -2.4355182872439814E-03_dp
      fit_coef(11, 2) = 1.1814205683589125E-03_dp
      fit_coef(12, 2) = -2.5780560125925923E-04_dp
      fit_coef(0, 3) = 8.5985451133225102E-02_dp
      fit_coef(1, 3) = 2.0835675405554938E-01_dp
      fit_coef(2, 3) = 2.6074125375875318E-01_dp
      fit_coef(3, 3) = 1.8534648507356308E-01_dp
      fit_coef(4, 3) = 5.3562405566639951E-02_dp
      fit_coef(5, 3) = -2.7120042620104229E-02_dp
      fit_coef(6, 3) = -2.9327190422043962E-02_dp
      fit_coef(7, 3) = -1.7334511896438729E-03_dp
      fit_coef(8, 3) = 7.0440054029393588E-03_dp
      fit_coef(9, 3) = 4.0411619935981176E-03_dp
      fit_coef(10, 3) = -5.4238824013001187E-03_dp
      fit_coef(11, 3) = 2.5591744878845188E-03_dp
      fit_coef(12, 3) = -5.4685946652433383E-04_dp
      fit_coef(0, 4) = 1.8459848405035967E-01_dp
      fit_coef(1, 4) = 4.0371818510223983E-01_dp
      fit_coef(2, 4) = 4.8557108030340618E-01_dp
      fit_coef(3, 4) = 3.3806374668738176E-01_dp
      fit_coef(4, 4) = 9.5858897409524330E-02_dp
      fit_coef(5, 4) = -5.0962592261426340E-02_dp
      fit_coef(6, 4) = -5.3586938181428435E-02_dp
      fit_coef(7, 4) = -2.2127455334925599E-03_dp
      fit_coef(8, 4) = 1.6831832652833605E-02_dp
      fit_coef(9, 4) = 3.0255279721402920E-03_dp
      fit_coef(10, 4) = -8.9664127863715425E-03_dp
      fit_coef(11, 4) = 5.2862341907787443E-03_dp
      fit_coef(12, 4) = -1.2820544869880796E-03_dp
      fit_coef(0, 5) = 3.6074777273994402E-01_dp
      fit_coef(1, 5) = 6.8428202854170017E-01_dp
      fit_coef(2, 5) = 7.7269339647894841E-01_dp
      fit_coef(3, 5) = 5.2037619710849736E-01_dp
      fit_coef(4, 5) = 1.5184339437662456E-01_dp
      fit_coef(5, 5) = -6.9927558231047474E-02_dp
      fit_coef(6, 5) = -8.7838357618212293E-02_dp
      fit_coef(7, 5) = -2.6261703222456140E-03_dp
      fit_coef(8, 5) = 7.4641639246942270E-04_dp
      fit_coef(9, 5) = 6.0318746256968377E-02_dp
      fit_coef(10, 5) = -6.1544498098345557E-02_dp
      fit_coef(11, 5) = 2.6700632837816051E-02_dp
      fit_coef(12, 5) = -4.6151870101819167E-03_dp
      fit_coef(0, 6) = 6.7343634143757625E-01_dp
      fit_coef(1, 6) = 1.0727120778542409E+00_dp
      fit_coef(2, 6) = 1.1129907008803974E+00_dp
      fit_coef(3, 6) = 7.0510865936118905E-01_dp
      fit_coef(4, 6) = 2.1018587810883321E-01_dp
      fit_coef(5, 6) = -6.7282163797883557E-02_dp
      fit_coef(6, 6) = -1.0027943290324723E-01_dp
      fit_coef(7, 6) = -3.5384540450618285E-02_dp
      fit_coef(8, 6) = 2.7298363384797509E-02_dp
      fit_coef(9, 6) = -1.7363812432330001E-02_dp
      fit_coef(10, 6) = 7.1528999902988041E-02_dp
      fit_coef(11, 6) = -6.3212218372623016E-02_dp
      fit_coef(12, 6) = 1.7361608127688759E-02_dp
      fit_coef(0, 7) = 1.2230817482271914E+00_dp
      fit_coef(1, 7) = 1.5911587218644585E+00_dp
      fit_coef(2, 7) = 1.4985802851816239E+00_dp
      fit_coef(3, 7) = 8.6838937807453775E-01_dp
      fit_coef(4, 7) = 2.5143747366687763E-01_dp
      fit_coef(5, 7) = -4.6498700536608512E-02_dp
      fit_coef(6, 7) = -7.6008644969091751E-02_dp
      fit_coef(7, 7) = -6.7454078584185465E-02_dp
      fit_coef(8, 7) = 6.4705093608466005E-02_dp
      fit_coef(9, 7) = -1.2642157251641012E-01_dp
      fit_coef(10, 7) = 1.5991737393743985E-01_dp
      fit_coef(11, 7) = -7.4991121735080518E-02_dp
      fit_coef(12, 7) = 1.0660115126825225E-02_dp
      fit_coef(0, 8) = 2.1806860083417416E+00_dp
      fit_coef(1, 8) = 2.2497241492451168E+00_dp
      fit_coef(2, 8) = 1.9203525211709502E+00_dp
      fit_coef(3, 8) = 9.9969616591345500E-01_dp
      fit_coef(4, 8) = 2.6834037626429591E-01_dp
      fit_coef(5, 8) = -1.4414826654982551E-02_dp
      fit_coef(6, 8) = -4.5077253602132578E-02_dp
      fit_coef(7, 8) = -2.2179502828276729E-02_dp
      fit_coef(8, 8) = -4.7039179966872090E-02_dp
      fit_coef(9, 8) = 2.8411736628542274E-02_dp
      fit_coef(10, 8) = -4.1123712302139126E-02_dp
      fit_coef(11, 8) = 5.9491575141153316E-02_dp
      fit_coef(12, 8) = -2.2126040786265226E-02_dp
      fit_coef(0, 9) = 3.8549835579963414E+00_dp
      fit_coef(1, 9) = 3.0378827588214725E+00_dp
      fit_coef(2, 9) = 2.3710751843734088E+00_dp
      fit_coef(3, 9) = 1.1032804253025836E+00_dp
      fit_coef(4, 9) = 2.6679454698903354E-01_dp
      fit_coef(5, 9) = 2.0134331862044223E-02_dp
      fit_coef(6, 9) = -3.9522591450155079E-04_dp
      fit_coef(7, 9) = 7.4771871672822615E-03_dp
      fit_coef(8, 9) = -8.2968039659936169E-02_dp
      fit_coef(9, 9) = 5.7925277454710461E-02_dp
      fit_coef(10, 9) = -6.7746771231294078E-02_dp
      fit_coef(11, 9) = 5.6397993380196895E-02_dp
      fit_coef(12, 9) = -1.4614807655446010E-02_dp
      fit_coef(0, 10) = 6.9382542461265118E+00_dp
      fit_coef(1, 10) = 3.9362780471139174E+00_dp
      fit_coef(2, 10) = 2.8610007550879040E+00_dp
      fit_coef(3, 10) = 1.1997673683708834E+00_dp
      fit_coef(4, 10) = 2.5581466127865921E-01_dp
      fit_coef(5, 10) = 5.3143801132552647E-02_dp
      fit_coef(6, 10) = 5.1417759773000900E-02_dp
      fit_coef(7, 10) = 1.7678646974295765E-02_dp
      fit_coef(8, 10) = -5.9075047650367432E-02_dp
      fit_coef(9, 10) = 3.3432194957975334E-03_dp
      fit_coef(10, 10) = 3.6827535376516749E-03_dp
      fit_coef(11, 10) = 1.1242355647481125E-03_dp
      fit_coef(12, 10) = 9.1270977158394349E-04_dp
      fit_coef(0, 11) = 1.4966733342014741E-02_dp
      fit_coef(1, 11) = 3.9453275484996823E-02_dp
      fit_coef(2, 11) = 5.1068715597025638E-02_dp
      fit_coef(3, 11) = 3.7329091135728701E-02_dp
      fit_coef(4, 11) = 1.1348488247113217E-02_dp
      fit_coef(5, 11) = -5.2730834502581489E-03_dp
      fit_coef(6, 11) = -6.1356690848579278E-03_dp
      fit_coef(7, 11) = -5.3965409284480723E-04_dp
      fit_coef(8, 11) = 1.4598866172702007E-03_dp
      fit_coef(9, 11) = 9.3119411997548570E-04_dp
      fit_coef(10, 11) = -1.2246113300432057E-03_dp
      fit_coef(11, 11) = 5.9768922566995949E-04_dp
      fit_coef(12, 11) = -1.3116685473996465E-04_dp
      fit_coef(0, 12) = 3.8347028757142743E-02_dp
      fit_coef(1, 12) = 9.5250205763913104E-02_dp
      fit_coef(2, 12) = 1.2014345228064902E-01_dp
      fit_coef(3, 12) = 8.5770121052551718E-02_dp
      fit_coef(4, 12) = 2.5026676932284449E-02_dp
      fit_coef(5, 12) = -1.2425272424457500E-02_dp
      fit_coef(6, 12) = -1.3711933166334719E-02_dp
      fit_coef(7, 12) = -9.4665591048253780E-04_dp
      fit_coef(8, 12) = 3.3417995277827857E-03_dp
      fit_coef(9, 12) = 1.9446347870976362E-03_dp
      fit_coef(10, 12) = -2.6369332046548466E-03_dp
      fit_coef(11, 12) = 1.2654540755567044E-03_dp
      fit_coef(12, 12) = -2.7319897896357627E-04_dp
      fit_coef(0, 13) = 7.2491464313657286E-02_dp
      fit_coef(1, 13) = 1.5872975379092427E-01_dp
      fit_coef(2, 13) = 1.9076646279106174E-01_dp
      fit_coef(3, 13) = 1.3193323939107945E-01_dp
      fit_coef(4, 13) = 3.6136743102133911E-02_dp
      fit_coef(5, 13) = -2.0313520615877108E-02_dp
      fit_coef(6, 13) = -1.9992174365017833E-02_dp
      fit_coef(7, 13) = -1.5921233024662168E-04_dp
      fit_coef(8, 13) = 5.1319901147546939E-03_dp
      fit_coef(9, 13) = 1.3342159647477779E-03_dp
      fit_coef(10, 13) = -2.2525420201435871E-03_dp
      fit_coef(11, 13) = 9.7269009421198751E-04_dp
      fit_coef(12, 13) = -2.0108283608299456E-04_dp
      fit_coef(0, 14) = 1.3008928188788060E-01_dp
      fit_coef(1, 14) = 2.3463488136486438E-01_dp
      fit_coef(2, 14) = 2.5759089406615687E-01_dp
      fit_coef(3, 14) = 1.7130271768667638E-01_dp
      fit_coef(4, 14) = 4.9321273379861225E-02_dp
      fit_coef(5, 14) = -2.5370622060839100E-02_dp
      fit_coef(6, 14) = -2.9809532415184900E-02_dp
      fit_coef(7, 14) = -1.5417020754737089E-03_dp
      fit_coef(8, 14) = 1.0746594499561295E-02_dp
      fit_coef(9, 14) = 9.4386979195149809E-03_dp
      fit_coef(10, 14) = -1.9447726689701503E-02_dp
      fit_coef(11, 14) = 1.2247210595444594E-02_dp
      fit_coef(12, 14) = -2.9200822587582783E-03_dp
      fit_coef(0, 15) = 2.3180186914323608E-01_dp
      fit_coef(1, 15) = 3.3040402783496398E-01_dp
      fit_coef(2, 15) = 3.1506693954203840E-01_dp
      fit_coef(3, 15) = 1.8838521756146479E-01_dp
      fit_coef(4, 15) = 6.0573038457003757E-02_dp
      fit_coef(5, 15) = -8.8248929607627029E-03_dp
      fit_coef(6, 15) = -3.2842705061116521E-02_dp
      fit_coef(7, 15) = -3.8281459849134332E-03_dp
      fit_coef(8, 15) = -3.5375767464361246E-02_dp
      fit_coef(9, 15) = 6.9092174166339268E-02_dp
      fit_coef(10, 15) = -3.1110614611766781E-02_dp
      fit_coef(11, 15) = -4.5781884002109160E-03_dp
      fit_coef(12, 15) = 4.8446541181232632E-03_dp
      fit_coef(0, 16) = 4.0986431246574251E-01_dp
      fit_coef(1, 16) = 4.5035445709902655E-01_dp
      fit_coef(2, 16) = 3.6429713371972472E-01_dp
      fit_coef(3, 16) = 1.7696011433627806E-01_dp
      fit_coef(4, 16) = 5.2443940811353028E-02_dp
      fit_coef(5, 16) = 1.2731147310719145E-02_dp
      fit_coef(6, 16) = 1.2109883658619791E-02_dp
      fit_coef(7, 16) = -5.9511768374893713E-02_dp
      fit_coef(8, 16) = 9.3282598295030036E-02_dp
      fit_coef(9, 16) = -2.1484669589298355E-01_dp
      fit_coef(10, 16) = 2.5063354943196953E-01_dp
      fit_coef(11, 16) = -1.3015721063956481E-01_dp
      fit_coef(12, 16) = 2.5163533574069471E-02_dp
      fit_coef(0, 17) = 7.1698706857622740E-01_dp
      fit_coef(1, 17) = 5.8839597496590568E-01_dp
      fit_coef(2, 17) = 4.0531879096886647E-01_dp
      fit_coef(3, 17) = 1.4784602258417653E-01_dp
      fit_coef(4, 17) = 2.8867088347120023E-02_dp
      fit_coef(5, 17) = 2.7968490127980788E-02_dp
      fit_coef(6, 17) = 2.8408348764650926E-02_dp
      fit_coef(7, 17) = 1.6630805867207835E-02_dp
      fit_coef(8, 17) = -6.5625993058929091E-02_dp
      fit_coef(9, 17) = 7.6969150282882134E-02_dp
      fit_coef(10, 17) = -1.3660437598545347E-01_dp
      fit_coef(11, 17) = 1.2266105788436464E-01_dp
      fit_coef(12, 17) = -3.7327615955604555E-02_dp
      fit_coef(0, 18) = 1.2469010615236285E+00_dp
      fit_coef(1, 18) = 7.2690110706329647E-01_dp
      fit_coef(2, 18) = 4.3679155042183249E-01_dp
      fit_coef(3, 18) = 1.1545454220850887E-01_dp
      fit_coef(4, 18) = 6.1457460689058442E-03_dp
      fit_coef(5, 18) = 3.4446450178461388E-02_dp
      fit_coef(6, 18) = 3.6803538637622653E-02_dp
      fit_coef(7, 18) = 4.9357342985472417E-02_dp
      fit_coef(8, 18) = -9.8393813847182057E-02_dp
      fit_coef(9, 18) = 1.3063016463318536E-01_dp
      fit_coef(10, 18) = -1.5004086198844574E-01_dp
      fit_coef(11, 18) = 7.8940423149777589E-02_dp
      fit_coef(12, 18) = -1.3100508732995670E-02_dp
      fit_coef(0, 19) = 2.2020857144301220E+00_dp
      fit_coef(1, 19) = 8.4496482669483042E-01_dp
      fit_coef(2, 19) = 4.6567133382072923E-01_dp
      fit_coef(3, 19) = 9.5046757527473213E-02_dp
      fit_coef(4, 19) = -7.6836222188503209E-03_dp
      fit_coef(5, 19) = 3.4125967755404704E-02_dp
      fit_coef(6, 19) = 5.0617003282733740E-02_dp
      fit_coef(7, 19) = 1.2669228101081176E-02_dp
      fit_coef(8, 19) = 1.4598367047566483E-02_dp
      fit_coef(9, 19) = -4.8146800081324943E-02_dp
      fit_coef(10, 19) = 6.4386938465572141E-02_dp
      fit_coef(11, 19) = -5.7857218864180007E-02_dp
      fit_coef(12, 19) = 1.9235315360114507E-02_dp
      fit_coef(0, 20) = 4.3050004435836158E+00_dp
      fit_coef(1, 20) = 9.5976454358976238E-01_dp
      fit_coef(2, 20) = 5.2630087854492191E-01_dp
      fit_coef(3, 20) = 1.0553767791618235E-01_dp
      fit_coef(4, 20) = -1.2890279709804050E-02_dp
      fit_coef(5, 20) = 3.1631245951617148E-02_dp
      fit_coef(6, 20) = 5.1846721167766188E-02_dp
      fit_coef(7, 20) = 1.3931333624011479E-02_dp
      fit_coef(8, 20) = 1.8289419861101722E-02_dp
      fit_coef(9, 20) = -3.6941691804435302E-02_dp
      fit_coef(10, 20) = 5.2037604092826852E-02_dp
      fit_coef(11, 20) = -3.6753728963893388E-02_dp
      fit_coef(12, 20) = 8.3623158625230273E-03_dp

   END SUBROUTINE get_fit_coef_k10_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k10_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 5.6502000000000000E+04_dp
      fit_coef(0, 1) = 6.8411577604394062E-05_dp
      fit_coef(1, 1) = 1.6187597727541771E-05_dp
      fit_coef(2, 1) = 4.5871906583389663E-03_dp
      fit_coef(3, 1) = -3.0179373289852002E-02_dp
      fit_coef(4, 1) = 2.0974369374411136E-01_dp
      fit_coef(5, 1) = -8.7935873765434858E-01_dp
      fit_coef(6, 1) = 2.5878525249126527E+00_dp
      fit_coef(7, 1) = -5.3179060500186415E+00_dp
      fit_coef(8, 1) = 7.5892153487240241E+00_dp
      fit_coef(9, 1) = -7.3662068168070789E+00_dp
      fit_coef(10, 1) = 4.6344768125153601E+00_dp
      fit_coef(11, 1) = -1.7032106343424485E+00_dp
      fit_coef(12, 1) = 2.7741991549383826E-01_dp
      fit_coef(0, 2) = 7.5275323624865811E-04_dp
      fit_coef(1, 2) = 9.9232481710528875E-05_dp
      fit_coef(2, 2) = 2.8524678297215784E-02_dp
      fit_coef(3, 2) = -1.9184050122553747E-01_dp
      fit_coef(4, 2) = 1.3126341229840768E+00_dp
      fit_coef(5, 2) = -5.5140690139707527E+00_dp
      fit_coef(6, 2) = 1.6223500597815761E+01_dp
      fit_coef(7, 2) = -3.3331549263061540E+01_dp
      fit_coef(8, 2) = 4.7557073952296840E+01_dp
      fit_coef(9, 2) = -4.6143966985722287E+01_dp
      fit_coef(10, 2) = 2.9018553303790554E+01_dp
      fit_coef(11, 2) = -1.0658983662222605E+01_dp
      fit_coef(12, 2) = 1.7351738260550609E+00_dp
      fit_coef(0, 3) = 3.9843606366604320E-03_dp
      fit_coef(1, 3) = 3.3764501333670899E-04_dp
      fit_coef(2, 3) = 9.9383477377645593E-02_dp
      fit_coef(3, 3) = -6.8774784004509004E-01_dp
      fit_coef(4, 3) = 4.6197752367269072E+00_dp
      fit_coef(5, 3) = -1.9617538089460556E+01_dp
      fit_coef(6, 3) = 5.7821554084121672E+01_dp
      fit_coef(7, 3) = -1.1866755198936438E+02_dp
      fit_coef(8, 3) = 1.6899851790071165E+02_dp
      fit_coef(9, 3) = -1.6365529688728941E+02_dp
      fit_coef(10, 3) = 1.0272818227085014E+02_dp
      fit_coef(11, 3) = -3.7670599569470930E+01_dp
      fit_coef(12, 3) = 6.1231647953439321E+00_dp
      fit_coef(0, 4) = 1.5378969066141557E-02_dp
      fit_coef(1, 4) = 9.1548502345901116E-04_dp
      fit_coef(2, 4) = 2.7372675256892554E-01_dp
      fit_coef(3, 4) = -1.9285608309876876E+00_dp
      fit_coef(4, 4) = 1.2803726769725465E+01_dp
      fit_coef(5, 4) = -5.4755620729323439E+01_dp
      fit_coef(6, 4) = 1.6149496474022160E+02_dp
      fit_coef(7, 4) = -3.3133212972051763E+02_dp
      fit_coef(8, 4) = 4.7168386050627021E+02_dp
      fit_coef(9, 4) = -4.5657552176306115E+02_dp
      fit_coef(10, 4) = 2.8644979600872671E+02_dp
      fit_coef(11, 4) = -1.0497888439910341E+02_dp
      fit_coef(12, 4) = 1.7052642506098298E+01_dp
      fit_coef(0, 5) = 4.9251180618047154E-02_dp
      fit_coef(1, 5) = 2.1434954065410114E-03_dp
      fit_coef(2, 5) = 6.4789716785013773E-01_dp
      fit_coef(3, 5) = -4.6201666914403452E+00_dp
      fit_coef(4, 5) = 3.0440812558042978E+01_dp
      fit_coef(5, 5) = -1.3079799887515426E+02_dp
      fit_coef(6, 5) = 3.8607450191762678E+02_dp
      fit_coef(7, 5) = -7.9209798040256953E+02_dp
      fit_coef(8, 5) = 1.1273243950495682E+03_dp
      fit_coef(9, 5) = -1.0907316178904814E+03_dp
      fit_coef(10, 5) = 6.8394210358544069E+02_dp
      fit_coef(11, 5) = -2.5050774766365097E+02_dp
      fit_coef(12, 5) = 4.0668458941697637E+01_dp
      fit_coef(0, 6) = 1.3896177117430439E-01_dp
      fit_coef(1, 6) = 4.4858622820344507E-03_dp
      fit_coef(2, 6) = 1.3666220193198582E+00_dp
      fit_coef(3, 6) = -9.8301126837720361E+00_dp
      fit_coef(4, 6) = 6.4421339807998407E+01_dp
      fit_coef(5, 6) = -2.7772758192096933E+02_dp
      fit_coef(6, 6) = 8.2039176215178543E+02_dp
      fit_coef(7, 6) = -1.6832836477463704E+03_dp
      fit_coef(8, 6) = 2.3951662567438020E+03_dp
      fit_coef(9, 6) = -2.3166358021581445E+03_dp
      fit_coef(10, 6) = 1.4520770947679368E+03_dp
      fit_coef(11, 6) = -5.3163751269775196E+02_dp
      fit_coef(12, 6) = 8.6273624073151339E+01_dp
      fit_coef(0, 7) = 3.5729488480195148E-01_dp
      fit_coef(1, 7) = 8.5259570236280655E-03_dp
      fit_coef(2, 7) = 2.6129406592915587E+00_dp
      fit_coef(3, 7) = -1.8917534003546230E+01_dp
      fit_coef(4, 7) = 1.2348529805386671E+02_dp
      fit_coef(5, 7) = -5.3366943392175176E+02_dp
      fit_coef(6, 7) = 1.5775492634594218E+03_dp
      fit_coef(7, 7) = -3.2371204369651055E+03_dp
      fit_coef(8, 7) = 4.6054333227458374E+03_dp
      fit_coef(9, 7) = -4.4533030686941875E+03_dp
      fit_coef(10, 7) = 2.7905356324814979E+03_dp
      fit_coef(11, 7) = -1.0213693901858213E+03_dp
      fit_coef(12, 7) = 1.6569764375256605E+02_dp
      fit_coef(0, 8) = 8.5576326888027310E-01_dp
      fit_coef(1, 8) = 1.4794510243735135E-02_dp
      fit_coef(2, 8) = 4.5547377508589388E+00_dp
      fit_coef(3, 8) = -3.3139938275450689E+01_dp
      fit_coef(4, 8) = 2.1567954334326745E+02_dp
      fit_coef(5, 8) = -9.3383264595716162E+02_dp
      fit_coef(6, 8) = 2.7621813661246701E+03_dp
      fit_coef(7, 8) = -5.6685491647807048E+03_dp
      fit_coef(8, 8) = 8.0637217679143905E+03_dp
      fit_coef(9, 8) = -7.7958512685636942E+03_dp
      fit_coef(10, 8) = 4.8839651073288114E+03_dp
      fit_coef(11, 8) = -1.7871820510264045E+03_dp
      fit_coef(12, 8) = 2.8987122244652051E+02_dp
      fit_coef(0, 9) = 1.9457453939982912E+00_dp
      fit_coef(1, 9) = 2.3412967111721320E-02_dp
      fit_coef(2, 9) = 7.2322818888009728E+00_dp
      fit_coef(3, 9) = -5.2813697786989756E+01_dp
      fit_coef(4, 9) = 3.4297595366907478E+02_dp
      fit_coef(5, 9) = -1.4869895685428958E+03_dp
      fit_coef(6, 9) = 4.4006178698965341E+03_dp
      fit_coef(7, 9) = -9.0318062878397541E+03_dp
      fit_coef(8, 9) = 1.2847105139427764E+04_dp
      fit_coef(9, 9) = -1.2418582634546168E+04_dp
      fit_coef(10, 9) = 7.7787564640697719E+03_dp
      fit_coef(11, 9) = -2.8459924565854521E+03_dp
      fit_coef(12, 9) = 4.6152903795316814E+02_dp
      fit_coef(0, 10) = 4.3456834390820589E+00_dp
      fit_coef(1, 10) = 3.3830817174057760E-02_dp
      fit_coef(2, 10) = 1.0473175894252773E+01_dp
      fit_coef(3, 10) = -7.6661420311861534E+01_dp
      fit_coef(4, 10) = 4.9715082005031553E+02_dp
      fit_coef(5, 10) = -2.1572931305206675E+03_dp
      fit_coef(6, 10) = 6.3865758669205397E+03_dp
      fit_coef(7, 10) = -1.3108704020006919E+04_dp
      fit_coef(8, 10) = 1.8645335870199138E+04_dp
      fit_coef(9, 10) = -1.8021767528763816E+04_dp
      fit_coef(10, 10) = 1.1287284370105295E+04_dp
      fit_coef(11, 10) = -4.1291992041477679E+03_dp
      fit_coef(12, 10) = 6.6955320372977019E+02_dp
      fit_coef(0, 11) = 2.1871793237239991E-04_dp
      fit_coef(1, 11) = 4.1996695354716680E-05_dp
      fit_coef(2, 11) = 1.1909044603489627E-02_dp
      fit_coef(3, 11) = -7.8370180707781900E-02_dp
      fit_coef(4, 11) = 5.4437980679476261E-01_dp
      fit_coef(5, 11) = -2.2731288629006801E+00_dp
      fit_coef(6, 11) = 6.6595389300613057E+00_dp
      fit_coef(7, 11) = -1.3637699933400350E+01_dp
      fit_coef(8, 11) = 1.9415873433968578E+01_dp
      fit_coef(9, 11) = -1.8815612032886889E+01_dp
      fit_coef(10, 11) = 1.1826000388483891E+01_dp
      fit_coef(11, 11) = -4.3434643915415521E+00_dp
      fit_coef(12, 11) = 7.0720998747701980E-01_dp
      fit_coef(0, 12) = 1.4372890424941905E-03_dp
      fit_coef(1, 12) = 1.3870980955592494E-04_dp
      fit_coef(2, 12) = 4.0623488722695424E-02_dp
      fit_coef(3, 12) = -2.7976182205365724E-01_dp
      fit_coef(4, 12) = 1.8859471674798522E+00_dp
      fit_coef(5, 12) = -7.9899664094381135E+00_dp
      fit_coef(6, 12) = 2.3590888671742746E+01_dp
      fit_coef(7, 12) = -4.8520694758407572E+01_dp
      fit_coef(8, 12) = 6.9210085374286393E+01_dp
      fit_coef(9, 12) = -6.7087759398956806E+01_dp
      fit_coef(10, 12) = 4.2134889802249333E+01_dp
      fit_coef(11, 12) = -1.5455367957037941E+01_dp
      fit_coef(12, 12) = 2.5125419036554626E+00_dp
      fit_coef(0, 13) = 5.8841380648918286E-03_dp
      fit_coef(1, 13) = 3.6831095166302676E-04_dp
      fit_coef(2, 13) = 1.1027978007247614E-01_dp
      fit_coef(3, 13) = -7.7842683175317762E-01_dp
      fit_coef(4, 13) = 5.1600304775061447E+00_dp
      fit_coef(5, 13) = -2.2078165304729261E+01_dp
      fit_coef(6, 13) = 6.5091858701369759E+01_dp
      fit_coef(7, 13) = -1.3345586681486435E+02_dp
      fit_coef(8, 13) = 1.8988777441894680E+02_dp
      fit_coef(9, 13) = -1.8375507399027742E+02_dp
      fit_coef(10, 13) = 1.1527676579024002E+02_dp
      fit_coef(11, 13) = -4.2248714108699794E+01_dp
      fit_coef(12, 13) = 6.8635207622576706E+00_dp
      fit_coef(0, 14) = 1.9128947829258891E-02_dp
      fit_coef(1, 14) = 8.3862411750444427E-04_dp
      fit_coef(2, 14) = 2.5428661374901718E-01_dp
      fit_coef(3, 14) = -1.8194680686167399E+00_dp
      fit_coef(4, 14) = 1.1961590588350870E+01_dp
      fit_coef(5, 14) = -5.1467768848618732E+01_dp
      fit_coef(6, 14) = 1.5192487004386140E+02_dp
      fit_coef(7, 14) = -3.1168857906892379E+02_dp
      fit_coef(8, 14) = 4.4359497200294106E+02_dp
      fit_coef(9, 14) = -4.2917874847834634E+02_dp
      fit_coef(10, 14) = 2.6909129512594421E+02_dp
      fit_coef(11, 14) = -9.8546840419864964E+01_dp
      fit_coef(12, 14) = 1.5995921552846303E+01_dp
      fit_coef(0, 15) = 5.3855099007751053E-02_dp
      fit_coef(1, 15) = 1.6951202894427357E-03_dp
      fit_coef(2, 15) = 5.1828197057808256E-01_dp
      fit_coef(3, 15) = -3.7423907453765741E+00_dp
      fit_coef(4, 15) = 2.4465983174552605E+01_dp
      fit_coef(5, 15) = -1.0563582904908276E+02_dp
      fit_coef(6, 15) = 3.1210473924282360E+02_dp
      fit_coef(7, 15) = -6.4036928643586384E+02_dp
      fit_coef(8, 15) = 9.1110483535702542E+02_dp
      fit_coef(9, 15) = -8.8110226224210760E+02_dp
      fit_coef(10, 15) = 5.5218119224018108E+02_dp
      fit_coef(11, 15) = -2.0212810880984674E+02_dp
      fit_coef(12, 15) = 3.2795085750230733E+01_dp
      fit_coef(0, 16) = 1.3708229195964663E-01_dp
      fit_coef(1, 16) = 3.0894965681102504E-03_dp
      fit_coef(2, 16) = 9.5049659850978763E-01_dp
      fit_coef(3, 16) = -6.9100459255840390E+00_dp
      fit_coef(4, 16) = 4.4990748234228938E+01_dp
      fit_coef(5, 16) = -1.9474683678177831E+02_dp
      fit_coef(6, 16) = 5.7588070574709627E+02_dp
      fit_coef(7, 16) = -1.1817334181055951E+03_dp
      fit_coef(8, 16) = 1.6810743266760862E+03_dp
      fit_coef(9, 16) = -1.6252829227456073E+03_dp
      fit_coef(10, 16) = 1.0182493640037462E+03_dp
      fit_coef(11, 16) = -3.7261997804392706E+02_dp
      fit_coef(12, 16) = 6.0438944235281710E+01_dp
      fit_coef(0, 17) = 3.2363428278830764E-01_dp
      fit_coef(1, 17) = 5.0854788099187554E-03_dp
      fit_coef(2, 17) = 1.5721133863814079E+00_dp
      fit_coef(3, 17) = -1.1489215181637956E+01_dp
      fit_coef(4, 17) = 7.4573534244656926E+01_dp
      fit_coef(5, 17) = -3.2342087496454349E+02_dp
      fit_coef(6, 17) = 9.5710734682470525E+02_dp
      fit_coef(7, 17) = -1.9642973854914649E+03_dp
      fit_coef(8, 17) = 2.7939950801758901E+03_dp
      fit_coef(9, 17) = -2.7007178237082067E+03_dp
      fit_coef(10, 17) = 1.6916202883905507E+03_dp
      fit_coef(11, 17) = -6.1888646191941052E+02_dp
      fit_coef(12, 17) = 1.0035990724886408E+02_dp
      fit_coef(0, 18) = 7.2227658357344648E-01_dp
      fit_coef(1, 18) = 7.4788858723735070E-03_dp
      fit_coef(2, 18) = 2.3204406383068310E+00_dp
      fit_coef(3, 18) = -1.7025391312911438E+01_dp
      fit_coef(4, 18) = 1.1025201945927834E+02_dp
      fit_coef(5, 18) = -4.7884421266123894E+02_dp
      fit_coef(6, 18) = 1.4179656218358291E+03_dp
      fit_coef(7, 18) = -2.9105197632238601E+03_dp
      fit_coef(8, 18) = 4.1395786832377917E+03_dp
      fit_coef(9, 18) = -4.0007691450093103E+03_dp
      fit_coef(10, 18) = 2.5054759345520029E+03_dp
      fit_coef(11, 18) = -9.1647339178452830E+02_dp
      fit_coef(12, 18) = 1.4859116561702967E+02_dp
      fit_coef(0, 19) = 1.5633146329155199E+00_dp
      fit_coef(1, 19) = 9.6549690016911830E-03_dp
      fit_coef(2, 19) = 3.0024364887949786E+00_dp
      fit_coef(3, 19) = -2.2083956048519831E+01_dp
      fit_coef(4, 19) = 1.4280565323321272E+02_dp
      fit_coef(5, 19) = -6.2078293004938610E+02_dp
      fit_coef(6, 19) = 1.8390901346312414E+03_dp
      fit_coef(7, 19) = -3.7753165401987344E+03_dp
      fit_coef(8, 19) = 5.3693456762596288E+03_dp
      fit_coef(9, 19) = -5.1888078370584299E+03_dp
      fit_coef(10, 19) = 3.2491199407176041E+03_dp
      fit_coef(11, 19) = -1.1883539144960864E+03_dp
      fit_coef(12, 19) = 1.9265088532272543E+02_dp
      fit_coef(0, 20) = 3.5765140128766371E+00_dp
      fit_coef(1, 20) = 1.1092017818854312E-02_dp
      fit_coef(2, 20) = 3.4506220465727910E+00_dp
      fit_coef(3, 20) = -2.5391055124716122E+01_dp
      fit_coef(4, 20) = 1.6415262153898755E+02_dp
      fit_coef(5, 20) = -7.1368338873477956E+02_dp
      fit_coef(6, 20) = 2.1144976599358906E+03_dp
      fit_coef(7, 20) = -4.3407867181744632E+03_dp
      fit_coef(8, 20) = 6.1735437097235817E+03_dp
      fit_coef(9, 20) = -5.9658726159362077E+03_dp
      fit_coef(10, 20) = 3.7356297419849939E+03_dp
      fit_coef(11, 20) = -1.3662663971818399E+03_dp
      fit_coef(12, 20) = 2.2148914181168524E+02_dp

   END SUBROUTINE get_fit_coef_k10_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k11_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 5.2171158533482745E-03_dp
      fit_coef(1, 1) = 1.3928278169011961E-02_dp
      fit_coef(2, 1) = 1.8150859480949404E-02_dp
      fit_coef(3, 1) = 1.3341207416860837E-02_dp
      fit_coef(4, 1) = 4.0847907407787400E-03_dp
      fit_coef(5, 1) = -1.8862057580171247E-03_dp
      fit_coef(6, 1) = -2.2005927300723840E-03_dp
      fit_coef(7, 1) = -2.1569598773101882E-04_dp
      fit_coef(8, 1) = 5.3374114714500355E-04_dp
      fit_coef(9, 1) = 3.6786248990713488E-04_dp
      fit_coef(10, 1) = -4.9716456588384794E-04_dp
      fit_coef(11, 1) = 2.5364063357271070E-04_dp
      fit_coef(12, 1) = -5.7560771366250216E-05_dp
      fit_coef(0, 2) = 2.8561321333081039E-02_dp
      fit_coef(1, 2) = 7.4577106385596906E-02_dp
      fit_coef(2, 2) = 9.6095992956011567E-02_dp
      fit_coef(3, 2) = 6.9921004219163585E-02_dp
      fit_coef(4, 2) = 2.1046006027983205E-02_dp
      fit_coef(5, 2) = -1.0007180795562516E-02_dp
      fit_coef(6, 2) = -1.1398359734540125E-02_dp
      fit_coef(7, 2) = -1.0344153575843251E-03_dp
      fit_coef(8, 2) = 2.8187351548754256E-03_dp
      fit_coef(9, 2) = 1.8266346026223723E-03_dp
      fit_coef(10, 2) = -2.5228015867621856E-03_dp
      fit_coef(11, 2) = 1.2855087989079601E-03_dp
      fit_coef(12, 2) = -2.9136730780610450E-04_dp
      fit_coef(0, 3) = 7.5565672474912859E-02_dp
      fit_coef(1, 3) = 1.8837838893012882E-01_dp
      fit_coef(2, 3) = 2.3788031159177950E-01_dp
      fit_coef(3, 3) = 1.7019649338979928E-01_dp
      fit_coef(4, 3) = 4.9702597603960480E-02_dp
      fit_coef(5, 3) = -2.4843028083216520E-02_dp
      fit_coef(6, 3) = -2.7107463499384012E-02_dp
      fit_coef(7, 3) = -2.0998400355075556E-03_dp
      fit_coef(8, 3) = 6.8537414053824698E-03_dp
      fit_coef(9, 3) = 4.1037972870635875E-03_dp
      fit_coef(10, 3) = -5.8343544319931773E-03_dp
      fit_coef(11, 3) = 2.9563811258876472E-03_dp
      fit_coef(12, 3) = -6.6659895709530450E-04_dp
      fit_coef(0, 4) = 1.5746184306040770E-01_dp
      fit_coef(1, 4) = 3.6307100829710298E-01_dp
      fit_coef(2, 4) = 4.4411067220455136E-01_dp
      fit_coef(3, 4) = 3.1156735023866838E-01_dp
      fit_coef(4, 4) = 8.8318376917750713E-02_dp
      fit_coef(5, 4) = -4.6905172790974697E-02_dp
      fit_coef(6, 4) = -4.8633790525923246E-02_dp
      fit_coef(7, 4) = -2.6374461756591546E-03_dp
      fit_coef(8, 4) = 1.4152819458865665E-02_dp
      fit_coef(9, 4) = 4.0713403471031064E-03_dp
      fit_coef(10, 4) = -8.1648840133065614E-03_dp
      fit_coef(11, 4) = 4.4172781177938061E-03_dp
      fit_coef(12, 4) = -1.0468849549327278E-03_dp
      fit_coef(0, 5) = 2.9628095459795201E-01_dp
      fit_coef(1, 5) = 6.1067548081671130E-01_dp
      fit_coef(2, 5) = 7.1073485839444905E-01_dp
      fit_coef(3, 5) = 4.8623515872109696E-01_dp
      fit_coef(4, 5) = 1.3873344982657515E-01_dp
      fit_coef(5, 5) = -7.1183381958692823E-02_dp
      fit_coef(6, 5) = -7.9018683634448778E-02_dp
      fit_coef(7, 5) = -4.1849522345858052E-03_dp
      fit_coef(8, 5) = 1.8669075914354103E-02_dp
      fit_coef(9, 5) = 2.5640549363637480E-02_dp
      fit_coef(10, 5) = -3.7238167104347006E-02_dp
      fit_coef(11, 5) = 2.0235521804447935E-02_dp
      fit_coef(12, 5) = -4.4633644085237409E-03_dp
      fit_coef(0, 6) = 5.3079425354928900E-01_dp
      fit_coef(1, 6) = 9.4909599974912540E-01_dp
      fit_coef(2, 6) = 1.0298933488572886E+00_dp
      fit_coef(3, 6) = 6.7464315041251843E-01_dp
      fit_coef(4, 6) = 1.9785397712598238E-01_dp
      fit_coef(5, 6) = -8.0277433825610117E-02_dp
      fit_coef(6, 6) = -1.0862148180963889E-01_dp
      fit_coef(7, 6) = -1.2190533137530357E-02_dp
      fit_coef(8, 6) = -4.3235690156727307E-03_dp
      fit_coef(9, 6) = 7.1406154351795634E-02_dp
      fit_coef(10, 6) = -4.7520764363596951E-02_dp
      fit_coef(11, 6) = 6.9200248215417976E-03_dp
      fit_coef(12, 6) = 1.7348958867823125E-03_dp
      fit_coef(0, 7) = 9.2513564519354430E-01_dp
      fit_coef(1, 7) = 1.3995257053484436E+00_dp
      fit_coef(2, 7) = 1.3941989426918036E+00_dp
      fit_coef(3, 7) = 8.5298558301770655E-01_dp
      fit_coef(4, 7) = 2.5029080019290989E-01_dp
      fit_coef(5, 7) = -6.9418886584503764E-02_dp
      fit_coef(6, 7) = -1.0048676489518746E-01_dp
      fit_coef(7, 7) = -6.5140068180983265E-02_dp
      fit_coef(8, 7) = 7.1961159826162302E-02_dp
      fit_coef(9, 7) = -1.0845933619226010E-01_dp
      fit_coef(10, 7) = 1.7018602738764302E-01_dp
      fit_coef(11, 7) = -1.0954176694571671E-01_dp
      fit_coef(12, 7) = 2.4968149530473070E-02_dp
      fit_coef(0, 8) = 1.5835046294778308E+00_dp
      fit_coef(1, 8) = 1.9791294738794438E+00_dp
      fit_coef(2, 8) = 1.7961980138940286E+00_dp
      fit_coef(3, 8) = 1.0051925409840725E+00_dp
      fit_coef(4, 8) = 2.8272786197343069E-01_dp
      fit_coef(5, 8) = -4.3458147169371041E-02_dp
      fit_coef(6, 8) = -7.1388547183826534E-02_dp
      fit_coef(7, 8) = -6.9105034008714997E-02_dp
      fit_coef(8, 8) = 4.9061162554567488E-02_dp
      fit_coef(9, 8) = -1.0568300110809625E-01_dp
      fit_coef(10, 8) = 1.1802618075924728E-01_dp
      fit_coef(11, 8) = -3.4740368047083123E-02_dp
      fit_coef(12, 8) = -2.0070688342111883E-03_dp
      fit_coef(0, 9) = 2.6780451266380028E+00_dp
      fit_coef(1, 9) = 2.6926034999618897E+00_dp
      fit_coef(2, 9) = 2.2273953918391149E+00_dp
      fit_coef(3, 9) = 1.1263680392434321E+00_dp
      fit_coef(4, 9) = 2.9395288063392006E-01_dp
      fit_coef(5, 9) = -8.9510678035254702E-03_dp
      fit_coef(6, 9) = -3.6852206543507286E-02_dp
      fit_coef(7, 9) = -2.0316323307418317E-02_dp
      fit_coef(8, 9) = -6.2348519142117052E-02_dp
      fit_coef(9, 9) = 5.0424313584454393E-02_dp
      fit_coef(10, 9) = -7.4016849780124558E-02_dp
      fit_coef(11, 9) = 8.1891720384978092E-02_dp
      fit_coef(12, 9) = -2.7202540553906623E-02_dp
      fit_coef(0, 10) = 4.5149548116067084E+00_dp
      fit_coef(1, 10) = 3.5279902941181556E+00_dp
      fit_coef(2, 10) = 2.6834197969535900E+00_dp
      fit_coef(3, 10) = 1.2234297318240759E+00_dp
      fit_coef(4, 10) = 2.9061691203379159E-01_dp
      fit_coef(5, 10) = 2.6998851164527955E-02_dp
      fit_coef(6, 10) = 1.1503151087891473E-02_dp
      fit_coef(7, 10) = -2.9263858448428439E-04_dp
      fit_coef(8, 10) = -7.4172023234369305E-02_dp
      fit_coef(9, 10) = 4.4507509039308939E-02_dp
      fit_coef(10, 10) = -6.0661296025507475E-02_dp
      fit_coef(11, 10) = 5.4037676664233973E-02_dp
      fit_coef(12, 10) = -1.3914979544265946E-02_dp
      fit_coef(0, 11) = 7.7831432490395818E+00_dp
      fit_coef(1, 11) = 4.4728899281150136E+00_dp
      fit_coef(2, 11) = 3.1793306487074551E+00_dp
      fit_coef(3, 11) = 1.3166282129264584E+00_dp
      fit_coef(4, 11) = 2.7989649807784495E-01_dp
      fit_coef(5, 11) = 6.1709167955964331E-02_dp
      fit_coef(6, 11) = 6.5030571914413512E-02_dp
      fit_coef(7, 11) = 5.8652910860924629E-03_dp
      fit_coef(8, 11) = -4.3479817133538171E-02_dp
      fit_coef(9, 11) = -1.5946703925855100E-02_dp
      fit_coef(10, 11) = 1.4840762855674984E-02_dp
      fit_coef(11, 11) = -1.6764580405545060E-03_dp
      fit_coef(12, 11) = 1.0361409280594497E-03_dp
      fit_coef(0, 12) = 1.3505623056792662E-02_dp
      fit_coef(1, 12) = 3.5881778901904388E-02_dp
      fit_coef(2, 12) = 4.6635585408641191E-02_dp
      fit_coef(3, 12) = 3.4196109564738703E-02_dp
      fit_coef(4, 12) = 1.0427820008519861E-02_dp
      fit_coef(5, 12) = -4.8496513706743740E-03_dp
      fit_coef(6, 12) = -5.6249874188544387E-03_dp
      fit_coef(7, 12) = -5.4152291061643394E-04_dp
      fit_coef(8, 12) = 1.3707635274310744E-03_dp
      fit_coef(9, 12) = 9.3081024123118999E-04_dp
      fit_coef(10, 12) = -1.2643999344426028E-03_dp
      fit_coef(11, 12) = 6.4488398114361249E-04_dp
      fit_coef(12, 12) = -1.4631445736861331E-04_dp
      fit_coef(0, 13) = 3.3953336933937019E-02_dp
      fit_coef(1, 13) = 8.6200605204820452E-02_dp
      fit_coef(2, 13) = 1.0956860661445041E-01_dp
      fit_coef(3, 13) = 7.8748352145278117E-02_dp
      fit_coef(4, 13) = 2.3201917648793789E-02_dp
      fit_coef(5, 13) = -1.1428692068356564E-02_dp
      fit_coef(6, 13) = -1.2646975774724427E-02_dp
      fit_coef(7, 13) = -1.0364501016809282E-03_dp
      fit_coef(8, 13) = 3.2018633266795151E-03_dp
      fit_coef(9, 13) = 1.9223995549073674E-03_dp
      fit_coef(10, 13) = -2.7320063801979958E-03_dp
      fit_coef(11, 13) = 1.3899117897301057E-03_dp
      fit_coef(12, 13) = -3.1441281031980347E-04_dp
      fit_coef(0, 14) = 6.1851572757860851E-02_dp
      fit_coef(1, 14) = 1.4267393289529001E-01_dp
      fit_coef(2, 14) = 1.7418539035135686E-01_dp
      fit_coef(3, 14) = 1.2141235983496969E-01_dp
      fit_coef(4, 14) = 3.3706874354447809E-02_dp
      fit_coef(5, 14) = -1.8360609527326543E-02_dp
      fit_coef(6, 14) = -1.8541339179191114E-02_dp
      fit_coef(7, 14) = -8.8529750610819111E-04_dp
      fit_coef(8, 14) = 4.7747853958221521E-03_dp
      fit_coef(9, 14) = 2.4689097936542491E-03_dp
      fit_coef(10, 14) = -3.6719910751692379E-03_dp
      fit_coef(11, 14) = 1.8252247936231616E-03_dp
      fit_coef(12, 14) = -4.0633639834620388E-04_dp
      fit_coef(0, 15) = 1.0551618877519962E-01_dp
      fit_coef(1, 15) = 2.0866982892328420E-01_dp
      fit_coef(2, 15) = 2.3758645562105263E-01_dp
      fit_coef(3, 15) = 1.6022711660180100E-01_dp
      fit_coef(4, 15) = 4.3995832787546506E-02_dp
      fit_coef(5, 15) = -2.5306890576195295E-02_dp
      fit_coef(6, 15) = -2.5168288399177011E-02_dp
      fit_coef(7, 15) = -9.3194538735970800E-04_dp
      fit_coef(8, 15) = 1.1886232926624151E-02_dp
      fit_coef(9, 15) = -3.1847717078905507E-03_dp
      fit_coef(10, 15) = -2.3213252827938664E-03_dp
      fit_coef(11, 15) = 2.4083825205789490E-03_dp
      fit_coef(12, 15) = -7.2529680267051182E-04_dp
      fit_coef(0, 16) = 1.7844852042063208E-01_dp
      fit_coef(1, 16) = 2.8956711047461353E-01_dp
      fit_coef(2, 16) = 2.9424320012524147E-01_dp
      fit_coef(3, 16) = 1.8566958335038000E-01_dp
      fit_coef(4, 16) = 5.6552619188388630E-02_dp
      fit_coef(5, 16) = -1.9536083908459775E-02_dp
      fit_coef(6, 16) = -3.5122564007741566E-02_dp
      fit_coef(7, 16) = -1.5160502263432738E-04_dp
      fit_coef(8, 16) = -1.8782648209524570E-02_dp
      fit_coef(9, 16) = 6.7898727329487735E-02_dp
      fit_coef(10, 16) = -6.6952173441897606E-02_dp
      fit_coef(11, 16) = 3.0244889278910881E-02_dp
      fit_coef(12, 16) = -5.4548300942684411E-03_dp
      fit_coef(0, 17) = 3.0103516162173788E-01_dp
      fit_coef(1, 17) = 3.9093832828366182E-01_dp
      fit_coef(2, 17) = 3.4287403235055175E-01_dp
      fit_coef(3, 17) = 1.8694079224619659E-01_dp
      fit_coef(4, 17) = 5.8759826809024841E-02_dp
      fit_coef(5, 17) = 1.8203291052277864E-03_dp
      fit_coef(6, 17) = -1.4847774478906702E-02_dp
      fit_coef(7, 17) = -3.0225726948560502E-02_dp
      fit_coef(8, 17) = 1.7580090339776970E-02_dp
      fit_coef(9, 17) = -6.3682901560076013E-02_dp
      fit_coef(10, 17) = 1.2826078376018740E-01_dp
      fit_coef(11, 17) = -9.3634524667592312E-02_dp
      fit_coef(12, 17) = 2.4194956180134675E-02_dp
      fit_coef(0, 18) = 5.0459047494440223E-01_dp
      fit_coef(1, 18) = 5.1293218387604667E-01_dp
      fit_coef(2, 18) = 3.8452848514240090E-01_dp
      fit_coef(3, 18) = 1.6703912013390276E-01_dp
      fit_coef(4, 18) = 4.3712588651811067E-02_dp
      fit_coef(5, 18) = 1.9051249953627954E-02_dp
      fit_coef(6, 18) = 2.5883902146128887E-02_dp
      fit_coef(7, 18) = -5.0770600600935537E-02_dp
      fit_coef(8, 18) = 7.7584684197698039E-02_dp
      fit_coef(9, 18) = -1.7926872234928357E-01_dp
      fit_coef(10, 18) = 1.6678298865222213E-01_dp
      fit_coef(11, 18) = -5.5645139244720254E-02_dp
      fit_coef(12, 18) = 2.7420024112079135E-03_dp
      fit_coef(0, 19) = 8.3992248227121313E-01_dp
      fit_coef(1, 19) = 6.4715555436372529E-01_dp
      fit_coef(2, 19) = 4.1796785083991977E-01_dp
      fit_coef(3, 19) = 1.3654955607272076E-01_dp
      fit_coef(4, 19) = 2.1070840652137535E-02_dp
      fit_coef(5, 19) = 3.1743033314222624E-02_dp
      fit_coef(6, 19) = 2.9784087155825711E-02_dp
      fit_coef(7, 19) = 4.1081836375772660E-02_dp
      fit_coef(8, 19) = -1.1439519939475226E-01_dp
      fit_coef(9, 19) = 1.6305361885018524E-01_dp
      fit_coef(10, 19) = -2.2346150468788073E-01_dp
      fit_coef(11, 19) = 1.6036296981154594E-01_dp
      fit_coef(12, 19) = -4.2164356684994284E-02_dp
      fit_coef(0, 20) = 1.3973915504793548E+00_dp
      fit_coef(1, 20) = 7.7764466691628698E-01_dp
      fit_coef(2, 20) = 4.4348607986397731E-01_dp
      fit_coef(3, 20) = 1.0694241439192023E-01_dp
      fit_coef(4, 20) = 2.6327294220292721E-03_dp
      fit_coef(5, 20) = 3.5966741485323687E-02_dp
      fit_coef(6, 20) = 4.1705374001031803E-02_dp
      fit_coef(7, 20) = 3.8907029803020722E-02_dp
      fit_coef(8, 20) = -6.8120058033489070E-02_dp
      fit_coef(9, 20) = 8.2786813906914350E-02_dp
      fit_coef(10, 20) = -9.1368844284995504E-02_dp
      fit_coef(11, 20) = 3.7752943929935098E-02_dp
      fit_coef(12, 20) = -1.9232870385022373E-03_dp
      fit_coef(0, 21) = 2.3760047811723362E+00_dp
      fit_coef(1, 21) = 8.8976406037305356E-01_dp
      fit_coef(2, 21) = 4.7037419059997421E-01_dp
      fit_coef(3, 21) = 9.0428943264851042E-02_dp
      fit_coef(4, 21) = -8.0951525285273343E-03_dp
      fit_coef(5, 21) = 3.5586564490451805E-02_dp
      fit_coef(6, 21) = 5.2837656508520334E-02_dp
      fit_coef(7, 21) = 6.4161552882951054E-03_dp
      fit_coef(8, 21) = 2.7554689120903888E-02_dp
      fit_coef(9, 21) = -6.3909420973179670E-02_dp
      fit_coef(10, 21) = 7.9771394321258737E-02_dp
      fit_coef(11, 21) = -6.5034308453826573E-02_dp
      fit_coef(12, 21) = 2.0094355264141887E-02_dp
      fit_coef(0, 22) = 4.5028432526721405E+00_dp
      fit_coef(1, 22) = 1.0114482591044458E+00_dp
      fit_coef(2, 22) = 5.3486119037036706E-01_dp
      fit_coef(3, 22) = 1.0317932221763447E-01_dp
      fit_coef(4, 22) = -1.2314024715229079E-02_dp
      fit_coef(5, 22) = 3.3690522473128488E-02_dp
      fit_coef(6, 22) = 5.3663513056946451E-02_dp
      fit_coef(7, 22) = 1.0785244432659307E-02_dp
      fit_coef(8, 22) = 2.2642841476459809E-02_dp
      fit_coef(9, 22) = -3.9218764804072916E-02_dp
      fit_coef(10, 22) = 5.2080967974292493E-02_dp
      fit_coef(11, 22) = -3.5217909102649027E-02_dp
      fit_coef(12, 22) = 7.5721478110433137E-03_dp

   END SUBROUTINE get_fit_coef_k11_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k11_mi(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.1533853671352736E+03_dp
      fit_coef(0, 1) = 7.3828894587718401E-04_dp
      fit_coef(1, 1) = 2.0775080037103563E-03_dp
      fit_coef(2, 1) = 2.9605504652653850E-03_dp
      fit_coef(3, 1) = 2.5239761597293463E-03_dp
      fit_coef(4, 1) = 1.1723403873630633E-03_dp
      fit_coef(5, 1) = 4.6491653525586301E-05_dp
      fit_coef(6, 1) = -2.5819644062415552E-04_dp
      fit_coef(7, 1) = -3.9831509458327039E-05_dp
      fit_coef(8, 1) = -4.0413534671757017E-05_dp
      fit_coef(9, 1) = 2.5700229767237757E-04_dp
      fit_coef(10, 1) = -2.8744580127669869E-04_dp
      fit_coef(11, 1) = 1.4858790649392180E-04_dp
      fit_coef(12, 1) = -3.1394148877043504E-05_dp
      fit_coef(0, 2) = 4.3620891975219853E-03_dp
      fit_coef(1, 2) = 1.1465965739592934E-02_dp
      fit_coef(2, 2) = 1.6004195686742342E-02_dp
      fit_coef(3, 2) = 1.3390829030018066E-02_dp
      fit_coef(4, 2) = 6.0723492355871950E-03_dp
      fit_coef(5, 2) = 1.7271227012111563E-04_dp
      fit_coef(6, 2) = -1.3295013600067558E-03_dp
      fit_coef(7, 2) = -2.0819503017025333E-04_dp
      fit_coef(8, 2) = -1.9485312280291003E-04_dp
      fit_coef(9, 2) = 1.3461436514495772E-03_dp
      fit_coef(10, 2) = -1.5381432319053760E-03_dp
      fit_coef(11, 2) = 8.0871429267979347E-04_dp
      fit_coef(12, 2) = -1.7336433357221408E-04_dp
      fit_coef(0, 3) = 1.3393965217340013E-02_dp
      fit_coef(1, 3) = 3.0758021534147448E-02_dp
      fit_coef(2, 3) = 4.0863508679567187E-02_dp
      fit_coef(3, 3) = 3.3137344564780624E-02_dp
      fit_coef(4, 3) = 1.4642382183959756E-02_dp
      fit_coef(5, 3) = 1.8546563333903454E-04_dp
      fit_coef(6, 3) = -3.2981507356429299E-03_dp
      fit_coef(7, 3) = -5.0371329478405473E-04_dp
      fit_coef(8, 3) = -3.3526812040834221E-04_dp
      fit_coef(9, 3) = 3.3638986135397790E-03_dp
      fit_coef(10, 3) = -4.0699164251672606E-03_dp
      fit_coef(11, 3) = 2.2132500698992790E-03_dp
      fit_coef(12, 3) = -4.8530598400547189E-04_dp
      fit_coef(0, 4) = 3.4178972479544684E-02_dp
      fit_coef(1, 4) = 6.5546031770227484E-02_dp
      fit_coef(2, 4) = 8.0042248964935656E-02_dp
      fit_coef(3, 4) = 6.0702076580615758E-02_dp
      fit_coef(4, 4) = 2.5719498592455581E-02_dp
      fit_coef(5, 4) = 6.3119790847243661E-04_dp
      fit_coef(6, 4) = -5.4349702628550432E-03_dp
      fit_coef(7, 4) = -1.2659291521318623E-03_dp
      fit_coef(8, 4) = -1.2216416836694759E-03_dp
      fit_coef(9, 4) = 6.6635331234771687E-03_dp
      fit_coef(10, 4) = -7.2400923612061709E-03_dp
      fit_coef(11, 4) = 3.7136157185037903E-03_dp
      fit_coef(12, 4) = -7.8916627439796453E-04_dp
      fit_coef(0, 5) = 8.0104328961022445E-02_dp
      fit_coef(1, 5) = 1.2586958566619899E-01_dp
      fit_coef(2, 5) = 1.3877416617876437E-01_dp
      fit_coef(3, 5) = 9.4435072367624009E-02_dp
      fit_coef(4, 5) = 3.5574728421505095E-02_dp
      fit_coef(5, 5) = 4.6709813565982229E-04_dp
      fit_coef(6, 5) = -5.5823066584693949E-03_dp
      fit_coef(7, 5) = -2.2508074433607930E-03_dp
      fit_coef(8, 5) = -4.8236909977057691E-04_dp
      fit_coef(9, 5) = 5.8248465053641401E-03_dp
      fit_coef(10, 5) = -7.4300192099394369E-03_dp
      fit_coef(11, 5) = 4.5059815374626708E-03_dp
      fit_coef(12, 5) = -1.1055416581292031E-03_dp
      fit_coef(0, 6) = 1.7731760423288304E-01_dp
      fit_coef(1, 6) = 2.2598655727257688E-01_dp
      fit_coef(2, 6) = 2.2430968436915458E-01_dp
      fit_coef(3, 6) = 1.3311060621110593E-01_dp
      fit_coef(4, 6) = 4.0488867018282172E-02_dp
      fit_coef(5, 6) = -2.1669989611095342E-03_dp
      fit_coef(6, 6) = -3.8187866091584140E-03_dp
      fit_coef(7, 6) = -2.1306428597222633E-03_dp
      fit_coef(8, 6) = 7.7683777527477076E-04_dp
      fit_coef(9, 6) = 4.9822370243334966E-03_dp
      fit_coef(10, 6) = -9.5726170222928900E-03_dp
      fit_coef(11, 6) = 6.8634964686467978E-03_dp
      fit_coef(12, 6) = -1.7999298724748951E-03_dp
      fit_coef(0, 7) = 3.7517258865761766E-01_dp
      fit_coef(1, 7) = 3.8333519695677720E-01_dp
      fit_coef(2, 7) = 3.4375707015157120E-01_dp
      fit_coef(3, 7) = 1.7423153407071543E-01_dp
      fit_coef(4, 7) = 3.7337456774139702E-02_dp
      fit_coef(5, 7) = -7.4660271581285850E-03_dp
      fit_coef(6, 7) = 5.8400107303188187E-04_dp
      fit_coef(7, 7) = -2.5244872448215655E-03_dp
      fit_coef(8, 7) = 4.2884385157431048E-03_dp
      fit_coef(9, 7) = 1.0453917543656458E-03_dp
      fit_coef(10, 7) = -9.2913505401427225E-03_dp
      fit_coef(11, 7) = 8.1349139258150831E-03_dp
      fit_coef(12, 7) = -2.3025720564300821E-03_dp
      fit_coef(0, 8) = 7.6460515920837613E-01_dp
      fit_coef(1, 8) = 6.1549142645707966E-01_dp
      fit_coef(2, 8) = 5.0182800509604020E-01_dp
      fit_coef(3, 8) = 2.1325320085467817E-01_dp
      fit_coef(4, 8) = 2.3360017105729571E-02_dp
      fit_coef(5, 8) = -1.4054384598267826E-02_dp
      fit_coef(6, 8) = 9.0661378414736218E-03_dp
      fit_coef(7, 8) = -4.7433494326856486E-03_dp
      fit_coef(8, 8) = 1.0160041902537189E-02_dp
      fit_coef(9, 8) = -6.9427342461860627E-03_dp
      fit_coef(10, 8) = -5.4907207434525202E-03_dp
      fit_coef(11, 8) = 8.0763791261621969E-03_dp
      fit_coef(12, 8) = -2.6484224072421108E-03_dp
      fit_coef(0, 9) = 1.5114021780363067E+00_dp
      fit_coef(1, 9) = 9.3368349741314893E-01_dp
      fit_coef(2, 9) = 6.9854171505894413E-01_dp
      fit_coef(3, 9) = 2.4498257864430506E-01_dp
      fit_coef(4, 9) = -2.7117345776915117E-03_dp
      fit_coef(5, 9) = -1.9684397407321690E-02_dp
      fit_coef(6, 9) = 2.2750102305246520E-02_dp
      fit_coef(7, 9) = -9.9596165022990610E-03_dp
      fit_coef(8, 9) = 1.8166084113835305E-02_dp
      fit_coef(9, 9) = -1.9043752557357005E-02_dp
      fit_coef(10, 9) = 2.2864257000369307E-03_dp
      fit_coef(11, 9) = 6.3623059846125004E-03_dp
      fit_coef(12, 9) = -2.7468174322941873E-03_dp
      fit_coef(0, 10) = 2.9272838679658890E+00_dp
      fit_coef(1, 10) = 1.3344095709119908E+00_dp
      fit_coef(2, 10) = 9.2872428481763736E-01_dp
      fit_coef(3, 10) = 2.6731240098857367E-01_dp
      fit_coef(4, 10) = -3.9323628602554263E-02_dp
      fit_coef(5, 10) = -2.2664729979783912E-02_dp
      fit_coef(6, 10) = 4.1418459538349613E-02_dp
      fit_coef(7, 10) = -1.8765822971615966E-02_dp
      fit_coef(8, 10) = 2.8023480739935883E-02_dp
      fit_coef(9, 10) = -3.4574149384271238E-02_dp
      fit_coef(10, 10) = 1.3734227321286851E-02_dp
      fit_coef(11, 10) = 2.9343888906998264E-03_dp
      fit_coef(12, 10) = -2.5433341475355697E-03_dp
      fit_coef(0, 11) = 5.7131838450036696E+00_dp
      fit_coef(1, 11) = 1.8005538920264632E+00_dp
      fit_coef(2, 11) = 1.1880536006356908E+00_dp
      fit_coef(3, 11) = 2.8559106021602204E-01_dp
      fit_coef(4, 11) = -8.2639115843313488E-02_dp
      fit_coef(5, 11) = -2.3761780940649874E-02_dp
      fit_coef(6, 11) = 6.3493180678645345E-02_dp
      fit_coef(7, 11) = -3.0273355255164173E-02_dp
      fit_coef(8, 11) = 3.9310329422917428E-02_dp
      fit_coef(9, 11) = -5.2371339180397693E-02_dp
      fit_coef(10, 11) = 2.7484059561085201E-02_dp
      fit_coef(11, 11) = -1.5309209297749292E-03_dp
      fit_coef(12, 11) = -2.1659478188002048E-03_dp
      fit_coef(0, 12) = 1.9426697312926121E-03_dp
      fit_coef(1, 12) = 5.3902355339823542E-03_dp
      fit_coef(2, 12) = 7.6460050852952361E-03_dp
      fit_coef(3, 12) = 6.4861449445603606E-03_dp
      fit_coef(4, 12) = 2.9975041609363725E-03_dp
      fit_coef(5, 12) = 1.1267866675204334E-04_dp
      fit_coef(6, 12) = -6.6052089732986108E-04_dp
      fit_coef(7, 12) = -1.0429593050340101E-04_dp
      fit_coef(8, 12) = -9.6642089460867778E-05_dp
      fit_coef(9, 12) = 6.5468723131356925E-04_dp
      fit_coef(10, 12) = -7.3909338726213167E-04_dp
      fit_coef(11, 12) = 3.8434400387220237E-04_dp
      fit_coef(12, 12) = -8.1587630115757833E-05_dp
      fit_coef(0, 13) = 5.6766583157203338E-03_dp
      fit_coef(1, 13) = 1.3746861086398525E-02_dp
      fit_coef(2, 13) = 1.8677286574891313E-02_dp
      fit_coef(3, 13) = 1.5318320571135682E-02_dp
      fit_coef(4, 13) = 6.7566786656071713E-03_dp
      fit_coef(5, 13) = 7.1919811741634156E-05_dp
      fit_coef(6, 13) = -1.4743272069088146E-03_dp
      fit_coef(7, 13) = -2.0331439016954778E-04_dp
      fit_coef(8, 13) = -1.9066551850151903E-04_dp
      fit_coef(9, 13) = 1.4614366852185527E-03_dp
      fit_coef(10, 13) = -1.7182804418019637E-03_dp
      fit_coef(11, 13) = 9.2385541265890814E-04_dp
      fit_coef(12, 13) = -2.0177461300662419E-04_dp
      fit_coef(0, 14) = 1.3396971882867550E-02_dp
      fit_coef(1, 14) = 2.5745664924886561E-02_dp
      fit_coef(2, 14) = 3.1409108695941872E-02_dp
      fit_coef(3, 14) = 2.3978671080145342E-02_dp
      fit_coef(4, 14) = 1.0268221460486631E-02_dp
      fit_coef(5, 14) = 1.2930798224983239E-04_dp
      fit_coef(6, 14) = -2.4415003996848300E-03_dp
      fit_coef(7, 14) = -4.3479841388968817E-04_dp
      fit_coef(8, 14) = -5.2954762444991550E-04_dp
      fit_coef(9, 14) = 3.3920573298805374E-03_dp
      fit_coef(10, 14) = -3.9601620536706839E-03_dp
      fit_coef(11, 14) = 2.1033319510085288E-03_dp
      fit_coef(12, 14) = -4.5128467093780031E-04_dp
      fit_coef(0, 15) = 3.0313612771117847E-02_dp
      fit_coef(1, 15) = 4.5511722981043855E-02_dp
      fit_coef(2, 15) = 4.7843577097949709E-02_dp
      fit_coef(3, 15) = 3.0848847485647288E-02_dp
      fit_coef(4, 15) = 1.1173391444445794E-02_dp
      fit_coef(5, 15) = 5.3673413927424196E-04_dp
      fit_coef(6, 15) = -1.2788112460830441E-03_dp
      fit_coef(7, 15) = -1.2106504102044236E-03_dp
      fit_coef(8, 15) = -1.3614299287656570E-04_dp
      fit_coef(9, 15) = 8.4193596534928113E-04_dp
      fit_coef(10, 15) = -4.9266376268823403E-04_dp
      fit_coef(11, 15) = 2.4264482966524924E-04_dp
      fit_coef(12, 15) = -8.0800915180974188E-05_dp
      fit_coef(0, 16) = 6.5727773851375626E-02_dp
      fit_coef(1, 16) = 7.7524625497591287E-02_dp
      fit_coef(2, 16) = 7.0866135281020856E-02_dp
      fit_coef(3, 16) = 3.6455648551271543E-02_dp
      fit_coef(4, 16) = 7.9253539051012834E-03_dp
      fit_coef(5, 16) = -1.2114748812317041E-03_dp
      fit_coef(6, 16) = 8.2015105160121449E-04_dp
      fit_coef(7, 16) = -3.0470170403903876E-04_dp
      fit_coef(8, 16) = 8.4675959307711432E-04_dp
      fit_coef(9, 16) = -6.9307175814985243E-04_dp
      fit_coef(10, 16) = -1.7074186347168978E-03_dp
      fit_coef(11, 16) = 2.1004252363773243E-03_dp
      fit_coef(12, 16) = -6.5954156111989743E-04_dp
      fit_coef(0, 17) = 1.3663204308698063E-01_dp
      fit_coef(1, 17) = 1.2566212501641080E-01_dp
      fit_coef(2, 17) = 1.0142604274144125E-01_dp
      fit_coef(3, 17) = 4.0492590765185132E-02_dp
      fit_coef(4, 17) = 1.3893599733611069E-03_dp
      fit_coef(5, 17) = -4.0829840345847287E-03_dp
      fit_coef(6, 17) = 2.8641937242071414E-03_dp
      fit_coef(7, 17) = 1.4034997283466152E-04_dp
      fit_coef(8, 17) = 2.2130109369625853E-03_dp
      fit_coef(9, 17) = -1.9780643203207446E-03_dp
      fit_coef(10, 17) = -1.3761639970619644E-03_dp
      fit_coef(11, 17) = 1.9654825472193074E-03_dp
      fit_coef(12, 17) = -6.0411821444792237E-04_dp
      fit_coef(0, 18) = 2.7371713350977883E-01_dp
      fit_coef(1, 18) = 1.9207380047459294E-01_dp
      fit_coef(2, 18) = 1.3829663189483482E-01_dp
      fit_coef(3, 18) = 4.0969959876443064E-02_dp
      fit_coef(4, 18) = -8.1744595653702960E-03_dp
      fit_coef(5, 18) = -6.2758002950377006E-03_dp
      fit_coef(6, 18) = 6.1933985486500821E-03_dp
      fit_coef(7, 18) = -1.1137257661482085E-03_dp
      fit_coef(8, 18) = 4.7491704854226244E-03_dp
      fit_coef(9, 18) = -5.9366942385150516E-03_dp
      fit_coef(10, 18) = 1.9770890986731604E-03_dp
      fit_coef(11, 18) = 6.2460690633304441E-04_dp
      fit_coef(12, 18) = -4.2628349475002297E-04_dp
      fit_coef(0, 19) = 5.3176214769076025E-01_dp
      fit_coef(1, 19) = 2.7427114346332437E-01_dp
      fit_coef(2, 19) = 1.7786670623161183E-01_dp
      fit_coef(3, 19) = 3.6141984545340726E-02_dp
      fit_coef(4, 19) = -2.0019481045233233E-02_dp
      fit_coef(5, 19) = -6.5029597037300265E-03_dp
      fit_coef(6, 19) = 1.0978881646568977E-02_dp
      fit_coef(7, 19) = -3.5326165348373130E-03_dp
      fit_coef(8, 19) = 6.9664331496620468E-03_dp
      fit_coef(9, 19) = -1.0072536637310938E-02_dp
      fit_coef(10, 19) = 5.7316495954974223E-03_dp
      fit_coef(11, 19) = -8.1951094019536695E-04_dp
      fit_coef(12, 19) = -2.4352440718963261E-04_dp
      fit_coef(0, 20) = 1.0110078977345236E+00_dp
      fit_coef(1, 20) = 3.6158811802096996E-01_dp
      fit_coef(2, 20) = 2.1462317218559554E-01_dp
      fit_coef(3, 20) = 2.6882407712661614E-02_dp
      fit_coef(4, 20) = -3.1860390881090057E-02_dp
      fit_coef(5, 20) = -4.4365197992393256E-03_dp
      fit_coef(6, 20) = 1.6356522282533974E-02_dp
      fit_coef(7, 20) = -7.0191772084550030E-03_dp
      fit_coef(8, 20) = 9.0001998656736130E-03_dp
      fit_coef(9, 20) = -1.4010401343075048E-02_dp
      fit_coef(10, 20) = 9.7773201744224356E-03_dp
      fit_coef(11, 20) = -2.6263089272930114E-03_dp
      fit_coef(12, 20) = 5.8079667692654174E-05_dp
      fit_coef(0, 21) = 1.9229902544572308E+00_dp
      fit_coef(1, 21) = 4.3639570130727606E-01_dp
      fit_coef(2, 21) = 2.4465083184076966E-01_dp
      fit_coef(3, 21) = 1.8632226021993891E-02_dp
      fit_coef(4, 21) = -4.0654823075276181E-02_dp
      fit_coef(5, 21) = -1.6438863848394413E-03_dp
      fit_coef(6, 21) = 2.0679654468894963E-02_dp
      fit_coef(7, 21) = -1.0433989825777933E-02_dp
      fit_coef(8, 21) = 1.0631019238151107E-02_dp
      fit_coef(9, 21) = -1.6823436574756146E-02_dp
      fit_coef(10, 21) = 1.2861656736955324E-02_dp
      fit_coef(11, 21) = -4.1063260520761953E-03_dp
      fit_coef(12, 21) = 3.2620599079153290E-04_dp
      fit_coef(0, 22) = 3.9876478956152801E+00_dp
      fit_coef(1, 22) = 4.9699973939632702E-01_dp
      fit_coef(2, 22) = 2.7726332870042114E-01_dp
      fit_coef(3, 22) = 2.0557461247107046E-02_dp
      fit_coef(4, 22) = -4.5589994159230353E-02_dp
      fit_coef(5, 22) = -1.1953594192753122E-03_dp
      fit_coef(6, 22) = 2.3317546922334211E-02_dp
      fit_coef(7, 22) = -1.2287704350934838E-02_dp
      fit_coef(8, 22) = 1.1987190724120066E-02_dp
      fit_coef(9, 22) = -1.8798325292619430E-02_dp
      fit_coef(10, 22) = 1.4466982948530036E-02_dp
      fit_coef(11, 22) = -4.6627474463343845E-03_dp
      fit_coef(12, 22) = 3.8343338955871208E-04_dp

   END SUBROUTINE get_fit_coef_k11_mi

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k11_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.0890000000000000E+05_dp
      fit_coef(0, 1) = 3.5476589277685659E-05_dp
      fit_coef(1, 1) = 1.4120924025631340E-05_dp
      fit_coef(2, 1) = 3.1076536118129984E-03_dp
      fit_coef(3, 1) = -2.0920445403620737E-02_dp
      fit_coef(4, 1) = 1.5503824842219388E-01_dp
      fit_coef(5, 1) = -6.5946634106067126E-01_dp
      fit_coef(6, 1) = 1.9567718308121806E+00_dp
      fit_coef(7, 1) = -4.0223034504541495E+00_dp
      fit_coef(8, 1) = 5.7073753184657754E+00_dp
      fit_coef(9, 1) = -5.4854200724345779E+00_dp
      fit_coef(10, 1) = 3.4084108118688103E+00_dp
      fit_coef(11, 1) = -1.2351133915217054E+00_dp
      fit_coef(12, 1) = 1.9818096546374844E-01_dp
      fit_coef(0, 2) = 3.9040597969308965E-04_dp
      fit_coef(1, 2) = 8.7958193593917847E-05_dp
      fit_coef(2, 2) = 1.9518122056059333E-02_dp
      fit_coef(3, 2) = -1.3498019724834617E-01_dp
      fit_coef(4, 2) = 9.8104943789195032E-01_dp
      fit_coef(5, 2) = -4.1896253590373966E+00_dp
      fit_coef(6, 2) = 1.2439540817050382E+01_dp
      fit_coef(7, 2) = -2.5590801859617393E+01_dp
      fit_coef(8, 2) = 3.6347431140247103E+01_dp
      fit_coef(9, 2) = -3.4967385668188783E+01_dp
      fit_coef(10, 2) = 2.1745312479785735E+01_dp
      fit_coef(11, 2) = -7.8852461430485148E+00_dp
      fit_coef(12, 2) = 1.2659121491127978E+00_dp
      fit_coef(0, 3) = 2.0665460609606825E-03_dp
      fit_coef(1, 3) = 3.0623598182770673E-04_dp
      fit_coef(2, 3) = 6.8781938169008475E-02_dp
      fit_coef(3, 3) = -4.9168954382575458E-01_dp
      fit_coef(4, 3) = 3.4930362601920923E+00_dp
      fit_coef(5, 3) = -1.5167014240164407E+01_dp
      fit_coef(6, 3) = 4.5295440057571241E+01_dp
      fit_coef(7, 3) = -9.3377006161805866E+01_dp
      fit_coef(8, 3) = 1.3271885253969418E+02_dp
      fit_coef(9, 3) = -1.2770077780864555E+02_dp
      fit_coef(10, 3) = 7.9411953123430422E+01_dp
      fit_coef(11, 3) = -2.8793930151510978E+01_dp
      fit_coef(12, 3) = 4.6222166428791578E+00_dp
      fit_coef(0, 4) = 7.9767528386571820E-03_dp
      fit_coef(1, 4) = 8.4910081412428783E-04_dp
      fit_coef(2, 4) = 1.9200935759540066E-01_dp
      fit_coef(3, 4) = -1.4007008225024304E+00_dp
      fit_coef(4, 4) = 9.8035873307907426E+00_dp
      fit_coef(5, 4) = -4.2999162330329014E+01_dp
      fit_coef(6, 4) = 1.2867236226296961E+02_dp
      fit_coef(7, 4) = -2.6546890298183149E+02_dp
      fit_coef(8, 4) = 3.7764422948289194E+02_dp
      fit_coef(9, 4) = -3.6371847327447421E+02_dp
      fit_coef(10, 4) = 2.2639866330975411E+02_dp
      fit_coef(11, 4) = -8.2161179736409451E+01_dp
      fit_coef(12, 4) = 1.3199047691946380E+01_dp
      fit_coef(0, 5) = 2.5545969787194346E-02_dp
      fit_coef(1, 5) = 2.0342070978811074E-03_dp
      fit_coef(2, 5) = 4.6216853849385225E-01_dp
      fit_coef(3, 5) = -3.4177493858793202E+00_dp
      fit_coef(4, 5) = 2.3691045823517857E+01_dp
      fit_coef(5, 5) = -1.0461938927709393E+02_dp
      fit_coef(6, 5) = 3.1370386858745627E+02_dp
      fit_coef(7, 5) = -6.4798973211308271E+02_dp
      fit_coef(8, 5) = 9.2266003511934855E+02_dp
      fit_coef(9, 5) = -8.8927522399842849E+02_dp
      fit_coef(10, 5) = 5.5383510055109002E+02_dp
      fit_coef(11, 5) = -2.0107166636215129E+02_dp
      fit_coef(12, 5) = 3.2311819146896454E+01_dp
      fit_coef(0, 6) = 7.2076697268762174E-02_dp
      fit_coef(1, 6) = 4.3697938638884411E-03_dp
      fit_coef(2, 6) = 9.9608072915738466E-01_dp
      fit_coef(3, 6) = -7.4383652629013977E+00_dp
      fit_coef(4, 6) = 5.1208554327170454E+01_dp
      fit_coef(5, 6) = -2.2721354633426219E+02_dp
      fit_coef(6, 6) = 6.8250267908846160E+02_dp
      fit_coef(7, 6) = -1.4111404970254437E+03_dp
      fit_coef(8, 6) = 2.0105465406328890E+03_dp
      fit_coef(9, 6) = -1.9386324350975158E+03_dp
      fit_coef(10, 6) = 1.2077495454590323E+03_dp
      fit_coef(11, 6) = -4.3858368013964258E+02_dp
      fit_coef(12, 6) = 7.0492945521981113E+01_dp
      fit_coef(0, 7) = 1.8528902440326048E-01_dp
      fit_coef(1, 7) = 8.5761877628560846E-03_dp
      fit_coef(2, 7) = 1.9596400785319805E+00_dp
      fit_coef(3, 7) = -1.4742195445863704E+01_dp
      fit_coef(4, 7) = 1.0097053276986436E+02_dp
      fit_coef(5, 7) = -4.4958224203673950E+02_dp
      fit_coef(6, 7) = 1.3525059945742282E+03_dp
      fit_coef(7, 7) = -2.7986628560219651E+03_dp
      fit_coef(8, 7) = 3.9893201483629923E+03_dp
      fit_coef(9, 7) = -3.8478121346267285E+03_dp
      fit_coef(10, 7) = 2.3976749678349697E+03_dp
      fit_coef(11, 7) = -8.7083742454840910E+02_dp
      fit_coef(12, 7) = 1.3998610570704449E+02_dp
      fit_coef(0, 8) = 4.4331726198419791E-01_dp
      fit_coef(1, 8) = 1.5511309832479217E-02_dp
      fit_coef(2, 8) = 3.5507741823714190E+00_dp
      fit_coef(3, 8) = -2.6865693016984693E+01_dp
      fit_coef(4, 8) = 1.8327560251896824E+02_dp
      fit_coef(5, 8) = -8.1823836532717723E+02_dp
      fit_coef(6, 8) = 2.4647768105667869E+03_dp
      fit_coef(7, 8) = -5.1036331419270782E+03_dp
      fit_coef(8, 8) = 7.2776337187971394E+03_dp
      fit_coef(9, 8) = -7.0211118090420277E+03_dp
      fit_coef(10, 8) = 4.3757325625709809E+03_dp
      fit_coef(11, 8) = -1.5894501984451088E+03_dp
      fit_coef(12, 8) = 2.5552396179932322E+02_dp
      fit_coef(0, 9) = 1.0023832039292297E+00_dp
      fit_coef(1, 9) = 2.5898302934575013E-02_dp
      fit_coef(2, 9) = 5.9367320140376343E+00_dp
      fit_coef(3, 9) = -4.5119163651256045E+01_dp
      fit_coef(4, 9) = 3.0685214719931554E+02_dp
      fit_coef(5, 9) = -1.3727570869476449E+03_dp
      fit_coef(6, 9) = 4.1397102479948817E+03_dp
      fit_coef(7, 9) = -8.5765932664076063E+03_dp
      fit_coef(8, 9) = 1.2233664390647617E+04_dp
      fit_coef(9, 9) = -1.1804562465563917E+04_dp
      fit_coef(10, 9) = 7.3577527032167036E+03_dp
      fit_coef(11, 9) = -2.6728608608039476E+03_dp
      fit_coef(12, 9) = 4.2972086937637658E+02_dp
      fit_coef(0, 10) = 2.1747295228003427E+00_dp
      fit_coef(1, 10) = 3.9825856965749605E-02_dp
      fit_coef(2, 10) = 9.1386038693541103E+00_dp
      fit_coef(3, 10) = -6.9684875108374484E+01_dp
      fit_coef(4, 10) = 4.7283781908541010E+02_dp
      fit_coef(5, 10) = -2.1185130397402854E+03_dp
      fit_coef(6, 10) = 6.3941973699549608E+03_dp
      fit_coef(7, 10) = -1.3253289025658554E+04_dp
      fit_coef(8, 10) = 1.8908933653764056E+04_dp
      fit_coef(9, 10) = -1.8248136599625253E+04_dp
      fit_coef(10, 10) = 1.1374965884570694E+04_dp
      fit_coef(11, 10) = -4.1324304505338878E+03_dp
      fit_coef(12, 10) = 6.6440472838840810E+02_dp
      fit_coef(0, 11) = 4.6737127719821085E+00_dp
      fit_coef(1, 11) = 5.6471682863233792E-02_dp
      fit_coef(2, 11) = 1.2966744420426899E+01_dp
      fit_coef(3, 11) = -9.9092986263943558E+01_dp
      fit_coef(4, 11) = 6.7136997890125622E+02_dp
      fit_coef(5, 11) = -3.0109895935607155E+03_dp
      fit_coef(6, 11) = 9.0933208665348502E+03_dp
      fit_coef(7, 11) = -1.8853514883165644E+04_dp
      fit_coef(8, 11) = 2.6903247532779991E+04_dp
      fit_coef(9, 11) = -2.5965395036313621E+04_dp
      fit_coef(10, 11) = 1.6186393937419709E+04_dp
      fit_coef(11, 11) = -5.8805907336782202E+03_dp
      fit_coef(12, 11) = 9.4549381336496526E+02_dp
      fit_coef(0, 12) = 1.1342816745586836E-04_dp
      fit_coef(1, 12) = 3.6576968324737305E-05_dp
      fit_coef(2, 12) = 8.1035507833661872E-03_dp
      fit_coef(3, 12) = -5.4754734634044767E-02_dp
      fit_coef(4, 12) = 4.0635660271826546E-01_dp
      fit_coef(5, 12) = -1.7252011567424357E+00_dp
      fit_coef(6, 12) = 5.0991785559775868E+00_dp
      fit_coef(7, 12) = -1.0443089131832354E+01_dp
      fit_coef(8, 12) = 1.4772607027520992E+01_dp
      fit_coef(9, 12) = -1.4163309607497387E+01_dp
      fit_coef(10, 12) = 8.7833832339673332E+00_dp
      fit_coef(11, 12) = -3.1779319684136769E+00_dp
      fit_coef(12, 12) = 5.0928471655731977E-01_dp
      fit_coef(0, 13) = 7.4546004903846121E-04_dp
      fit_coef(1, 13) = 1.2519052310895946E-04_dp
      fit_coef(2, 13) = 2.7990742019097205E-02_dp
      fit_coef(3, 13) = -1.9878830068107795E-01_dp
      fit_coef(4, 13) = 1.4175816998714001E+00_dp
      fit_coef(5, 13) = -6.1266447654072609E+00_dp
      fit_coef(6, 13) = 1.8316162032568503E+01_dp
      fit_coef(7, 13) = -3.7851604364859433E+01_dp
      fit_coef(8, 13) = 5.3922831469585589E+01_dp
      fit_coef(9, 13) = -5.1975301110470006E+01_dp
      fit_coef(10, 13) = 3.2361327142046470E+01_dp
      fit_coef(11, 13) = -1.1743599901510185E+01_dp
      fit_coef(12, 13) = 1.8861805928420110E+00_dp
      fit_coef(0, 14) = 3.0519693802336465E-03_dp
      fit_coef(1, 14) = 3.4087975510808114E-04_dp
      fit_coef(2, 14) = 7.7167269766524146E-02_dp
      fit_coef(3, 14) = -5.6411267648401553E-01_dp
      fit_coef(4, 14) = 3.9421610586910765E+00_dp
      fit_coef(5, 14) = -1.7305564712941042E+01_dp
      fit_coef(6, 14) = 5.1769587328365290E+01_dp
      fit_coef(7, 14) = -1.0670184531108517E+02_dp
      fit_coef(8, 14) = 1.5162660340890758E+02_dp
      fit_coef(9, 14) = -1.4590770106048981E+02_dp
      fit_coef(10, 14) = 9.0764737150408479E+01_dp
      fit_coef(11, 14) = -3.2925715085693440E+01_dp
      fit_coef(12, 14) = 5.2881884424153274E+00_dp
      fit_coef(0, 15) = 9.9219657650822073E-03_dp
      fit_coef(1, 15) = 7.9669846000388241E-04_dp
      fit_coef(2, 15) = 1.8124096946957050E-01_dp
      fit_coef(3, 15) = -1.3450166921157356E+00_dp
      fit_coef(4, 15) = 9.2993195019520840E+00_dp
      fit_coef(5, 15) = -4.1141390542063078E+01_dp
      fit_coef(6, 15) = 1.2339197373801343E+02_dp
      fit_coef(7, 15) = -2.5492312646946198E+02_dp
      fit_coef(8, 15) = 3.6308859133318100E+02_dp
      fit_coef(9, 15) = -3.5006957300050141E+02_dp
      fit_coef(10, 15) = 2.1808828321735942E+02_dp
      fit_coef(11, 15) = -7.9197426579451133E+01_dp
      fit_coef(12, 15) = 1.2729290836050698E+01_dp
      fit_coef(0, 16) = 2.7934155736110411E-02_dp
      fit_coef(1, 16) = 1.6577256632979880E-03_dp
      fit_coef(2, 16) = 3.7842445866661195E-01_dp
      fit_coef(3, 16) = -2.8373077002350295E+00_dp
      fit_coef(4, 16) = 1.9477723226227106E+01_dp
      fit_coef(5, 16) = -8.6598640916868774E+01_dp
      fit_coef(6, 16) = 2.6025925292176282E+02_dp
      fit_coef(7, 16) = -5.3828360062495688E+02_dp
      fit_coef(8, 16) = 7.6712698683802557E+02_dp
      fit_coef(9, 16) = -7.3983293131025414E+02_dp
      fit_coef(10, 16) = 4.6097553662229262E+02_dp
      fit_coef(11, 16) = -1.6741742993909008E+02_dp
      fit_coef(12, 16) = 2.6911040347084413E+01_dp
      fit_coef(0, 17) = 7.1098068080604329E-02_dp
      fit_coef(1, 17) = 3.1344747909274207E-03_dp
      fit_coef(2, 17) = 7.1727119916633186E-01_dp
      fit_coef(3, 17) = -5.4188366020725613E+00_dp
      fit_coef(4, 17) = 3.7004431896986794E+01_dp
      fit_coef(5, 17) = -1.6510685856094676E+02_dp
      fit_coef(6, 17) = 4.9706150634055501E+02_dp
      fit_coef(7, 17) = -1.0289510132785595E+03_dp
      fit_coef(8, 17) = 1.4670886352439184E+03_dp
      fit_coef(9, 17) = -1.4153076846774175E+03_dp
      fit_coef(10, 17) = 8.8203709163857422E+02_dp
      fit_coef(11, 17) = -3.2039085416412871E+02_dp
      fit_coef(12, 17) = 5.1506838572455884E+01_dp
      fit_coef(0, 18) = 1.6776001286461931E-01_dp
      fit_coef(1, 18) = 5.4256243545789283E-03_dp
      fit_coef(2, 18) = 1.2438301964496239E+00_dp
      fit_coef(3, 18) = -9.4525365959341521E+00_dp
      fit_coef(4, 18) = 6.4287530580546374E+01_dp
      fit_coef(5, 18) = -2.8761451534054112E+02_dp
      fit_coef(6, 18) = 8.6715581047564865E+02_dp
      fit_coef(7, 18) = -1.7963912098780600E+03_dp
      fit_coef(8, 18) = 2.5623050325486606E+03_dp
      fit_coef(9, 18) = -2.4724333743568150E+03_dp
      fit_coef(10, 18) = 1.5410800293083496E+03_dp
      fit_coef(11, 18) = -5.5983897666404175E+02_dp
      fit_coef(12, 18) = 9.0007748058412503E+01_dp
      fit_coef(0, 19) = 3.7318901847591823E-01_dp
      fit_coef(1, 19) = 8.5741894158780207E-03_dp
      fit_coef(2, 19) = 1.9684022818381477E+00_dp
      fit_coef(3, 19) = -1.5029124759115470E+01_dp
      fit_coef(4, 19) = 1.0188619253386351E+02_dp
      fit_coef(5, 19) = -4.5678606319996516E+02_dp
      fit_coef(6, 19) = 1.3789467638499050E+03_dp
      fit_coef(7, 19) = -2.8584140258154093E+03_dp
      fit_coef(8, 19) = 4.0784585124208406E+03_dp
      fit_coef(9, 19) = -3.9361263103870356E+03_dp
      fit_coef(10, 19) = 2.4536814783734530E+03_dp
      fit_coef(11, 19) = -8.9143170846733369E+02_dp
      fit_coef(12, 19) = 1.4332667121452164E+02_dp
      fit_coef(0, 20) = 7.9411999812197465E-01_dp
      fit_coef(1, 20) = 1.2222970002666296E-02_dp
      fit_coef(2, 20) = 2.8089778868006130E+00_dp
      fit_coef(3, 20) = -2.1523724383534741E+01_dp
      fit_coef(4, 20) = 1.4555880546046066E+02_dp
      fit_coef(5, 20) = -6.5361423760505409E+02_dp
      fit_coef(6, 20) = 1.9751646032271344E+03_dp
      fit_coef(7, 20) = -4.0964475425550763E+03_dp
      fit_coef(8, 20) = 5.8464616672416569E+03_dp
      fit_coef(9, 20) = -5.6432337456166642E+03_dp
      fit_coef(10, 20) = 3.5181407559250561E+03_dp
      fit_coef(11, 20) = -1.2782197709057280E+03_dp
      fit_coef(12, 20) = 2.0552243167477678E+02_dp
      fit_coef(0, 21) = 1.6552401388491180E+00_dp
      fit_coef(1, 21) = 1.5475484985637764E-02_dp
      fit_coef(2, 21) = 3.5586577312092178E+00_dp
      fit_coef(3, 21) = -2.7328423530790474E+01_dp
      fit_coef(4, 21) = 1.8453595318851117E+02_dp
      fit_coef(5, 21) = -8.2943501604533219E+02_dp
      fit_coef(6, 21) = 2.5081783317335353E+03_dp
      fit_coef(7, 21) = -5.2037296899464609E+03_dp
      fit_coef(8, 21) = 7.4280891345363261E+03_dp
      fit_coef(9, 21) = -7.1705455707166202E+03_dp
      fit_coef(10, 21) = 4.4705361605949456E+03_dp
      fit_coef(11, 21) = -1.6242951789399854E+03_dp
      fit_coef(12, 21) = 2.6117201512280747E+02_dp
      fit_coef(0, 22) = 3.6819715912926032E+00_dp
      fit_coef(1, 22) = 1.7724765507747240E-02_dp
      fit_coef(2, 22) = 4.0762275279679621E+00_dp
      fit_coef(3, 22) = -3.1313218015726726E+01_dp
      fit_coef(4, 22) = 2.1139684580404636E+02_dp
      fit_coef(5, 22) = -9.5029615455661258E+02_dp
      fit_coef(6, 22) = 2.8739893385614673E+03_dp
      fit_coef(7, 22) = -5.9630507500813519E+03_dp
      fit_coef(8, 22) = 8.5122490258366070E+03_dp
      fit_coef(9, 22) = -8.2172415036707516E+03_dp
      fit_coef(10, 22) = 5.1231461949378099E+03_dp
      fit_coef(11, 22) = -1.8614171279369891E+03_dp
      fit_coef(12, 22) = 2.9929960798773016E+02_dp

   END SUBROUTINE get_fit_coef_k11_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k12_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 4.7621173842280597E-03_dp
      fit_coef(1, 1) = 1.2779124179361450E-02_dp
      fit_coef(2, 1) = 1.6703410753614655E-02_dp
      fit_coef(3, 1) = 1.2304881409550545E-02_dp
      fit_coef(4, 1) = 3.7749685318705707E-03_dp
      fit_coef(5, 1) = -1.7434526435310602E-03_dp
      fit_coef(6, 1) = -2.0397823012582900E-03_dp
      fit_coef(7, 1) = -1.9667552585414541E-04_dp
      fit_coef(8, 1) = 4.9733991832847789E-04_dp
      fit_coef(9, 1) = 3.3438174942184307E-04_dp
      fit_coef(10, 1) = -4.5000255284931892E-04_dp
      fit_coef(11, 1) = 2.2568471250347118E-04_dp
      fit_coef(12, 1) = -5.0451245801427364E-05_dp
      fit_coef(0, 2) = 2.5892231310748540E-02_dp
      fit_coef(1, 2) = 6.8272349594299786E-02_dp
      fit_coef(2, 2) = 8.8380776497241700E-02_dp
      fit_coef(3, 2) = 6.4547463586032264E-02_dp
      fit_coef(4, 2) = 1.9512355251434594E-02_dp
      fit_coef(5, 2) = -9.2451673534127178E-03_dp
      fit_coef(6, 2) = -1.0592986518352070E-02_dp
      fit_coef(7, 2) = -9.5129033192623013E-04_dp
      fit_coef(8, 2) = 2.6222990485926024E-03_dp
      fit_coef(9, 2) = 1.6749335092584598E-03_dp
      fit_coef(10, 2) = -2.2956689993050795E-03_dp
      fit_coef(11, 2) = 1.1496135531876962E-03_dp
      fit_coef(12, 2) = -2.5655891388629779E-04_dp
      fit_coef(0, 3) = 6.7567184232783226E-02_dp
      fit_coef(1, 3) = 1.7186061352294887E-01_dp
      fit_coef(2, 3) = 2.1865383385874831E-01_dp
      fit_coef(3, 3) = 1.5733674301194298E-01_dp
      fit_coef(4, 3) = 4.6341817240479470E-02_dp
      fit_coef(5, 3) = -2.2920074974133117E-02_dp
      fit_coef(6, 3) = -2.5337146772633146E-02_dp
      fit_coef(7, 3) = -1.9861858355440082E-03_dp
      fit_coef(8, 3) = 6.4137367343312197E-03_dp
      fit_coef(9, 3) = 3.7881539374578659E-03_dp
      fit_coef(10, 3) = -5.3502220389785920E-03_dp
      fit_coef(11, 3) = 2.6697483796138425E-03_dp
      fit_coef(12, 3) = -5.9317689227171935E-04_dp
      fit_coef(0, 4) = 1.3773723556451323E-01_dp
      fit_coef(1, 4) = 3.2976491191003016E-01_dp
      fit_coef(2, 4) = 4.0862257137690161E-01_dp
      fit_coef(3, 4) = 2.8864406700326678E-01_dp
      fit_coef(4, 4) = 8.2368391636866115E-02_dp
      fit_coef(5, 4) = -4.3117011034440302E-02_dp
      fit_coef(6, 4) = -4.5352842098712653E-02_dp
      fit_coef(7, 4) = -2.5963365652366087E-03_dp
      fit_coef(8, 4) = 1.2026849031951965E-02_dp
      fit_coef(9, 4) = 5.5083059863503094E-03_dp
      fit_coef(10, 4) = -8.5515536696519658E-03_dp
      fit_coef(11, 4) = 4.2807103276345489E-03_dp
      fit_coef(12, 4) = -9.5472120253031029E-04_dp
      fit_coef(0, 5) = 2.5166827573794959E-01_dp
      fit_coef(1, 5) = 5.5126096550816395E-01_dp
      fit_coef(2, 5) = 6.5643230451689671E-01_dp
      fit_coef(3, 5) = 4.5390238114819420E-01_dp
      fit_coef(4, 5) = 1.2801802686296021E-01_dp
      fit_coef(5, 5) = -6.8438409988043361E-02_dp
      fit_coef(6, 5) = -7.1938494310622186E-02_dp
      fit_coef(7, 5) = -3.7818617497551153E-03_dp
      fit_coef(8, 5) = 2.2624574898233631E-02_dp
      fit_coef(9, 5) = 7.8302540970164088E-03_dp
      fit_coef(10, 5) = -1.7136424214551770E-02_dp
      fit_coef(11, 5) = 9.9150152558899352E-03_dp
      fit_coef(12, 5) = -2.3421367209343525E-03_dp
      fit_coef(0, 6) = 4.3589892832098903E-01_dp
      fit_coef(1, 6) = 8.5015779653011903E-01_dp
      fit_coef(2, 6) = 9.5612932461007116E-01_dp
      fit_coef(3, 6) = 6.4033701313877844E-01_dp
      fit_coef(4, 6) = 1.8409747302118074E-01_dp
      fit_coef(5, 6) = -8.6870444790686202E-02_dp
      fit_coef(6, 6) = -1.0566638795756129E-01_dp
      fit_coef(7, 6) = -4.6671972769434591E-03_dp
      fit_coef(8, 6) = 3.7389422606596834E-03_dp
      fit_coef(9, 6) = 6.9425181522878909E-02_dp
      fit_coef(10, 6) = -7.1984629461529587E-02_dp
      fit_coef(11, 6) = 3.1646309655079689E-02_dp
      fit_coef(12, 6) = -5.6112243049874127E-03_dp
      fit_coef(0, 7) = 7.3351111106633227E-01_dp
      fit_coef(1, 7) = 1.2445631654341636E+00_dp
      fit_coef(2, 7) = 1.3007750122826698E+00_dp
      fit_coef(3, 7) = 8.2747315487381601E-01_dp
      fit_coef(4, 7) = 2.4170033735886270E-01_dp
      fit_coef(5, 7) = -8.6239790315945988E-02_dp
      fit_coef(6, 7) = -1.2071103017063417E-01_dp
      fit_coef(7, 7) = -3.2654987483749273E-02_dp
      fit_coef(8, 7) = 1.8998066989690098E-02_dp
      fit_coef(9, 7) = 1.4613127895351082E-02_dp
      fit_coef(10, 7) = 4.0069718709076961E-02_dp
      fit_coef(11, 7) = -5.0009415078673951E-02_dp
      fit_coef(12, 7) = 1.5360303269937313E-02_dp
      fit_coef(0, 8) = 1.2123305472983175E+00_dp
      fit_coef(1, 8) = 1.7528922824057094E+00_dp
      fit_coef(2, 8) = 1.6839617270590670E+00_dp
      fit_coef(3, 8) = 9.9693839370582804E-01_dp
      fit_coef(4, 8) = 2.8670772642350745E-01_dp
      fit_coef(5, 8) = -6.8905744965899102E-02_dp
      fit_coef(6, 8) = -9.8654503643879704E-02_dp
      fit_coef(7, 8) = -8.3833733804509133E-02_dp
      fit_coef(8, 8) = 9.6718797026550299E-02_dp
      fit_coef(9, 8) = -1.6741156687903433E-01_dp
      fit_coef(10, 8) = 2.2372278410722946E-01_dp
      fit_coef(11, 8) = -1.2395029372556110E-01_dp
      fit_coef(12, 8) = 2.4118562336799924E-02_dp
      fit_coef(0, 9) = 1.9789993531060730E+00_dp
      fit_coef(1, 9) = 2.3876726133812909E+00_dp
      fit_coef(2, 9) = 2.0983662040251230E+00_dp
      fit_coef(3, 9) = 1.1384504548954626E+00_dp
      fit_coef(4, 9) = 3.1166074671817184E-01_dp
      fit_coef(5, 9) = -3.8607395105363457E-02_dp
      fit_coef(6, 9) = -6.8592782585259429E-02_dp
      fit_coef(7, 9) = -5.9599143650154687E-02_dp
      fit_coef(8, 9) = 1.5691989539307449E-02_dp
      fit_coef(9, 9) = -6.2725923481489740E-02_dp
      fit_coef(10, 9) = 6.0324082207683011E-02_dp
      fit_coef(11, 9) = 8.7062543323561756E-03_dp
      fit_coef(12, 9) = -1.3873325190108787E-02_dp
      fit_coef(0, 10) = 3.2058934682668765E+00_dp
      fit_coef(1, 10) = 3.1496633677047572E+00_dp
      fit_coef(2, 10) = 2.5365793238247520E+00_dp
      fit_coef(3, 10) = 1.2508497217707608E+00_dp
      fit_coef(4, 10) = 3.1872392758788848E-01_dp
      fit_coef(5, 10) = -2.0600061181166207E-03_dp
      fit_coef(6, 10) = -3.1041845624322214E-02_dp
      fit_coef(7, 10) = -1.2321496532646232E-02_dp
      fit_coef(8, 10) = -8.3121266244633590E-02_dp
      fit_coef(9, 10) = 7.0964635310385341E-02_dp
      fit_coef(10, 10) = -9.6237663309511903E-02_dp
      fit_coef(11, 10) = 9.4035522256218335E-02_dp
      fit_coef(12, 10) = -2.9123579370706108E-02_dp
      fit_coef(0, 11) = 5.1968506439191477E+00_dp
      fit_coef(1, 11) = 4.0268304189360205E+00_dp
      fit_coef(2, 11) = 2.9970061410634696E+00_dp
      fit_coef(3, 11) = 1.3425771299145335E+00_dp
      fit_coef(4, 11) = 3.1423213151024637E-01_dp
      fit_coef(5, 11) = 3.5373029368206500E-02_dp
      fit_coef(6, 11) = 1.9061738352272550E-02_dp
      fit_coef(7, 11) = 2.6409588889696359E-03_dp
      fit_coef(8, 11) = -8.0417831524322025E-02_dp
      fit_coef(9, 11) = 4.2275831208402226E-02_dp
      fit_coef(10, 11) = -5.5418627131360906E-02_dp
      fit_coef(11, 11) = 4.9224692889938519E-02_dp
      fit_coef(12, 11) = -1.2095050837781235E-02_dp
      fit_coef(0, 12) = 8.6390361849437944E+00_dp
      fit_coef(1, 12) = 5.0136507853028878E+00_dp
      fit_coef(2, 12) = 3.4986136939910475E+00_dp
      fit_coef(3, 12) = 1.4331799180952967E+00_dp
      fit_coef(4, 12) = 3.0383064834623563E-01_dp
      fit_coef(5, 12) = 7.1865533368803650E-02_dp
      fit_coef(6, 12) = 7.2720154429931433E-02_dp
      fit_coef(7, 12) = 9.0385929014071293E-03_dp
      fit_coef(8, 12) = -5.0215078840340335E-02_dp
      fit_coef(9, 12) = -1.6285869510121675E-02_dp
      fit_coef(10, 12) = 1.8810174768165218E-02_dp
      fit_coef(11, 12) = -4.9947183805966527E-03_dp
      fit_coef(12, 12) = 2.0646826965866082E-03_dp
      fit_coef(0, 13) = 1.2309045987914696E-02_dp
      fit_coef(1, 13) = 3.2903305474272386E-02_dp
      fit_coef(2, 13) = 4.2910210824443497E-02_dp
      fit_coef(3, 13) = 3.1546310458948443E-02_dp
      fit_coef(4, 13) = 9.6441395261976804E-03_dp
      fit_coef(5, 13) = -4.4817476203074733E-03_dp
      fit_coef(6, 13) = -5.2171109042265578E-03_dp
      fit_coef(7, 13) = -4.9480341958360466E-04_dp
      fit_coef(8, 13) = 1.2767786392278010E-03_dp
      fit_coef(9, 13) = 8.4783300412091171E-04_dp
      fit_coef(10, 13) = -1.1459374668372458E-03_dp
      fit_coef(11, 13) = 5.7451555318326379E-04_dp
      fit_coef(12, 13) = -1.2838651639065425E-04_dp
      fit_coef(0, 14) = 3.0520606650124064E-02_dp
      fit_coef(1, 14) = 7.8713058722533771E-02_dp
      fit_coef(2, 14) = 1.0070665977979393E-01_dp
      fit_coef(3, 14) = 7.2776550328791670E-02_dp
      fit_coef(4, 14) = 2.1600668660132882E-02_dp
      fit_coef(5, 14) = -1.0555564054529464E-02_dp
      fit_coef(6, 14) = -1.1793835693071972E-02_dp
      fit_coef(7, 14) = -9.6297146505358621E-04_dp
      fit_coef(8, 14) = 2.9732450457490918E-03_dp
      fit_coef(9, 14) = 1.7819794140024070E-03_dp
      fit_coef(10, 14) = -2.5012216511924486E-03_dp
      fit_coef(11, 14) = 1.2501076056785740E-03_dp
      fit_coef(12, 14) = -2.7830976323044436E-04_dp
      fit_coef(0, 15) = 5.4104615066817303E-02_dp
      fit_coef(1, 15) = 1.2949473893963137E-01_dp
      fit_coef(2, 15) = 1.6007411910969396E-01_dp
      fit_coef(3, 15) = 1.1246986879057709E-01_dp
      fit_coef(4, 15) = 3.1687096013120898E-02_dp
      fit_coef(5, 15) = -1.6843831926749306E-02_dp
      fit_coef(6, 15) = -1.7499787344842720E-02_dp
      fit_coef(7, 15) = -1.0003535746501569E-03_dp
      fit_coef(8, 15) = 4.5492268626677436E-03_dp
      fit_coef(9, 15) = 2.3968134835377171E-03_dp
      fit_coef(10, 15) = -3.5437690918133768E-03_dp
      fit_coef(11, 15) = 1.7505580597451528E-03_dp
      fit_coef(12, 15) = -3.8471361244187636E-04_dp
      fit_coef(0, 16) = 8.8713962386069239E-02_dp
      fit_coef(1, 16) = 1.8782619870062861E-01_dp
      fit_coef(2, 16) = 2.1958966435206734E-01_dp
      fit_coef(3, 16) = 1.4951615375763633E-01_dp
      fit_coef(4, 16) = 4.0416134448872314E-02_dp
      fit_coef(5, 16) = -2.3556358508589833E-02_dp
      fit_coef(6, 16) = -2.2705784163903314E-02_dp
      fit_coef(7, 16) = -2.4239170766668074E-04_dp
      fit_coef(8, 16) = 7.8042834525613135E-03_dp
      fit_coef(9, 16) = -4.9113258408041134E-04_dp
      fit_coef(10, 16) = -1.9479282277976330E-03_dp
      fit_coef(11, 16) = 1.2196593353575873E-03_dp
      fit_coef(12, 16) = -3.1694142020982456E-04_dp
      fit_coef(0, 17) = 1.4348763712172333E-01_dp
      fit_coef(1, 17) = 2.5746461675937976E-01_dp
      fit_coef(2, 17) = 2.7495919513020706E-01_dp
      fit_coef(3, 17) = 1.7887972591601542E-01_dp
      fit_coef(4, 17) = 5.1311785540495640E-02_dp
      fit_coef(5, 17) = -2.4925938907640350E-02_dp
      fit_coef(6, 17) = -3.1142444913188724E-02_dp
      fit_coef(7, 17) = -2.1097798874644825E-03_dp
      fit_coef(8, 17) = 6.5586673591109508E-03_dp
      fit_coef(9, 17) = 2.0478515204581191E-02_dp
      fit_coef(10, 17) = -2.9823701948863770E-02_dp
      fit_coef(11, 17) = 1.6872186230723108E-02_dp
      fit_coef(12, 17) = -3.7548105037558685E-03_dp
      fit_coef(0, 18) = 2.3202358256485947E-01_dp
      fit_coef(1, 18) = 3.4341973939523573E-01_dp
      fit_coef(2, 18) = 3.2325120093458798E-01_dp
      fit_coef(3, 18) = 1.9035817858464060E-01_dp
      fit_coef(4, 18) = 5.9225378877332230E-02_dp
      fit_coef(5, 18) = -9.4096717658205004E-03_dp
      fit_coef(6, 18) = -3.1716501677007125E-02_dp
      fit_coef(7, 18) = -3.4342725566089858E-03_dp
      fit_coef(8, 18) = -3.6117207531663649E-02_dp
      fit_coef(9, 18) = 7.0489947329726826E-02_dp
      fit_coef(10, 18) = -3.4005775663369868E-02_dp
      fit_coef(11, 18) = -2.1420043382904190E-03_dp
      fit_coef(12, 18) = 4.1297341956319060E-03_dp
      fit_coef(0, 19) = 3.7428328098642150E-01_dp
      fit_coef(1, 19) = 4.4856859975928659E-01_dp
      fit_coef(2, 19) = 3.6501433131086519E-01_dp
      fit_coef(3, 19) = 1.8076014515005295E-01_dp
      fit_coef(4, 19) = 5.3394895281397445E-02_dp
      fit_coef(5, 19) = 9.7271953924831494E-03_dp
      fit_coef(6, 19) = 6.5024998766183614E-03_dp
      fit_coef(7, 19) = -5.6807667997414525E-02_dp
      fit_coef(8, 19) = 8.5060380392462415E-02_dp
      fit_coef(9, 19) = -2.0020496460057130E-01_dp
      fit_coef(10, 19) = 2.5161185176172424E-01_dp
      fit_coef(11, 19) = -1.4235287380012143E-01_dp
      fit_coef(12, 19) = 3.0493092100125858E-02_dp
      fit_coef(0, 20) = 6.0072363160772879E-01_dp
      fit_coef(1, 20) = 5.7028422329589301E-01_dp
      fit_coef(2, 20) = 4.0013007703573727E-01_dp
      fit_coef(3, 20) = 1.5645364840008563E-01_dp
      fit_coef(4, 20) = 3.5325876772107849E-02_dp
      fit_coef(5, 20) = 2.4573063136551693E-02_dp
      fit_coef(6, 20) = 3.0002223471025873E-02_dp
      fit_coef(7, 20) = -1.9881304587976378E-02_dp
      fit_coef(8, 20) = 1.0115048553542312E-02_dp
      fit_coef(9, 20) = -5.4992791578093891E-02_dp
      fit_coef(10, 20) = 7.6777053071169160E-03_dp
      fit_coef(11, 20) = 4.6072145627107347E-02_dp
      fit_coef(12, 20) = -2.2183512517221421E-02_dp
      fit_coef(0, 21) = 9.6043920592225018E-01_dp
      fit_coef(1, 21) = 6.9944535349071602E-01_dp
      fit_coef(2, 21) = 4.2734572555968520E-01_dp
      fit_coef(3, 21) = 1.2636666009694217E-01_dp
      fit_coef(4, 21) = 1.5064634695905944E-02_dp
      fit_coef(5, 21) = 3.4716238870572014E-02_dp
      fit_coef(6, 21) = 3.1317370170347261E-02_dp
      fit_coef(7, 21) = 5.3593201219119903E-02_dp
      fit_coef(8, 21) = -1.3301585636639335E-01_dp
      fit_coef(9, 21) = 1.9267279137188889E-01_dp
      fit_coef(10, 21) = -2.4206546149930569E-01_dp
      fit_coef(11, 21) = 1.5617787264354296E-01_dp
      fit_coef(12, 21) = -3.7036502375757153E-02_dp
      fit_coef(0, 22) = 1.5411563268918131E+00_dp
      fit_coef(1, 22) = 8.2232573092114225E-01_dp
      fit_coef(2, 22) = 4.4862113308179175E-01_dp
      fit_coef(3, 22) = 9.9831855190677568E-02_dp
      fit_coef(4, 22) = 2.0555253555204072E-04_dp
      fit_coef(5, 22) = 3.7473053552323943E-02_dp
      fit_coef(6, 22) = 4.4849690487304791E-02_dp
      fit_coef(7, 22) = 3.0975652422523699E-02_dp
      fit_coef(8, 22) = -4.3734383490695965E-02_dp
      fit_coef(9, 22) = 4.2410222686447063E-02_dp
      fit_coef(10, 22) = -4.1372075255841384E-02_dp
      fit_coef(11, 22) = 3.9840636595015544E-03_dp
      fit_coef(12, 22) = 6.7429310166998080E-03_dp
      fit_coef(0, 23) = 2.5398045384820884E+00_dp
      fit_coef(1, 23) = 9.2972371375931584E-01_dp
      fit_coef(2, 23) = 4.7467813354768001E-01_dp
      fit_coef(3, 23) = 8.6733603768484421E-02_dp
      fit_coef(4, 23) = -8.2527901024651990E-03_dp
      fit_coef(5, 23) = 3.7146900277764211E-02_dp
      fit_coef(6, 23) = 5.3317290655855602E-02_dp
      fit_coef(7, 23) = 4.9026643386959707E-03_dp
      fit_coef(8, 23) = 3.1627402676112949E-02_dp
      fit_coef(9, 23) = -6.8971662666575367E-02_dp
      fit_coef(10, 23) = 8.6156199337311154E-02_dp
      fit_coef(11, 23) = -6.7821612909471307E-02_dp
      fit_coef(12, 23) = 2.0108151024555578E-02_dp
      fit_coef(0, 24) = 4.6893655829446335E+00_dp
      fit_coef(1, 24) = 1.0584609550444721E+00_dp
      fit_coef(2, 24) = 5.4301054603506849E-01_dp
      fit_coef(3, 24) = 1.0134661869424059E-01_dp
      fit_coef(4, 24) = -1.1773405982811336E-02_dp
      fit_coef(5, 24) = 3.5814305836818772E-02_dp
      fit_coef(6, 24) = 5.3896955561631099E-02_dp
      fit_coef(7, 24) = 1.1794245491552634E-02_dp
      fit_coef(8, 24) = 2.0345857593909816E-02_dp
      fit_coef(9, 24) = -3.4998963970591324E-02_dp
      fit_coef(10, 24) = 4.8444694673587890E-02_dp
      fit_coef(11, 24) = -3.2817041946018727E-02_dp
      fit_coef(12, 24) = 6.8092769866170090E-03_dp

   END SUBROUTINE get_fit_coef_k12_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k12_mi(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.1533853671352736E+03_dp
      fit_coef(0, 1) = 6.6326907649095958E-04_dp
      fit_coef(1, 1) = 1.8894159031673920E-03_dp
      fit_coef(2, 1) = 2.7089355355061177E-03_dp
      fit_coef(3, 1) = 2.3224757764888754E-03_dp
      fit_coef(4, 1) = 1.0863117203975862E-03_dp
      fit_coef(5, 1) = 4.5151540703808010E-05_dp
      fit_coef(6, 1) = -2.3515149052107272E-04_dp
      fit_coef(7, 1) = -8.2929268974049080E-05_dp
      fit_coef(8, 1) = 1.1614919956202285E-04_dp
      fit_coef(9, 1) = -2.5284724218147597E-05_dp
      fit_coef(10, 1) = -1.3858258705778020E-05_dp
      fit_coef(11, 1) = 8.0009405004549032E-06_dp
      fit_coef(12, 1) = -1.0748701295750969E-06_dp
      fit_coef(0, 2) = 3.8227101929591994E-03_dp
      fit_coef(1, 2) = 1.0364280264194323E-02_dp
      fit_coef(2, 2) = 1.4590446568112637E-02_dp
      fit_coef(3, 2) = 1.2309386146930817E-02_dp
      fit_coef(4, 2) = 5.6468571531384185E-03_dp
      fit_coef(5, 2) = 1.8293346115620324E-04_dp
      fit_coef(6, 2) = -1.2260374384726667E-03_dp
      fit_coef(7, 2) = -4.2679063675895478E-04_dp
      fit_coef(8, 2) = 6.4562440386666691E-04_dp
      fit_coef(9, 2) = -2.0425912250358465E-04_dp
      fit_coef(10, 2) = -1.0718250746829363E-05_dp
      fit_coef(11, 2) = 1.3450401493497087E-05_dp
      fit_coef(12, 2) = 3.6903867124824597E-08_dp
      fit_coef(0, 3) = 1.1208738658104969E-02_dp
      fit_coef(1, 3) = 2.7433467450107527E-02_dp
      fit_coef(2, 3) = 3.7122781203738041E-02_dp
      fit_coef(3, 3) = 3.0505000476205493E-02_dp
      fit_coef(4, 3) = 1.3613669155986628E-02_dp
      fit_coef(5, 3) = 2.0805962290796248E-04_dp
      fit_coef(6, 3) = -3.0074089272987834E-03_dp
      fit_coef(7, 3) = -1.0066712315667174E-03_dp
      fit_coef(8, 3) = 1.8121892538119837E-03_dp
      fit_coef(9, 3) = -8.6670466868856387E-04_dp
      fit_coef(10, 3) = 2.3031588351142570E-04_dp
      fit_coef(11, 3) = -6.1244719299360145E-05_dp
      fit_coef(12, 3) = 1.5021644001300220E-05_dp
      fit_coef(0, 4) = 2.7060109838528964E-02_dp
      fit_coef(1, 4) = 5.7181276021668774E-02_dp
      fit_coef(2, 4) = 7.2163197996693393E-02_dp
      fit_coef(3, 4) = 5.6401492360043291E-02_dp
      fit_coef(4, 4) = 2.4397958649191483E-02_dp
      fit_coef(5, 4) = 4.1542953141707141E-04_dp
      fit_coef(6, 4) = -5.4082092785117776E-03_dp
      fit_coef(7, 4) = -1.9248439551360079E-03_dp
      fit_coef(8, 4) = 2.9200980290588186E-03_dp
      fit_coef(9, 4) = -9.7952619465276037E-04_dp
      fit_coef(10, 4) = 1.6789215463763337E-04_dp
      fit_coef(11, 4) = -1.4147604636971736E-04_dp
      fit_coef(12, 4) = 6.0667969636576076E-05_dp
      fit_coef(0, 5) = 6.0162354447464560E-02_dp
      fit_coef(1, 5) = 1.0720429058942169E-01_dp
      fit_coef(2, 5) = 1.2343759565688829E-01_dp
      fit_coef(3, 5) = 8.8531585067086430E-02_dp
      fit_coef(4, 5) = 3.5364287351100372E-02_dp
      fit_coef(5, 5) = 6.0629507679625032E-04_dp
      fit_coef(6, 5) = -6.6115157196993245E-03_dp
      fit_coef(7, 5) = -3.1579068142273776E-03_dp
      fit_coef(8, 5) = 5.0366482632659849E-03_dp
      fit_coef(9, 5) = -5.0447339830031095E-03_dp
      fit_coef(10, 5) = 4.5014753449279963E-03_dp
      fit_coef(11, 5) = -2.3546672147255598E-03_dp
      fit_coef(12, 5) = 5.0845127538551175E-04_dp
      fit_coef(0, 6) = 1.2715764725238982E-01_dp
      fit_coef(1, 6) = 1.8881057094422343E-01_dp
      fit_coef(2, 6) = 1.9672779228400020E-01_dp
      fit_coef(3, 6) = 1.2577665875739202E-01_dp
      fit_coef(4, 6) = 4.3061470368101699E-02_dp
      fit_coef(5, 6) = -9.4453422287111893E-04_dp
      fit_coef(6, 6) = -6.1901712508776350E-03_dp
      fit_coef(7, 6) = -2.9123796915168040E-03_dp
      fit_coef(8, 6) = 7.5406628009243720E-03_dp
      fit_coef(9, 6) = -1.0878759759659724E-02_dp
      fit_coef(10, 6) = 9.0511431231051973E-03_dp
      fit_coef(11, 6) = -3.9841133010950968E-03_dp
      fit_coef(12, 6) = 7.4476527248770687E-04_dp
      fit_coef(0, 7) = 2.5849300626981397E-01_dp
      fit_coef(1, 7) = 3.1665344505552956E-01_dp
      fit_coef(2, 7) = 2.9840462942608115E-01_dp
      fit_coef(3, 7) = 1.6655149775226027E-01_dp
      fit_coef(4, 7) = 4.4635943287943819E-02_dp
      fit_coef(5, 7) = -5.0840338248153878E-03_dp
      fit_coef(6, 7) = -4.1935958573150935E-03_dp
      fit_coef(7, 7) = -1.2238069004985115E-03_dp
      fit_coef(8, 7) = 1.0355993093508651E-02_dp
      fit_coef(9, 7) = -1.8364631667690871E-02_dp
      fit_coef(10, 7) = 1.5416014403986910E-02_dp
      fit_coef(11, 7) = -6.8470052029501608E-03_dp
      fit_coef(12, 7) = 1.3249810031683059E-03_dp
      fit_coef(0, 8) = 5.0864125871117427E-01_dp
      fit_coef(1, 8) = 5.0748696881076827E-01_dp
      fit_coef(2, 8) = 4.3378556875785329E-01_dp
      fit_coef(3, 8) = 2.0761241946020209E-01_dp
      fit_coef(4, 8) = 3.7494223843576979E-02_dp
      fit_coef(5, 8) = -1.1199453066603072E-02_dp
      fit_coef(6, 8) = 2.6211861457654248E-04_dp
      fit_coef(7, 8) = 1.0181856045809953E-03_dp
      fit_coef(8, 8) = 1.3617506208780007E-02_dp
      fit_coef(9, 8) = -2.9083172474373836E-02_dp
      fit_coef(10, 8) = 2.6093940387631132E-02_dp
      fit_coef(11, 8) = -1.2187538547457615E-02_dp
      fit_coef(12, 8) = 2.4442839760494672E-03_dp
      fit_coef(0, 9) = 9.7351269086421022E-01_dp
      fit_coef(1, 9) = 7.7682804453441712E-01_dp
      fit_coef(2, 9) = 6.0537341804741274E-01_dp
      fit_coef(3, 9) = 2.4417103093525630E-01_dp
      fit_coef(4, 9) = 1.9643095133629673E-02_dp
      fit_coef(5, 9) = -1.7552027242727909E-02_dp
      fit_coef(6, 9) = 8.2268681071634708E-03_dp
      fit_coef(7, 9) = 3.3800170464644732E-03_dp
      fit_coef(8, 9) = 1.5951520448433745E-02_dp
      fit_coef(9, 9) = -4.2082858702406162E-02_dp
      fit_coef(10, 9) = 4.0832752772189854E-02_dp
      fit_coef(11, 9) = -1.9914993594640506E-02_dp
      fit_coef(12, 9) = 4.0828025408897019E-03_dp
      fit_coef(0, 10) = 1.8218150409063079E+00_dp
      fit_coef(1, 10) = 1.1328940612144838E+00_dp
      fit_coef(2, 10) = 8.1145165005204600E-01_dp
      fit_coef(3, 10) = 2.7169894714361381E-01_dp
      fit_coef(4, 10) = -9.2391365403536683E-03_dp
      fit_coef(5, 10) = -2.1969808519198136E-02_dp
      fit_coef(6, 10) = 2.0246181533686430E-02_dp
      fit_coef(7, 10) = 5.2443671348299629E-03_dp
      fit_coef(8, 10) = 1.6294437997198600E-02_dp
      fit_coef(9, 10) = -5.6113566142827496E-02_dp
      fit_coef(10, 10) = 5.9108716717708169E-02_dp
      fit_coef(11, 10) = -3.0017787615102216E-02_dp
      fit_coef(12, 10) = 6.2851355640176588E-03_dp
      fit_coef(0, 11) = 3.3638892281286350E+00_dp
      fit_coef(1, 11) = 1.5701436592256290E+00_dp
      fit_coef(2, 11) = 1.0468794820440899E+00_dp
      fit_coef(3, 11) = 2.8932244841622851E-01_dp
      fit_coef(4, 11) = -4.7147330695069850E-02_dp
      fit_coef(5, 11) = -2.3231305463748575E-02_dp
      fit_coef(6, 11) = 3.5814349205865105E-02_dp
      fit_coef(7, 11) = 6.1681105190574941E-03_dp
      fit_coef(8, 11) = 1.4380388978420344E-02_dp
      fit_coef(9, 11) = -7.0072767054186372E-02_dp
      fit_coef(10, 11) = 7.9925293980236525E-02_dp
      fit_coef(11, 11) = -4.2106160715496485E-02_dp
      fit_coef(12, 11) = 8.9934277746186089E-03_dp
      fit_coef(0, 12) = 6.2958157177735696E+00_dp
      fit_coef(1, 12) = 2.0736278555613925E+00_dp
      fit_coef(2, 12) = 1.3100027104859728E+00_dp
      fit_coef(3, 12) = 3.0309841684283917E-01_dp
      fit_coef(4, 12) = -9.0874734045887973E-02_dp
      fit_coef(5, 12) = -2.2511990033847892E-02_dp
      fit_coef(6, 12) = 5.3762327238007911E-02_dp
      fit_coef(7, 12) = 6.4301275560599072E-03_dp
      fit_coef(8, 12) = 1.1130561807738493E-02_dp
      fit_coef(9, 12) = -8.4343774419336590E-02_dp
      fit_coef(10, 12) = 1.0282039809053546E-01_dp
      fit_coef(11, 12) = -5.5735942843285882E-02_dp
      fit_coef(12, 12) = 1.2089953955292771E-02_dp
      fit_coef(0, 13) = 1.7361176995326563E-03_dp
      fit_coef(1, 13) = 4.8950464700643088E-03_dp
      fit_coef(2, 13) = 6.9887889000087690E-03_dp
      fit_coef(3, 13) = 5.9681592147681978E-03_dp
      fit_coef(4, 13) = 2.7797694969593904E-03_dp
      fit_coef(5, 13) = 1.0974079237570183E-04_dp
      fit_coef(6, 13) = -6.0276261358583726E-04_dp
      fit_coef(7, 13) = -2.1159566479273911E-04_dp
      fit_coef(8, 13) = 3.0294248744461766E-04_dp
      fit_coef(9, 13) = -7.4399992147188479E-05_dp
      fit_coef(10, 13) = -2.7468045558828552E-05_dp
      fit_coef(11, 13) = 1.6865141251068427E-05_dp
      fit_coef(12, 13) = -2.0377510281004492E-06_dp
      fit_coef(0, 14) = 4.8355092096145566E-03_dp
      fit_coef(1, 14) = 1.2333932959289709E-02_dp
      fit_coef(2, 14) = 1.6972825457174837E-02_dp
      fit_coef(3, 14) = 1.4060909191254103E-02_dp
      fit_coef(4, 14) = 6.2955683700635671E-03_dp
      fit_coef(5, 14) = 1.2177301196491723E-04_dp
      fit_coef(6, 14) = -1.3677107176344808E-03_dp
      fit_coef(7, 14) = -4.6095049437820380E-04_dp
      fit_coef(8, 14) = 7.6937933826695876E-04_dp
      fit_coef(9, 14) = -3.3164200459717756E-04_dp
      fit_coef(10, 14) = 8.1724364865420408E-05_dp
      fit_coef(11, 14) = -2.9312146170772606E-05_dp
      fit_coef(12, 14) = 8.9970999683901846E-06_dp
      fit_coef(0, 15) = 1.0612378818881555E-02_dp
      fit_coef(1, 15) = 2.2459336506101473E-02_dp
      fit_coef(2, 15) = 2.8371129792457078E-02_dp
      fit_coef(3, 15) = 2.2261852601842439E-02_dp
      fit_coef(4, 15) = 9.5999731643326893E-03_dp
      fit_coef(5, 15) = -1.8558563529812452E-05_dp
      fit_coef(6, 15) = -2.2357421491240670E-03_dp
      fit_coef(7, 15) = -7.4423534393391617E-04_dp
      fit_coef(8, 15) = 1.5057122702193376E-03_dp
      fit_coef(9, 15) = -7.0297537707829122E-04_dp
      fit_coef(10, 15) = 7.1959722147797053E-05_dp
      fit_coef(11, 15) = 4.6639800520364974E-05_dp
      fit_coef(12, 15) = -1.2660883153173911E-05_dp
      fit_coef(0, 16) = 2.2502651094608313E-02_dp
      fit_coef(1, 16) = 3.8302915711040345E-02_dp
      fit_coef(2, 16) = 4.2327776536488830E-02_dp
      fit_coef(3, 16) = 2.9242971747304312E-02_dp
      fit_coef(4, 16) = 1.1483493368870591E-02_dp
      fit_coef(5, 16) = 4.4372686263633332E-04_dp
      fit_coef(6, 16) = -2.1409365324177557E-03_dp
      fit_coef(7, 16) = -1.1499475056719786E-03_dp
      fit_coef(8, 16) = 9.9305143157658372E-04_dp
      fit_coef(9, 16) = -8.6357250324013148E-04_dp
      fit_coef(10, 16) = 1.4327807191970657E-03_dp
      fit_coef(11, 16) = -1.0670863396023819E-03_dp
      fit_coef(12, 16) = 2.8011891218022886E-04_dp
      fit_coef(0, 17) = 4.6379852772014826E-02_dp
      fit_coef(1, 17) = 6.3632440293261802E-02_dp
      fit_coef(2, 17) = 6.1214298124389097E-02_dp
      fit_coef(3, 17) = 3.4856680478514791E-02_dp
      fit_coef(4, 17) = 9.8511055688256591E-03_dp
      fit_coef(5, 17) = -4.1408424159966774E-04_dp
      fit_coef(6, 17) = -2.7116403636660304E-04_dp
      fit_coef(7, 17) = -8.3345412989501344E-04_dp
      fit_coef(8, 17) = 2.8388169463227923E-03_dp
      fit_coef(9, 17) = -6.0896330652201345E-03_dp
      fit_coef(10, 17) = 5.4431369422386086E-03_dp
      fit_coef(11, 17) = -2.2497263945134662E-03_dp
      fit_coef(12, 17) = 3.6399313731776849E-04_dp
      fit_coef(0, 18) = 9.2521202823536430E-02_dp
      fit_coef(1, 18) = 1.0203973195050642E-01_dp
      fit_coef(2, 18) = 8.6453931343429755E-02_dp
      fit_coef(3, 18) = 3.9391412232875435E-02_dp
      fit_coef(4, 18) = 5.0815288567403350E-03_dp
      fit_coef(5, 18) = -2.8315021323427666E-03_dp
      fit_coef(6, 18) = 1.1044933352897194E-03_dp
      fit_coef(7, 18) = 1.2026713128526101E-03_dp
      fit_coef(8, 18) = 2.3687414607481347E-03_dp
      fit_coef(9, 18) = -5.9882091197473600E-03_dp
      fit_coef(10, 18) = 4.5247728538859366E-03_dp
      fit_coef(11, 18) = -1.7071153745195653E-03_dp
      fit_coef(12, 18) = 2.9758468692744693E-04_dp
      fit_coef(0, 19) = 1.7895123701529778E-01_dp
      fit_coef(1, 19) = 1.5649170004661517E-01_dp
      fit_coef(2, 19) = 1.1781626034461101E-01_dp
      fit_coef(3, 19) = 4.1616478405439042E-02_dp
      fit_coef(4, 19) = -2.3717162520277809E-03_dp
      fit_coef(5, 19) = -5.3386440882515086E-03_dp
      fit_coef(6, 19) = 3.0469549759779353E-03_dp
      fit_coef(7, 19) = 2.0334444707206275E-03_dp
      fit_coef(8, 19) = 3.2479855073872382E-03_dp
      fit_coef(9, 19) = -9.1986286972097365E-03_dp
      fit_coef(10, 19) = 8.5242332274741671E-03_dp
      fit_coef(11, 19) = -4.1329263555985006E-03_dp
      fit_coef(12, 19) = 8.6616145809027139E-04_dp
      fit_coef(0, 20) = 3.3684611079928783E-01_dp
      fit_coef(1, 20) = 2.2786610373174537E-01_dp
      fit_coef(2, 20) = 1.5339071016684352E-01_dp
      fit_coef(3, 20) = 3.9677854839472344E-02_dp
      fit_coef(4, 20) = -1.2264564209904332E-02_dp
      fit_coef(5, 20) = -6.5927521643919379E-03_dp
      fit_coef(6, 20) = 6.0512561434243684E-03_dp
      fit_coef(7, 20) = 2.3902542529467478E-03_dp
      fit_coef(8, 20) = 3.0195138462315793E-03_dp
      fit_coef(9, 20) = -1.2050135993152236E-02_dp
      fit_coef(10, 20) = 1.2762387029635221E-02_dp
      fit_coef(11, 20) = -6.5239286773536758E-03_dp
      fit_coef(12, 20) = 1.3694468177831019E-03_dp
      fit_coef(0, 21) = 6.2004160201378677E-01_dp
      fit_coef(1, 21) = 3.1231452836704038E-01_dp
      fit_coef(2, 21) = 1.8949289978076236E-01_dp
      fit_coef(3, 21) = 3.2623505491414392E-02_dp
      fit_coef(4, 21) = -2.3505999675060181E-02_dp
      fit_coef(5, 21) = -5.7260158711334771E-03_dp
      fit_coef(6, 21) = 9.9807237729862951E-03_dp
      fit_coef(7, 21) = 2.2222797839610000E-03_dp
      fit_coef(8, 21) = 1.4664130053888486E-03_dp
      fit_coef(9, 21) = -1.3740779304315174E-02_dp
      fit_coef(10, 21) = 1.6640549191351679E-02_dp
      fit_coef(11, 21) = -8.9416960852337218E-03_dp
      fit_coef(12, 21) = 1.9187092435646925E-03_dp
      fit_coef(0, 22) = 1.1255181435117496E+00_dp
      fit_coef(1, 22) = 3.9889498776431231E-01_dp
      fit_coef(2, 22) = 2.2171584629075095E-01_dp
      fit_coef(3, 22) = 2.2225826529528901E-02_dp
      fit_coef(4, 22) = -3.3907556115484502E-02_dp
      fit_coef(5, 22) = -2.8764770890914323E-03_dp
      fit_coef(6, 22) = 1.3968131777336038E-02_dp
      fit_coef(7, 22) = 1.4214928484319333E-03_dp
      fit_coef(8, 22) = -8.3453566281261021E-04_dp
      fit_coef(9, 22) = -1.4109958929870586E-02_dp
      fit_coef(10, 22) = 1.9726438913848971E-02_dp
      fit_coef(11, 22) = -1.1193142315748616E-02_dp
      fit_coef(12, 22) = 2.4755702682764356E-03_dp
      fit_coef(0, 23) = 2.0598610404318460E+00_dp
      fit_coef(1, 23) = 4.7243372474599793E-01_dp
      fit_coef(2, 23) = 2.4849210820527318E-01_dp
      fit_coef(3, 23) = 1.4006813172765054E-02_dp
      fit_coef(4, 23) = -4.1308949825231950E-02_dp
      fit_coef(5, 23) = 1.5690842161415170E-04_dp
      fit_coef(6, 23) = 1.6942861921432139E-02_dp
      fit_coef(7, 23) = 4.7152452086483685E-04_dp
      fit_coef(8, 23) = -2.8341187428894861E-03_dp
      fit_coef(9, 23) = -1.3846753524951650E-02_dp
      fit_coef(10, 23) = 2.1790170949455349E-02_dp
      fit_coef(11, 23) = -1.2886556590999076E-02_dp
      fit_coef(12, 23) = 2.9162146561646620E-03_dp
      fit_coef(0, 24) = 4.1434250198174896E+00_dp
      fit_coef(1, 24) = 5.3740865227175982E-01_dp
      fit_coef(2, 24) = 2.8197846608440591E-01_dp
      fit_coef(3, 24) = 1.6092456027878326E-02_dp
      fit_coef(4, 24) = -4.6040710751708042E-02_dp
      fit_coef(5, 24) = 6.5042394570127837E-04_dp
      fit_coef(6, 24) = 1.8932659149369899E-02_dp
      fit_coef(7, 24) = 2.4082473716184700E-04_dp
      fit_coef(8, 24) = -3.3484787514316723E-03_dp
      fit_coef(9, 24) = -1.5284105573369426E-02_dp
      fit_coef(10, 24) = 2.4407280609082902E-02_dp
      fit_coef(11, 24) = -1.4513713351654377E-02_dp
      fit_coef(12, 24) = 3.2970008560632581E-03_dp

   END SUBROUTINE get_fit_coef_k12_mi

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k12_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.0890000000000000E+05_dp
      fit_coef(0, 1) = 2.3129808933247259E-05_dp
      fit_coef(1, 1) = 1.3229791040971861E-04_dp
      fit_coef(2, 1) = 8.9078224043542527E-04_dp
      fit_coef(3, 1) = 3.5906725339575205E-04_dp
      fit_coef(4, 1) = 1.4513983588735334E-02_dp
      fit_coef(5, 1) = -5.2233430638158315E-02_dp
      fit_coef(6, 1) = 1.6441632828707831E-01_dp
      fit_coef(7, 1) = -3.6058840741824150E-01_dp
      fit_coef(8, 1) = 5.3428676974850686E-01_dp
      fit_coef(9, 1) = -5.2731672647478189E-01_dp
      fit_coef(10, 1) = 3.3265094311193522E-01_dp
      fit_coef(11, 1) = -1.2145741283688967E-01_dp
      fit_coef(12, 1) = 1.9537487536713263E-02_dp
      fit_coef(0, 2) = 2.3455247909237882E-04_dp
      fit_coef(1, 2) = 8.2160398869974953E-04_dp
      fit_coef(2, 2) = 5.3319821419640272E-03_dp
      fit_coef(3, 2) = 6.4575814983609142E-04_dp
      fit_coef(4, 2) = 8.7822372114111175E-02_dp
      fit_coef(5, 2) = -3.2631392232557216E-01_dp
      fit_coef(6, 2) = 1.0208260875409669E+00_dp
      fit_coef(7, 2) = -2.2300026314422441E+00_dp
      fit_coef(8, 2) = 3.3057374126568657E+00_dp
      fit_coef(9, 2) = -3.2713127430468321E+00_dp
      fit_coef(10, 2) = 2.0709853451512608E+00_dp
      fit_coef(11, 2) = -7.5901118153106462E-01_dp
      fit_coef(12, 2) = 1.2254511188136619E-01_dp
      fit_coef(0, 3) = 1.1944107046876868E-03_dp
      fit_coef(1, 3) = 2.8490802943664383E-03_dp
      fit_coef(2, 3) = 1.7507904217475573E-02_dp
      fit_coef(3, 3) = -6.4522716917725272E-03_dp
      fit_coef(4, 3) = 2.8224130041588374E-01_dp
      fit_coef(5, 3) = -1.1757205798515102E+00_dp
      fit_coef(6, 3) = 3.7822941249117266E+00_dp
      fit_coef(7, 3) = -8.2843383024642065E+00_dp
      fit_coef(8, 3) = 1.2238481082417001E+01_dp
      fit_coef(9, 3) = -1.2055852436649243E+01_dp
      fit_coef(10, 3) = 7.6007782878436343E+00_dp
      fit_coef(11, 3) = -2.7765219925860354E+00_dp
      fit_coef(12, 3) = 4.4718773429739378E-01_dp
      fit_coef(0, 4) = 4.5120683484607405E-03_dp
      fit_coef(1, 4) = 7.9095288831978021E-03_dp
      fit_coef(2, 4) = 4.6697868261753268E-02_dp
      fit_coef(3, 4) = -3.3544522740173804E-02_dp
      fit_coef(4, 4) = 7.2425137774094894E-01_dp
      fit_coef(5, 4) = -3.2633003302313521E+00_dp
      fit_coef(6, 4) = 1.0450673822567490E+01_dp
      fit_coef(7, 4) = -2.2594031368355793E+01_dp
      fit_coef(8, 4) = 3.3109462208134374E+01_dp
      fit_coef(9, 4) = -3.2506292757603717E+01_dp
      fit_coef(10, 4) = 2.0487989626895530E+01_dp
      fit_coef(11, 4) = -7.4948539272232253E+00_dp
      fit_coef(12, 4) = 1.2099196036738793E+00_dp
      fit_coef(0, 5) = 1.4248952254027877E-02_dp
      fit_coef(1, 5) = 1.9074798942523639E-02_dp
      fit_coef(2, 5) = 1.0943167990343619E-01_dp
      fit_coef(3, 5) = -1.0671670298403642E-01_dp
      fit_coef(4, 5) = 1.6526275275600326E+00_dp
      fit_coef(5, 5) = -7.8775178049382024E+00_dp
      fit_coef(6, 5) = 2.5305646502994179E+01_dp
      fit_coef(7, 5) = -5.4666053461159009E+01_dp
      fit_coef(8, 5) = 8.0167059237163343E+01_dp
      fit_coef(9, 5) = -7.8825322292091528E+01_dp
      fit_coef(10, 5) = 4.9764390489009351E+01_dp
      fit_coef(11, 5) = -1.8233568839424560E+01_dp
      fit_coef(12, 5) = 2.9478228744193333E+00_dp
      fit_coef(0, 6) = 3.9801596054429601E-02_dp
      fit_coef(1, 6) = 4.1447021761176796E-02_dp
      fit_coef(2, 6) = 2.3270073234920230E-01_dp
      fit_coef(3, 6) = -2.7292288542272475E-01_dp
      fit_coef(4, 6) = 3.4506136511101664E+00_dp
      fit_coef(5, 6) = -1.7124854042561445E+01_dp
      fit_coef(6, 6) = 5.5335626634127287E+01_dp
      fit_coef(7, 6) = -1.1966893675071402E+02_dp
      fit_coef(8, 6) = 1.7554731067692325E+02_dp
      fit_coef(9, 6) = -1.7262672563793396E+02_dp
      fit_coef(10, 6) = 1.0899761858660371E+02_dp
      fit_coef(11, 6) = -3.9946406626328411E+01_dp
      fit_coef(12, 6) = 6.4605650754876507E+00_dp
      fit_coef(0, 7) = 1.0153848301405891E-01_dp
      fit_coef(1, 7) = 8.2744575115798746E-02_dp
      fit_coef(2, 7) = 4.5675281550800451E-01_dp
      fit_coef(3, 7) = -6.0755960983591895E-01_dp
      fit_coef(4, 7) = 6.6840425945230333E+00_dp
      fit_coef(5, 7) = -3.4181153764390729E+01_dp
      fit_coef(6, 7) = 1.1120663926638191E+02_dp
      fit_coef(7, 7) = -2.4092641084144972E+02_dp
      fit_coef(8, 7) = 3.5355143552413938E+02_dp
      fit_coef(9, 7) = -3.4765931274190069E+02_dp
      fit_coef(10, 7) = 2.1950666813842750E+02_dp
      fit_coef(11, 7) = -8.0451830347244339E+01_dp
      fit_coef(12, 7) = 1.3013672905602384E+01_dp
      fit_coef(0, 8) = 2.4141648240665145E-01_dp
      fit_coef(1, 8) = 1.5335953037887173E-01_dp
      fit_coef(2, 8) = 8.3504832050626976E-01_dp
      fit_coef(3, 8) = -1.2178157194978934E+00_dp
      fit_coef(4, 8) = 1.2099585014816819E+01_dp
      fit_coef(5, 8) = -6.3313617826158982E+01_dp
      fit_coef(6, 8) = 2.0741722770815861E+02_dp
      fit_coef(7, 8) = -4.5029790616194964E+02_dp
      fit_coef(8, 8) = 6.6110515192856462E+02_dp
      fit_coef(9, 8) = -6.5007509234254940E+02_dp
      fit_coef(10, 8) = 4.1040648800886731E+02_dp
      fit_coef(11, 8) = -1.5041026682712433E+02_dp
      fit_coef(12, 8) = 2.4330155922976292E+01_dp
      fit_coef(0, 9) = 5.4247999033308114E-01_dp
      fit_coef(1, 9) = 2.6499445034705121E-01_dp
      fit_coef(2, 9) = 1.4268709589831525E+00_dp
      fit_coef(3, 9) = -2.2313886617958429E+00_dp
      fit_coef(4, 9) = 2.0518890229005105E+01_dp
      fit_coef(5, 9) = -1.0930947078396949E+02_dp
      fit_coef(6, 9) = 3.6045833696584418E+02_dp
      fit_coef(7, 9) = -7.8423527312908436E+02_dp
      fit_coef(8, 9) = 1.1520122421721351E+03_dp
      fit_coef(9, 9) = -1.1328174992310501E+03_dp
      fit_coef(10, 9) = 7.1509445979181294E+02_dp
      fit_coef(11, 9) = -2.6204807849728383E+02_dp
      fit_coef(12, 9) = 4.2385855620433802E+01_dp
      fit_coef(0, 10) = 1.1649235984334310E+00_dp
      fit_coef(1, 10) = 4.2669860434201973E-01_dp
      fit_coef(2, 10) = 2.2769218254498678E+00_dp
      fit_coef(3, 10) = -3.7556598199553419E+00_dp
      fit_coef(4, 10) = 3.2554151177465044E+01_dp
      fit_coef(5, 10) = -1.7584595871324311E+02_dp
      fit_coef(6, 10) = 5.8329724780519552E+02_dp
      fit_coef(7, 10) = -1.2717006032060506E+03_dp
      fit_coef(8, 10) = 1.8691949694429402E+03_dp
      fit_coef(9, 10) = -1.8381765181945018E+03_dp
      fit_coef(10, 10) = 1.1602510916037054E+03_dp
      fit_coef(11, 10) = -4.2512958735964088E+02_dp
      fit_coef(12, 10) = 6.8758003199532453E+01_dp
      fit_coef(0, 11) = 2.4212426991212057E+00_dp
      fit_coef(1, 11) = 6.3843606877537262E-01_dp
      fit_coef(2, 11) = 3.3832944732267647E+00_dp
      fit_coef(3, 11) = -5.8032061767073824E+00_dp
      fit_coef(4, 11) = 4.8167947114080192E+01_dp
      fit_coef(5, 11) = -2.6286759432487781E+02_dp
      fit_coef(6, 11) = 8.7620572915548814E+02_dp
      fit_coef(7, 11) = -1.9137694452215110E+03_dp
      fit_coef(8, 11) = 2.8145274374678247E+03_dp
      fit_coef(9, 11) = -2.7680919699076530E+03_dp
      fit_coef(10, 11) = 1.7471091763340671E+03_dp
      fit_coef(11, 11) = -6.4010183684923209E+02_dp
      fit_coef(12, 11) = 1.0351773205889742E+02_dp
      fit_coef(0, 12) = 5.0200579894610868E+00_dp
      fit_coef(1, 12) = 8.8885387561015605E-01_dp
      fit_coef(2, 12) = 4.6883210461508416E+00_dp
      fit_coef(3, 12) = -8.2512193578220483E+00_dp
      fit_coef(4, 12) = 6.6560605008919623E+01_dp
      fit_coef(5, 12) = -3.6572207879590201E+02_dp
      fit_coef(6, 12) = 1.2232199700982205E+03_dp
      fit_coef(7, 12) = -2.6752268847851810E+03_dp
      fit_coef(8, 12) = 3.9360821093174841E+03_dp
      fit_coef(9, 12) = -3.8714869910206562E+03_dp
      fit_coef(10, 12) = 2.4434517488097472E+03_dp
      fit_coef(11, 12) = -8.9517007084930697E+02_dp
      fit_coef(12, 12) = 1.4475878519184542E+02_dp
      fit_coef(0, 13) = 7.1480283160792286E-05_dp
      fit_coef(1, 13) = 3.4522475460442413E-04_dp
      fit_coef(2, 13) = 2.3289358018025357E-03_dp
      fit_coef(3, 13) = 8.6564456187697154E-04_dp
      fit_coef(4, 13) = 3.9444759289228873E-02_dp
      fit_coef(5, 13) = -1.4607809825480644E-01_dp
      fit_coef(6, 13) = 4.5801849703119135E-01_dp
      fit_coef(7, 13) = -9.9551749039632609E-01_dp
      fit_coef(8, 13) = 1.4633991701280600E+00_dp
      fit_coef(9, 13) = -1.4351206048684466E+00_dp
      fit_coef(10, 13) = 9.0075138053758774E-01_dp
      fit_coef(11, 13) = -3.2755040336927682E-01_dp
      fit_coef(12, 13) = 5.2516019203549216E-02_dp
      fit_coef(0, 14) = 4.3512255060467424E-04_dp
      fit_coef(1, 14) = 1.1619032009340215E-03_dp
      fit_coef(2, 14) = 7.2125439844170202E-03_dp
      fit_coef(3, 14) = -1.5997165186082474E-03_dp
      fit_coef(4, 14) = 1.1580962215615681E-01_dp
      fit_coef(5, 14) = -4.5588347118168826E-01_dp
      fit_coef(6, 14) = 1.4668513671508041E+00_dp
      fit_coef(7, 14) = -3.2658979984403209E+00_dp
      fit_coef(8, 14) = 4.9101588767499029E+00_dp
      fit_coef(9, 14) = -4.9103950358935151E+00_dp
      fit_coef(10, 14) = 3.1335239775801118E+00_dp
      fit_coef(11, 14) = -1.1555673477193260E+00_dp
      fit_coef(12, 14) = 1.8749662399467892E-01_dp
      fit_coef(0, 15) = 1.7310205801983041E-03_dp
      fit_coef(1, 15) = 3.1669171624321819E-03_dp
      fit_coef(2, 15) = 1.8644195870254832E-02_dp
      fit_coef(3, 15) = -1.4397262238866691E-02_dp
      fit_coef(4, 15) = 2.9169977531123992E-01_dp
      fit_coef(5, 15) = -1.3313481914883316E+00_dp
      fit_coef(6, 15) = 4.2756251478512226E+00_dp
      fit_coef(7, 15) = -9.1958972503690166E+00_dp
      fit_coef(8, 15) = 1.3362193435862070E+01_dp
      fit_coef(9, 15) = -1.3004184796302695E+01_dp
      fit_coef(10, 15) = 8.1326698622082620E+00_dp
      fit_coef(11, 15) = -2.9559537831500120E+00_dp
      fit_coef(12, 15) = 4.7473220658578025E-01_dp
      fit_coef(0, 16) = 5.5360416394316771E-03_dp
      fit_coef(1, 16) = 7.4544925165475080E-03_dp
      fit_coef(2, 16) = 4.2442503620852397E-02_dp
      fit_coef(3, 16) = -4.4215606043728234E-02_dp
      fit_coef(4, 16) = 6.3386466456182333E-01_dp
      fit_coef(5, 16) = -3.0715020425513488E+00_dp
      fit_coef(6, 16) = 9.8361982985015253E+00_dp
      fit_coef(7, 16) = -2.1180667612485266E+01_dp
      fit_coef(8, 16) = 3.1043088235838166E+01_dp
      fit_coef(9, 16) = -3.0559569250108336E+01_dp
      fit_coef(10, 16) = 1.9330505471141983E+01_dp
      fit_coef(11, 16) = -7.0975069725987057E+00_dp
      fit_coef(12, 16) = 1.1497097136190815E+00_dp
      fit_coef(0, 17) = 1.5414916387788672E-02_dp
      fit_coef(1, 17) = 1.5715678851517083E-02_dp
      fit_coef(2, 17) = 8.7468870899845605E-02_dp
      fit_coef(3, 17) = -1.0975140358097349E-01_dp
      fit_coef(4, 17) = 1.2852966018847800E+00_dp
      fit_coef(5, 17) = -6.4984894179776092E+00_dp
      fit_coef(6, 17) = 2.1000940295849389E+01_dp
      fit_coef(7, 17) = -4.5396783297216288E+01_dp
      fit_coef(8, 17) = 6.6608738269027413E+01_dp
      fit_coef(9, 17) = -6.5525352583359691E+01_dp
      fit_coef(10, 17) = 4.1386945120960981E+01_dp
      fit_coef(11, 17) = -1.5171963717870183E+01_dp
      fit_coef(12, 17) = 2.4543585537353447E+00_dp
      fit_coef(0, 18) = 3.8916898079420920E-02_dp
      fit_coef(1, 18) = 3.0312084028360865E-02_dp
      fit_coef(2, 18) = 1.6576784827897006E-01_dp
      fit_coef(3, 18) = -2.3516121046359606E-01_dp
      fit_coef(4, 18) = 2.4050118552826443E+00_dp
      fit_coef(5, 18) = -1.2525237294764612E+01_dp
      fit_coef(6, 18) = 4.0835140514931439E+01_dp
      fit_coef(7, 18) = -8.8493887923230830E+01_dp
      fit_coef(8, 18) = 1.2986858749858015E+02_dp
      fit_coef(9, 18) = -1.2770555910039072E+02_dp
      fit_coef(10, 18) = 8.0636021175329191E+01_dp
      fit_coef(11, 18) = -2.9558075550791788E+01_dp
      fit_coef(12, 18) = 4.7822089064411317E+00_dp
      fit_coef(0, 19) = 9.1240612545516228E-02_dp
      fit_coef(1, 19) = 5.4034493382966092E-02_dp
      fit_coef(2, 19) = 2.9135349130073102E-01_dp
      fit_coef(3, 19) = -4.5212301023387952E-01_dp
      fit_coef(4, 19) = 4.1878151793085188E+00_dp
      fit_coef(5, 19) = -2.2306599564702282E+01_dp
      fit_coef(6, 19) = 7.3348603292660968E+01_dp
      fit_coef(7, 19) = -1.5938432150171050E+02_dp
      fit_coef(8, 19) = 2.3403743136945869E+02_dp
      fit_coef(9, 19) = -2.3012589740957122E+02_dp
      fit_coef(10, 19) = 1.4528076150095299E+02_dp
      fit_coef(11, 19) = -5.3246238298596822E+01_dp
      fit_coef(12, 19) = 8.6138613725490618E+00_dp
      fit_coef(0, 20) = 2.0179304564503545E-01_dp
      fit_coef(1, 20) = 8.9210315804133342E-02_dp
      fit_coef(2, 20) = 4.7545204864558721E-01_dp
      fit_coef(3, 20) = -7.9032122760761392E-01_dp
      fit_coef(4, 20) = 6.7848792396969912E+00_dp
      fit_coef(5, 20) = -3.6782770745377228E+01_dp
      fit_coef(6, 20) = 1.2191045579708324E+02_dp
      fit_coef(7, 20) = -2.6562927550077575E+02_dp
      fit_coef(8, 20) = 3.9032415807076654E+02_dp
      fit_coef(9, 20) = -3.8381139592009845E+02_dp
      fit_coef(10, 20) = 2.4226291326011409E+02_dp
      fit_coef(11, 20) = -8.8773520901146185E+01_dp
      fit_coef(12, 20) = 1.4358956500859366E+01_dp
      fit_coef(0, 21) = 4.2593952773320098E-01_dp
      fit_coef(1, 21) = 1.3568653118311338E-01_dp
      fit_coef(2, 21) = 7.1623822743042354E-01_dp
      fit_coef(3, 21) = -1.2560727864754864E+00_dp
      fit_coef(4, 21) = 1.0163959042284175E+01_dp
      fit_coef(5, 21) = -5.5870348063658447E+01_dp
      fit_coef(6, 21) = 1.8650175155799977E+02_dp
      fit_coef(7, 21) = -4.0744671452434335E+02_dp
      fit_coef(8, 21) = 5.9919650717200875E+02_dp
      fit_coef(9, 21) = -5.8926562705630113E+02_dp
      fit_coef(10, 21) = 3.7190276855784708E+02_dp
      fit_coef(11, 21) = -1.3625501508902701E+02_dp
      fit_coef(12, 21) = 2.2035594424138168E+01_dp
      fit_coef(0, 22) = 8.6820953800091549E-01_dp
      fit_coef(1, 22) = 1.8780376015732225E-01_dp
      fit_coef(2, 22) = 9.8395493031048986E-01_dp
      fit_coef(3, 22) = -1.7957997667101941E+00_dp
      fit_coef(4, 22) = 1.3905769981465895E+01_dp
      fit_coef(5, 22) = -7.7228503961529029E+01_dp
      fit_coef(6, 22) = 2.5935228694562625E+02_dp
      fit_coef(7, 22) = -5.6796297694678856E+02_dp
      fit_coef(8, 22) = 8.3593005414822915E+02_dp
      fit_coef(9, 22) = -8.2221938772784824E+02_dp
      fit_coef(10, 22) = 5.1889603350223433E+02_dp
      fit_coef(11, 22) = -1.9008494460566484E+02_dp
      fit_coef(12, 22) = 3.0737084830973959E+01_dp
      fit_coef(0, 23) = 1.7481856703690022E+00_dp
      fit_coef(1, 23) = 2.3346689762517508E-01_dp
      fit_coef(2, 23) = 1.2175326565887574E+00_dp
      fit_coef(3, 23) = -2.2759160332795787E+00_dp
      fit_coef(4, 23) = 1.7165218972872630E+01_dp
      fit_coef(5, 23) = -9.5910398735102618E+01_dp
      fit_coef(6, 23) = 3.2337485666569609E+02_dp
      fit_coef(7, 23) = -7.0937150908633851E+02_dp
      fit_coef(8, 23) = 1.0447103125794226E+03_dp
      fit_coef(9, 23) = -1.0277535716241930E+03_dp
      fit_coef(10, 23) = 6.4860244114600755E+02_dp
      fit_coef(11, 23) = -2.3758314620598588E+02_dp
      fit_coef(12, 23) = 3.8414439806414237E+01_dp
      fit_coef(0, 24) = 3.7882878590945706E+00_dp
      fit_coef(1, 24) = 2.6672643511994515E-01_dp
      fit_coef(2, 24) = 1.3901597873378935E+00_dp
      fit_coef(3, 24) = -2.6062582764354025E+00_dp
      fit_coef(4, 24) = 1.9594239409959744E+01_dp
      fit_coef(5, 24) = -1.0955141547937794E+02_dp
      fit_coef(6, 24) = 3.6959856539692777E+02_dp
      fit_coef(7, 24) = -8.1102106918375466E+02_dp
      fit_coef(8, 24) = 1.1945701938386849E+03_dp
      fit_coef(9, 24) = -1.1752376758852338E+03_dp
      fit_coef(10, 24) = 7.4168531728428854E+02_dp
      fit_coef(11, 24) = -2.7167776208361096E+02_dp
      fit_coef(12, 24) = 4.3926613983789913E+01_dp

   END SUBROUTINE get_fit_coef_k12_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k13_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 4.3806959251380109E-03_dp
      fit_coef(1, 1) = 1.1805095681240434E-02_dp
      fit_coef(2, 1) = 1.5469511812292528E-02_dp
      fit_coef(3, 1) = 1.1417683853936199E-02_dp
      fit_coef(4, 1) = 3.5089172870450503E-03_dp
      fit_coef(5, 1) = -1.6212006877380839E-03_dp
      fit_coef(6, 1) = -1.8991808677819510E-03_dp
      fit_coef(7, 1) = -1.8247535956842715E-04_dp
      fit_coef(8, 1) = 4.6296227124327700E-04_dp
      fit_coef(9, 1) = 3.1714051554687601E-04_dp
      fit_coef(10, 1) = -4.2555794580904834E-04_dp
      fit_coef(11, 1) = 2.1260104216928313E-04_dp
      fit_coef(12, 1) = -4.7302059547381066E-05_dp
      fit_coef(0, 2) = 2.3696032787686165E-02_dp
      fit_coef(1, 2) = 6.2953090652548949E-02_dp
      fit_coef(2, 2) = 8.1812508977763682E-02_dp
      fit_coef(3, 2) = 5.9935099274848913E-02_dp
      fit_coef(4, 2) = 1.8183698651738394E-02_dp
      fit_coef(5, 2) = -8.5926086043909885E-03_dp
      fit_coef(6, 2) = -9.8822906424356552E-03_dp
      fit_coef(7, 2) = -8.9212086720690355E-04_dp
      fit_coef(8, 2) = 2.4413793189226370E-03_dp
      fit_coef(9, 2) = 1.5997746303520928E-03_dp
      fit_coef(10, 2) = -2.1816355097527334E-03_dp
      fit_coef(11, 2) = 1.0890520189793466E-03_dp
      fit_coef(12, 2) = -2.4201086567176617E-04_dp
      fit_coef(0, 3) = 6.1206818631275799E-02_dp
      fit_coef(1, 3) = 1.5800300818436122E-01_dp
      fit_coef(2, 3) = 2.0228324543170234E-01_dp
      fit_coef(3, 3) = 1.4627045924471369E-01_dp
      fit_coef(4, 3) = 4.3383951728004834E-02_dp
      fit_coef(5, 3) = -2.1290904984850527E-02_dp
      fit_coef(6, 3) = -2.3732647813102206E-02_dp
      fit_coef(7, 3) = -1.9047854738674234E-03_dp
      fit_coef(8, 3) = 5.9890949427493156E-03_dp
      fit_coef(9, 3) = 3.6484212780910821E-03_dp
      fit_coef(10, 3) = -5.1182607229190253E-03_dp
      fit_coef(11, 3) = 2.5506650132351419E-03_dp
      fit_coef(12, 3) = -5.6515120961717785E-04_dp
      fit_coef(0, 4) = 1.2273976907574904E-01_dp
      fit_coef(1, 4) = 3.0202610531899604E-01_dp
      fit_coef(2, 4) = 3.7810599512667076E-01_dp
      fit_coef(3, 4) = 2.6877539512343529E-01_dp
      fit_coef(4, 4) = 7.7369990738750394E-02_dp
      fit_coef(5, 4) = -3.9931781197404183E-02_dp
      fit_coef(6, 4) = -4.2604215569189550E-02_dp
      fit_coef(7, 4) = -2.7223623128689494E-03_dp
      fit_coef(8, 4) = 1.1010550845635398E-02_dp
      fit_coef(9, 4) = 5.9806779753877508E-03_dp
      fit_coef(10, 4) = -8.7651481516437067E-03_dp
      fit_coef(11, 4) = 4.3537113139725109E-03_dp
      fit_coef(12, 4) = -9.6064235434420908E-04_dp
      fit_coef(0, 5) = 2.1922884952130442E-01_dp
      fit_coef(1, 5) = 5.0239277267534854E-01_dp
      fit_coef(2, 5) = 6.0878683083947316E-01_dp
      fit_coef(3, 5) = 4.2445143636733651E-01_dp
      fit_coef(4, 5) = 1.1949507494490993E-01_dp
      fit_coef(5, 5) = -6.4343269144918072E-02_dp
      fit_coef(6, 5) = -6.6581789191562657E-02_dp
      fit_coef(7, 5) = -3.2621357858438767E-03_dp
      fit_coef(8, 5) = 2.0013098730936332E-02_dp
      fit_coef(9, 5) = 5.6869475340117050E-03_dp
      fit_coef(10, 5) = -1.2059838120484294E-02_dp
      fit_coef(11, 5) = 6.5584888368600807E-03_dp
      fit_coef(12, 5) = -1.5261431175938014E-03_dp
      fit_coef(0, 6) = 3.6944915773408882E-01_dp
      fit_coef(1, 6) = 7.6979797793437166E-01_dp
      fit_coef(2, 6) = 8.9033669564805873E-01_dp
      fit_coef(3, 6) = 6.0532992150591614E-01_dp
      fit_coef(4, 6) = 1.7130397954922461E-01_dp
      fit_coef(5, 6) = -8.8192182985475279E-02_dp
      fit_coef(6, 6) = -9.7966070392074389E-02_dp
      fit_coef(7, 6) = -4.8436829553671259E-03_dp
      fit_coef(8, 6) = 2.0833002694907763E-02_dp
      fit_coef(9, 6) = 3.5063121112200099E-02_dp
      fit_coef(10, 6) = -4.7438965646280404E-02_dp
      fit_coef(11, 6) = 2.4703796048967528E-02_dp
      fit_coef(12, 6) = -5.2577610515280489E-03_dp
      fit_coef(0, 7) = 6.0348926202620468E-01_dp
      fit_coef(1, 7) = 1.1189609399540368E+00_dp
      fit_coef(2, 7) = 1.2166554041301936E+00_dp
      fit_coef(3, 7) = 7.9564460112103685E-01_dp
      fit_coef(4, 7) = 2.2942989193522059E-01_dp
      fit_coef(5, 7) = -9.7714358227207448E-02_dp
      fit_coef(6, 7) = -1.2746711453160517E-01_dp
      fit_coef(7, 7) = -1.1047808976133488E-02_dp
      fit_coef(8, 7) = -6.1793196597745570E-03_dp
      fit_coef(9, 7) = 8.7267129421221348E-02_dp
      fit_coef(10, 7) = -6.4166653350011871E-02_dp
      fit_coef(11, 7) = 1.4747313101772477E-02_dp
      fit_coef(12, 7) = 2.3977301657934416E-04_dp
      fit_coef(0, 8) = 9.6768723967411674E-01_dp
      fit_coef(1, 8) = 1.5668848161581139E+00_dp
      fit_coef(2, 8) = 1.5818658236819807E+00_dp
      fit_coef(3, 8) = 9.7733015664802736E-01_dp
      fit_coef(4, 8) = 2.8255279451553161E-01_dp
      fit_coef(5, 8) = -8.9391803088038346E-02_dp
      fit_coef(6, 8) = -1.2421016153669295E-01_dp
      fit_coef(7, 8) = -6.1021516900378492E-02_dp
      fit_coef(8, 8) = 6.2677519441489124E-02_dp
      fit_coef(9, 8) = -8.0250399433389449E-02_dp
      fit_coef(10, 8) = 1.5369427952945630E-01_dp
      fit_coef(11, 8) = -1.1092684421625504E-01_dp
      fit_coef(12, 8) = 2.7549923051188644E-02_dp
      fit_coef(0, 9) = 1.5325340819017725E+00_dp
      fit_coef(1, 9) = 2.1287996982316222E+00_dp
      fit_coef(2, 9) = 1.9799152367840203E+00_dp
      fit_coef(3, 9) = 1.1370787116857501E+00_dp
      fit_coef(4, 9) = 3.2008697898446903E-01_dp
      fit_coef(5, 9) = -6.6640990523988475E-02_dp
      fit_coef(6, 9) = -9.5466135919243425E-02_dp
      fit_coef(7, 9) = -9.2042888056781805E-02_dp
      fit_coef(8, 9) = 9.7360896464160462E-02_dp
      fit_coef(9, 9) = -1.8127317675763444E-01_dp
      fit_coef(10, 9) = 2.2214526244390012E-01_dp
      fit_coef(11, 9) = -1.0494335632769292E-01_dp
      fit_coef(12, 9) = 1.5481275818011413E-02_dp
      fit_coef(0, 10) = 2.4061866675798735E+00_dp
      fit_coef(1, 10) = 2.8132877393110070E+00_dp
      fit_coef(2, 10) = 2.4037430753636428E+00_dp
      fit_coef(3, 10) = 1.2687163022362127E+00_dp
      fit_coef(4, 10) = 3.3899108614443485E-01_dp
      fit_coef(5, 10) = -3.2980562806816459E-02_dp
      fit_coef(6, 10) = -6.4495629707069638E-02_dp
      fit_coef(7, 10) = -4.9196501013761162E-02_dp
      fit_coef(8, 10) = -1.8973498191644747E-02_dp
      fit_coef(9, 10) = -1.6127741943174356E-02_dp
      fit_coef(10, 10) = 1.0778655156935451E-03_dp
      fit_coef(11, 10) = 4.9661103678829296E-02_dp
      fit_coef(12, 10) = -2.4208904460265530E-02_dp
      fit_coef(0, 11) = 3.7607303248227026E+00_dp
      fit_coef(1, 11) = 3.6183981611913807E+00_dp
      fit_coef(2, 11) = 2.8473208041198181E+00_dp
      fit_coef(3, 11) = 1.3735956585457292E+00_dp
      fit_coef(4, 11) = 3.4303721265244597E-01_dp
      fit_coef(5, 11) = 5.1365790260944994E-03_dp
      fit_coef(6, 11) = -2.3683792375462247E-02_dp
      fit_coef(7, 11) = -7.7157824160531199E-03_dp
      fit_coef(8, 11) = -9.7425777898074148E-02_dp
      fit_coef(9, 11) = 8.3518836722114576E-02_dp
      fit_coef(10, 11) = -1.1044052160199369E-01_dp
      fit_coef(11, 11) = 1.0133376735634138E-01_dp
      fit_coef(12, 11) = -2.9906470245636543E-02_dp
      fit_coef(0, 12) = 5.8979582300825903E+00_dp
      fit_coef(1, 12) = 4.5329131993702125E+00_dp
      fit_coef(2, 12) = 3.3115813786557839E+00_dp
      fit_coef(3, 12) = 1.4609520855498697E+00_dp
      fit_coef(4, 12) = 3.3783414333423767E-01_dp
      fit_coef(5, 12) = 4.3826302981445429E-02_dp
      fit_coef(6, 12) = 2.7923917468353154E-02_dp
      fit_coef(7, 12) = 2.7776299343854542E-03_dp
      fit_coef(8, 12) = -8.3036915110656417E-02_dp
      fit_coef(9, 12) = 3.7742211424302625E-02_dp
      fit_coef(10, 12) = -4.9885068616704199E-02_dp
      fit_coef(11, 12) = 4.5143513558608316E-02_dp
      fit_coef(12, 12) = -1.0663317459135431E-02_dp
      fit_coef(0, 13) = 9.5045784867374294E+00_dp
      fit_coef(1, 13) = 5.5579146894547398E+00_dp
      fit_coef(2, 13) = 3.8186762276220274E+00_dp
      fit_coef(3, 13) = 1.5494988454410836E+00_dp
      fit_coef(4, 13) = 3.2780283598581783E-01_dp
      fit_coef(5, 13) = 8.1835332761182517E-02_dp
      fit_coef(6, 13) = 8.1960790462737673E-02_dp
      fit_coef(7, 13) = 8.9017910412122984E-03_dp
      fit_coef(8, 13) = -5.2756529738479566E-02_dp
      fit_coef(9, 13) = -1.9132069822996558E-02_dp
      fit_coef(10, 13) = 2.2550262525116004E-02_dp
      fit_coef(11, 13) = -7.1174631455651191E-03_dp
      fit_coef(12, 13) = 2.6153386964695999E-03_dp
      fit_coef(0, 14) = 1.1310131998968936E-02_dp
      fit_coef(1, 14) = 3.0381882917892972E-02_dp
      fit_coef(2, 14) = 3.9735561262063152E-02_dp
      fit_coef(3, 14) = 2.9276417146478988E-02_dp
      fit_coef(4, 14) = 8.9698683531803033E-03_dp
      fit_coef(5, 14) = -4.1668138620200067E-03_dp
      fit_coef(6, 14) = -4.8597312635549313E-03_dp
      fit_coef(7, 14) = -4.6021605223065968E-04_dp
      fit_coef(8, 14) = 1.1885176206365188E-03_dp
      fit_coef(9, 14) = 8.0548775761279056E-04_dp
      fit_coef(10, 14) = -1.0849907036440899E-03_dp
      fit_coef(11, 14) = 5.4194768198807949E-04_dp
      fit_coef(12, 14) = -1.2054873767366686E-04_dp
      fit_coef(0, 15) = 2.7754333434518740E-02_dp
      fit_coef(1, 15) = 7.2425221320615574E-02_dp
      fit_coef(2, 15) = 9.3171796368255483E-02_dp
      fit_coef(3, 15) = 6.7635165400501615E-02_dp
      fit_coef(4, 15) = 2.0195124261742557E-02_dp
      fit_coef(5, 15) = -9.8067712377566344E-03_dp
      fit_coef(6, 15) = -1.1030330619304156E-02_dp
      fit_coef(7, 15) = -9.1656688079229311E-04_dp
      fit_coef(8, 15) = 2.7689224656143080E-03_dp
      fit_coef(9, 15) = 1.7173968423888279E-03_dp
      fit_coef(10, 15) = -2.3913024611684291E-03_dp
      fit_coef(11, 15) = 1.1924633353927081E-03_dp
      fit_coef(12, 15) = -2.6453134445543388E-04_dp
      fit_coef(0, 16) = 4.8205145888088087E-02_dp
      fit_coef(1, 16) = 1.1851851999761799E-01_dp
      fit_coef(2, 16) = 1.4800321928791688E-01_dp
      fit_coef(3, 16) = 1.0475875377187537E-01_dp
      fit_coef(4, 16) = 2.9891263684815553E-02_dp
      fit_coef(5, 16) = -1.5617974138582616E-02_dp
      fit_coef(6, 16) = -1.6518323399901058E-02_dp
      fit_coef(7, 16) = -1.0386542240169942E-03_dp
      fit_coef(8, 16) = 4.3031804541452875E-03_dp
      fit_coef(9, 16) = 2.3296942214223468E-03_dp
      fit_coef(10, 16) = -3.4336179918521620E-03_dp
      fit_coef(11, 16) = 1.7042723642156064E-03_dp
      fit_coef(12, 16) = -3.7515429160435323E-04_dp
      fit_coef(0, 17) = 7.6637000090091376E-02_dp
      fit_coef(1, 17) = 1.7073623076407599E-01_dp
      fit_coef(2, 17) = 2.0360358978335744E-01_dp
      fit_coef(3, 17) = 1.3981725357818653E-01_dp
      fit_coef(4, 17) = 3.7911755072739847E-02_dp
      fit_coef(5, 17) = -2.1730826974108042E-02_dp
      fit_coef(6, 17) = -2.1165853445559209E-02_dp
      fit_coef(7, 17) = -4.3396413995129904E-04_dp
      fit_coef(8, 17) = 5.9819715180725796E-03_dp
      fit_coef(9, 17) = 1.7458780018996711E-03_dp
      fit_coef(10, 17) = -3.3770780241614367E-03_dp
      fit_coef(11, 17) = 1.6967821161902551E-03_dp
      fit_coef(12, 17) = -3.7810404338462888E-04_dp
      fit_coef(0, 18) = 1.1941431468383547E-01_dp
      fit_coef(1, 18) = 2.3176044822900177E-01_dp
      fit_coef(2, 18) = 2.5706191693734187E-01_dp
      fit_coef(3, 18) = 1.7029046402780762E-01_dp
      fit_coef(4, 18) = 4.6762674262981342E-02_dp
      fit_coef(5, 18) = -2.6181781674092602E-02_dp
      fit_coef(6, 18) = -2.7357912314188080E-02_dp
      fit_coef(7, 18) = -1.4374470171699591E-03_dp
      fit_coef(8, 18) = 1.2160163836047608E-02_dp
      fit_coef(9, 18) = 7.2851082934719002E-04_dp
      fit_coef(10, 18) = -7.7359567257111982E-03_dp
      fit_coef(11, 18) = 5.4157372236747318E-03_dp
      fit_coef(12, 18) = -1.3750847828525607E-03_dp
      fit_coef(0, 19) = 1.8596241285720252E-01_dp
      fit_coef(1, 19) = 3.0553364801134103E-01_dp
      fit_coef(2, 19) = 3.0496827359436179E-01_dp
      fit_coef(3, 19) = 1.8871536236922606E-01_dp
      fit_coef(4, 19) = 5.6375245286351912E-02_dp
      fit_coef(5, 19) = -1.8700350814395448E-02_dp
      fit_coef(6, 19) = -3.4527217392530681E-02_dp
      fit_coef(7, 19) = -6.3872694617221477E-05_dp
      fit_coef(8, 19) = -2.2205235124562510E-02_dp
      fit_coef(9, 19) = 7.0816767471594697E-02_dp
      fit_coef(10, 19) = -6.6120085860146768E-02_dp
      fit_coef(11, 19) = 2.8289249488746906E-02_dp
      fit_coef(12, 19) = -4.7843132337791906E-03_dp
      fit_coef(0, 20) = 2.8971434006808466E-01_dp
      fit_coef(1, 20) = 3.9569789839598735E-01_dp
      fit_coef(2, 20) = 3.4670622356516184E-01_dp
      fit_coef(3, 20) = 1.8876868004229172E-01_dp
      fit_coef(4, 20) = 5.7788519919382296E-02_dp
      fit_coef(5, 20) = 5.4887363554269302E-05_dp
      fit_coef(6, 20) = -1.7484482268241279E-02_dp
      fit_coef(7, 20) = -2.3812061294207602E-02_dp
      fit_coef(8, 20) = 3.6673498428854141E-03_dp
      fit_coef(9, 20) = -3.6033615474019322E-02_dp
      fit_coef(10, 20) = 9.9463288944652603E-02_dp
      fit_coef(11, 20) = -7.9311225375931868E-02_dp
      fit_coef(12, 20) = 2.1453780314514743E-02_dp
      fit_coef(0, 21) = 4.5021808525283785E-01_dp
      fit_coef(1, 21) = 5.0277275639459118E-01_dp
      fit_coef(2, 21) = 3.8272861573261618E-01_dp
      fit_coef(3, 21) = 1.7235807446302881E-01_dp
      fit_coef(4, 21) = 4.6572363862520262E-02_dp
      fit_coef(5, 21) = 1.5881707525432769E-02_dp
      fit_coef(6, 21) = 2.2025322210179048E-02_dp
      fit_coef(7, 21) = -6.2308984089511119E-02_dp
      fit_coef(8, 21) = 1.0309535046393066E-01_dp
      fit_coef(9, 21) = -2.2967091111578883E-01_dp
      fit_coef(10, 21) = 2.4313367642185343E-01_dp
      fit_coef(11, 21) = -1.1137814444194188E-01_dp
      fit_coef(12, 21) = 1.7887958008893388E-02_dp
      fit_coef(0, 22) = 6.9710873795862005E-01_dp
      fit_coef(1, 22) = 6.2256041357045133E-01_dp
      fit_coef(2, 22) = 4.1216050290755374E-01_dp
      fit_coef(3, 22) = 1.4610725607133404E-01_dp
      fit_coef(4, 22) = 2.8009636014951746E-02_dp
      fit_coef(5, 22) = 2.9171197185645873E-02_dp
      fit_coef(6, 22) = 3.0389459142154723E-02_dp
      fit_coef(7, 22) = 1.2605729292193118E-02_dp
      fit_coef(8, 22) = -6.0485512980266130E-02_dp
      fit_coef(9, 22) = 7.0433726753735382E-02_dp
      fit_coef(10, 22) = -1.3375453829923470E-01_dp
      fit_coef(11, 22) = 1.2456403150558767E-01_dp
      fit_coef(12, 22) = -3.8728626230412060E-02_dp
      fit_coef(0, 23) = 1.0779765313876117E+00_dp
      fit_coef(1, 23) = 7.4609279169328480E-01_dp
      fit_coef(2, 23) = 4.3439321318574559E-01_dp
      fit_coef(3, 23) = 1.1738611843507005E-01_dp
      fit_coef(4, 23) = 1.0530309563700853E-02_dp
      fit_coef(5, 23) = 3.6897191893439256E-02_dp
      fit_coef(6, 23) = 3.4044088146621569E-02_dp
      fit_coef(7, 23) = 5.5410057699209495E-02_dp
      fit_coef(8, 23) = -1.2792045350588788E-01_dp
      fit_coef(9, 23) = 1.8262592506787398E-01_dp
      fit_coef(10, 23) = -2.1867054150540904E-01_dp
      fit_coef(11, 23) = 1.2980623399937813E-01_dp
      fit_coef(12, 23) = -2.7516040284682932E-02_dp
      fit_coef(0, 24) = 1.6785411350231978E+00_dp
      fit_coef(1, 24) = 8.6201990862233702E-01_dp
      fit_coef(2, 24) = 4.5273971353495590E-01_dp
      fit_coef(3, 24) = 9.3905664009676981E-02_dp
      fit_coef(4, 24) = -1.4517367291959843E-03_dp
      fit_coef(5, 24) = 3.8726538172099936E-02_dp
      fit_coef(6, 24) = 4.7579808176484301E-02_dp
      fit_coef(7, 24) = 2.3323849168230236E-02_dp
      fit_coef(8, 24) = -2.2238667006060801E-02_dp
      fit_coef(9, 24) = 8.3019689769436942E-03_dp
      fit_coef(10, 24) = -8.7391379765216180E-04_dp
      fit_coef(11, 24) = -2.1892298772107112E-02_dp
      fit_coef(12, 24) = 1.2948019279576950E-02_dp
      fit_coef(0, 25) = 2.6946568599370697E+00_dp
      fit_coef(1, 25) = 9.6576339244071330E-01_dp
      fit_coef(2, 25) = 4.7873508727946390E-01_dp
      fit_coef(3, 25) = 8.3762183149891992E-02_dp
      fit_coef(4, 25) = -8.2300319329272915E-03_dp
      fit_coef(5, 25) = 3.8469555958125901E-02_dp
      fit_coef(6, 25) = 5.3859007054592130E-02_dp
      fit_coef(7, 25) = 3.4301845918780008E-03_dp
      fit_coef(8, 25) = 3.4783546363171880E-02_dp
      fit_coef(9, 25) = -7.1844324491727757E-02_dp
      fit_coef(10, 25) = 8.9320542108673212E-02_dp
      fit_coef(11, 25) = -6.8466536447383028E-02_dp
      fit_coef(12, 25) = 1.9631577406052880E-02_dp
      fit_coef(0, 26) = 4.8659851910149703E+00_dp
      fit_coef(1, 26) = 1.1016222687932846E+00_dp
      fit_coef(2, 26) = 5.5082141072412572E-01_dp
      fit_coef(3, 26) = 9.9924055011494373E-02_dp
      fit_coef(4, 26) = -1.1242348443643024E-02_dp
      fit_coef(5, 26) = 3.7638525869382183E-02_dp
      fit_coef(6, 26) = 5.4458588514965843E-02_dp
      fit_coef(7, 26) = 1.1975138580003340E-02_dp
      fit_coef(8, 26) = 1.9387184067595927E-02_dp
      fit_coef(9, 26) = -3.2120751380981701E-02_dp
      fit_coef(10, 26) = 4.5502034364108004E-02_dp
      fit_coef(11, 26) = -3.0627416778773675E-02_dp
      fit_coef(12, 26) = 6.0865046197086294E-03_dp

   END SUBROUTINE get_fit_coef_k13_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k13_mi(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.1693251854461644E+03_dp
      fit_coef(0, 1) = 5.9571247085285603E-04_dp
      fit_coef(1, 1) = 1.7200118003776515E-03_dp
      fit_coef(2, 1) = 2.4894492611009689E-03_dp
      fit_coef(3, 1) = 2.1537471215151561E-03_dp
      fit_coef(4, 1) = 1.0174139924949508E-03_dp
      fit_coef(5, 1) = 4.6357691642043961E-05_dp
      fit_coef(6, 1) = -2.2922455590426252E-04_dp
      fit_coef(7, 1) = -4.5027856096292727E-05_dp
      fit_coef(8, 1) = 3.2368140787386204E-05_dp
      fit_coef(9, 1) = 7.8015031488317424E-05_dp
      fit_coef(10, 1) = -9.1410611282450338E-05_dp
      fit_coef(11, 1) = 4.0455359969903499E-05_dp
      fit_coef(12, 1) = -6.8554946194150892E-06_dp
      fit_coef(0, 2) = 3.3749978736579119E-03_dp
      fit_coef(1, 2) = 9.3865303071796868E-03_dp
      fit_coef(2, 2) = 1.3370050944546014E-02_dp
      fit_coef(3, 2) = 1.1408817834538747E-02_dp
      fit_coef(4, 2) = 5.3018599368705056E-03_dp
      fit_coef(5, 2) = 1.9895208136664570E-04_dp
      fit_coef(6, 2) = -1.2013959538153690E-03_dp
      fit_coef(7, 2) = -2.1522942347644450E-04_dp
      fit_coef(8, 2) = 1.6227159652722640E-04_dp
      fit_coef(9, 2) = 4.0937759938008522E-04_dp
      fit_coef(10, 2) = -4.8203337203911678E-04_dp
      fit_coef(11, 2) = 2.1503525076000546E-04_dp
      fit_coef(12, 2) = -3.6767133200206401E-05_dp
      fit_coef(0, 3) = 9.5695858129820257E-03_dp
      fit_coef(1, 3) = 2.4594246941685691E-02_dp
      fit_coef(2, 3) = 3.3909912181295111E-02_dp
      fit_coef(3, 3) = 2.8275272197259928E-02_dp
      fit_coef(4, 3) = 1.2793089131376839E-02_dp
      fit_coef(5, 3) = 2.8137691410500335E-04_dp
      fit_coef(6, 3) = -2.9348509030753785E-03_dp
      fit_coef(7, 3) = -4.3154495298106305E-04_dp
      fit_coef(8, 3) = 4.1400164345487900E-04_dp
      fit_coef(9, 3) = 9.1497967805617892E-04_dp
      fit_coef(10, 3) = -1.1243533488522371E-03_dp
      fit_coef(11, 3) = 5.1331094158916764E-04_dp
      fit_coef(12, 3) = -8.9543140456745687E-05_dp
      fit_coef(0, 4) = 2.2102062421629307E-02_dp
      fit_coef(1, 4) = 5.0388239734107464E-02_dp
      fit_coef(2, 4) = 6.5597474400466413E-02_dp
      fit_coef(3, 4) = 5.2593837264420067E-02_dp
      fit_coef(4, 4) = 2.3138245186338603E-02_dp
      fit_coef(5, 4) = 3.5363359503456162E-04_dp
      fit_coef(6, 4) = -5.4308630565349861E-03_dp
      fit_coef(7, 4) = -7.1690573641355025E-04_dp
      fit_coef(8, 4) = 4.5525716627285349E-04_dp
      fit_coef(9, 4) = 2.3355643424946655E-03_dp
      fit_coef(10, 4) = -2.6681070004726240E-03_dp
      fit_coef(11, 4) = 1.2090694842315898E-03_dp
      fit_coef(12, 4) = -2.1156836663214308E-04_dp
      fit_coef(0, 5) = 4.6953984288651295E-02_dp
      fit_coef(1, 5) = 9.2513506061771736E-02_dp
      fit_coef(2, 5) = 1.1114827337176617E-01_dp
      fit_coef(3, 5) = 8.3233678613110429E-02_dp
      fit_coef(4, 5) = 3.4608571546062436E-02_dp
      fit_coef(5, 5) = 6.2863207858912962E-04_dp
      fit_coef(6, 5) = -7.4478035609513958E-03_dp
      fit_coef(7, 5) = -1.1406622073555459E-03_dp
      fit_coef(8, 5) = 6.1142108219932139E-05_dp
      fit_coef(9, 5) = 2.8713082420645268E-03_dp
      fit_coef(10, 5) = -2.3152373063920332E-03_dp
      fit_coef(11, 5) = 6.9898905114079072E-04_dp
      fit_coef(12, 5) = -6.0162701114697736E-05_dp
      fit_coef(0, 6) = 9.5182645385563286E-02_dp
      fit_coef(1, 6) = 1.5985393660033628E-01_dp
      fit_coef(2, 6) = 1.7511073340807867E-01_dp
      fit_coef(3, 6) = 1.1909736973485001E-01_dp
      fit_coef(4, 6) = 4.4162492084356325E-02_dp
      fit_coef(5, 6) = -1.7481218668970254E-04_dp
      fit_coef(6, 6) = -7.8382486311276943E-03_dp
      fit_coef(7, 6) = -6.5498122365884547E-04_dp
      fit_coef(8, 6) = 7.2110889562901611E-05_dp
      fit_coef(9, 6) = 9.5978488993906188E-04_dp
      fit_coef(10, 6) = -6.1518613580109974E-04_dp
      fit_coef(11, 6) = 2.4316501239763139E-04_dp
      fit_coef(12, 6) = -5.3529912569972472E-05_dp
      fit_coef(0, 7) = 1.8645165311079259E-01_dp
      fit_coef(1, 7) = 2.6436160380537138E-01_dp
      fit_coef(2, 7) = 2.6302083875944771E-01_dp
      fit_coef(3, 7) = 1.5910774510414938E-01_dp
      fit_coef(4, 7) = 4.9073653325355709E-02_dp
      fit_coef(5, 7) = -3.2795187276607025E-03_dp
      fit_coef(6, 7) = -6.8810055747994857E-03_dp
      fit_coef(7, 7) = 2.0021507448155459E-03_dp
      fit_coef(8, 7) = -1.4027814134567758E-03_dp
      fit_coef(9, 7) = 6.9098043550271603E-04_dp
      fit_coef(10, 7) = -1.0637655745511332E-03_dp
      fit_coef(11, 7) = 9.2536193146607883E-04_dp
      fit_coef(12, 7) = -2.4940185896462312E-04_dp
      fit_coef(0, 8) = 3.5501666895531042E-01_dp
      fit_coef(1, 8) = 4.2062390191703830E-01_dp
      fit_coef(2, 8) = 3.8014762337012359E-01_dp
      fit_coef(3, 8) = 2.0109387264368472E-01_dp
      fit_coef(4, 8) = 4.6954945508850141E-02_dp
      fit_coef(5, 8) = -8.7244019551801829E-03_dp
      fit_coef(6, 8) = -4.0487338235444901E-03_dp
      fit_coef(7, 8) = 5.9512521905139611E-03_dp
      fit_coef(8, 8) = -3.7498978634030503E-03_dp
      fit_coef(9, 8) = 2.1099763933573197E-04_dp
      fit_coef(10, 8) = -4.3311783252141623E-04_dp
      fit_coef(11, 8) = 7.0199947379022344E-04_dp
      fit_coef(12, 8) = -2.0527979943177016E-04_dp
      fit_coef(0, 9) = 6.5960294052628099E-01_dp
      fit_coef(1, 9) = 6.4434244249393713E-01_dp
      fit_coef(2, 9) = 5.3018928176744540E-01_dp
      fit_coef(3, 9) = 2.4130639214143917E-01_dp
      fit_coef(4, 9) = 3.5653141239669908E-02_dp
      fit_coef(5, 9) = -1.5407439897727408E-02_dp
      fit_coef(6, 9) = 1.6157905391085510E-03_dp
      fit_coef(7, 9) = 1.0708866394811738E-02_dp
      fit_coef(8, 9) = -8.0705877921707705E-03_dp
      fit_coef(9, 9) = 1.9585451901869654E-04_dp
      fit_coef(10, 9) = 1.1606320581042068E-03_dp
      fit_coef(11, 9) = -2.3231700639256246E-04_dp
      fit_coef(12, 9) = -1.9173218066136157E-05_dp
      fit_coef(0, 10) = 1.1999125471871359E+00_dp
      fit_coef(1, 10) = 9.4887754086202536E-01_dp
      fit_coef(2, 10) = 7.1392148306788861E-01_dp
      fit_coef(3, 10) = 2.7502418030970321E-01_dp
      fit_coef(4, 10) = 1.3864789499781904E-02_dp
      fit_coef(5, 10) = -2.1424737017769674E-02_dp
      fit_coef(6, 10) = 1.0990606216211212E-02_dp
      fit_coef(7, 10) = 1.5754376025894169E-02_dp
      fit_coef(8, 10) = -1.5720284237712618E-02_dp
      fit_coef(9, 10) = 2.0987333422688541E-03_dp
      fit_coef(10, 10) = 2.9263737858274525E-03_dp
      fit_coef(11, 10) = -1.6119558060080142E-03_dp
      fit_coef(12, 10) = 2.7700177106853099E-04_dp
      fit_coef(0, 11) = 2.1464589921277430E+00_dp
      fit_coef(1, 11) = 1.3398008082135153E+00_dp
      fit_coef(2, 11) = 9.2855004111736639E-01_dp
      fit_coef(3, 11) = 2.9850369229739804E-01_dp
      fit_coef(4, 11) = -1.8024501913163959E-02_dp
      fit_coef(5, 11) = -2.4792994220460485E-02_dp
      fit_coef(6, 11) = 2.4309270311746085E-02_dp
      fit_coef(7, 11) = 2.0376050200980579E-02_dp
      fit_coef(8, 11) = -2.7413278630826550E-02_dp
      fit_coef(9, 11) = 7.1577405110793838E-03_dp
      fit_coef(10, 11) = 4.1376588773677659E-03_dp
      fit_coef(11, 11) = -3.2837333351403799E-03_dp
      fit_coef(12, 11) = 6.9015069153665373E-04_dp
      fit_coef(0, 12) = 3.8079327315676679E+00_dp
      fit_coef(1, 12) = 1.8102136985271529E+00_dp
      fit_coef(2, 12) = 1.1693768610866579E+00_dp
      fit_coef(3, 12) = 3.1192154530912264E-01_dp
      fit_coef(4, 12) = -5.7801420614916721E-02_dp
      fit_coef(5, 12) = -2.4725077881519446E-02_dp
      fit_coef(6, 12) = 4.0887008920651363E-02_dp
      fit_coef(7, 12) = 2.4169421898620948E-02_dp
      fit_coef(8, 12) = -4.2916907464272722E-02_dp
      fit_coef(9, 12) = 1.5836520381657740E-02_dp
      fit_coef(10, 12) = 4.2153911444231572E-03_dp
      fit_coef(11, 12) = -5.0069163439141369E-03_dp
      fit_coef(12, 12) = 1.1850222139382642E-03_dp
      fit_coef(0, 13) = 6.8767988424609685E+00_dp
      fit_coef(1, 13) = 2.3477041904049458E+00_dp
      fit_coef(2, 13) = 1.4371470680965179E+00_dp
      fit_coef(3, 13) = 3.2177727032276354E-01_dp
      fit_coef(4, 13) = -1.0285806633893342E-01_dp
      fit_coef(5, 13) = -2.2677975661685234E-02_dp
      fit_coef(6, 13) = 5.9703880119426990E-02_dp
      fit_coef(7, 13) = 2.7586476735575311E-02_dp
      fit_coef(8, 13) = -6.1182232779440265E-02_dp
      fit_coef(9, 13) = 2.7103387955214536E-02_dp
      fit_coef(10, 13) = 3.4260870398828006E-03_dp
      fit_coef(11, 13) = -6.7215676437135359E-03_dp
      fit_coef(12, 13) = 1.7275038459318248E-03_dp
      fit_coef(0, 14) = 1.5536180830079045E-03_dp
      fit_coef(1, 14) = 4.4506140115349617E-03_dp
      fit_coef(2, 14) = 6.4178048727277781E-03_dp
      fit_coef(3, 14) = 5.5342602456232345E-03_dp
      fit_coef(4, 14) = 2.6046871795835742E-03_dp
      fit_coef(5, 14) = 1.1368344431687511E-04_dp
      fit_coef(6, 14) = -5.8781627121643181E-04_dp
      fit_coef(7, 14) = -1.1285527376227048E-04_dp
      fit_coef(8, 14) = 8.2159134979200638E-05_dp
      fit_coef(9, 14) = 1.9994801604583968E-04_dp
      fit_coef(10, 14) = -2.3463016476463693E-04_dp
      fit_coef(11, 14) = 1.0405614708917490E-04_dp
      fit_coef(12, 14) = -1.7674175132111243E-05_dp
      fit_coef(0, 15) = 4.1834098874360473E-03_dp
      fit_coef(1, 15) = 1.1102877368401057E-02_dp
      fit_coef(2, 15) = 1.5508145773336595E-02_dp
      fit_coef(3, 15) = 1.3020239177249622E-02_dp
      fit_coef(4, 15) = 5.9265544529938010E-03_dp
      fit_coef(5, 15) = 1.6028851543557779E-04_dp
      fit_coef(6, 15) = -1.3496676957738614E-03_dp
      fit_coef(7, 15) = -2.1194238008936958E-04_dp
      fit_coef(8, 15) = 1.6858962619466225E-04_dp
      fit_coef(9, 15) = 4.6274241631720910E-04_dp
      fit_coef(10, 15) = -5.4510567016319549E-04_dp
      fit_coef(11, 15) = 2.4471543070832161E-04_dp
      fit_coef(12, 15) = -4.2146500465828684E-05_dp
      fit_coef(0, 16) = 8.6718931437848486E-03_dp
      fit_coef(1, 16) = 1.9798273846496010E-02_dp
      fit_coef(2, 16) = 2.5805909237192746E-02_dp
      fit_coef(3, 16) = 2.0709762497138948E-02_dp
      fit_coef(4, 16) = 9.0267828044289580E-03_dp
      fit_coef(5, 16) = -8.2030242003155534E-06_dp
      fit_coef(6, 16) = -2.1466923201882496E-03_dp
      fit_coef(7, 16) = -2.3839418666152536E-04_dp
      fit_coef(8, 16) = 3.9713557969496511E-04_dp
      fit_coef(9, 16) = 5.4418611978137037E-04_dp
      fit_coef(10, 16) = -7.9800986012929951E-04_dp
      fit_coef(11, 16) = 3.9897082718071900E-04_dp
      fit_coef(12, 16) = -7.5484271874357470E-05_dp
      fit_coef(0, 17) = 1.7358581030722377E-02_dp
      fit_coef(1, 17) = 3.2741569114568887E-02_dp
      fit_coef(2, 17) = 3.8005765925330802E-02_dp
      fit_coef(3, 17) = 2.7705507640870534E-02_dp
      fit_coef(4, 17) = 1.1375410029678385E-02_dp
      fit_coef(5, 17) = 2.5826247104955189E-04_dp
      fit_coef(6, 17) = -2.6224716007508839E-03_dp
      fit_coef(7, 17) = -3.1373038771318260E-04_dp
      fit_coef(8, 17) = -6.4482138467522902E-04_dp
      fit_coef(9, 17) = 2.3188017066251647E-03_dp
      fit_coef(10, 17) = -1.8974564925171979E-03_dp
      fit_coef(11, 17) = 6.6385240884256846E-04_dp
      fit_coef(12, 17) = -8.2217281900096977E-05_dp
      fit_coef(0, 18) = 3.4131291193668846E-02_dp
      fit_coef(1, 18) = 5.2991311639389929E-02_dp
      fit_coef(2, 18) = 5.3873497163302024E-02_dp
      fit_coef(3, 18) = 3.3395673340616244E-02_dp
      fit_coef(4, 18) = 1.1018355439494237E-02_dp
      fit_coef(5, 18) = 1.9602178372730269E-05_dp
      fit_coef(6, 18) = -1.1688901927899766E-03_dp
      fit_coef(7, 18) = -4.0116811266019766E-04_dp
      fit_coef(8, 18) = 3.2072961899086736E-04_dp
      fit_coef(9, 18) = -1.9201145783602021E-03_dp
      fit_coef(10, 18) = 2.5996845486259772E-03_dp
      fit_coef(11, 18) = -1.3613851440031023E-03_dp
      fit_coef(12, 18) = 2.6284577266591210E-04_dp
      fit_coef(0, 19) = 6.5506726065987148E-02_dp
      fit_coef(1, 19) = 8.3686618695792928E-02_dp
      fit_coef(2, 19) = 7.4994836901648418E-02_dp
      fit_coef(3, 19) = 3.8175402517462660E-02_dp
      fit_coef(4, 19) = 7.6429339054575234E-03_dp
      fit_coef(5, 19) = -1.8419060705385749E-03_dp
      fit_coef(6, 19) = 2.7308690208177329E-04_dp
      fit_coef(7, 19) = 1.6627660280614655E-03_dp
      fit_coef(8, 19) = -8.3731649392740858E-04_dp
      fit_coef(9, 19) = -8.4805699699116934E-04_dp
      fit_coef(10, 19) = 1.7342012245697370E-05_dp
      fit_coef(11, 19) = 6.1766339244338580E-04_dp
      fit_coef(12, 19) = -2.3138634064074576E-04_dp
      fit_coef(0, 20) = 1.2258176989295384E-01_dp
      fit_coef(1, 20) = 1.2778305763600906E-01_dp
      fit_coef(2, 20) = 1.0172491446523661E-01_dp
      fit_coef(3, 20) = 4.1491871675887949E-02_dp
      fit_coef(4, 20) = 1.7844232528150912E-03_dp
      fit_coef(5, 20) = -4.3718842072007049E-03_dp
      fit_coef(6, 20) = 1.7506621259053797E-03_dp
      fit_coef(7, 20) = 3.4297516381714774E-03_dp
      fit_coef(8, 20) = -1.8677852767092760E-03_dp
      fit_coef(9, 20) = -2.5672855213340693E-04_dp
      fit_coef(10, 20) = -1.3122155242537307E-04_dp
      fit_coef(11, 20) = 3.0827036993633989E-04_dp
      fit_coef(12, 20) = -7.5397667904590979E-05_dp
      fit_coef(0, 21) = 2.2407341908539469E-01_dp
      fit_coef(1, 21) = 1.8741517537369246E-01_dp
      fit_coef(2, 21) = 1.3316931550316194E-01_dp
      fit_coef(3, 21) = 4.1846304068269342E-02_dp
      fit_coef(4, 21) = -6.3923892927403468E-03_dp
      fit_coef(5, 21) = -6.3290486135616480E-03_dp
      fit_coef(6, 21) = 4.0785090266698392E-03_dp
      fit_coef(7, 21) = 4.4129649982871406E-03_dp
      fit_coef(8, 21) = -3.0701955863932482E-03_dp
      fit_coef(9, 21) = -5.0155854005366068E-04_dp
      fit_coef(10, 21) = 1.2677059967610719E-03_dp
      fit_coef(11, 21) = -6.5728561506278344E-04_dp
      fit_coef(12, 21) = 1.3236733395452350E-04_dp
      fit_coef(0, 22) = 4.0125849549125764E-01_dp
      fit_coef(1, 22) = 2.6234441978160788E-01_dp
      fit_coef(2, 22) = 1.6705808644309805E-01_dp
      fit_coef(3, 22) = 3.7759081921101723E-02_dp
      fit_coef(4, 22) = -1.6444808051945672E-02_dp
      fit_coef(5, 22) = -6.7086735514358875E-03_dp
      fit_coef(6, 22) = 7.4079840573032128E-03_dp
      fit_coef(7, 22) = 5.0153636055701280E-03_dp
      fit_coef(8, 22) = -5.8034207254719751E-03_dp
      fit_coef(9, 22) = 7.1487044000835261E-04_dp
      fit_coef(10, 22) = 1.7935699565438628E-03_dp
      fit_coef(11, 22) = -1.1756448775298205E-03_dp
      fit_coef(12, 22) = 2.3788942753303003E-04_dp
      fit_coef(0, 23) = 7.0683650900241923E-01_dp
      fit_coef(1, 23) = 3.4775269742355053E-01_dp
      fit_coef(2, 23) = 1.9991809868316918E-01_dp
      fit_coef(3, 23) = 2.8998979248009617E-02_dp
      fit_coef(4, 23) = -2.7072247921548855E-02_dp
      fit_coef(5, 23) = -4.9720416004814215E-03_dp
      fit_coef(6, 23) = 1.1386666180451437E-02_dp
      fit_coef(7, 23) = 4.9476484479521700E-03_dp
      fit_coef(8, 23) = -9.6210622467224182E-03_dp
      fit_coef(9, 23) = 3.3069098482828903E-03_dp
      fit_coef(10, 23) = 1.6164601539057236E-03_dp
      fit_coef(11, 23) = -1.5619532900798120E-03_dp
      fit_coef(12, 23) = 3.5694695591620168E-04_dp
      fit_coef(0, 24) = 1.2349305749929025E+00_dp
      fit_coef(1, 24) = 4.3289986648386525E-01_dp
      fit_coef(2, 24) = 2.2845857393019400E-01_dp
      fit_coef(3, 24) = 1.7946385370009432E-02_dp
      fit_coef(4, 24) = -3.6314829777361717E-02_dp
      fit_coef(5, 24) = -1.6120574652244938E-03_dp
      fit_coef(6, 24) = 1.5127251697612880E-02_dp
      fit_coef(7, 24) = 4.2179377246450477E-03_dp
      fit_coef(8, 24) = -1.3727271503534512E-02_dp
      fit_coef(9, 24) = 6.9063396049124318E-03_dp
      fit_coef(10, 24) = 6.9366312151399967E-04_dp
      fit_coef(11, 24) = -1.7353304432909385E-03_dp
      fit_coef(12, 24) = 4.6246138244891298E-04_dp
      fit_coef(0, 25) = 2.1884594395954422E+00_dp
      fit_coef(1, 25) = 5.0512490214011563E-01_dp
      fit_coef(2, 25) = 2.5294331999124109E-01_dp
      fit_coef(3, 25) = 9.9552977533875799E-03_dp
      fit_coef(4, 25) = -4.2741646955111831E-02_dp
      fit_coef(5, 25) = 1.4996753556095149E-03_dp
      fit_coef(6, 25) = 1.7828818970966834E-02_dp
      fit_coef(7, 25) = 3.4400974982679682E-03_dp
      fit_coef(8, 25) = -1.7077573476609768E-02_dp
      fit_coef(9, 25) = 1.0264977475261727E-02_dp
      fit_coef(10, 25) = -4.9655867540557698E-04_dp
      fit_coef(11, 25) = -1.6960953173867819E-03_dp
      fit_coef(12, 25) = 5.2021831684559222E-04_dp
      fit_coef(0, 26) = 4.2896812115374434E+00_dp
      fit_coef(1, 26) = 5.7441806238555526E-01_dp
      fit_coef(2, 26) = 2.8764315098250082E-01_dp
      fit_coef(3, 26) = 1.2205511996262715E-02_dp
      fit_coef(4, 26) = -4.7500863719757744E-02_dp
      fit_coef(5, 26) = 1.9853897004354611E-03_dp
      fit_coef(6, 26) = 1.9866889885557001E-02_dp
      fit_coef(7, 26) = 3.6718675673585712E-03_dp
      fit_coef(8, 26) = -1.9405833778650398E-02_dp
      fit_coef(9, 26) = 1.1968919361045174E-02_dp
      fit_coef(10, 26) = -8.6190027805330149E-04_dp
      fit_coef(11, 26) = -1.8010240017088701E-03_dp
      fit_coef(12, 26) = 5.7196032777695038E-04_dp

   END SUBROUTINE get_fit_coef_k13_mi

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k13_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 3.7370000000000000E+05_dp
      fit_coef(0, 1) = 1.0328373040080946E-05_dp
      fit_coef(1, 1) = 1.0231661765573473E-05_dp
      fit_coef(2, 1) = 1.4274421628573503E-03_dp
      fit_coef(3, 1) = -9.2218261005801595E-03_dp
      fit_coef(4, 1) = 7.9453452208337924E-02_dp
      fit_coef(5, 1) = -3.4395924855303062E-01_dp
      fit_coef(6, 1) = 1.0503011873064374E+00_dp
      fit_coef(7, 1) = -2.2107082564421301E+00_dp
      fit_coef(8, 1) = 3.1904625274566456E+00_dp
      fit_coef(9, 1) = -3.1005928635963018E+00_dp
      fit_coef(10, 1) = 1.9388301331438162E+00_dp
      fit_coef(11, 1) = -7.0437721043920587E-01_dp
      fit_coef(12, 1) = 1.1297211307192120E-01_dp
      fit_coef(0, 2) = 1.1370944761350369E-04_dp
      fit_coef(1, 2) = 6.5461398441333884E-05_dp
      fit_coef(2, 2) = 9.0857762128942112E-03_dp
      fit_coef(3, 2) = -6.1037659291954947E-02_dp
      fit_coef(4, 2) = 5.0701553071696137E-01_dp
      fit_coef(5, 2) = -2.2094585548616172E+00_dp
      fit_coef(6, 2) = 6.7357245459730972E+00_dp
      fit_coef(7, 2) = -1.4152000667009908E+01_dp
      fit_coef(8, 2) = 2.0402127724921442E+01_dp
      fit_coef(9, 2) = -1.9817215032938513E+01_dp
      fit_coef(10, 2) = 1.2389779868005123E+01_dp
      fit_coef(11, 2) = -4.5013495035293962E+00_dp
      fit_coef(12, 2) = 7.2205686510589684E-01_dp
      fit_coef(0, 3) = 6.0201360833921547E-04_dp
      fit_coef(1, 3) = 2.3713142412530957E-04_dp
      fit_coef(2, 3) = 3.2408075864146658E-02_dp
      fit_coef(3, 3) = -2.2633769436834197E-01_dp
      fit_coef(4, 3) = 1.7944183579800783E+00_dp
      fit_coef(5, 3) = -8.0322858690413170E+00_dp
      fit_coef(6, 3) = 2.4758947246186221E+01_dp
      fit_coef(7, 3) = -5.2297882977111676E+01_dp
      fit_coef(8, 3) = 7.5599300021495196E+01_dp
      fit_coef(9, 3) = -7.3529871053182731E+01_dp
      fit_coef(10, 3) = 4.5999098353208382E+01_dp
      fit_coef(11, 3) = -1.6715853539005085E+01_dp
      fit_coef(12, 3) = 2.6814682156594003E+00_dp
      fit_coef(0, 4) = 2.3239931463444295E-03_dp
      fit_coef(1, 4) = 6.7915483564538726E-04_dp
      fit_coef(2, 4) = 9.2065143976354222E-02_dp
      fit_coef(3, 4) = -6.5920735047309242E-01_dp
      fit_coef(4, 4) = 5.0701238208447936E+00_dp
      fit_coef(5, 4) = -2.3078149685475569E+01_dp
      fit_coef(6, 4) = 7.1252230915973797E+01_dp
      fit_coef(7, 4) = -1.5026332268406557E+02_dp
      fit_coef(8, 4) = 2.1682358395745959E+02_dp
      fit_coef(9, 4) = -2.1059821196596170E+02_dp
      fit_coef(10, 4) = 1.3162443597110655E+02_dp
      fit_coef(11, 4) = -4.7803836464227409E+01_dp
      fit_coef(12, 4) = 7.6658373170675951E+00_dp
      fit_coef(0, 5) = 7.4432405183934578E-03_dp
      fit_coef(1, 5) = 1.6828152409194407E-03_dp
      fit_coef(2, 5) = 2.2645498811467052E-01_dp
      fit_coef(3, 5) = -1.6470414947742027E+00_dp
      fit_coef(4, 5) = 1.2406213814972343E+01_dp
      fit_coef(5, 5) = -5.7052515747077990E+01_dp
      fit_coef(6, 5) = 1.7648148382559137E+02_dp
      fit_coef(7, 5) = -3.7241542082575245E+02_dp
      fit_coef(8, 5) = 5.3766501171349421E+02_dp
      fit_coef(9, 5) = -5.2248797832997366E+02_dp
      fit_coef(10, 5) = 3.2669873907570917E+02_dp
      fit_coef(11, 5) = -1.1869478309408113E+02_dp
      fit_coef(12, 5) = 1.9039603955629481E+01_dp
      fit_coef(0, 6) = 2.1001669881294191E-02_dp
      fit_coef(1, 6) = 3.7493611984017344E-03_dp
      fit_coef(2, 6) = 5.0174562015489310E-01_dp
      fit_coef(3, 6) = -3.6903350557230272E+00_dp
      fit_coef(4, 6) = 2.7380167524745204E+01_dp
      fit_coef(5, 6) = -1.2680883651700998E+02_dp
      fit_coef(6, 6) = 3.9315723208637161E+02_dp
      fit_coef(7, 6) = -8.3061013874245850E+02_dp
      fit_coef(8, 6) = 1.2000077154718190E+03_dp
      fit_coef(9, 6) = -1.1666285380811528E+03_dp
      fit_coef(10, 6) = 7.2965528810276930E+02_dp
      fit_coef(11, 6) = -2.6514064826477511E+02_dp
      fit_coef(12, 6) = 4.2535848707897181E+01_dp
      fit_coef(0, 7) = 5.3988578814712126E-02_dp
      fit_coef(1, 7) = 7.6787040162035062E-03_dp
      fit_coef(2, 7) = 1.0229613253807799E+00_dp
      fit_coef(3, 7) = -7.5874025191500545E+00_dp
      fit_coef(4, 7) = 5.5644821860203081E+01_dp
      fit_coef(5, 7) = -2.5903592714657515E+02_dp
      fit_coef(6, 7) = 8.0491520561775883E+02_dp
      fit_coef(7, 7) = -1.7025740847831944E+03_dp
      fit_coef(8, 7) = 2.4614597158125571E+03_dp
      fit_coef(9, 7) = -2.3939670515607563E+03_dp
      fit_coef(10, 7) = 1.4976634167179416E+03_dp
      fit_coef(11, 7) = -5.4431186678214578E+02_dp
      fit_coef(12, 7) = 8.7333538775326986E+01_dp
      fit_coef(0, 8) = 1.2913416270542125E-01_dp
      fit_coef(1, 8) = 1.4637308467888006E-02_dp
      fit_coef(2, 8) = 1.9427142860381650E+00_dp
      fit_coef(3, 8) = -1.4504186559531496E+01_dp
      fit_coef(4, 8) = 1.0539283356781200E+02_dp
      fit_coef(5, 8) = -4.9250973009160060E+02_dp
      fit_coef(6, 8) = 1.5335961739523964E+03_dp
      fit_coef(7, 8) = -3.2477036056116422E+03_dp
      fit_coef(8, 8) = 4.6984634605620395E+03_dp
      fit_coef(9, 8) = -4.5714722547915217E+03_dp
      fit_coef(10, 8) = 2.8606360597296798E+03_dp
      fit_coef(11, 8) = -1.0398487637572459E+03_dp
      fit_coef(12, 8) = 1.6686131315440755E+02_dp
      fit_coef(0, 9) = 2.9151438803227880E-01_dp
      fit_coef(1, 9) = 2.6141059594544370E-02_dp
      fit_coef(2, 9) = 3.4585779341218914E+00_dp
      fit_coef(3, 9) = -2.5957371800387559E+01_dp
      fit_coef(4, 9) = 1.8719902735916295E+02_dp
      fit_coef(5, 9) = -8.7738902954943887E+02_dp
      fit_coef(6, 9) = 2.7372194924276941E+03_dp
      fit_coef(7, 9) = -5.8030086010828445E+03_dp
      fit_coef(8, 9) = 8.4006232674596267E+03_dp
      fit_coef(9, 9) = -8.1767388699954290E+03_dp
      fit_coef(10, 9) = 5.1179187652543951E+03_dp
      fit_coef(11, 9) = -1.8606812820509087E+03_dp
      fit_coef(12, 9) = 2.9861178044260612E+02_dp
      fit_coef(0, 10) = 6.2761071501945276E-01_dp
      fit_coef(1, 10) = 4.3829697389904873E-02_dp
      fit_coef(2, 10) = 5.7833019860251502E+00_dp
      fit_coef(3, 10) = -4.3589104975069581E+01_dp
      fit_coef(4, 10) = 3.1240972247660102E+02_dp
      fit_coef(5, 10) = -1.4676188400419255E+03_dp
      fit_coef(6, 10) = 4.5862715204341948E+03_dp
      fit_coef(7, 10) = -9.7328856974721111E+03_dp
      fit_coef(8, 10) = 1.4098068171267152E+04_dp
      fit_coef(9, 10) = -1.3727352353131695E+04_dp
      fit_coef(10, 10) = 8.5941131147490141E+03_dp
      fit_coef(11, 10) = -3.1249775864789372E+03_dp
      fit_coef(12, 10) = 5.0156651089860543E+02_dp
      fit_coef(0, 11) = 1.3001830139334465E+00_dp
      fit_coef(1, 11) = 6.8889203557683909E-02_dp
      fit_coef(2, 11) = 9.0695201480430416E+00_dp
      fit_coef(3, 11) = -6.8588523830477001E+01_dp
      fit_coef(4, 11) = 4.8911175818178219E+02_dp
      fit_coef(5, 11) = -2.3017939006101778E+03_dp
      fit_coef(6, 11) = 7.2034170829290761E+03_dp
      fit_coef(7, 11) = -1.5300571255050740E+04_dp
      fit_coef(8, 11) = 2.2174747250984612E+04_dp
      fit_coef(9, 11) = -2.1598799766875782E+04_dp
      fit_coef(10, 11) = 1.3524962200868669E+04_dp
      fit_coef(11, 11) = -4.9186216136818812E+03_dp
      fit_coef(12, 11) = 7.8952838092425554E+02_dp
      fit_coef(0, 12) = 2.6215772490899263E+00_dp
      fit_coef(1, 12) = 1.0119187877095515E-01_dp
      fit_coef(2, 12) = 1.3298790450785793E+01_dp
      fit_coef(3, 12) = -1.0082948191003250E+02_dp
      fit_coef(4, 12) = 7.1624051882694800E+02_dp
      fit_coef(5, 12) = -3.3750615698310812E+03_dp
      fit_coef(6, 12) = 1.0574350973373779E+04_dp
      fit_coef(7, 12) = -2.2477084129138388E+04_dp
      fit_coef(8, 12) = 3.2590012647996798E+04_dp
      fit_coef(9, 12) = -3.1752355473368389E+04_dp
      fit_coef(10, 12) = 1.9886579242865319E+04_dp
      fit_coef(11, 12) = -7.2330170769743845E+03_dp
      fit_coef(12, 12) = 1.1611278007055782E+03_dp
      fit_coef(0, 13) = 5.2970688415539806E+00_dp
      fit_coef(1, 13) = 1.3913815821587258E-01_dp
      fit_coef(2, 13) = 1.8263424659829006E+01_dp
      fit_coef(3, 13) = -1.3870952790580026E+02_dp
      fit_coef(4, 13) = 9.8271245072947818E+02_dp
      fit_coef(5, 13) = -4.6347205083298186E+03_dp
      fit_coef(6, 13) = 1.4532606648412351E+04_dp
      fit_coef(7, 13) = -3.0906783800217447E+04_dp
      fit_coef(8, 13) = 4.4826644713674119E+04_dp
      fit_coef(9, 13) = -4.3683152384513414E+04_dp
      fit_coef(10, 13) = 2.7362388150846749E+04_dp
      fit_coef(11, 13) = -9.9529300467619232E+03_dp
      fit_coef(12, 13) = 1.5978560147630030E+03_dp
      fit_coef(0, 14) = 3.3029079169995746E-05_dp
      fit_coef(1, 14) = 2.6643361318081685E-05_dp
      fit_coef(2, 14) = 3.7453854295855472E-03_dp
      fit_coef(3, 14) = -2.4461573064643052E-02_dp
      fit_coef(4, 14) = 2.1192677272898949E-01_dp
      fit_coef(5, 14) = -9.2298813642259792E-01_dp
      fit_coef(6, 14) = 2.8202954491441887E+00_dp
      fit_coef(7, 14) = -5.9307601899579083E+00_dp
      fit_coef(8, 14) = 8.5478350235106983E+00_dp
      fit_coef(9, 14) = -8.2958924766363573E+00_dp
      fit_coef(10, 14) = 5.1810926707806537E+00_dp
      fit_coef(11, 14) = -1.8802314403008740E+00_dp
      fit_coef(12, 14) = 3.0127397290756769E-01_dp
      fit_coef(0, 15) = 2.1715120936896851E-04_dp
      fit_coef(1, 15) = 9.6158838084655870E-05_dp
      fit_coef(2, 15) = 1.3102183838061880E-02_dp
      fit_coef(3, 15) = -9.0442571305710165E-02_dp
      fit_coef(4, 15) = 7.2073640737136146E-01_dp
      fit_coef(5, 15) = -3.1779329684848223E+00_dp
      fit_coef(6, 15) = 9.7478690828954182E+00_dp
      fit_coef(7, 15) = -2.0582241665854877E+01_dp
      fit_coef(8, 15) = 2.9794991676852089E+01_dp
      fit_coef(9, 15) = -2.9037927118237850E+01_dp
      fit_coef(10, 15) = 1.8203528955433185E+01_dp
      fit_coef(11, 15) = -6.6278898854648514E+00_dp
      fit_coef(12, 15) = 1.0650411372773607E+00_dp
      fit_coef(0, 16) = 8.8916557085975359E-04_dp
      fit_coef(1, 16) = 2.7147003112010167E-04_dp
      fit_coef(2, 16) = 3.6886711972600503E-02_dp
      fit_coef(3, 16) = -2.6522390438189170E-01_dp
      fit_coef(4, 16) = 2.0383361789339589E+00_dp
      fit_coef(5, 16) = -9.3159312869424529E+00_dp
      fit_coef(6, 16) = 2.8836230573430583E+01_dp
      fit_coef(7, 16) = -6.0862901561188608E+01_dp
      fit_coef(8, 16) = 8.7796437530891168E+01_dp
      fit_coef(9, 16) = -8.5203527246658652E+01_dp
      fit_coef(10, 16) = 5.3197256397040526E+01_dp
      fit_coef(11, 16) = -1.9300241563707537E+01_dp
      fit_coef(12, 16) = 3.0920034178994409E+00_dp
      fit_coef(0, 17) = 2.8909222386173440E-03_dp
      fit_coef(1, 17) = 6.5986499705299726E-04_dp
      fit_coef(2, 17) = 8.8661109522972403E-02_dp
      fit_coef(3, 17) = -6.4720090249728535E-01_dp
      fit_coef(4, 17) = 4.8509444516645317E+00_dp
      fit_coef(5, 17) = -2.2365482960360204E+01_dp
      fit_coef(6, 17) = 6.9142722370552491E+01_dp
      fit_coef(7, 17) = -1.4576807130419496E+02_dp
      fit_coef(8, 17) = 2.1031616948122442E+02_dp
      fit_coef(9, 17) = -2.0432151206915759E+02_dp
      fit_coef(10, 17) = 1.2775135126432421E+02_dp
      fit_coef(11, 17) = -4.6418220663487915E+01_dp
      fit_coef(12, 17) = 7.4470095221066224E+00_dp
      fit_coef(0, 18) = 8.1395124276193651E-03_dp
      fit_coef(1, 18) = 1.4295080919928922E-03_dp
      fit_coef(2, 18) = 1.9095038913922552E-01_dp
      fit_coef(3, 18) = -1.4103276835369478E+00_dp
      fit_coef(4, 18) = 1.0406506856630042E+01_dp
      fit_coef(5, 18) = -4.8337660891823411E+01_dp
      fit_coef(6, 18) = 1.4991827465085015E+02_dp
      fit_coef(7, 18) = -3.1676057422657112E+02_dp
      fit_coef(8, 18) = 4.5769739592136420E+02_dp
      fit_coef(9, 18) = -4.4502674091606121E+02_dp
      fit_coef(10, 18) = 2.7836644169317543E+02_dp
      fit_coef(11, 18) = -1.0115946637317019E+02_dp
      fit_coef(12, 18) = 1.6229500712924121E+01_dp
      fit_coef(0, 19) = 2.0717076029345759E-02_dp
      fit_coef(1, 19) = 2.8356072593395381E-03_dp
      fit_coef(2, 19) = 3.7697026585720034E-01_dp
      fit_coef(3, 19) = -2.8080135576835485E+00_dp
      fit_coef(4, 19) = 2.0474949195105673E+01_dp
      fit_coef(5, 19) = -9.5582095227962725E+01_dp
      fit_coef(6, 19) = 2.9724241292709326E+02_dp
      fit_coef(7, 19) = -6.2897377141618983E+02_dp
      fit_coef(8, 19) = 9.0952762860883638E+02_dp
      fit_coef(9, 19) = -8.8470183512709411E+02_dp
      fit_coef(10, 19) = 5.5351047071161520E+02_dp
      fit_coef(11, 19) = -2.0117820626571816E+02_dp
      fit_coef(12, 19) = 3.2279767274864902E+01_dp
      fit_coef(0, 20) = 4.8876460483642845E-02_dp
      fit_coef(1, 20) = 5.2174808763164821E-03_dp
      fit_coef(2, 20) = 6.9090200030553239E-01_dp
      fit_coef(3, 20) = -5.1804295291616596E+00_dp
      fit_coef(4, 20) = 3.7420371100773927E+01_dp
      fit_coef(5, 20) = -1.7533661861799825E+02_dp
      fit_coef(6, 20) = 5.4658265990827033E+02_dp
      fit_coef(7, 20) = -1.1581864472600616E+03_dp
      fit_coef(8, 20) = 1.6761108303901358E+03_dp
      fit_coef(9, 20) = -1.6311266679849095E+03_dp
      fit_coef(10, 20) = 1.0208160076430613E+03_dp
      fit_coef(11, 20) = -3.7110037560401321E+02_dp
      fit_coef(12, 20) = 5.9552981560922035E+01_dp
      fit_coef(0, 21) = 1.0862470653291458E-01_dp
      fit_coef(1, 21) = 8.9557210597009419E-03_dp
      fit_coef(2, 21) = 1.1819761188507145E+00_dp
      fit_coef(3, 21) = -8.9092575248149508E+00_dp
      fit_coef(4, 21) = 6.3862158642253107E+01_dp
      fit_coef(5, 21) = -3.0008394292824050E+02_dp
      fit_coef(6, 21) = 9.3746815308631687E+02_dp
      fit_coef(7, 21) = -1.9889858954576416E+03_dp
      fit_coef(8, 21) = 2.8805735949038526E+03_dp
      fit_coef(9, 21) = -2.8045287776281298E+03_dp
      fit_coef(10, 21) = 1.7556741628884376E+03_dp
      fit_coef(11, 21) = -6.3836601287630617E+02_dp
      fit_coef(12, 21) = 1.0245612268571801E+02_dp
      fit_coef(0, 22) = 2.3000575937054105E-01_dp
      fit_coef(1, 22) = 1.4334607293667794E-02_dp
      fit_coef(2, 22) = 1.8864545891903439E+00_dp
      fit_coef(3, 22) = -1.4280352788842627E+01_dp
      fit_coef(4, 22) = 1.0170803299679834E+02_dp
      fit_coef(5, 22) = -4.7898216088621695E+02_dp
      fit_coef(6, 22) = 1.4991893689445949E+03_dp
      fit_coef(7, 22) = -3.1844797903553772E+03_dp
      fit_coef(8, 22) = 4.6151825676230128E+03_dp
      fit_coef(9, 22) = -4.4952731839222251E+03_dp
      fit_coef(10, 22) = 2.8148739884496081E+03_dp
      fit_coef(11, 22) = -1.0236785820690926E+03_dp
      fit_coef(12, 22) = 1.6431836962570071E+02_dp
      fit_coef(0, 23) = 4.6825429870000890E-01_dp
      fit_coef(1, 23) = 2.1256443046165384E-02_dp
      fit_coef(2, 23) = 2.7905151165222755E+00_dp
      fit_coef(3, 23) = -2.1197457436160700E+01_dp
      fit_coef(4, 23) = 1.5017192005334016E+02_dp
      fit_coef(5, 23) = -7.0843773023664141E+02_dp
      fit_coef(6, 23) = 2.2210314991220407E+03_dp
      fit_coef(7, 23) = -4.7227286038286647E+03_dp
      fit_coef(8, 23) = 6.8489146048021994E+03_dp
      fit_coef(9, 23) = -6.6736306428096368E+03_dp
      fit_coef(10, 23) = 4.1800022473986091E+03_dp
      fit_coef(11, 23) = -1.5203921664836910E+03_dp
      fit_coef(12, 23) = 2.4407860279057451E+02_dp
      fit_coef(0, 24) = 9.2615015167881409E-01_dp
      fit_coef(1, 24) = 2.8851029357722529E-02_dp
      fit_coef(2, 24) = 3.7800720586882299E+00_dp
      fit_coef(3, 24) = -2.8790145745264784E+01_dp
      fit_coef(4, 24) = 2.0311830055185206E+02_dp
      fit_coef(5, 24) = -9.5941022550967875E+02_dp
      fit_coef(6, 24) = 3.0118661422758710E+03_dp
      fit_coef(7, 24) = -6.4099639194822603E+03_dp
      fit_coef(8, 24) = 9.3008315617510980E+03_dp
      fit_coef(9, 24) = -9.0659183663546519E+03_dp
      fit_coef(10, 24) = 5.6796728484798878E+03_dp
      fit_coef(11, 24) = -2.0661781228629306E+03_dp
      fit_coef(12, 24) = 3.3173195728147749E+02_dp
      fit_coef(0, 25) = 1.8197406882188483E+00_dp
      fit_coef(1, 25) = 3.5439117125219140E-02_dp
      fit_coef(2, 25) = 4.6374741775348376E+00_dp
      fit_coef(3, 25) = -3.5376154233726169E+01_dp
      fit_coef(4, 25) = 2.4894862911423743E+02_dp
      fit_coef(5, 25) = -1.1767196794716788E+03_dp
      fit_coef(6, 25) = 3.6971703023290502E+03_dp
      fit_coef(7, 25) = -7.8729723077281906E+03_dp
      fit_coef(8, 25) = 1.1427813147275856E+04_dp
      fit_coef(9, 25) = -1.1141769814002466E+04_dp
      fit_coef(10, 25) = 6.9812302395124789E+03_dp
      fit_coef(11, 25) = -2.5399272897012484E+03_dp
      fit_coef(12, 25) = 4.0782349862907824E+02_dp
      fit_coef(0, 26) = 3.8699714083339569E+00_dp
      fit_coef(1, 26) = 4.0430697875251408E-02_dp
      fit_coef(2, 26) = 5.2898404083820747E+00_dp
      fit_coef(3, 26) = -4.0359075619765512E+01_dp
      fit_coef(4, 26) = 2.8393314652682147E+02_dp
      fit_coef(5, 26) = -1.3421543359055597E+03_dp
      fit_coef(6, 26) = 4.2173994053467277E+03_dp
      fit_coef(7, 26) = -8.9814986542291372E+03_dp
      fit_coef(8, 26) = 1.3037566274562520E+04_dp
      fit_coef(9, 26) = -1.2711681344703467E+04_dp
      fit_coef(10, 26) = 7.9650982603751236E+03_dp
      fit_coef(11, 26) = -2.8979280791052265E+03_dp
      fit_coef(12, 26) = 4.6531132823912742E+02_dp

   END SUBROUTINE get_fit_coef_k13_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k14_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 4.0562245241676015E-03_dp
      fit_coef(1, 1) = 1.0969040206552460E-02_dp
      fit_coef(2, 1) = 1.4405191275649688E-02_dp
      fit_coef(3, 1) = 1.0649635686233358E-02_dp
      fit_coef(4, 1) = 3.2778746107621219E-03_dp
      fit_coef(5, 1) = -1.5147809086133006E-03_dp
      fit_coef(6, 1) = -1.7768378644379805E-03_dp
      fit_coef(7, 1) = -1.7040745050157837E-04_dp
      fit_coef(8, 1) = 4.3464032657914310E-04_dp
      fit_coef(9, 1) = 2.9772666014833774E-04_dp
      fit_coef(10, 1) = -3.9919353211016935E-04_dp
      fit_coef(11, 1) = 1.9832458267342404E-04_dp
      fit_coef(12, 1) = -4.3867764193471502E-05_dp
      fit_coef(0, 2) = 2.1853974229517521E-02_dp
      fit_coef(1, 2) = 5.8406086298247388E-02_dp
      fit_coef(2, 2) = 7.6152918249589960E-02_dp
      fit_coef(3, 2) = 5.5933948505163415E-02_dp
      fit_coef(4, 2) = 1.7021801576352837E-02_dp
      fit_coef(5, 2) = -8.0247279189159457E-03_dp
      fit_coef(6, 2) = -9.2606257683723776E-03_dp
      fit_coef(7, 2) = -8.4039046464910933E-04_dp
      fit_coef(8, 2) = 2.2920465602860878E-03_dp
      fit_coef(9, 2) = 1.5096114542523550E-03_dp
      fit_coef(10, 2) = -2.0536902896596300E-03_dp
      fit_coef(11, 2) = 1.0198430679936156E-03_dp
      fit_coef(12, 2) = -2.2535568541917208E-04_dp
      fit_coef(0, 3) = 5.6010245407281063E-02_dp
      fit_coef(1, 3) = 1.4622279778697206E-01_dp
      fit_coef(2, 3) = 1.8818537023914636E-01_dp
      fit_coef(3, 3) = 1.3664268008941158E-01_dp
      fit_coef(4, 3) = 4.0762241157569830E-02_dp
      fit_coef(5, 3) = -1.9875823087429162E-02_dp
      fit_coef(6, 3) = -2.2308862126762034E-02_dp
      fit_coef(7, 3) = -1.8247566176568901E-03_dp
      fit_coef(8, 3) = 5.6284318752787164E-03_dp
      fit_coef(9, 3) = 3.4697564866730192E-03_dp
      fit_coef(10, 3) = -4.8438240685035059E-03_dp
      fit_coef(11, 3) = 2.4030731708201193E-03_dp
      fit_coef(12, 3) = -5.2978518261630908E-04_dp
      fit_coef(0, 4) = 1.1092328992892893E-01_dp
      fit_coef(1, 4) = 2.7859680338891873E-01_dp
      fit_coef(2, 4) = 3.5169352998000031E-01_dp
      fit_coef(3, 4) = 2.5142478390455381E-01_dp
      fit_coef(4, 4) = 7.2990020531424102E-02_dp
      fit_coef(5, 4) = -3.7221206950039602E-02_dp
      fit_coef(6, 4) = -4.0210444834685513E-02_dp
      fit_coef(7, 4) = -2.7577403917926737E-03_dp
      fit_coef(8, 4) = 1.0362436751946298E-02_dp
      fit_coef(9, 4) = 5.8637315949354504E-03_dp
      fit_coef(10, 4) = -8.4785224237194489E-03_dp
      fit_coef(11, 4) = 4.1965838668479400E-03_dp
      fit_coef(12, 4) = -9.2121396633372532E-04_dp
      fit_coef(0, 5) = 1.9465685621238632E-01_dp
      fit_coef(1, 5) = 4.6152338108270852E-01_dp
      fit_coef(2, 5) = 5.6693470576329441E-01_dp
      fit_coef(3, 5) = 3.9807380880317028E-01_dp
      fit_coef(4, 5) = 1.1253953929433713E-01_dp
      fit_coef(5, 5) = -6.0165003179139975E-02_dp
      fit_coef(6, 5) = -6.2492236340615785E-02_dp
      fit_coef(7, 5) = -3.1516877242773393E-03_dp
      fit_coef(8, 5) = 1.7371050989592705E-02_dp
      fit_coef(9, 5) = 6.8885779131374857E-03_dp
      fit_coef(10, 5) = -1.1712499657108166E-02_dp
      fit_coef(11, 5) = 5.9558759087716033E-03_dp
      fit_coef(12, 5) = -1.3302080306339819E-03_dp
      fit_coef(0, 6) = 3.2087033276941335E-01_dp
      fit_coef(1, 6) = 7.0343427128207259E-01_dp
      fit_coef(2, 6) = 8.3155424687895485E-01_dp
      fit_coef(3, 6) = 5.7158103074002242E-01_dp
      fit_coef(4, 6) = 1.6033482878158128E-01_dp
      fit_coef(5, 6) = -8.5970624138830110E-02_dp
      fit_coef(6, 6) = -9.0709492543919121E-02_dp
      fit_coef(7, 6) = -4.6906764752169330E-03_dp
      fit_coef(8, 6) = 2.7054261624413843E-02_dp
      fit_coef(9, 6) = 1.4112534091895116E-02_dp
      fit_coef(10, 6) = -2.5828598206741545E-02_dp
      fit_coef(11, 6) = 1.4364571586594627E-02_dp
      fit_coef(12, 6) = -3.2730781292866435E-03_dp
      fit_coef(0, 7) = 5.1123719677822987E-01_dp
      fit_coef(1, 7) = 1.0160745695891158E+00_dp
      fit_coef(2, 7) = 1.1406376033091816E+00_dp
      fit_coef(3, 7) = 7.6072431018364839E-01_dp
      fit_coef(4, 7) = 2.1623548132157380E-01_dp
      fit_coef(5, 7) = -1.0385662211042419E-01_dp
      fit_coef(6, 7) = -1.2398076049296762E-01_dp
      fit_coef(7, 7) = -5.5258832288643421E-03_dp
      fit_coef(8, 7) = 5.4749667731452062E-03_dp
      fit_coef(9, 7) = 7.8756125030946533E-02_dp
      fit_coef(10, 7) = -8.1450501282583770E-02_dp
      fit_coef(11, 7) = 3.5582766973436604E-02_dp
      fit_coef(12, 7) = -6.2851848529127960E-03_dp
      fit_coef(0, 8) = 7.9862454256130755E-01_dp
      fit_coef(1, 8) = 1.4140917795956391E+00_dp
      fit_coef(2, 8) = 1.4887439734322301E+00_dp
      fit_coef(3, 8) = 9.4959856751799310E-01_dp
      fit_coef(4, 8) = 2.7310813416061092E-01_dp
      fit_coef(5, 8) = -1.0471167639057280E-01_dp
      fit_coef(6, 8) = -1.4117735675446755E-01_dp
      fit_coef(7, 8) = -2.8720200153060829E-02_dp
      fit_coef(8, 8) = 9.9988658384949844E-03_dp
      fit_coef(9, 8) = 4.4994288100147821E-02_dp
      fit_coef(10, 8) = 9.8296298107774047E-03_dp
      fit_coef(11, 8) = -3.6510317397394888E-02_dp
      fit_coef(12, 8) = 1.3048706811534378E-02_dp
      fit_coef(0, 9) = 1.2318540903074773E+00_dp
      fit_coef(1, 9) = 1.9126535284122799E+00_dp
      fit_coef(2, 9) = 1.8706375498146159E+00_dp
      fit_coef(3, 9) = 1.1237046593985318E+00_dp
      fit_coef(4, 9) = 3.2033236887832522E-01_dp
      fit_coef(5, 9) = -9.0264012269780411E-02_dp
      fit_coef(6, 9) = -1.2294629287092433E-01_dp
      fit_coef(7, 9) = -8.7215860683106217E-02_dp
      fit_coef(8, 9) = 1.0323515038061323E-01_dp
      fit_coef(9, 9) = -1.6774909880232289E-01_dp
      fit_coef(10, 9) = 2.4446069536217355E-01_dp
      fit_coef(11, 9) = -1.4978114145464044E-01_dp
      fit_coef(12, 9) = 3.2845034249025538E-02_dp
      fit_coef(0, 10) = 1.8833267667765345E+00_dp
      fit_coef(1, 10) = 2.5238922534191865E+00_dp
      fit_coef(2, 10) = 2.2803847971690017E+00_dp
      fit_coef(3, 10) = 1.2737526428801398E+00_dp
      fit_coef(4, 10) = 3.5105123224520401E-01_dp
      fit_coef(5, 10) = -6.2857188281637041E-02_dp
      fit_coef(6, 10) = -9.2210310288645597E-02_dp
      fit_coef(7, 10) = -9.0722133116613465E-02_dp
      fit_coef(8, 10) = 7.7837218665255559E-02_dp
      fit_coef(9, 10) = -1.6015129533752295E-01_dp
      fit_coef(10, 10) = 1.8283995908504222E-01_dp
      fit_coef(11, 10) = -6.6325970667744996E-02_dp
      fit_coef(12, 10) = 3.0820539987197176E-03_dp
      fit_coef(0, 11) = 2.8620729078671019E+00_dp
      fit_coef(1, 11) = 3.2531972126791606E+00_dp
      fit_coef(2, 11) = 2.7114279585314924E+00_dp
      fit_coef(3, 11) = 1.3964902193137803E+00_dp
      fit_coef(4, 11) = 3.6517639119142786E-01_dp
      fit_coef(5, 11) = -2.6597116040773625E-02_dp
      fit_coef(6, 11) = -5.9694762140569908E-02_dp
      fit_coef(7, 11) = -3.8758341163376515E-02_dp
      fit_coef(8, 11) = -5.1332218924274670E-02_dp
      fit_coef(9, 11) = 2.5767378796536407E-02_dp
      fit_coef(10, 11) = -4.9870158292428322E-02_dp
      fit_coef(11, 11) = 8.2792987390825273E-02_dp
      fit_coef(12, 11) = -3.1947335725901589E-02_dp
      fit_coef(0, 12) = 4.3396108413641770E+00_dp
      fit_coef(1, 12) = 4.0968755914648138E+00_dp
      fit_coef(2, 12) = 3.1592385203472810E+00_dp
      fit_coef(3, 12) = 1.4949577953660931E+00_dp
      fit_coef(4, 12) = 3.6703288840113468E-01_dp
      fit_coef(5, 12) = 1.2827934117818336E-02_dp
      fit_coef(6, 12) = -1.5966885082492444E-02_dp
      fit_coef(7, 12) = -4.1769471548348015E-03_dp
      fit_coef(8, 12) = -1.0829887623228632E-01_dp
      fit_coef(9, 12) = 9.0063525913958845E-02_dp
      fit_coef(10, 12) = -1.1718852204865197E-01_dp
      fit_coef(11, 12) = 1.0379976023591850E-01_dp
      fit_coef(12, 12) = -2.9513628601409953E-02_dp
      fit_coef(0, 13) = 6.6160518775236428E+00_dp
      fit_coef(1, 13) = 5.0450926097390525E+00_dp
      fit_coef(2, 13) = 3.6269712810913406E+00_dp
      fit_coef(3, 13) = 1.5787210925267479E+00_dp
      fit_coef(4, 13) = 3.6140765183946116E-01_dp
      fit_coef(5, 13) = 5.2660795733846105E-02_dp
      fit_coef(6, 13) = 3.6679818035805827E-02_dp
      fit_coef(7, 13) = 3.0817238462290870E-03_dp
      fit_coef(8, 13) = -8.5513616434734735E-02_dp
      fit_coef(9, 13) = 3.2480437301242795E-02_dp
      fit_coef(10, 13) = -4.3097395900725401E-02_dp
      fit_coef(11, 13) = 4.0285499498477921E-02_dp
      fit_coef(12, 13) = -9.0795169799425725E-03_dp
      fit_coef(0, 14) = 1.0378663082522372E+01_dp
      fit_coef(1, 14) = 6.1051748119732530E+00_dp
      fit_coef(2, 14) = 4.1393929927065685E+00_dp
      fit_coef(3, 14) = 1.6656362128427478E+00_dp
      fit_coef(4, 14) = 3.5175946384380213E-01_dp
      fit_coef(5, 14) = 9.2053766494033781E-02_dp
      fit_coef(6, 14) = 9.0977018487087655E-02_dp
      fit_coef(7, 14) = 9.2641513651703887E-03_dp
      fit_coef(8, 14) = -5.5705669679808022E-02_dp
      fit_coef(9, 14) = -2.2165364299652109E-02_dp
      fit_coef(10, 14) = 2.7050501330544607E-02_dp
      fit_coef(11, 14) = -9.8544712706888447E-03_dp
      fit_coef(12, 14) = 3.3393463606462780E-03_dp
      fit_coef(0, 15) = 1.0463061146787979E-02_dp
      fit_coef(1, 15) = 2.8219919737545576E-02_dp
      fit_coef(2, 15) = 3.6997970832918815E-02_dp
      fit_coef(3, 15) = 2.7310442581655728E-02_dp
      fit_coef(4, 15) = 8.3833807447771629E-03_dp
      fit_coef(5, 15) = -3.8927429936348136E-03_dp
      fit_coef(6, 15) = -4.5483815262270738E-03_dp
      fit_coef(7, 15) = -4.3064267744574074E-04_dp
      fit_coef(8, 15) = 1.1157764240345824E-03_dp
      fit_coef(9, 15) = 7.5712564333480544E-04_dp
      fit_coef(10, 15) = -1.0186473443121876E-03_dp
      fit_coef(11, 15) = 5.0602854051267157E-04_dp
      fit_coef(12, 15) = -1.1190592947238425E-04_dp
      fit_coef(0, 16) = 2.5471251770259992E-02_dp
      fit_coef(1, 16) = 6.7073734523031128E-02_dp
      fit_coef(2, 16) = 8.6686390148806164E-02_dp
      fit_coef(3, 16) = 6.3163726138437945E-02_dp
      fit_coef(4, 16) = 1.8954077181435283E-02_dp
      fit_coef(5, 16) = -9.1546883000010616E-03_dp
      fit_coef(6, 16) = -1.0357597642743745E-02_dp
      fit_coef(7, 16) = -8.7389112118697681E-04_dp
      fit_coef(8, 16) = 2.5999493392185090E-03_dp
      fit_coef(9, 16) = 1.6312506180299385E-03_dp
      fit_coef(10, 16) = -2.2608110574895506E-03_dp
      fit_coef(11, 16) = 1.1220090537728696E-03_dp
      fit_coef(12, 16) = -2.4758294844337260E-04_dp
      fit_coef(0, 17) = 4.3551694119749178E-02_dp
      fit_coef(1, 17) = 1.0925569190527253E-01_dp
      fit_coef(2, 17) = 1.3759390563727888E-01_dp
      fit_coef(3, 17) = 9.8026280776137426E-02_dp
      fit_coef(4, 17) = 2.8264281628926104E-02_dp
      fit_coef(5, 17) = -1.4570182357127179E-02_dp
      fit_coef(6, 17) = -1.5618085831617121E-02_dp
      fit_coef(7, 17) = -1.0386286860505113E-03_dp
      fit_coef(8, 17) = 4.0623648670812191E-03_dp
      fit_coef(9, 17) = 2.2393521384181845E-03_dp
      fit_coef(10, 17) = -3.2758370140128020E-03_dp
      fit_coef(11, 17) = 1.6218936224125669E-03_dp
      fit_coef(12, 17) = -3.5585220637753326E-04_dp
      fit_coef(0, 18) = 6.7585710312796060E-02_dp
      fit_coef(1, 18) = 1.5648129031717098E-01_dp
      fit_coef(2, 18) = 1.8950121838341139E-01_dp
      fit_coef(3, 18) = 1.3119730563134613E-01_dp
      fit_coef(4, 18) = 3.5948846183710491E-02_dp
      fit_coef(5, 18) = -2.0164654606078806E-02_dp
      fit_coef(6, 18) = -2.0055622446027831E-02_dp
      fit_coef(7, 18) = -6.9996789003826023E-04_dp
      fit_coef(8, 18) = 5.4059447836637305E-03_dp
      fit_coef(9, 18) = 2.3753238258316855E-03_dp
      fit_coef(10, 18) = -3.8252579297800944E-03_dp
      fit_coef(11, 18) = 1.8835729864486405E-03_dp
      fit_coef(12, 18) = -4.0977704401502611E-04_dp
      fit_coef(0, 19) = 1.0211824263133877E-01_dp
      fit_coef(1, 19) = 2.1075731039578760E-01_dp
      fit_coef(2, 19) = 2.4059034712033525E-01_dp
      fit_coef(3, 19) = 1.6134884680771011E-01_dp
      fit_coef(4, 19) = 4.3350649255430605E-02_dp
      fit_coef(5, 19) = -2.5429945851545656E-02_dp
      fit_coef(6, 19) = -2.4808194209897237E-02_dp
      fit_coef(7, 19) = -4.5717313943595194E-04_dp
      fit_coef(8, 19) = 9.7364264196386652E-03_dp
      fit_coef(9, 19) = -9.7009624494576196E-04_dp
      fit_coef(10, 19) = -2.9110931341643355E-03_dp
      fit_coef(11, 19) = 2.0554438206790630E-03_dp
      fit_coef(12, 19) = -5.3799541981056737E-04_dp
      fit_coef(0, 20) = 1.5387007707293665E-01_dp
      fit_coef(1, 20) = 2.7502702312847982E-01_dp
      fit_coef(2, 20) = 2.8778273930986864E-01_dp
      fit_coef(3, 20) = 1.8376769806731383E-01_dp
      fit_coef(4, 20) = 5.2375278020144772E-02_dp
      fit_coef(5, 20) = -2.4188675683434760E-02_dp
      fit_coef(6, 20) = -3.1878082611254147E-02_dp
      fit_coef(7, 20) = -2.0440498451142863E-03_dp
      fit_coef(8, 20) = 1.7743467113242066E-03_dp
      fit_coef(9, 20) = 2.9993028484423083E-02_dp
      fit_coef(10, 20) = -3.7493397222034638E-02_dp
      fit_coef(11, 20) = 1.9771570359339716E-02_dp
      fit_coef(12, 20) = -4.1700973428850313E-03_dp
      fit_coef(0, 21) = 2.3226494166758732E-01_dp
      fit_coef(1, 21) = 3.5273736189295163E-01_dp
      fit_coef(2, 21) = 3.2946524583838299E-01_dp
      fit_coef(3, 21) = 1.9174149746173408E-01_dp
      fit_coef(4, 21) = 5.8141887677596829E-02_dp
      fit_coef(5, 21) = -9.7825016932211791E-03_dp
      fit_coef(6, 21) = -3.0972124419450422E-02_dp
      fit_coef(7, 21) = -2.8877801446827007E-03_dp
      fit_coef(8, 21) = -3.6764584475810547E-02_dp
      fit_coef(9, 21) = 7.0959072404975718E-02_dp
      fit_coef(10, 21) = -3.5178291993696188E-02_dp
      fit_coef(11, 21) = -9.8088379200671617E-04_dp
      fit_coef(12, 21) = 3.7618284525752117E-03_dp
      fit_coef(0, 22) = 3.5052467022994122E-01_dp
      fit_coef(1, 22) = 4.4588641535984419E-01_dp
      fit_coef(2, 22) = 3.6589831022304731E-01_dp
      fit_coef(3, 22) = 1.8352935077226612E-01_dp
      fit_coef(4, 22) = 5.3697669457182857E-02_dp
      fit_coef(5, 22) = 7.4149644860502300E-03_dp
      fit_coef(6, 22) = 1.2817157120917000E-03_dp
      fit_coef(7, 22) = -5.0033726005623343E-02_dp
      fit_coef(8, 22) = 6.8349818056791897E-02_dp
      fit_coef(9, 22) = -1.7000421175238506E-01_dp
      fit_coef(10, 22) = 2.3036729897077399E-01_dp
      fit_coef(11, 22) = -1.3857914112066141E-01_dp
      fit_coef(12, 22) = 3.1478942190343064E-02_dp
      fit_coef(0, 23) = 5.2779750339667775E-01_dp
      fit_coef(1, 23) = 5.5332389225129053E-01_dp
      fit_coef(2, 23) = 3.9690230471520443E-01_dp
      fit_coef(3, 23) = 1.6312883734215308E-01_dp
      fit_coef(4, 23) = 3.9552650126779464E-02_dp
      fit_coef(5, 23) = 2.1230303801123101E-02_dp
      fit_coef(6, 23) = 2.9571154052233411E-02_dp
      fit_coef(7, 23) = -4.5494947806276824E-02_dp
      fit_coef(8, 23) = 6.6833147701713960E-02_dp
      fit_coef(9, 23) = -1.5861734029387040E-01_dp
      fit_coef(10, 23) = 1.3362195160908977E-01_dp
      fit_coef(11, 23) = -3.0781412640455877E-02_dp
      fit_coef(12, 23) = -4.0950711211332350E-03_dp
      fit_coef(0, 24) = 7.9293469396592242E-01_dp
      fit_coef(1, 24) = 6.7009858113053966E-01_dp
      fit_coef(2, 24) = 4.2145481861873185E-01_dp
      fit_coef(3, 24) = 1.3642048218640376E-01_dp
      fit_coef(4, 24) = 2.1922752252294477E-02_dp
      fit_coef(5, 24) = 3.2874744499942446E-02_dp
      fit_coef(6, 24) = 3.0227410874256264E-02_dp
      fit_coef(7, 24) = 3.7127938956607313E-02_dp
      fit_coef(8, 24) = -1.1080272293058828E-01_dp
      fit_coef(9, 24) = 1.5696651468543765E-01_dp
      fit_coef(10, 24) = -2.2158444445106978E-01_dp
      fit_coef(11, 24) = 1.6507639589038958E-01_dp
      fit_coef(12, 24) = -4.4931540554666242E-02_dp
      fit_coef(0, 25) = 1.1922412385770986E+00_dp
      fit_coef(1, 25) = 7.8786890944540700E-01_dp
      fit_coef(2, 25) = 4.3978600997964423E-01_dp
      fit_coef(3, 25) = 1.0955812452827124E-01_dp
      fit_coef(4, 25) = 7.1243000175333282E-03_dp
      fit_coef(5, 25) = 3.8586652562380884E-02_dp
      fit_coef(6, 25) = 3.7152220150416283E-02_dp
      fit_coef(7, 25) = 5.1558439446380777E-02_dp
      fit_coef(8, 25) = -1.1083564016687726E-01_dp
      fit_coef(9, 25) = 1.5290802796167138E-01_dp
      fit_coef(10, 25) = -1.7624313473286510E-01_dp
      fit_coef(11, 25) = 9.5327501416994614E-02_dp
      fit_coef(12, 25) = -1.7028005023958050E-02_dp
      fit_coef(0, 26) = 1.8099472486487456E+00_dp
      fit_coef(1, 26) = 8.9758639724651268E-01_dp
      fit_coef(2, 26) = 4.5619407404374801E-01_dp
      fit_coef(3, 26) = 8.8963048208584819E-02_dp
      fit_coef(4, 26) = -2.5842727308316919E-03_dp
      fit_coef(5, 26) = 3.9852736978965124E-02_dp
      fit_coef(6, 26) = 4.9637846366422017E-02_dp
      fit_coef(7, 26) = 1.7104716626654703E-02_dp
      fit_coef(8, 26) = -5.1404130029074147E-03_dp
      fit_coef(9, 26) = -1.8271352476782278E-02_dp
      fit_coef(10, 26) = 3.0152636564646153E-02_dp
      fit_coef(11, 26) = -4.0913238884136026E-02_dp
      fit_coef(12, 26) = 1.7209300973137127E-02_dp
      fit_coef(0, 27) = 2.8415520058389325E+00_dp
      fit_coef(1, 27) = 9.9858013711246008E-01_dp
      fit_coef(2, 27) = 4.8263553345863242E-01_dp
      fit_coef(3, 27) = 8.1362436826961093E-02_dp
      fit_coef(4, 27) = -8.1050564154871527E-03_dp
      fit_coef(5, 27) = 3.9674059022610446E-02_dp
      fit_coef(6, 27) = 5.4176205880144122E-02_dp
      fit_coef(7, 27) = 2.6408258437705550E-03_dp
      fit_coef(8, 27) = 3.6418581121829920E-02_dp
      fit_coef(9, 27) = -7.2590963148842197E-02_dp
      fit_coef(10, 27) = 9.0190833590845856E-02_dp
      fit_coef(11, 27) = -6.7850000466319185E-02_dp
      fit_coef(12, 27) = 1.8925208943418582E-02_dp
      fit_coef(0, 28) = 5.0338617688312999E+00_dp
      fit_coef(1, 28) = 1.1415572014123594E+00_dp
      fit_coef(2, 28) = 5.5834539482924794E-01_dp
      fit_coef(3, 28) = 9.8824807093161326E-02_dp
      fit_coef(4, 28) = -1.0740041647079411E-02_dp
      fit_coef(5, 28) = 3.9309450063690091E-02_dp
      fit_coef(6, 28) = 5.4929203884089291E-02_dp
      fit_coef(7, 28) = 1.2322580840140561E-02_dp
      fit_coef(8, 28) = 1.8323929549560497E-02_dp
      fit_coef(9, 28) = -2.9383377972038723E-02_dp
      fit_coef(10, 28) = 4.2820470862085666E-02_dp
      fit_coef(11, 28) = -2.8686857835350920E-02_dp
      fit_coef(12, 28) = 5.4565805278581868E-03_dp

   END SUBROUTINE get_fit_coef_k14_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k14_mi(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.1722808202408223E+03_dp
      fit_coef(0, 1) = 5.4467895546803416E-04_dp
      fit_coef(1, 1) = 1.5857771854287006E-03_dp
      fit_coef(2, 1) = 2.3073428933440593E-03_dp
      fit_coef(3, 1) = 2.0056780241026078E-03_dp
      fit_coef(4, 1) = 9.5252806530394737E-04_dp
      fit_coef(5, 1) = 4.5225480425719547E-05_dp
      fit_coef(6, 1) = -2.1217629613097376E-04_dp
      fit_coef(7, 1) = -5.0216895355491796E-05_dp
      fit_coef(8, 1) = 3.8287705248810476E-05_dp
      fit_coef(9, 1) = 7.2535948309621214E-05_dp
      fit_coef(10, 1) = -9.2252835912099736E-05_dp
      fit_coef(11, 1) = 4.3983534472290816E-05_dp
      fit_coef(12, 1) = -8.2237204784032551E-06_dp
      fit_coef(0, 2) = 3.0464535721029173E-03_dp
      fit_coef(1, 2) = 8.6149108928923676E-03_dp
      fit_coef(2, 2) = 1.2362165768730629E-02_dp
      fit_coef(3, 2) = 1.0619459122609453E-02_dp
      fit_coef(4, 2) = 4.9727394573446972E-03_dp
      fit_coef(5, 2) = 2.0040047900640227E-04_dp
      fit_coef(6, 2) = -1.1121250175071501E-03_dp
      fit_coef(7, 2) = -2.5253190910672047E-04_dp
      fit_coef(8, 2) = 2.0786942223541925E-04_dp
      fit_coef(9, 2) = 3.6647453679821519E-04_dp
      fit_coef(10, 2) = -4.7725200991211036E-04_dp
      fit_coef(11, 2) = 2.3017003059152598E-04_dp
      fit_coef(12, 2) = -4.3440212796964786E-05_dp
      fit_coef(0, 3) = 8.4179094658504636E-03_dp
      fit_coef(1, 3) = 2.2382112298387765E-02_dp
      fit_coef(2, 3) = 3.1254689785869226E-02_dp
      fit_coef(3, 3) = 2.6305836769925709E-02_dp
      fit_coef(4, 3) = 1.2019230418388028E-02_dp
      fit_coef(5, 3) = 3.2210723317093755E-04_dp
      fit_coef(6, 3) = -2.7056207855165560E-03_dp
      fit_coef(7, 3) = -5.6040022077270704E-04_dp
      fit_coef(8, 3) = 5.4121141114160339E-04_dp
      fit_coef(9, 3) = 8.1004167044932199E-04_dp
      fit_coef(10, 3) = -1.1112740690897719E-03_dp
      fit_coef(11, 3) = 5.4779892766252693E-04_dp
      fit_coef(12, 3) = -1.0512838684553561E-04_dp
      fit_coef(0, 4) = 1.8743574617056585E-02_dp
      fit_coef(1, 4) = 4.5210340778486291E-02_dp
      fit_coef(2, 4) = 6.0224002984973317E-02_dp
      fit_coef(3, 4) = 4.9094876694850402E-02_dp
      fit_coef(4, 4) = 2.1831501314762986E-02_dp
      fit_coef(5, 4) = 3.2839489812930983E-04_dp
      fit_coef(6, 4) = -5.0144416787406901E-03_dp
      fit_coef(7, 4) = -9.7843840294365748E-04_dp
      fit_coef(8, 4) = 1.0342060306701268E-03_dp
      fit_coef(9, 4) = 1.5910763103707255E-03_dp
      fit_coef(10, 4) = -2.2573252364893278E-03_dp
      fit_coef(11, 4) = 1.1445941503911742E-03_dp
      fit_coef(12, 4) = -2.2506139447365695E-04_dp
      fit_coef(0, 5) = 3.8243937154880572E-02_dp
      fit_coef(1, 5) = 8.1467829047918827E-02_dp
      fit_coef(2, 5) = 1.0129786520932144E-01_dp
      fit_coef(3, 5) = 7.8216152968887814E-02_dp
      fit_coef(4, 5) = 3.3350941393280653E-02_dp
      fit_coef(5, 5) = 5.2823533556460426E-04_dp
      fit_coef(6, 5) = -7.3749832466795265E-03_dp
      fit_coef(7, 5) = -1.5646554141616070E-03_dp
      fit_coef(8, 5) = 9.1282587135457208E-04_dp
      fit_coef(9, 5) = 2.8481769249969957E-03_dp
      fit_coef(10, 5) = -3.1292839864729078E-03_dp
      fit_coef(11, 5) = 1.3433122010827140E-03_dp
      fit_coef(12, 5) = -2.2591769996675308E-04_dp
      fit_coef(0, 6) = 7.4571737250878262E-02_dp
      fit_coef(1, 6) = 1.3811565890497141E-01_dp
      fit_coef(2, 6) = 1.5798315551672285E-01_dp
      fit_coef(3, 6) = 1.1260753901599996E-01_dp
      fit_coef(4, 6) = 4.4100222118730312E-02_dp
      fit_coef(5, 6) = 2.1773862891778066E-04_dp
      fit_coef(6, 6) = -8.3087822491922161E-03_dp
      fit_coef(7, 6) = -2.1043028241427880E-03_dp
      fit_coef(8, 6) = 1.8792681898943972E-03_dp
      fit_coef(9, 6) = 3.1606237510210035E-04_dp
      fit_coef(10, 6) = -5.4242874432157808E-04_dp
      fit_coef(11, 6) = 2.2538605364452114E-04_dp
      fit_coef(12, 6) = -5.0219363398529674E-05_dp
      fit_coef(0, 7) = 1.4098571554579828E-01_dp
      fit_coef(1, 7) = 2.2478771441370507E-01_dp
      fit_coef(2, 7) = 2.3494329165976915E-01_dp
      fit_coef(3, 7) = 1.5140646411833683E-01_dp
      fit_coef(4, 7) = 5.1480514391188670E-02_dp
      fit_coef(5, 7) = -1.9136093640078735E-03_dp
      fit_coef(6, 7) = -7.8257777325216969E-03_dp
      fit_coef(7, 7) = -9.9807084346081507E-04_dp
      fit_coef(8, 7) = 2.2963464957662357E-03_dp
      fit_coef(9, 7) = -1.9447302259471774E-03_dp
      fit_coef(10, 7) = 4.3017240890955453E-04_dp
      fit_coef(11, 7) = 4.7144332953281682E-04_dp
      fit_coef(12, 7) = -2.1997791857150025E-04_dp
      fit_coef(0, 8) = 2.5999858418452071E-01_dp
      fit_coef(1, 8) = 3.5368499241454637E-01_dp
      fit_coef(2, 8) = 3.3704821274394664E-01_dp
      fit_coef(3, 8) = 1.9318903274447125E-01_dp
      fit_coef(4, 8) = 5.3280944838378592E-02_dp
      fit_coef(5, 8) = -6.3861100195045095E-03_dp
      fit_coef(6, 8) = -5.8003810200407740E-03_dp
      fit_coef(7, 8) = 1.4075796032845642E-03_dp
      fit_coef(8, 8) = 2.4434092283670862E-03_dp
      fit_coef(9, 8) = -4.5616600603175456E-03_dp
      fit_coef(10, 8) = 1.8926356543785563E-03_dp
      fit_coef(11, 8) = 2.2201716818935917E-04_dp
      fit_coef(12, 8) = -2.2032045330927993E-04_dp
      fit_coef(0, 9) = 4.6927925194228570E-01_dp
      fit_coef(1, 9) = 5.3896976247930295E-01_dp
      fit_coef(2, 9) = 4.6838489154168161E-01_dp
      fit_coef(3, 9) = 2.3526229136436330E-01_dp
      fit_coef(4, 9) = 4.7454952056688732E-02_dp
      fit_coef(5, 9) = -1.2685981031239016E-02_dp
      fit_coef(6, 9) = -1.4744220873732451E-03_dp
      fit_coef(7, 9) = 4.2905695407292228E-03_dp
      fit_coef(8, 9) = 2.4711991371951058E-03_dp
      fit_coef(9, 9) = -8.6550441967466144E-03_dp
      fit_coef(10, 9) = 5.5989569069497076E-03_dp
      fit_coef(11, 9) = -1.3463966560621509E-03_dp
      fit_coef(12, 9) = 5.6917396181052615E-05_dp
      fit_coef(0, 10) = 8.3113745246756776E-01_dp
      fit_coef(1, 10) = 7.9501654600703975E-01_dp
      fit_coef(2, 10) = 6.3118071249392171E-01_dp
      fit_coef(3, 10) = 2.7364530371013518E-01_dp
      fit_coef(4, 10) = 3.2333401581906181E-02_dp
      fit_coef(5, 10) = -1.9422043063657940E-02_dp
      fit_coef(6, 10) = 6.0392759901066601E-03_dp
      fit_coef(7, 10) = 7.1646103557214402E-03_dp
      fit_coef(8, 10) = 1.3209371123133807E-03_dp
      fit_coef(9, 10) = -1.3329493616005855E-02_dp
      fit_coef(10, 10) = 1.1118391120775293E-02_dp
      fit_coef(11, 10) = -3.9802162789608605E-03_dp
      fit_coef(12, 10) = 5.4198727348417862E-04_dp
      fit_coef(0, 11) = 1.4482023274879392E+00_dp
      fit_coef(1, 11) = 1.1330048924955463E+00_dp
      fit_coef(2, 11) = 8.2485522761522434E-01_dp
      fit_coef(3, 11) = 3.0394925510152571E-01_dp
      fit_coef(4, 11) = 7.2638847029286640E-03_dp
      fit_coef(5, 11) = -2.4671422910119499E-02_dp
      fit_coef(6, 11) = 1.7366475605224755E-02_dp
      fit_coef(7, 11) = 9.3584526548286841E-03_dp
      fit_coef(8, 11) = -1.8731252836985278E-03_dp
      fit_coef(9, 11) = -1.7532663892185343E-02_dp
      fit_coef(10, 11) = 1.7972283371407362E-02_dp
      fit_coef(11, 11) = -7.6180802914050756E-03_dp
      fit_coef(12, 11) = 1.2561263859800901E-03_dp
      fit_coef(0, 12) = 2.4920360801417596E+00_dp
      fit_coef(1, 12) = 1.5562052067140830E+00_dp
      fit_coef(2, 12) = 1.0460115631746023E+00_dp
      fit_coef(3, 12) = 3.2329508963929127E-01_dp
      fit_coef(4, 12) = -2.6843427897215052E-02_dp
      fit_coef(5, 12) = -2.6728379775572533E-02_dp
      fit_coef(6, 12) = 3.2440057594826148E-02_dp
      fit_coef(7, 12) = 1.0203339477948091E-02_dp
      fit_coef(8, 12) = -7.5116389543029823E-03_dp
      fit_coef(9, 12) = -2.0251551518009794E-02_dp
      fit_coef(10, 12) = 2.5495522854135564E-02_dp
      fit_coef(11, 12) = -1.2104045379629679E-02_dp
      fit_coef(12, 12) = 2.2002184184182066E-03_dp
      fit_coef(0, 13) = 4.2695725351678870E+00_dp
      fit_coef(1, 13) = 2.0569142521947663E+00_dp
      fit_coef(2, 13) = 1.2906701222000276E+00_dp
      fit_coef(3, 13) = 3.3274287972587252E-01_dp
      fit_coef(4, 13) = -6.7760240413352227E-02_dp
      fit_coef(5, 13) = -2.5210323779882499E-02_dp
      fit_coef(6, 13) = 5.0470840558373135E-02_dp
      fit_coef(7, 13) = 9.5160539695092316E-03_dp
      fit_coef(8, 13) = -1.5360600121122914E-02_dp
      fit_coef(9, 13) = -2.1060084216392275E-02_dp
      fit_coef(10, 13) = 3.3061137892312159E-02_dp
      fit_coef(11, 13) = -1.7135171327692631E-02_dp
      fit_coef(12, 13) = 3.3198097857474905E-03_dp
      fit_coef(0, 14) = 7.4706366451203330E+00_dp
      fit_coef(1, 14) = 2.6255475168075386E+00_dp
      fit_coef(2, 14) = 1.5617856992687866E+00_dp
      fit_coef(3, 14) = 3.3898821868503676E-01_dp
      fit_coef(4, 14) = -1.1349423621812627E-01_dp
      fit_coef(5, 14) = -2.1752415544934678E-02_dp
      fit_coef(6, 14) = 7.0659941513993574E-02_dp
      fit_coef(7, 14) = 7.9261305614328585E-03_dp
      fit_coef(8, 14) = -2.4678202093857821E-02_dp
      fit_coef(9, 14) = -2.0791238512333610E-02_dp
      fit_coef(10, 14) = 4.0834605352408383E-02_dp
      fit_coef(11, 14) = -2.2601414883140436E-02_dp
      fit_coef(12, 14) = 4.5678225526306271E-03_dp
      fit_coef(0, 15) = 1.4166144548081244E-03_dp
      fit_coef(1, 15) = 4.0988329863504650E-03_dp
      fit_coef(2, 15) = 5.9448187504450257E-03_dp
      fit_coef(3, 15) = 5.1533064863913264E-03_dp
      fit_coef(4, 15) = 2.4394418491756548E-03_dp
      fit_coef(5, 15) = 1.1163214246871172E-04_dp
      fit_coef(6, 15) = -5.4395417380967090E-04_dp
      fit_coef(7, 15) = -1.2744617101115058E-04_dp
      fit_coef(8, 15) = 9.9020248580942340E-05_dp
      fit_coef(9, 15) = 1.8429376434903194E-04_dp
      fit_coef(10, 15) = -2.3572941138030507E-04_dp
      fit_coef(11, 15) = 1.1270257688728590E-04_dp
      fit_coef(12, 15) = -2.1120306082703953E-05_dp
      fit_coef(0, 16) = 3.7180060821504178E-03_dp
      fit_coef(1, 16) = 1.0136250443655282E-02_dp
      fit_coef(2, 16) = 1.4300959605183702E-02_dp
      fit_coef(3, 16) = 1.2111006697402966E-02_dp
      fit_coef(4, 16) = 5.5715231524265116E-03_dp
      fit_coef(5, 16) = 1.7441962760050931E-04_dp
      fit_coef(6, 16) = -1.2509495746884976E-03_dp
      fit_coef(7, 16) = -2.6932600471792030E-04_dp
      fit_coef(8, 16) = 2.4202856508010202E-04_dp
      fit_coef(9, 16) = 3.9544082532111499E-04_dp
      fit_coef(10, 16) = -5.2934481355649469E-04_dp
      fit_coef(11, 16) = 2.5870684688921450E-04_dp
      fit_coef(12, 16) = -4.9352748645601205E-05_dp
      fit_coef(0, 17) = 7.3573792852102415E-03_dp
      fit_coef(1, 17) = 1.7769053739351975E-02_dp
      fit_coef(2, 17) = 2.3686656634361428E-02_dp
      fit_coef(3, 17) = 1.9285714809028087E-02_dp
      fit_coef(4, 17) = 8.4872745881098346E-03_dp
      fit_coef(5, 17) = 3.4696339797671155E-05_dp
      fit_coef(6, 17) = -1.9464426054708341E-03_dp
      fit_coef(7, 17) = -3.4174035665184455E-04_dp
      fit_coef(8, 17) = 4.7554753790861559E-04_dp
      fit_coef(9, 17) = 4.2925599708472817E-04_dp
      fit_coef(10, 17) = -7.0959124228041825E-04_dp
      fit_coef(11, 17) = 3.7293802063667668E-04_dp
      fit_coef(12, 17) = -7.4872640261499496E-05_dp
      fit_coef(0, 18) = 1.3977177746865988E-02_dp
      fit_coef(1, 18) = 2.8606023538180633E-02_dp
      fit_coef(2, 18) = 3.4571612966699722E-02_dp
      fit_coef(3, 18) = 2.6149166759159050E-02_dp
      fit_coef(4, 18) = 1.0984913521363774E-02_dp
      fit_coef(5, 18) = 7.7055787642742803E-05_dp
      fit_coef(6, 18) = -2.5922755331103330E-03_dp
      fit_coef(7, 18) = -4.8489977776779012E-04_dp
      fit_coef(8, 18) = 1.5912901418092180E-04_dp
      fit_coef(9, 18) = 1.5881868365226204E-03_dp
      fit_coef(10, 18) = -1.8151464311114490E-03_dp
      fit_coef(11, 18) = 8.5592821993186205E-04_dp
      fit_coef(12, 18) = -1.5937383356627851E-04_dp
      fit_coef(0, 19) = 2.6264101814506675E-02_dp
      fit_coef(1, 19) = 4.5069161923578932E-02_dp
      fit_coef(2, 19) = 4.8169189713700107E-02_dp
      fit_coef(3, 19) = 3.1904492882075079E-02_dp
      fit_coef(4, 19) = 1.1603766261612461E-02_dp
      fit_coef(5, 19) = 1.8449908333587032E-04_dp
      fit_coef(6, 19) = -1.7877726122463406E-03_dp
      fit_coef(7, 19) = -7.5764206916605059E-04_dp
      fit_coef(8, 19) = 2.8343231172872694E-04_dp
      fit_coef(9, 19) = -5.9914897854961122E-04_dp
      fit_coef(10, 19) = 1.3517055951379382E-03_dp
      fit_coef(11, 19) = -9.4046015820634316E-04_dp
      fit_coef(12, 19) = 2.2350581726021887E-04_dp
      fit_coef(0, 20) = 4.8542207101173400E-02_dp
      fit_coef(1, 20) = 6.9832820033889687E-02_dp
      fit_coef(2, 20) = 6.5995740079255238E-02_dp
      fit_coef(3, 20) = 3.6752596227584124E-02_dp
      fit_coef(4, 20) = 9.4571244340008204E-03_dp
      fit_coef(5, 20) = -1.0484402051494153E-03_dp
      fit_coef(6, 20) = -1.5752390773495777E-04_dp
      fit_coef(7, 20) = 1.0876735653687303E-04_dp
      fit_coef(8, 20) = 1.0073470936970389E-03_dp
      fit_coef(9, 20) = -3.0335098186363291E-03_dp
      fit_coef(10, 20) = 2.1855940575857963E-03_dp
      fit_coef(11, 20) = -4.5502298284867644E-04_dp
      fit_coef(12, 20) = -3.4790305788945764E-05_dp
      fit_coef(0, 21) = 8.7942456610255598E-02_dp
      fit_coef(1, 21) = 1.0555324383884353E-01_dp
      fit_coef(2, 21) = 8.8752093517154698E-02_dp
      fit_coef(3, 21) = 4.0620455348294761E-02_dp
      fit_coef(4, 21) = 4.9418858666300082E-03_dp
      fit_coef(5, 21) = -3.3032510647944919E-03_dp
      fit_coef(6, 21) = 1.1546343697416441E-03_dp
      fit_coef(7, 21) = 1.9549734572505361E-03_dp
      fit_coef(8, 21) = 5.6213884018718063E-05_dp
      fit_coef(9, 21) = -1.9404012109824945E-03_dp
      fit_coef(10, 21) = 5.3863092062487572E-04_dp
      fit_coef(11, 21) = 3.8099987718174182E-04_dp
      fit_coef(12, 21) = -1.6063218378014254E-04_dp
      fit_coef(0, 22) = 1.5620880999122425E-01_dp
      fit_coef(1, 22) = 1.5462937480065686E-01_dp
      fit_coef(2, 22) = 1.1615728166593994E-01_dp
      fit_coef(3, 22) = 4.2496601438037490E-02_dp
      fit_coef(4, 22) = -1.6862132065151886E-03_dp
      fit_coef(5, 22) = -5.5527712863935753E-03_dp
      fit_coef(6, 22) = 3.0288231280584822E-03_dp
      fit_coef(7, 22) = 2.7239270206899919E-03_dp
      fit_coef(8, 22) = 2.3923207067007299E-04_dp
      fit_coef(9, 22) = -3.4467494238847617E-03_dp
      fit_coef(10, 22) = 2.6143123206984532E-03_dp
      fit_coef(11, 22) = -9.5623661147149601E-04_dp
      fit_coef(12, 22) = 1.5884529743136178E-04_dp
      fit_coef(0, 23) = 2.7248883355523357E-01_dp
      fit_coef(1, 23) = 2.1839288076669880E-01_dp
      fit_coef(2, 23) = 1.4690508945871661E-01_dp
      fit_coef(3, 23) = 4.0972405769053262E-02_dp
      fit_coef(4, 23) = -1.0255836491498508E-02_dp
      fit_coef(5, 23) = -6.8103039861294205E-03_dp
      fit_coef(6, 23) = 5.7774256001662866E-03_dp
      fit_coef(7, 23) = 2.9716307657906005E-03_dp
      fit_coef(8, 23) = -3.8415214985403008E-04_dp
      fit_coef(9, 23) = -4.5597005257117969E-03_dp
      fit_coef(10, 23) = 4.6944833426410494E-03_dp
      fit_coef(11, 23) = -2.1231572186145029E-03_dp
      fit_coef(12, 23) = 3.8196379187575797E-04_dp
      fit_coef(0, 24) = 4.6789642997790454E-01_dp
      fit_coef(1, 24) = 2.9564830559342764E-01_dp
      fit_coef(2, 24) = 1.7859318172912267E-01_dp
      fit_coef(3, 24) = 3.5029171279950352E-02_dp
      fit_coef(4, 24) = -2.0106998533844407E-02_dp
      fit_coef(5, 24) = -6.3290051255165099E-03_dp
      fit_coef(6, 24) = 9.3627910136797072E-03_dp
      fit_coef(7, 24) = 2.6599389408438201E-03_dp
      fit_coef(8, 24) = -2.0615332759583711E-03_dp
      fit_coef(9, 24) = -4.6173560793899586E-03_dp
      fit_coef(10, 24) = 6.2747096995400117E-03_dp
      fit_coef(11, 24) = -3.1434465968595239E-03_dp
      fit_coef(12, 24) = 5.9448478612450029E-04_dp
      fit_coef(0, 25) = 7.9390404143601678E-01_dp
      fit_coef(1, 25) = 3.8094350576553354E-01_dp
      fit_coef(2, 25) = 2.0815634677459935E-01_dp
      fit_coef(3, 25) = 2.5054416108850112E-02_dp
      fit_coef(4, 25) = -2.9877418297843340E-02_dp
      fit_coef(5, 25) = -3.8648566065354233E-03_dp
      fit_coef(6, 25) = 1.3281740014574807E-02_dp
      fit_coef(7, 25) = 1.6073039898392629E-03_dp
      fit_coef(8, 25) = -4.3975022763734321E-03_dp
      fit_coef(9, 25) = -3.6092937145507167E-03_dp
      fit_coef(10, 25) = 7.3177671091434133E-03_dp
      fit_coef(11, 25) = -4.1069025396153367E-03_dp
      fit_coef(12, 25) = 8.3442522049861077E-04_dp
      fit_coef(0, 26) = 1.3421299372909123E+00_dp
      fit_coef(1, 26) = 4.6409313440446603E-01_dp
      fit_coef(2, 26) = 2.3339441406477554E-01_dp
      fit_coef(3, 26) = 1.3785053308326188E-02_dp
      fit_coef(4, 26) = -3.7941211960007708E-02_dp
      fit_coef(5, 26) = -1.7132484132012456E-04_dp
      fit_coef(6, 26) = 1.6722419910236191E-02_dp
      fit_coef(7, 26) = 4.5287606479249242E-05_dp
      fit_coef(8, 26) = -6.8396430449082519E-03_dp
      fit_coef(9, 26) = -1.7433447611092249E-03_dp
      fit_coef(10, 26) = 7.6094558147861901E-03_dp
      fit_coef(11, 26) = -4.8080400030221532E-03_dp
      fit_coef(12, 26) = 1.0435040492483968E-03_dp
      fit_coef(0, 27) = 2.3127647897779702E+00_dp
      fit_coef(1, 27) = 5.3494805569190185E-01_dp
      fit_coef(2, 27) = 2.5603491110440707E-01_dp
      fit_coef(3, 27) = 6.2124469115487774E-03_dp
      fit_coef(4, 27) = -4.3500571617098645E-02_dp
      fit_coef(5, 27) = 2.9267102310614004E-03_dp
      fit_coef(6, 27) = 1.9176839114778386E-02_dp
      fit_coef(7, 27) = -1.3114239457504276E-03_dp
      fit_coef(8, 27) = -8.7196471933934649E-03_dp
      fit_coef(9, 27) = -1.4236987601663990E-05_dp
      fit_coef(10, 27) = 7.5377946011018363E-03_dp
      fit_coef(11, 27) = -5.2272346094548762E-03_dp
      fit_coef(12, 27) = 1.1857786912012739E-03_dp
      fit_coef(0, 28) = 4.4310556767314351E+00_dp
      fit_coef(1, 28) = 6.0851379611857248E-01_dp
      fit_coef(2, 28) = 2.9194514151733547E-01_dp
      fit_coef(3, 28) = 8.6001943701959644E-03_dp
      fit_coef(4, 28) = -4.8295972382001021E-02_dp
      fit_coef(5, 28) = 3.4189484756573089E-03_dp
      fit_coef(6, 28) = 2.1364281048084820E-02_dp
      fit_coef(7, 28) = -1.6460148108292097E-03_dp
      fit_coef(8, 28) = -9.8467342034871362E-03_dp
      fit_coef(9, 28) = 1.9090555066223413E-04_dp
      fit_coef(10, 28) = 8.3032333199491645E-03_dp
      fit_coef(11, 28) = -5.8101678643582262E-03_dp
      fit_coef(12, 28) = 1.3244384992621098E-03_dp

   END SUBROUTINE get_fit_coef_k14_mi

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k14_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 9.4708858956711737E+01_dp
      U_b = 6.6910000000000000E+05_dp
      fit_coef(0, 1) = 5.7644508135195325E-06_dp
      fit_coef(1, 1) = 8.8059899546543633E-06_dp
      fit_coef(2, 1) = 9.8246367093007174E-04_dp
      fit_coef(3, 1) = -5.7370459929060305E-03_dp
      fit_coef(4, 1) = 5.4909691175756163E-02_dp
      fit_coef(5, 1) = -2.3201875112604708E-01_dp
      fit_coef(6, 1) = 7.0724150950356890E-01_dp
      fit_coef(7, 1) = -1.4925604266418715E+00_dp
      fit_coef(8, 1) = 2.1603279291789366E+00_dp
      fit_coef(9, 1) = -2.1043657951755335E+00_dp
      fit_coef(10, 1) = 1.3179096181212082E+00_dp
      fit_coef(11, 1) = -4.7916040323369041E-01_dp
      fit_coef(12, 1) = 7.6853944770310323E-02_dp
      fit_coef(0, 2) = 6.3482217934141788E-05_dp
      fit_coef(1, 2) = 5.6715849587820812E-05_dp
      fit_coef(2, 2) = 6.2902935744076748E-03_dp
      fit_coef(3, 2) = -3.8825646732371183E-02_dp
      fit_coef(4, 2) = 3.5215367034308359E-01_dp
      fit_coef(5, 2) = -1.5068239503633283E+00_dp
      fit_coef(6, 2) = 4.5854281378525714E+00_dp
      fit_coef(7, 2) = -9.6460393022334507E+00_dp
      fit_coef(8, 2) = 1.3925086681528807E+01_dp
      fit_coef(9, 2) = -1.3538171863016146E+01_dp
      fit_coef(10, 2) = 8.4668043138765512E+00_dp
      fit_coef(11, 2) = -3.0752413149585256E+00_dp
      fit_coef(12, 2) = 4.9288840241270471E-01_dp
      fit_coef(0, 3) = 3.3613588297925877E-04_dp
      fit_coef(1, 3) = 2.0810718441537009E-04_dp
      fit_coef(2, 3) = 2.2496352429605227E-02_dp
      fit_coef(3, 3) = -1.4508680125471160E-01_dp
      fit_coef(4, 3) = 1.2240878955617038E+00_dp
      fit_coef(5, 3) = -5.4269306653344804E+00_dp
      fit_coef(6, 3) = 1.6783138844933998E+01_dp
      fit_coef(7, 3) = -3.5628143752717627E+01_dp
      fit_coef(8, 3) = 5.1720421721162715E+01_dp
      fit_coef(9, 3) = -5.0458114469889409E+01_dp
      fit_coef(10, 3) = 3.1625635386312588E+01_dp
      fit_coef(11, 3) = -1.1502746293672999E+01_dp
      fit_coef(12, 3) = 1.8452512921088244E+00_dp
      fit_coef(0, 4) = 1.2977058344000711E-03_dp
      fit_coef(1, 4) = 6.0051493044492425E-04_dp
      fit_coef(2, 4) = 6.4379175573630529E-02_dp
      fit_coef(3, 4) = -4.2944616668718744E-01_dp
      fit_coef(4, 4) = 3.4645304257752727E+00_dp
      fit_coef(5, 4) = -1.5790291716256178E+01_dp
      fit_coef(6, 4) = 4.9078910409345937E+01_dp
      fit_coef(7, 4) = -1.0406732095355643E+02_dp
      fit_coef(8, 4) = 1.5069319179330122E+02_dp
      fit_coef(9, 4) = -1.4664743032705482E+02_dp
      fit_coef(10, 4) = 9.1717900592693681E+01_dp
      fit_coef(11, 4) = -3.3301586691674665E+01_dp
      fit_coef(12, 4) = 5.3348395481918729E+00_dp
      fit_coef(0, 5) = 4.1564666328203293E-03_dp
      fit_coef(1, 5) = 1.5024272029496644E-03_dp
      fit_coef(2, 5) = 1.5961690430876282E-01_dp
      fit_coef(3, 5) = -1.0856048657785689E+00_dp
      fit_coef(4, 5) = 8.4819759793848739E+00_dp
      fit_coef(5, 5) = -3.9236853866944621E+01_dp
      fit_coef(6, 5) = 1.2225270389761864E+02_dp
      fit_coef(7, 5) = -2.5932080520941867E+02_dp
      fit_coef(8, 5) = 3.7560378396472731E+02_dp
      fit_coef(9, 5) = -3.6564399683386677E+02_dp
      fit_coef(10, 5) = 2.2876919730710458E+02_dp
      fit_coef(11, 5) = -8.3091860520998807E+01_dp
      fit_coef(12, 5) = 1.3315159702194723E+01_dp
      fit_coef(0, 6) = 1.1728185806559138E-02_dp
      fit_coef(1, 6) = 3.3854671547705814E-03_dp
      fit_coef(2, 6) = 3.5708079041567914E-01_dp
      fit_coef(3, 6) = -2.4616715243590788E+00_dp
      fit_coef(4, 6) = 1.8784735712246597E+01_dp
      fit_coef(5, 6) = -8.7758450553219319E+01_dp
      fit_coef(6, 6) = 2.7430420638340945E+02_dp
      fit_coef(7, 6) = -5.8284545441781790E+02_dp
      fit_coef(8, 6) = 8.4515131830014900E+02_dp
      fit_coef(9, 6) = -8.2334556633746092E+02_dp
      fit_coef(10, 6) = 5.1537716566398308E+02_dp
      fit_coef(11, 6) = -1.8724712324686436E+02_dp
      fit_coef(12, 6) = 3.0011313266344054E+01_dp
      fit_coef(0, 7) = 3.0150116177462289E-02_dp
      fit_coef(1, 7) = 7.0271320856721093E-03_dp
      fit_coef(2, 7) = 7.3674268137292198E-01_dp
      fit_coef(3, 7) = -5.1296976804148287E+00_dp
      fit_coef(4, 7) = 3.8429636567626652E+01_dp
      fit_coef(5, 7) = -1.8076150868535993E+02_dp
      fit_coef(6, 7) = 5.6676422159052879E+02_dp
      fit_coef(7, 7) = -1.2064911489924480E+03_dp
      fit_coef(8, 7) = 1.7514726329556534E+03_dp
      fit_coef(9, 7) = -1.7075068700281913E+03_dp
      fit_coef(10, 7) = 1.0693144394521294E+03_dp
      fit_coef(11, 7) = -3.8862358061084939E+02_dp
      fit_coef(12, 7) = 6.2300629367049794E+01_dp
      fit_coef(0, 8) = 7.2115669678109895E-02_dp
      fit_coef(1, 8) = 1.3616506829774023E-02_dp
      fit_coef(2, 8) = 1.4203294247117737E+00_dp
      fit_coef(3, 8) = -9.9648479573798241E+00_dp
      fit_coef(4, 8) = 7.3547304516631314E+01_dp
      fit_coef(5, 8) = -3.4763570685939845E+02_dp
      fit_coef(6, 8) = 1.0931366005010798E+03_dp
      fit_coef(7, 8) = -2.3311962044304441E+03_dp
      fit_coef(8, 8) = 3.3880559926822129E+03_dp
      fit_coef(9, 8) = -3.3054142191968185E+03_dp
      fit_coef(10, 8) = 2.0709975619174938E+03_dp
      fit_coef(11, 8) = -7.5292079976886191E+02_dp
      fit_coef(12, 8) = 1.2073055461067925E+02_dp
      fit_coef(0, 9) = 1.6278592399460171E-01_dp
      fit_coef(1, 9) = 2.4823462702639722E-02_dp
      fit_coef(2, 9) = 2.5779844651162875E+00_dp
      fit_coef(3, 9) = -1.8195793409897270E+01_dp
      fit_coef(4, 9) = 1.3264400733833980E+02_dp
      fit_coef(5, 9) = -6.2921784847028903E+02_dp
      fit_coef(6, 9) = 1.9837274033715200E+03_dp
      fit_coef(7, 9) = -4.2376252580133068E+03_dp
      fit_coef(8, 9) = 6.1654716911349278E+03_dp
      fit_coef(9, 9) = -6.0193635625448478E+03_dp
      fit_coef(10, 9) = 3.7732320110600213E+03_dp
      fit_coef(11, 9) = -1.3722389889916640E+03_dp
      fit_coef(12, 9) = 2.2009155940297515E+02_dp
      fit_coef(0, 10) = 3.5032689539851275E-01_dp
      fit_coef(1, 10) = 4.2737381729215401E-02_dp
      fit_coef(2, 10) = 4.4215419106194087E+00_dp
      fit_coef(3, 10) = -3.1358951004261929E+01_dp
      fit_coef(4, 10) = 2.2622710383931096E+02_dp
      fit_coef(5, 10) = -1.0760101708193520E+03_dp
      fit_coef(6, 10) = 3.4001822298289899E+03_dp
      fit_coef(7, 10) = -7.2748054353872876E+03_dp
      fit_coef(8, 10) = 1.0595181287227540E+04_dp
      fit_coef(9, 10) = -1.0351115266133067E+04_dp
      fit_coef(10, 10) = 6.4916005963247544E+03_dp
      fit_coef(11, 10) = -2.3616229307375006E+03_dp
      fit_coef(12, 10) = 3.7886792111335615E+02_dp
      fit_coef(0, 11) = 7.2447713579898398E-01_dp
      fit_coef(1, 11) = 6.9519327565965100E-02_dp
      fit_coef(2, 11) = 7.1687589104301583E+00_dp
      fit_coef(3, 11) = -5.1041793387291833E+01_dp
      fit_coef(4, 11) = 3.6499013684520196E+02_dp
      fit_coef(5, 11) = -1.7394623733769292E+03_dp
      fit_coef(6, 11) = 5.5078086466424802E+03_dp
      fit_coef(7, 11) = -1.1800797526080492E+04_dp
      fit_coef(8, 11) = 1.7203067323847517E+04_dp
      fit_coef(9, 11) = -1.6817396268950357E+04_dp
      fit_coef(10, 11) = 1.0551469125979243E+04_dp
      fit_coef(11, 11) = -3.8397827922054616E+03_dp
      fit_coef(12, 11) = 6.1614478300189171E+02_dp
      fit_coef(0, 12) = 1.4502467390425808E+00_dp
      fit_coef(1, 12) = 1.0661107753602478E-01_dp
      fit_coef(2, 12) = 1.0963185897291746E+01_dp
      fit_coef(3, 12) = -7.8300692460550920E+01_dp
      fit_coef(4, 12) = 5.5584223144558530E+02_dp
      fit_coef(5, 12) = -2.6528572306059414E+03_dp
      fit_coef(6, 12) = 8.4143371492021706E+03_dp
      fit_coef(7, 12) = -1.8050418749368942E+04_dp
      fit_coef(8, 12) = 2.6335573301758028E+04_dp
      fit_coef(9, 12) = -2.5759713754102195E+04_dp
      fit_coef(10, 12) = 1.6168369713625771E+04_dp
      fit_coef(11, 12) = -5.8854889464635708E+03_dp
      fit_coef(12, 12) = 9.4460223869839535E+02_dp
      fit_coef(0, 13) = 2.8395717023886369E+00_dp
      fit_coef(1, 13) = 1.5366537352589704E-01_dp
      fit_coef(2, 13) = 1.5767280021363437E+01_dp
      fit_coef(3, 13) = -1.1287621706718255E+02_dp
      fit_coef(4, 13) = 7.9673183153700211E+02_dp
      fit_coef(5, 13) = -3.8063687624389750E+03_dp
      fit_coef(6, 13) = 1.2089346236379877E+04_dp
      fit_coef(7, 13) = -2.5959938226844133E+04_dp
      fit_coef(8, 13) = 3.7901426581098414E+04_dp
      fit_coef(9, 13) = -3.7090020172809200E+04_dp
      fit_coef(10, 13) = 2.3287606333954678E+04_dp
      fit_coef(11, 13) = -8.4789859745242975E+03_dp
      fit_coef(12, 13) = 1.3610885265396084E+03_dp
      fit_coef(0, 14) = 5.5945431735079278E+00_dp
      fit_coef(1, 14) = 2.0856787151109613E-01_dp
      fit_coef(2, 14) = 2.1367891445332752E+01_dp
      fit_coef(3, 14) = -1.5321402782127791E+02_dp
      fit_coef(4, 14) = 1.0771822980689246E+03_dp
      fit_coef(5, 14) = -5.1495462635748599E+03_dp
      fit_coef(6, 14) = 1.6370767484953010E+04_dp
      fit_coef(7, 14) = -3.5178430776828289E+04_dp
      fit_coef(8, 14) = 5.1385383840132439E+04_dp
      fit_coef(9, 14) = -5.0302167209417865E+04_dp
      fit_coef(10, 14) = 3.1590534978260897E+04_dp
      fit_coef(11, 14) = -1.1504033487227012E+04_dp
      fit_coef(12, 14) = 1.8469177824866308E+03_dp
      fit_coef(0, 15) = 1.8436235820019048E-05_dp
      fit_coef(1, 15) = 2.3107371865340239E-05_dp
      fit_coef(2, 15) = 2.5779879028853276E-03_dp
      fit_coef(3, 15) = -1.5220660379968103E-02_dp
      fit_coef(4, 15) = 1.4674238877197790E-01_dp
      fit_coef(5, 15) = -6.2720720593806745E-01_dp
      fit_coef(6, 15) = 1.9196496962860754E+00_dp
      fit_coef(7, 15) = -4.0563204420089054E+00_dp
      fit_coef(8, 15) = 5.8719543662268974E+00_dp
      fit_coef(9, 15) = -5.7179931070241752E+00_dp
      fit_coef(10, 15) = 3.5792099623773965E+00_dp
      fit_coef(11, 15) = -1.3005810138247331E+00_dp
      fit_coef(12, 15) = 2.0848698369596702E-01_dp
      fit_coef(0, 16) = 1.2124348449172285E-04_dp
      fit_coef(1, 16) = 8.3813269418456530E-05_dp
      fit_coef(2, 16) = 9.0897102638036447E-03_dp
      fit_coef(3, 16) = -5.7911668892308865E-02_dp
      fit_coef(4, 16) = 4.9294344453578842E-01_dp
      fit_coef(5, 16) = -2.1314596717186927E+00_dp
      fit_coef(6, 16) = 6.5094911565878411E+00_dp
      fit_coef(7, 16) = -1.3738226450588236E+01_dp
      fit_coef(8, 16) = 1.9897297045097488E+01_dp
      fit_coef(9, 16) = -1.9402396684702406E+01_dp
      fit_coef(10, 16) = 1.2166199138851475E+01_dp
      fit_coef(11, 16) = -4.4288923635919755E+00_dp
      fit_coef(12, 16) = 7.1121697870477651E-01_dp
      fit_coef(0, 17) = 4.9649968043970563E-04_dp
      fit_coef(1, 17) = 2.4013320731623524E-04_dp
      fit_coef(2, 17) = 2.5742527096276644E-02_dp
      fit_coef(3, 17) = -1.7243873414770514E-01_dp
      fit_coef(4, 17) = 1.3895762899290705E+00_dp
      fit_coef(5, 17) = -6.3761453454497428E+00_dp
      fit_coef(6, 17) = 1.9936814180089332E+01_dp
      fit_coef(7, 17) = -4.2440835139836899E+01_dp
      fit_coef(8, 17) = 6.1588254511673313E+01_dp
      fit_coef(9, 17) = -5.9991544128378500E+01_dp
      fit_coef(10, 17) = 3.7528872717238954E+01_dp
      fit_coef(11, 17) = -1.3623706398668176E+01_dp
      fit_coef(12, 17) = 2.1816169197074027E+00_dp
      fit_coef(0, 18) = 1.6143519027172607E-03_dp
      fit_coef(1, 18) = 5.8859169571603552E-04_dp
      fit_coef(2, 18) = 6.2478735725773553E-02_dp
      fit_coef(3, 18) = -4.2711342899465948E-01_dp
      fit_coef(4, 18) = 3.3147135916470036E+00_dp
      fit_coef(5, 18) = -1.5401373294982442E+01_dp
      fit_coef(6, 18) = 4.7954953793067048E+01_dp
      fit_coef(7, 18) = -1.0153978302120807E+02_dp
      fit_coef(8, 18) = 1.4682909652330889E+02_dp
      fit_coef(9, 18) = -1.4276036864090156E+02_dp
      fit_coef(10, 18) = 8.9245840751766565E+01_dp
      fit_coef(11, 18) = -3.2398300584951542E+01_dp
      fit_coef(12, 18) = 5.1900908130178340E+00_dp
      fit_coef(0, 19) = 4.5454454347221614E-03_dp
      fit_coef(1, 19) = 1.2919916988711280E-03_dp
      fit_coef(2, 19) = 1.3597866025632654E-01_dp
      fit_coef(3, 19) = -9.4213191119816730E-01_dp
      fit_coef(4, 19) = 7.1313517958816419E+00_dp
      fit_coef(5, 19) = -3.3458248656849491E+01_dp
      fit_coef(6, 19) = 1.0462115349553649E+02_dp
      fit_coef(7, 19) = -2.2230013365840125E+02_dp
      fit_coef(8, 19) = 3.2237995602287663E+02_dp
      fit_coef(9, 19) = -3.1411722628668201E+02_dp
      fit_coef(10, 19) = 1.9665741705687799E+02_dp
      fit_coef(11, 19) = -7.1459843200825901E+01_dp
      fit_coef(12, 19) = 1.1454531272981891E+01_dp
      fit_coef(0, 20) = 1.1569593836045351E-02_dp
      fit_coef(1, 20) = 2.6018407470241437E-03_dp
      fit_coef(2, 20) = 2.7208902236045929E-01_dp
      fit_coef(3, 20) = -1.9039272645212668E+00_dp
      fit_coef(4, 20) = 1.4140811044678967E+01_dp
      fit_coef(5, 20) = -6.6772929853793315E+01_dp
      fit_coef(6, 20) = 2.0958553375037602E+02_dp
      fit_coef(7, 20) = -4.4639362838727067E+02_dp
      fit_coef(8, 20) = 6.4825457307516922E+02_dp
      fit_coef(9, 20) = -6.3211066495078410E+02_dp
      fit_coef(10, 20) = 3.9589944611561242E+02_dp
      fit_coef(11, 20) = -1.4389110171301371E+02_dp
      fit_coef(12, 20) = 2.3067990020945253E+01_dp
      fit_coef(0, 21) = 2.7295736154530334E-02_dp
      fit_coef(1, 21) = 4.8782968254850083E-03_dp
      fit_coef(2, 21) = 5.0739031138062951E-01_dp
      fit_coef(3, 21) = -3.5769970240701987E+00_dp
      fit_coef(4, 21) = 2.6164736431475280E+01_dp
      fit_coef(5, 21) = -1.2410237515226034E+02_dp
      fit_coef(6, 21) = 3.9082320961974159E+02_dp
      fit_coef(7, 21) = -8.3417851759708878E+02_dp
      fit_coef(8, 21) = 1.2129937424978116E+03_dp
      fit_coef(9, 21) = -1.1837877342560605E+03_dp
      fit_coef(10, 21) = 7.4185106985661810E+02_dp
      fit_coef(11, 21) = -2.6974056955563321E+02_dp
      fit_coef(12, 21) = 4.3256996209012456E+01_dp
      fit_coef(0, 22) = 6.0660781563425317E-02_dp
      fit_coef(1, 22) = 8.5794040217948859E-03_dp
      fit_coef(2, 22) = 8.8821388130041634E-01_dp
      fit_coef(3, 22) = -6.2985610532830956E+00_dp
      fit_coef(4, 22) = 4.5492645494685874E+01_dp
      fit_coef(5, 22) = -2.1647261584955578E+02_dp
      fit_coef(6, 22) = 6.8370044474800454E+02_dp
      fit_coef(7, 22) = -1.4621335320610413E+03_dp
      fit_coef(8, 22) = 2.1287754889346256E+03_dp
      fit_coef(9, 22) = -2.0792410286940094E+03_dp
      fit_coef(10, 22) = 1.3037501245395274E+03_dp
      fit_coef(11, 22) = -4.7423981275141949E+02_dp
      fit_coef(12, 22) = 7.6073607209164223E+01_dp
      fit_coef(0, 23) = 1.2841217227424742E-01_dp
      fit_coef(1, 23) = 1.4188057337444604E-02_dp
      fit_coef(2, 23) = 1.4629508408620973E+00_dp
      fit_coef(3, 23) = -1.0423327628190115E+01_dp
      fit_coef(4, 23) = 7.4479996989370790E+01_dp
      fit_coef(5, 23) = -3.5524319107213654E+02_dp
      fit_coef(6, 23) = 1.1248481685847394E+03_dp
      fit_coef(7, 23) = -2.4098075097537289E+03_dp
      fit_coef(8, 23) = 3.5126267297370300E+03_dp
      fit_coef(9, 23) = -3.4335802646233810E+03_dp
      fit_coef(10, 23) = 2.1541281893481446E+03_dp
      fit_coef(11, 23) = -7.8386608937552410E+02_dp
      fit_coef(12, 23) = 1.2577649211754711E+02_dp
      fit_coef(0, 24) = 2.6110218583026157E-01_dp
      fit_coef(1, 24) = 2.2013210715212957E-02_dp
      fit_coef(2, 24) = 2.2618021820392888E+00_dp
      fit_coef(3, 24) = -1.6177325764415730E+01_dp
      fit_coef(4, 24) = 1.1453533241330872E+02_dp
      fit_coef(5, 24) = -5.4723833721116853E+02_dp
      fit_coef(6, 24) = 1.7366429607985401E+03_dp
      fit_coef(7, 24) = -3.7264526069498952E+03_dp
      fit_coef(8, 24) = 5.4376949325389724E+03_dp
      fit_coef(9, 24) = -5.3192362305534161E+03_dp
      fit_coef(10, 24) = 3.3388433630606642E+03_dp
      fit_coef(11, 24) = -1.2154181459321210E+03_dp
      fit_coef(12, 24) = 1.9507456139845061E+02_dp
      fit_coef(0, 25) = 5.1367200823804149E-01_dp
      fit_coef(1, 25) = 3.1811354680997252E-02_dp
      fit_coef(2, 25) = 3.2585980335370284E+00_dp
      fit_coef(3, 25) = -2.3379379060808798E+01_dp
      fit_coef(4, 25) = 1.6424152452481010E+02_dp
      fit_coef(5, 25) = -7.8569152121234890E+02_dp
      fit_coef(6, 25) = 2.4981080171697954E+03_dp
      fit_coef(7, 25) = -5.3680206518529076E+03_dp
      fit_coef(8, 25) = 7.8407527495523063E+03_dp
      fit_coef(9, 25) = -7.6750999224540392E+03_dp
      fit_coef(10, 25) = 4.8198767192823361E+03_dp
      fit_coef(11, 25) = -1.7551476130824467E+03_dp
      fit_coef(12, 25) = 2.8177229346472944E+02_dp
      fit_coef(0, 26) = 9.8707480111071055E-01_dp
      fit_coef(1, 26) = 4.2323366844794610E-02_dp
      fit_coef(2, 26) = 4.3248620683488594E+00_dp
      fit_coef(3, 26) = -3.1101768723183273E+01_dp
      fit_coef(4, 26) = 2.1715876476583264E+02_dp
      fit_coef(5, 26) = -1.0396442548695900E+03_dp
      fit_coef(6, 26) = 3.3104927069027003E+03_dp
      fit_coef(7, 26) = -7.1219751261209803E+03_dp
      fit_coef(8, 26) = 1.0411115847803250E+04_dp
      fit_coef(9, 26) = -1.0196927796670321E+04_dp
      fit_coef(10, 26) = 6.4061259733442739E+03_dp
      fit_coef(11, 26) = -2.3334570806260454E+03_dp
      fit_coef(12, 26) = 3.7469580771614551E+02_dp
      fit_coef(0, 27) = 1.8940435475450061E+00_dp
      fit_coef(1, 27) = 5.1360926179740848E-02_dp
      fit_coef(2, 27) = 5.2404454262865769E+00_dp
      fit_coef(3, 27) = -3.7737331118803013E+01_dp
      fit_coef(4, 27) = 2.6250426447967351E+02_dp
      fit_coef(5, 27) = -1.2571969272125693E+03_dp
      fit_coef(6, 27) = 4.0068677681386607E+03_dp
      fit_coef(7, 27) = -8.6264349397777878E+03_dp
      fit_coef(8, 27) = 1.2616976771304368E+04_dp
      fit_coef(9, 27) = -1.2361959020546034E+04_dp
      fit_coef(10, 27) = 7.7683322759890580E+03_dp
      fit_coef(11, 27) = -2.8301914681580724E+03_dp
      fit_coef(12, 27) = 4.5452476900845880E+02_dp
      fit_coef(0, 28) = 3.9546758238329791E+00_dp
      fit_coef(1, 28) = 5.8523148761978797E-02_dp
      fit_coef(2, 28) = 5.9702778386920183E+00_dp
      fit_coef(3, 28) = -4.2997287866308362E+01_dp
      fit_coef(4, 28) = 2.9898597699426614E+02_dp
      fit_coef(5, 28) = -1.4318903472750299E+03_dp
      fit_coef(6, 28) = 4.5640259927402913E+03_dp
      fit_coef(7, 28) = -9.8267478597299942E+03_dp
      fit_coef(8, 28) = 1.4373446674375520E+04_dp
      fit_coef(9, 28) = -1.4083574897310356E+04_dp
      fit_coef(10, 28) = 8.8505067497423279E+03_dp
      fit_coef(11, 28) = -3.2245357392636670E+03_dp
      fit_coef(12, 28) = 5.1786594922938832E+02_dp

   END SUBROUTINE get_fit_coef_k14_lo

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k15_sh(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 2.0000000000000000E+01_dp
      U_b = 1.0000000000000000E+02_dp
      fit_coef(0, 1) = 3.7767636658328917E-03_dp
      fit_coef(1, 1) = 8.1211740094389931E-03_dp
      fit_coef(2, 1) = 8.4713080613072721E-03_dp
      fit_coef(3, 1) = 4.9722625039685525E-03_dp
      fit_coef(4, 1) = 1.2149725765861214E-03_dp
      fit_coef(5, 1) = -4.4528795735146096E-04_dp
      fit_coef(6, 1) = -4.1397080234777437E-04_dp
      fit_coef(7, 1) = -3.3426994813408613E-05_dp
      fit_coef(8, 1) = 6.8763616446013435E-05_dp
      fit_coef(9, 1) = 2.6807372741525578E-05_dp
      fit_coef(10, 1) = -2.8714180956136556E-05_dp
      fit_coef(11, 1) = 9.6272015928166124E-06_dp
      fit_coef(12, 1) = -1.2744413786197213E-06_dp
      fit_coef(0, 2) = 2.0284783343273730E-02_dp
      fit_coef(1, 2) = 4.3187845019499065E-02_dp
      fit_coef(2, 2) = 4.4768123513328789E-02_dp
      fit_coef(3, 2) = 2.6126720474645704E-02_dp
      fit_coef(4, 2) = 6.3201115983603541E-03_dp
      fit_coef(5, 2) = -2.3578759036246063E-03_dp
      fit_coef(6, 2) = -2.1606459505341194E-03_dp
      fit_coef(7, 2) = -1.6580374853951958E-04_dp
      fit_coef(8, 2) = 3.6074887455186000E-04_dp
      fit_coef(9, 2) = 1.3781750573407730E-04_dp
      fit_coef(10, 2) = -1.4977975033955343E-04_dp
      fit_coef(11, 2) = 5.0490744097337301E-05_dp
      fit_coef(12, 2) = -6.7504814416040953E-06_dp
      fit_coef(0, 3) = 5.1673260780627950E-02_dp
      fit_coef(1, 3) = 1.0789295508104972E-01_dp
      fit_coef(2, 3) = 1.1057351883667553E-01_dp
      fit_coef(3, 3) = 6.3878269616968825E-02_dp
      fit_coef(4, 3) = 1.5180871067702629E-02_dp
      fit_coef(5, 3) = -5.8375073829049018E-03_dp
      fit_coef(6, 3) = -5.2188107945168550E-03_dp
      fit_coef(7, 3) = -3.6420812718932399E-04_dp
      fit_coef(8, 3) = 8.7905419776698594E-04_dp
      fit_coef(9, 3) = 3.2415893028994188E-04_dp
      fit_coef(10, 3) = -3.6139234804234916E-04_dp
      fit_coef(11, 3) = 1.2287995757801328E-04_dp
      fit_coef(12, 3) = -1.6675660997752087E-05_dp
      fit_coef(0, 4) = 1.0134621267461434E-01_dp
      fit_coef(1, 4) = 2.0498960034213551E-01_dp
      fit_coef(2, 4) = 2.0657812668218858E-01_dp
      fit_coef(3, 4) = 1.1767379558549464E-01_dp
      fit_coef(4, 4) = 2.7288178413390367E-02_dp
      fit_coef(5, 4) = -1.0924148745488705E-02_dp
      fit_coef(6, 4) = -9.4432740413890777E-03_dp
      fit_coef(7, 4) = -5.6769752225927254E-04_dp
      fit_coef(8, 4) = 1.6082317938650511E-03_dp
      fit_coef(9, 4) = 5.6496352365318649E-04_dp
      fit_coef(10, 4) = -6.5204939503955980E-04_dp
      fit_coef(11, 4) = 2.2378537926777479E-04_dp
      fit_coef(12, 4) = -3.0821862755447693E-05_dp
      fit_coef(0, 5) = 1.7540646948915822E-01_dp
      fit_coef(1, 5) = 3.3841107332002812E-01_dp
      fit_coef(2, 5) = 3.3317753537976630E-01_dp
      fit_coef(3, 5) = 1.8664213256865744E-01_dp
      fit_coef(4, 5) = 4.2122973931643735E-02_dp
      fit_coef(5, 5) = -1.7653496250856538E-02_dp
      fit_coef(6, 5) = -1.4677820339754412E-02_dp
      fit_coef(7, 5) = -6.7172992100736459E-04_dp
      fit_coef(8, 5) = 2.5463509524145020E-03_dp
      fit_coef(9, 5) = 7.9689495762103719E-04_dp
      fit_coef(10, 5) = -9.9176850710568391E-04_dp
      fit_coef(11, 5) = 3.5040533908774847E-04_dp
      fit_coef(12, 5) = -5.0509200990651491E-05_dp
      fit_coef(0, 6) = 2.8403858696731338E-01_dp
      fit_coef(1, 6) = 5.1354379131383487E-01_dp
      fit_coef(2, 6) = 4.8963667336603639E-01_dp
      fit_coef(3, 6) = 2.6912774693427816E-01_dp
      fit_coef(4, 6) = 5.9724553513996380E-02_dp
      fit_coef(5, 6) = -2.5703789609353513E-02_dp
      fit_coef(6, 6) = -2.1114838016090463E-02_dp
      fit_coef(7, 6) = -7.0228514229756346E-04_dp
      fit_coef(8, 6) = 3.8581969229036879E-03_dp
      fit_coef(9, 6) = 1.2194767346868753E-03_dp
      fit_coef(10, 6) = -1.8021363838819021E-03_dp
      fit_coef(11, 6) = 7.5507648660773613E-04_dp
      fit_coef(12, 6) = -1.3075406318307264E-04_dp
      fit_coef(0, 7) = 4.4325614880929831E-01_dp
      fit_coef(1, 7) = 7.3781666079304498E-01_dp
      fit_coef(2, 7) = 6.7367679967887084E-01_dp
      fit_coef(3, 7) = 3.6132751402347024E-01_dp
      fit_coef(4, 7) = 8.0499595178773156E-02_dp
      fit_coef(5, 7) = -3.3029532312549834E-02_dp
      fit_coef(6, 7) = -2.8790475820923794E-02_dp
      fit_coef(7, 7) = -1.7795611862018949E-03_dp
      fit_coef(8, 7) = 5.1379466412685825E-03_dp
      fit_coef(9, 7) = 3.1353774828616044E-03_dp
      fit_coef(10, 7) = -3.3341189829340133E-03_dp
      fit_coef(11, 7) = 9.5780853222126921E-04_dp
      fit_coef(12, 7) = -4.7564006414347488E-05_dp
      fit_coef(0, 8) = 6.7707705919780792E-01_dp
      fit_coef(1, 8) = 1.0208946341910239E+00_dp
      fit_coef(2, 8) = 8.8225292727415439E-01_dp
      fit_coef(3, 8) = 4.5683990983568840E-01_dp
      fit_coef(4, 8) = 1.0315634645104643E-01_dp
      fit_coef(5, 8) = -3.6302702410040397E-02_dp
      fit_coef(6, 8) = -3.5178064788216744E-02_dp
      fit_coef(7, 8) = -5.4614672692391321E-03_dp
      fit_coef(8, 8) = 5.9947422484209678E-03_dp
      fit_coef(9, 8) = 2.7277018733018010E-03_dp
      fit_coef(10, 8) = 9.2005661112863006E-04_dp
      fit_coef(11, 8) = -2.7491424308858183E-03_dp
      fit_coef(12, 8) = 9.3417034494363172E-04_dp
      fit_coef(0, 9) = 1.0205594171746224E+00_dp
      fit_coef(1, 9) = 1.3737393203947723E+00_dp
      fit_coef(2, 9) = 1.1123746840237188E+00_dp
      fit_coef(3, 9) = 5.4861957447712506E-01_dp
      fit_coef(4, 9) = 1.2420055508662378E-01_dp
      fit_coef(5, 9) = -3.4071192988448372E-02_dp
      fit_coef(6, 9) = -3.6752879111842843E-02_dp
      fit_coef(7, 9) = -9.8874781786285233E-03_dp
      fit_coef(8, 9) = 5.0195565707678762E-03_dp
      fit_coef(9, 9) = -7.6508018734261539E-04_dp
      fit_coef(10, 9) = 7.1843501967062814E-03_dp
      fit_coef(11, 9) = -4.9459397204347057E-03_dp
      fit_coef(12, 9) = 8.9782856131910235E-04_dp
      fit_coef(0, 10) = 1.5244258156390766E+00_dp
      fit_coef(1, 10) = 1.8066993429217002E+00_dp
      fit_coef(2, 10) = 1.3609242838483797E+00_dp
      fit_coef(3, 10) = 6.3115404377443018E-01_dp
      fit_coef(4, 10) = 1.4030025341356947E-01_dp
      fit_coef(5, 10) = -2.7329038046514077E-02_dp
      fit_coef(6, 10) = -3.2676873494217773E-02_dp
      fit_coef(7, 10) = -1.1415330583361451E-02_dp
      fit_coef(8, 10) = -4.9009191168826780E-06_dp
      fit_coef(9, 10) = 1.3306667731994193E-04_dp
      fit_coef(10, 10) = 3.7437558621727469E-03_dp
      fit_coef(11, 10) = 4.8589108656243562E-04_dp
      fit_coef(12, 10) = -1.0469727509729387E-03_dp
      fit_coef(0, 11) = 2.2624235303155098E+00_dp
      fit_coef(1, 11) = 2.3271837597685869E+00_dp
      fit_coef(2, 11) = 1.6242573707008741E+00_dp
      fit_coef(3, 11) = 7.0131249760739722E-01_dp
      fit_coef(4, 11) = 1.5012768567411283E-01_dp
      fit_coef(5, 11) = -1.7685800296799118E-02_dp
      fit_coef(6, 11) = -2.4133853082992485E-02_dp
      fit_coef(7, 11) = -1.0124353385573369E-02_dp
      fit_coef(8, 11) = -5.4053355656654414E-03_dp
      fit_coef(9, 11) = 2.0601323692594885E-03_dp
      fit_coef(10, 11) = -1.9168010184212245E-03_dp
      fit_coef(11, 11) = 5.2635221561462978E-03_dp
      fit_coef(12, 11) = -2.0650943061245367E-03_dp
      fit_coef(0, 12) = 3.3440185138398908E+00_dp
      fit_coef(1, 12) = 2.9374869303015658E+00_dp
      fit_coef(2, 12) = 1.8986970103463769E+00_dp
      fit_coef(3, 12) = 7.5852460527436116E-01_dp
      fit_coef(4, 12) = 1.5431087315065109E-01_dp
      fit_coef(5, 12) = -6.2637714258166168E-03_dp
      fit_coef(6, 12) = -1.2867661430609648E-02_dp
      fit_coef(7, 12) = -7.7297513116621266E-03_dp
      fit_coef(8, 12) = -7.9684475195776500E-03_dp
      fit_coef(9, 12) = 7.0191640613392089E-04_dp
      fit_coef(10, 12) = -2.3205524009468916E-03_dp
      fit_coef(11, 12) = 4.7115703272364142E-03_dp
      fit_coef(12, 12) = -1.3939780575372204E-03_dp
      fit_coef(0, 13) = 4.9400406076482417E+00_dp
      fit_coef(1, 13) = 3.6338865240810336E+00_dp
      fit_coef(2, 13) = 2.1823364183861980E+00_dp
      fit_coef(3, 13) = 8.0486953703551811E-01_dp
      fit_coef(4, 13) = 1.5444706443796333E-01_dp
      fit_coef(5, 13) = 6.2245841342081663E-03_dp
      fit_coef(6, 13) = -4.5752093352805732E-04_dp
      fit_coef(7, 13) = -4.9804550111796943E-03_dp
      fit_coef(8, 13) = -8.4778512610873564E-03_dp
      fit_coef(9, 13) = -2.0214510332599842E-03_dp
      fit_coef(10, 13) = 2.4374313233334989E-04_dp
      fit_coef(11, 13) = 1.9195383981977568E-03_dp
      fit_coef(12, 13) = -3.7228313108121918E-04_dp
      fit_coef(0, 14) = 7.3492824088289943E+00_dp
      fit_coef(1, 14) = 4.4099530097240827E+00_dp
      fit_coef(2, 14) = 2.4783620449758916E+00_dp
      fit_coef(3, 14) = 8.4513384310546469E-01_dp
      fit_coef(4, 14) = 1.5211923704319483E-01_dp
      fit_coef(5, 14) = 1.9150817795865041E-02_dp
      fit_coef(6, 14) = 1.2175087669328011E-02_dp
      fit_coef(7, 14) = -1.9876602775734512E-03_dp
      fit_coef(8, 14) = -8.3964504059626105E-03_dp
      fit_coef(9, 14) = -3.7517484520995536E-03_dp
      fit_coef(10, 14) = 2.5098290732908429E-03_dp
      fit_coef(11, 14) = -2.2174814827198228E-04_dp
      fit_coef(12, 14) = 1.7762122632969973E-04_dp
      fit_coef(0, 15) = 1.1260372250288922E+01_dp
      fit_coef(1, 15) = 5.2761367256188780E+00_dp
      fit_coef(2, 15) = 2.8037057941841899E+00_dp
      fit_coef(3, 15) = 8.8780102116211801E-01_dp
      fit_coef(4, 15) = 1.4843982592096647E-01_dp
      fit_coef(5, 15) = 3.2027815482417170E-02_dp
      fit_coef(6, 15) = 2.5031908188376327E-02_dp
      fit_coef(7, 15) = 1.3905318224525349E-03_dp
      fit_coef(8, 15) = -8.2768851858847706E-03_dp
      fit_coef(9, 15) = -4.1010656916587729E-03_dp
      fit_coef(10, 15) = 3.7595505954445777E-03_dp
      fit_coef(11, 15) = -1.2076253125456195E-03_dp
      fit_coef(12, 15) = 2.7539917666371817E-04_dp
      fit_coef(0, 16) = 9.7353035224829500E-03_dp
      fit_coef(1, 16) = 2.0886976749800733E-02_dp
      fit_coef(2, 16) = 2.1755665531631045E-02_dp
      fit_coef(3, 16) = 1.2752382400370278E-02_dp
      fit_coef(4, 16) = 3.1086335160659716E-03_dp
      fit_coef(5, 16) = -1.1441703553915466E-03_dp
      fit_coef(6, 16) = -1.0600438281376386E-03_dp
      fit_coef(7, 16) = -8.4585212227881147E-05_dp
      fit_coef(8, 16) = 1.7629742631278603E-04_dp
      fit_coef(9, 16) = 6.8399637396386302E-05_dp
      fit_coef(10, 16) = -7.3517685522867944E-05_dp
      fit_coef(11, 16) = 2.4681600354643785E-05_dp
      fit_coef(12, 16) = -3.2754157510708014E-06_dp
      fit_coef(0, 17) = 2.3550935618358620E-02_dp
      fit_coef(1, 17) = 4.9522562600574506E-02_dp
      fit_coef(2, 17) = 5.0940261314895138E-02_dp
      fit_coef(3, 17) = 2.9520151355165270E-02_dp
      fit_coef(4, 17) = 7.0526431743143485E-03_dp
      fit_coef(5, 17) = -2.6886627324933282E-03_dp
      fit_coef(6, 17) = -2.4209379054065576E-03_dp
      fit_coef(7, 17) = -1.7381016208850244E-04_dp
      fit_coef(8, 17) = 4.0675731507373121E-04_dp
      fit_coef(9, 17) = 1.5152633434160025E-04_dp
      fit_coef(10, 17) = -1.6769671036439022E-04_dp
      fit_coef(11, 17) = 5.6896776687533419E-05_dp
      fit_coef(12, 17) = -7.6947045683246862E-06_dp
      fit_coef(0, 18) = 3.9777578377675531E-02_dp
      fit_coef(1, 18) = 8.0346278051875858E-02_dp
      fit_coef(2, 18) = 8.0791839568255266E-02_dp
      fit_coef(3, 18) = 4.5890106592349850E-02_dp
      fit_coef(4, 18) = 1.0581783818106148E-02_dp
      fit_coef(5, 18) = -4.2772167663816868E-03_dp
      fit_coef(6, 18) = -3.6711065186785822E-03_dp
      fit_coef(7, 18) = -2.1295067951996251E-04_dp
      fit_coef(8, 18) = 6.2748074852499950E-04_dp
      fit_coef(9, 18) = 2.1773104607889315E-04_dp
      fit_coef(10, 18) = -2.5378512352868573E-04_dp
      fit_coef(11, 18) = 8.7481966065534113E-05_dp
      fit_coef(12, 18) = -1.2141070019913610E-05_dp
      fit_coef(0, 19) = 6.0562131081817058E-02_dp
      fit_coef(1, 19) = 1.1449873938476936E-01_dp
      fit_coef(2, 19) = 1.1129944008292268E-01_dp
      fit_coef(3, 19) = 6.1560590847477989E-02_dp
      fit_coef(4, 19) = 1.3534891528073901E-02_dp
      fit_coef(5, 19) = -5.9033543536541871E-03_dp
      fit_coef(6, 19) = -4.7450366895801120E-03_dp
      fit_coef(7, 19) = -1.7764449860279051E-04_dp
      fit_coef(8, 19) = 8.2604161016353668E-04_dp
      fit_coef(9, 19) = 2.5651846839713858E-04_dp
      fit_coef(10, 19) = -3.2185438255892302E-04_dp
      fit_coef(11, 19) = 1.1225233677891509E-04_dp
      fit_coef(12, 19) = -1.5831921058481165E-05_dp
      fit_coef(0, 20) = 8.9226696029307800E-02_dp
      fit_coef(1, 20) = 1.5322995807536874E-01_dp
      fit_coef(2, 20) = 1.4178721914749884E-01_dp
      fit_coef(3, 20) = 7.6128391563680162E-02_dp
      fit_coef(4, 20) = 1.6136643158648808E-02_dp
      fit_coef(5, 20) = -7.5330493292977262E-03_dp
      fit_coef(6, 20) = -5.7529238965828419E-03_dp
      fit_coef(7, 20) = -1.1707762217534482E-05_dp
      fit_coef(8, 20) = 1.0894518853190246E-03_dp
      fit_coef(9, 20) = 1.8904025635640227E-04_dp
      fit_coef(10, 20) = -3.8503303956792745E-04_dp
      fit_coef(11, 20) = 1.6892465284271444E-04_dp
      fit_coef(12, 20) = -3.1503784546561519E-05_dp
      fit_coef(0, 21) = 1.3067133728179103E-01_dp
      fit_coef(1, 21) = 1.9826805579875501E-01_dp
      fit_coef(2, 21) = 1.7073403778885513E-01_dp
      fit_coef(3, 21) = 8.8223081054500926E-02_dp
      fit_coef(4, 21) = 1.9182613501886863E-02_dp
      fit_coef(5, 21) = -8.2278999968844989E-03_dp
      fit_coef(6, 21) = -7.1938427008996124E-03_dp
      fit_coef(7, 21) = -2.3630795339020511E-04_dp
      fit_coef(8, 21) = 1.4678849246082021E-03_dp
      fit_coef(9, 21) = 1.0431379287022221E-03_dp
      fit_coef(10, 21) = -1.5357712204035895E-03_dp
      fit_coef(11, 21) = 7.1343365682177630E-04_dp
      fit_coef(12, 21) = -1.2573122208893727E-04_dp
      fit_coef(0, 22) = 1.9171789962238633E-01_dp
      fit_coef(1, 22) = 2.5190501212777222E-01_dp
      fit_coef(2, 22) = 1.9682158960117535E-01_dp
      fit_coef(3, 22) = 9.5068036046292548E-02_dp
      fit_coef(4, 22) = 2.2164921662634080E-02_dp
      fit_coef(5, 22) = -5.7919264320229974E-03_dp
      fit_coef(6, 22) = -7.7045178800915467E-03_dp
      fit_coef(7, 22) = -2.2724175760734976E-03_dp
      fit_coef(8, 22) = 1.0118193318505186E-03_dp
      fit_coef(9, 22) = 2.0169812523148983E-03_dp
      fit_coef(10, 22) = -9.2577594680840871E-06_dp
      fit_coef(11, 22) = -1.3014429778002868E-03_dp
      fit_coef(12, 22) = 5.2474898728580737E-04_dp
      fit_coef(0, 23) = 2.8169225019506716E-01_dp
      fit_coef(1, 23) = 3.1611118986493281E-01_dp
      fit_coef(2, 23) = 2.1983832131612513E-01_dp
      fit_coef(3, 23) = 9.4742330833316193E-02_dp
      fit_coef(4, 23) = 2.2513794210120210E-02_dp
      fit_coef(5, 23) = -5.0847893755973886E-04_dp
      fit_coef(6, 23) = -4.3926536784528831E-03_dp
      fit_coef(7, 23) = -4.7592264209054850E-03_dp
      fit_coef(8, 23) = 5.4373390930081153E-04_dp
      fit_coef(9, 23) = -3.2404033623214612E-03_dp
      fit_coef(10, 23) = 7.9533405554355376E-03_dp
      fit_coef(11, 23) = -5.0278147899489754E-03_dp
      fit_coef(12, 23) = 1.0211564000482632E-03_dp
      fit_coef(0, 24) = 4.1360870667704458E-01_dp
      fit_coef(1, 24) = 3.9134279925186166E-01_dp
      fit_coef(2, 24) = 2.3989520225411959E-01_dp
      fit_coef(3, 24) = 8.7887165884511398E-02_dp
      fit_coef(4, 24) = 1.8986346126332081E-02_dp
      fit_coef(5, 24) = 4.7658422103979536E-03_dp
      fit_coef(6, 24) = 1.3492788485409764E-03_dp
      fit_coef(7, 24) = -3.3357808918414787E-03_dp
      fit_coef(8, 24) = -3.1088669483680088E-03_dp
      fit_coef(9, 24) = -1.8511419844152608E-03_dp
      fit_coef(10, 24) = 1.7736202675366283E-03_dp
      fit_coef(11, 24) = 2.1114330344624466E-03_dp
      fit_coef(12, 24) = -1.3226303250553409E-03_dp
      fit_coef(0, 25) = 6.0621802446427508E-01_dp
      fit_coef(1, 25) = 4.7585064174951042E-01_dp
      fit_coef(2, 25) = 2.5658915757225498E-01_dp
      fit_coef(3, 25) = 7.6658354342588539E-02_dp
      fit_coef(4, 25) = 1.2995988704122139E-02_dp
      fit_coef(5, 25) = 8.4318408777111463E-03_dp
      fit_coef(6, 25) = 6.5820665925556696E-03_dp
      fit_coef(7, 25) = 1.8845523065198023E-04_dp
      fit_coef(8, 25) = -6.1170321098400040E-03_dp
      fit_coef(9, 25) = 2.6963328854067922E-03_dp
      fit_coef(10, 25) = -6.7464598603185540E-03_dp
      fit_coef(11, 25) = 6.8191427136023554E-03_dp
      fit_coef(12, 25) = -1.9275464800632919E-03_dp
      fit_coef(0, 26) = 8.8764155033206149E-01_dp
      fit_coef(1, 26) = 5.6552140990150657E-01_dp
      fit_coef(2, 26) = 2.6943393742267280E-01_dp
      fit_coef(3, 26) = 6.3561620296985116E-02_dp
      fit_coef(4, 26) = 6.7821546029915858E-03_dp
      fit_coef(5, 26) = 1.0671885911308030E-02_dp
      fit_coef(6, 26) = 1.0202617186318210E-02_dp
      fit_coef(7, 26) = 2.0623230180061023E-03_dp
      fit_coef(8, 26) = -4.0999080350069333E-03_dp
      fit_coef(9, 26) = 3.3757569259654344E-04_dp
      fit_coef(10, 26) = -3.2168362434274113E-03_dp
      fit_coef(11, 26) = 1.9561977675803606E-03_dp
      fit_coef(12, 26) = -1.0127854359087637E-05_dp
      fit_coef(0, 27) = 1.3031121163620227E+00_dp
      fit_coef(1, 27) = 6.5442627957957644E-01_dp
      fit_coef(2, 27) = 2.7907544789512745E-01_dp
      fit_coef(3, 27) = 5.1209194718773846E-02_dp
      fit_coef(4, 27) = 1.8680237447576146E-03_dp
      fit_coef(5, 27) = 1.2060225663826541E-02_dp
      fit_coef(6, 27) = 1.2057587702851754E-02_dp
      fit_coef(7, 27) = 2.6179284359982300E-03_dp
      fit_coef(8, 27) = -1.2404837715336242E-03_dp
      fit_coef(9, 27) = -2.5831650098316326E-03_dp
      fit_coef(10, 27) = 1.8086181493242913E-03_dp
      fit_coef(11, 27) = -2.3645920962560019E-03_dp
      fit_coef(12, 27) = 1.0746225263111922E-03_dp
      fit_coef(0, 28) = 1.9357875429827360E+00_dp
      fit_coef(1, 28) = 7.3707936760259829E-01_dp
      fit_coef(2, 28) = 2.8863411383043175E-01_dp
      fit_coef(3, 28) = 4.2272421989736957E-02_dp
      fit_coef(4, 28) = -1.3236760602961552E-03_dp
      fit_coef(5, 28) = 1.2811295424819810E-02_dp
      fit_coef(6, 28) = 1.2606618885406514E-02_dp
      fit_coef(7, 28) = 2.8686367951981208E-03_dp
      fit_coef(8, 28) = -4.2949161661710804E-05_dp
      fit_coef(9, 28) = -2.4420917937843715E-03_dp
      fit_coef(10, 28) = 2.7346339142362607E-03_dp
      fit_coef(11, 28) = -2.7192205830173727E-03_dp
      fit_coef(12, 28) = 8.3334329375208966E-04_dp
      fit_coef(0, 29) = 2.9813305217729247E+00_dp
      fit_coef(1, 29) = 8.1556707900662084E-01_dp
      fit_coef(2, 29) = 3.0574250813204679E-01_dp
      fit_coef(3, 29) = 3.9576605747462988E-02_dp
      fit_coef(4, 29) = -3.1516191139201481E-03_dp
      fit_coef(5, 29) = 1.2947548079823706E-02_dp
      fit_coef(6, 29) = 1.2648514261782698E-02_dp
      fit_coef(7, 29) = 3.1343334418898880E-03_dp
      fit_coef(8, 29) = 1.2305058520847782E-04_dp
      fit_coef(9, 29) = -9.8020686398801802E-04_dp
      fit_coef(10, 29) = 1.7095022610177765E-03_dp
      fit_coef(11, 29) = -1.5105150611297968E-03_dp
      fit_coef(12, 29) = 2.8195377606949853E-04_dp
      fit_coef(0, 30) = 5.1939564822259765E+00_dp
      fit_coef(1, 30) = 9.3452037903374208E-01_dp
      fit_coef(2, 30) = 3.5551447159259092E-01_dp
      fit_coef(3, 30) = 4.8826933534516080E-02_dp
      fit_coef(4, 30) = -4.0672239419310479E-03_dp
      fit_coef(5, 30) = 1.2876916535847991E-02_dp
      fit_coef(6, 30) = 1.3347608970308627E-02_dp
      fit_coef(7, 30) = 3.7439290651003427E-03_dp
      fit_coef(8, 30) = 1.3798569341874731E-04_dp
      fit_coef(9, 30) = 1.7467286573212931E-04_dp
      fit_coef(10, 30) = 9.2180740758796703E-04_dp
      fit_coef(11, 30) = -5.6981792340525463E-04_dp
      fit_coef(12, 30) = -3.5378278861525159E-05_dp

   END SUBROUTINE get_fit_coef_k15_sh

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k15_mi(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 1.0000000000000000E+02_dp
      U_b = 1.1738988511192356E+03_dp
      fit_coef(0, 1) = 5.0212735786207190E-04_dp
      fit_coef(1, 1) = 1.4715598436560518E-03_dp
      fit_coef(2, 1) = 2.1503966610512443E-03_dp
      fit_coef(3, 1) = 1.8763892029773810E-03_dp
      fit_coef(4, 1) = 8.9509317267188646E-04_dp
      fit_coef(5, 1) = 4.4252413515572752E-05_dp
      fit_coef(6, 1) = -1.9929035057866966E-04_dp
      fit_coef(7, 1) = -5.0485804346307633E-05_dp
      fit_coef(8, 1) = 4.3968375574232470E-05_dp
      fit_coef(9, 1) = 5.5154425679931483E-05_dp
      fit_coef(10, 1) = -7.3825809523082769E-05_dp
      fit_coef(11, 1) = 3.4418190236750298E-05_dp
      fit_coef(12, 1) = -6.1654048128061428E-06_dp
      fit_coef(0, 2) = 2.7813466099976822E-03_dp
      fit_coef(1, 2) = 7.9639349588492489E-03_dp
      fit_coef(2, 2) = 1.1498418447714263E-02_dp
      fit_coef(3, 2) = 9.9309054844448665E-03_dp
      fit_coef(4, 2) = 4.6793464572880547E-03_dp
      fit_coef(5, 2) = 2.0134378481332076E-04_dp
      fit_coef(6, 2) = -1.0461830621068451E-03_dp
      fit_coef(7, 2) = -2.5547396359220325E-04_dp
      fit_coef(8, 2) = 2.3777699827707836E-04_dp
      fit_coef(9, 2) = 2.7502067371517350E-04_dp
      fit_coef(10, 2) = -3.7882475853089022E-04_dp
      fit_coef(11, 2) = 1.7838098098998129E-04_dp
      fit_coef(12, 2) = -3.2173474626321743E-05_dp
      fit_coef(0, 3) = 7.5356589508119715E-03_dp
      fit_coef(1, 3) = 2.0547011268364480E-02_dp
      fit_coef(2, 3) = 2.8986632022747129E-02_dp
      fit_coef(3, 3) = 2.4585742824888547E-02_dp
      fit_coef(4, 3) = 1.1332324461894762E-02_dp
      fit_coef(5, 3) = 3.5501376638767630E-04_dp
      fit_coef(6, 3) = -2.5488497267062661E-03_dp
      fit_coef(7, 3) = -5.7872409900516482E-04_dp
      fit_coef(8, 3) = 6.0733088515269604E-04_dp
      fit_coef(9, 3) = 6.0352400231591973E-04_dp
      fit_coef(10, 3) = -8.7889897017369044E-04_dp
      fit_coef(11, 3) = 4.2107155921709870E-04_dp
      fit_coef(12, 3) = -7.6791853388456923E-05_dp
      fit_coef(0, 4) = 1.6293654129500710E-02_dp
      fit_coef(1, 4) = 4.1034183058244827E-02_dp
      fit_coef(2, 4) = 5.5669340450143955E-02_dp
      fit_coef(3, 4) = 4.5958232995003026E-02_dp
      fit_coef(4, 4) = 2.0631055086044563E-02_dp
      fit_coef(5, 4) = 3.6554527705558160E-04_dp
      fit_coef(6, 4) = -4.7143452832445311E-03_dp
      fit_coef(7, 4) = -9.9264691382352382E-04_dp
      fit_coef(8, 4) = 1.2322009627532241E-03_dp
      fit_coef(9, 4) = 9.8185260435618648E-04_dp
      fit_coef(10, 4) = -1.5959947544011608E-03_dp
      fit_coef(11, 4) = 7.9825857947420098E-04_dp
      fit_coef(12, 4) = -1.5065823877361689E-04_dp
      fit_coef(0, 5) = 3.2116880043710744E-02_dp
      fit_coef(1, 5) = 7.2781942537050448E-02_dp
      fit_coef(2, 5) = 9.3125631841889350E-02_dp
      fit_coef(3, 5) = 7.3591074240215937E-02_dp
      fit_coef(4, 5) = 3.1935152140736278E-02_dp
      fit_coef(5, 5) = 4.5831309250835472E-04_dp
      fit_coef(6, 5) = -7.2579368033050861E-03_dp
      fit_coef(7, 5) = -1.5242456005729315E-03_dp
      fit_coef(8, 5) = 1.3880778676825922E-03_dp
      fit_coef(9, 5) = 2.2376670079413724E-03_dp
      fit_coef(10, 5) = -2.8229103749250520E-03_dp
      fit_coef(11, 5) = 1.2683701333535493E-03_dp
      fit_coef(12, 5) = -2.1655233738863679E-04_dp
      fit_coef(0, 6) = 6.0477196624071472E-02_dp
      fit_coef(1, 6) = 1.2125095461882505E-01_dp
      fit_coef(2, 6) = 1.4404594415379110E-01_dp
      fit_coef(3, 6) = 1.0653676582778267E-01_dp
      fit_coef(4, 6) = 4.3372431118558187E-02_dp
      fit_coef(5, 6) = 4.1288606887333705E-04_dp
      fit_coef(6, 6) = -8.8579280488838095E-03_dp
      fit_coef(7, 6) = -2.2179852026296800E-03_dp
      fit_coef(8, 6) = 2.0958278133874500E-03_dp
      fit_coef(9, 6) = 6.9556540802001902E-04_dp
      fit_coef(10, 6) = -7.7943766158600528E-04_dp
      fit_coef(11, 6) = 1.1446740869183900E-04_dp
      fit_coef(12, 6) = 2.7632689805845764E-05_dp
      fit_coef(0, 7) = 1.1064399319041990E-01_dp
      fit_coef(1, 7) = 1.9418438863145457E-01_dp
      fit_coef(2, 7) = 2.1228094292888233E-01_dp
      fit_coef(3, 7) = 1.4399775783971039E-01_dp
      fit_coef(4, 7) = 5.2541209241390302E-02_dp
      fit_coef(5, 7) = -9.2537655514377752E-04_dp
      fit_coef(6, 7) = -9.0102109734644890E-03_dp
      fit_coef(7, 7) = -1.8871157279751268E-03_dp
      fit_coef(8, 7) = 3.1309099680042760E-03_dp
      fit_coef(9, 7) = -2.5004066194694607E-03_dp
      fit_coef(10, 7) = 1.5656511931140620E-03_dp
      fit_coef(11, 7) = -4.6135496029970502E-04_dp
      fit_coef(12, 7) = 3.2764384769675960E-05_dp
      fit_coef(0, 8) = 1.9798263302865093E-01_dp
      fit_coef(1, 8) = 3.0162627763842398E-01_dp
      fit_coef(2, 8) = 3.0217808443088900E-01_dp
      fit_coef(3, 8) = 1.8500065299157431E-01_dp
      fit_coef(4, 8) = 5.7340548130963523E-02_dp
      fit_coef(5, 8) = -4.3765501944754721E-03_dp
      fit_coef(6, 8) = -7.8869300696153017E-03_dp
      fit_coef(7, 8) = 1.3835188624175719E-04_dp
      fit_coef(8, 8) = 3.3750918634805581E-03_dp
      fit_coef(9, 8) = -5.0244120026273721E-03_dp
      fit_coef(10, 8) = 2.8097352620888730E-03_dp
      fit_coef(11, 8) = -4.8749664373366714E-04_dp
      fit_coef(12, 8) = -3.4446077687337908E-05_dp
      fit_coef(0, 9) = 3.4762598832467834E-01_dp
      fit_coef(1, 9) = 4.5573634191690915E-01_dp
      fit_coef(2, 9) = 4.1780164140817255E-01_dp
      fit_coef(3, 9) = 2.2768957404533277E-01_dp
      fit_coef(4, 9) = 5.5879258865385507E-02_dp
      fit_coef(5, 9) = -9.9081094692702303E-03_dp
      fit_coef(6, 9) = -5.0068528611206149E-03_dp
      fit_coef(7, 9) = 3.0160271964011468E-03_dp
      fit_coef(8, 9) = 3.6222265118699710E-03_dp
      fit_coef(9, 9) = -8.8397638670103182E-03_dp
      fit_coef(10, 9) = 5.9259674015836239E-03_dp
      fit_coef(11, 9) = -1.7478385304278991E-03_dp
      fit_coef(12, 9) = 1.9988838059233530E-04_dp
      fit_coef(0, 10) = 6.0023140973611755E-01_dp
      fit_coef(1, 10) = 6.6999181383267858E-01_dp
      fit_coef(2, 10) = 5.6210771610729648E-01_dp
      fit_coef(3, 10) = 2.6900784986520249E-01_dp
      fit_coef(4, 10) = 4.6435418682624782E-02_dp
      fit_coef(5, 10) = -1.6642885693320342E-02_dp
      fit_coef(6, 10) = 3.8014327565738168E-04_dp
      fit_coef(7, 10) = 6.2610253146593775E-03_dp
      fit_coef(8, 10) = 3.2146802606417453E-03_dp
      fit_coef(9, 10) = -1.3711754453203115E-02_dp
      fit_coef(10, 10) = 1.1099834233825179E-02_dp
      fit_coef(11, 10) = -4.2111321446896338E-03_dp
      fit_coef(12, 10) = 6.7715733841587470E-04_dp
      fit_coef(0, 11) = 1.0210655880657196E+00_dp
      fit_coef(1, 11) = 9.5730192228839006E-01_dp
      fit_coef(2, 11) = 7.3607353846396317E-01_dp
      fit_coef(3, 11) = 3.0496319982420517E-01_dp
      fit_coef(4, 11) = 2.7825470684213698E-02_dp
      fit_coef(5, 11) = -2.3036487360505412E-02_dp
      fit_coef(6, 11) = 9.0002346429021024E-03_dp
      fit_coef(7, 11) = 9.4071699260299454E-03_dp
      fit_coef(8, 11) = 1.1246072531079006E-03_dp
      fit_coef(9, 11) = -1.8600382884449367E-02_dp
      fit_coef(10, 11) = 1.7805450549655721E-02_dp
      fit_coef(11, 11) = -7.6910226916182545E-03_dp
      fit_coef(12, 11) = 1.3712121437792250E-03_dp
      fit_coef(0, 12) = 1.7149561293339615E+00_dp
      fit_coef(1, 12) = 1.3267744256882403E+00_dp
      fit_coef(2, 12) = 9.3812895165305032E-01_dp
      fit_coef(3, 12) = 3.3164771999234482E-01_dp
      fit_coef(4, 12) = -4.9178162808578741E-05_dp
      fit_coef(5, 12) = -2.7249300178168019E-02_dp
      fit_coef(6, 12) = 2.1205458082480019E-02_dp
      fit_coef(7, 12) = 1.1807990512194026E-02_dp
      fit_coef(8, 12) = -3.3661868556830038E-03_dp
      fit_coef(9, 12) = -2.2430236127638940E-02_dp
      fit_coef(10, 12) = 2.5486188136087912E-02_dp
      fit_coef(11, 12) = -1.2102581018395429E-02_dp
      fit_coef(12, 12) = 2.3004381707239473E-03_dp
      fit_coef(0, 13) = 2.8539608561967360E+00_dp
      fit_coef(1, 13) = 1.7797688368810891E+00_dp
      fit_coef(2, 13) = 1.1646443061164560E+00_dp
      fit_coef(3, 13) = 3.4701820688332419E-01_dp
      fit_coef(4, 13) = -3.5929878683565182E-02_dp
      fit_coef(5, 13) = -2.7887790242899075E-02_dp
      fit_coef(6, 13) = 3.6702245498739489E-02_dp
      fit_coef(7, 13) = 1.2908675248958822E-02_dp
      fit_coef(8, 13) = -1.0507576982805505E-02_dp
      fit_coef(9, 13) = -2.4256099909163161E-02_dp
      fit_coef(10, 13) = 3.3404768253646827E-02_dp
      fit_coef(11, 13) = -1.7212687576065637E-02_dp
      fit_coef(12, 13) = 3.4428516812536531E-03_dp
      fit_coef(0, 14) = 4.7432936777438393E+00_dp
      fit_coef(1, 14) = 2.3082382103412251E+00_dp
      fit_coef(2, 14) = 1.4124938166206138E+00_dp
      fit_coef(3, 14) = 3.5283333252051396E-01_dp
      fit_coef(4, 14) = -7.7696363994240056E-02_dp
      fit_coef(5, 14) = -2.4914689387708116E-02_dp
      fit_coef(6, 14) = 5.4698262510218663E-02_dp
      fit_coef(7, 14) = 1.2640616949106807E-02_dp
      fit_coef(8, 14) = -1.9941323208816249E-02_dp
      fit_coef(9, 14) = -2.3878392941469071E-02_dp
      fit_coef(10, 14) = 4.1055603040417131E-02_dp
      fit_coef(11, 14) = -2.2734874408993563E-02_dp
      fit_coef(12, 14) = 4.7400588393736497E-03_dp
      fit_coef(0, 15) = 8.0710112216342775E+00_dp
      fit_coef(1, 15) = 2.9055595870824633E+00_dp
      fit_coef(2, 15) = 1.6866060964430423E+00_dp
      fit_coef(3, 15) = 3.5581815059695143E-01_dp
      fit_coef(4, 15) = -1.2393245100043256E-01_dp
      fit_coef(5, 15) = -2.0058383552301023E-02_dp
      fit_coef(6, 15) = 7.4640422271172679E-02_dp
      fit_coef(7, 15) = 1.1638820236828815E-02_dp
      fit_coef(8, 15) = -3.0928127527058424E-02_dp
      fit_coef(9, 15) = -2.2331476973823049E-02_dp
      fit_coef(10, 15) = 4.8851141745373840E-02_dp
      fit_coef(11, 15) = -2.8681346469043273E-02_dp
      fit_coef(12, 15) = 6.1686347014912515E-03_dp
      fit_coef(0, 16) = 1.3032066476244642E-03_dp
      fit_coef(1, 16) = 3.8001657920981851E-03_dp
      fit_coef(2, 16) = 5.5378064823689123E-03_dp
      fit_coef(3, 16) = 4.8206577899485006E-03_dp
      fit_coef(4, 16) = 2.2930144045863621E-03_dp
      fit_coef(5, 16) = 1.0986261495124486E-04_dp
      fit_coef(6, 16) = -5.1106328320878528E-04_dp
      fit_coef(7, 16) = -1.2833639292522648E-04_dp
      fit_coef(8, 16) = 1.1356979226852716E-04_dp
      fit_coef(9, 16) = 1.3973942184720303E-04_dp
      fit_coef(10, 16) = -1.8830360701542777E-04_dp
      fit_coef(11, 16) = 8.7997162572425538E-05_dp
      fit_coef(12, 16) = -1.5788992466518323E-05_dp
      fit_coef(0, 17) = 3.3544869799214628E-03_dp
      fit_coef(1, 17) = 9.3282083632741895E-03_dp
      fit_coef(2, 17) = 1.3271370703308300E-02_dp
      fit_coef(3, 17) = 1.1319835500098886E-02_dp
      fit_coef(4, 17) = 5.2524178309963383E-03_dp
      fit_coef(5, 17) = 1.8435087221054173E-04_dp
      fit_coef(6, 17) = -1.1797404048471125E-03_dp
      fit_coef(7, 17) = -2.7516373364939583E-04_dp
      fit_coef(8, 17) = 2.7734414812233385E-04_dp
      fit_coef(9, 17) = 2.9073439175898594E-04_dp
      fit_coef(10, 17) = -4.1550528469613449E-04_dp
      fit_coef(11, 17) = 1.9807893770328183E-04_dp
      fit_coef(12, 17) = -3.6025500379600989E-05_dp
      fit_coef(0, 18) = 6.3980511048867201E-03_dp
      fit_coef(1, 18) = 1.6128813328333793E-02_dp
      fit_coef(2, 18) = 2.1881033183460945E-02_dp
      fit_coef(3, 18) = 1.8021417293275506E-02_dp
      fit_coef(4, 18) = 8.0169003498435417E-03_dp
      fit_coef(5, 18) = 8.7835720128096682E-05_dp
      fit_coef(6, 18) = -1.8255106100472426E-03_dp
      fit_coef(7, 18) = -3.5831804749509092E-04_dp
      fit_coef(8, 18) = 4.8285099003593436E-04_dp
      fit_coef(9, 18) = 3.2465264639707240E-04_dp
      fit_coef(10, 18) = -5.5521076066674020E-04_dp
      fit_coef(11, 18) = 2.7719905231653717E-04_dp
      fit_coef(12, 18) = -5.1820236344809990E-05_dp
      fit_coef(0, 19) = 1.1614376264492066E-02_dp
      fit_coef(1, 19) = 2.5394323464202898E-02_dp
      fit_coef(2, 19) = 3.1733081124309684E-02_dp
      fit_coef(3, 19) = 2.4647015498881791E-02_dp
      fit_coef(4, 19) = 1.0500322632906088E-02_dp
      fit_coef(5, 19) = -1.2365877291301032E-05_dp
      fit_coef(6, 19) = -2.4904038430234622E-03_dp
      fit_coef(7, 19) = -4.8516369798151703E-04_dp
      fit_coef(8, 19) = 6.1030277836899595E-04_dp
      fit_coef(9, 19) = 7.4056213117855896E-04_dp
      fit_coef(10, 19) = -1.1330531422737624E-03_dp
      fit_coef(11, 19) = 5.8028785292689455E-04_dp
      fit_coef(12, 19) = -1.1292938481673382E-04_dp
      fit_coef(0, 20) = 2.0921907643278354E-02_dp
      fit_coef(1, 20) = 3.9009588349526253E-02_dp
      fit_coef(2, 20) = 4.3627398050961781E-02_dp
      fit_coef(3, 20) = 3.0443034676393681E-02_dp
      fit_coef(4, 20) = 1.1774589070303968E-02_dp
      fit_coef(5, 20) = 1.7741503302715373E-04_dp
      fit_coef(6, 20) = -2.3187359164888160E-03_dp
      fit_coef(7, 20) = -6.1637368205092694E-04_dp
      fit_coef(8, 20) = -5.3221714109514493E-05_dp
      fit_coef(9, 20) = 8.1942339340144803E-04_dp
      fit_coef(10, 20) = -2.2354982629820732E-04_dp
      fit_coef(11, 20) = -2.0305492301268794E-04_dp
      fit_coef(12, 20) = 9.4193850816527892E-05_dp
      fit_coef(0, 21) = 3.7306315930286466E-02_dp
      fit_coef(1, 21) = 5.9218420616905736E-02_dp
      fit_coef(2, 21) = 5.8868131527264646E-02_dp
      fit_coef(3, 21) = 3.5324612409341802E-02_dp
      fit_coef(4, 21) = 1.0684265753661753E-02_dp
      fit_coef(5, 21) = -4.9371106930299301E-04_dp
      fit_coef(6, 21) = -8.1292928450316614E-04_dp
      fit_coef(7, 21) = -5.1619725924012061E-04_dp
      fit_coef(8, 21) = 1.3804573307375741E-03_dp
      fit_coef(9, 21) = -3.4663308659499629E-03_dp
      fit_coef(10, 21) = 3.3745183791109305E-03_dp
      fit_coef(11, 21) = -1.4291379899595120E-03_dp
      fit_coef(12, 21) = 2.2454897381750800E-04_dp
      fit_coef(0, 22) = 6.5528592554674012E-02_dp
      fit_coef(1, 22) = 8.8324441729528691E-02_dp
      fit_coef(2, 22) = 7.8339185211631426E-02_dp
      fit_coef(3, 22) = 3.9448552741652813E-02_dp
      fit_coef(4, 22) = 7.3083551980621516E-03_dp
      fit_coef(5, 22) = -2.3247267042822940E-03_dp
      fit_coef(6, 22) = 4.5294755854382923E-04_dp
      fit_coef(7, 22) = 1.2969150143151570E-03_dp
      fit_coef(8, 22) = 4.6541067893071160E-04_dp
      fit_coef(9, 22) = -2.5032569665157803E-03_dp
      fit_coef(10, 22) = 1.1684561968462262E-03_dp
      fit_coef(11, 22) = 1.5370128957771531E-04_dp
      fit_coef(12, 22) = -1.4620664545141759E-04_dp
      fit_coef(0, 23) = 1.1324428651781197E-01_dp
      fit_coef(1, 23) = 1.2859806173609753E-01_dp
      fit_coef(2, 23) = 1.0214609707305304E-01_dp
      fit_coef(3, 23) = 4.2253398459926365E-02_dp
      fit_coef(4, 23) = 1.9762536371264882E-03_dp
      fit_coef(5, 23) = -4.5714996217850632E-03_dp
      fit_coef(6, 23) = 1.8868500721209416E-03_dp
      fit_coef(7, 23) = 2.5695589529597719E-03_dp
      fit_coef(8, 23) = 2.3850425890078301E-04_dp
      fit_coef(9, 23) = -3.0200640644662698E-03_dp
      fit_coef(10, 23) = 1.9503667616499671E-03_dp
      fit_coef(11, 23) = -5.5257674937112412E-04_dp
      fit_coef(12, 23) = 7.8113401323294423E-05_dp
      fit_coef(0, 24) = 1.9267245469427555E-01_dp
      fit_coef(1, 24) = 1.8189882938781768E-01_dp
      fit_coef(2, 24) = 1.2960901455552834E-01_dp
      fit_coef(3, 24) = 4.2611091038268031E-02_dp
      fit_coef(4, 24) = -5.1926218101019609E-03_dp
      fit_coef(5, 24) = -6.3425612842204224E-03_dp
      fit_coef(6, 24) = 4.0010538098000379E-03_dp
      fit_coef(7, 24) = 3.1214642561276981E-03_dp
      fit_coef(8, 24) = 1.0563824093820107E-04_dp
      fit_coef(9, 24) = -4.5110524003026762E-03_dp
      fit_coef(10, 24) = 4.2402114420073420E-03_dp
      fit_coef(11, 24) = -1.9163626754183359E-03_dp
      fit_coef(12, 24) = 3.6969348039378312E-04_dp
      fit_coef(0, 25) = 3.2317976642569807E-01_dp
      fit_coef(1, 25) = 2.4881850229416846E-01_dp
      fit_coef(2, 25) = 1.5917736470901570E-01_dp
      fit_coef(3, 25) = 3.9343328717200927E-02_dp
      fit_coef(4, 25) = -1.3902519808734082E-02_dp
      fit_coef(5, 25) = -6.8640926105066162E-03_dp
      fit_coef(6, 25) = 6.9006071006695079E-03_dp
      fit_coef(7, 25) = 3.2927740290179990E-03_dp
      fit_coef(8, 25) = -1.0983902176266030E-03_dp
      fit_coef(9, 25) = -5.0555007528979004E-03_dp
      fit_coef(10, 25) = 6.0190910606148324E-03_dp
      fit_coef(11, 25) = -2.9821080457428841E-03_dp
      fit_coef(12, 25) = 5.8221139171923042E-04_dp
      fit_coef(0, 26) = 5.3554089995575782E-01_dp
      fit_coef(1, 26) = 3.2739155636224210E-01_dp
      fit_coef(2, 26) = 1.8848689258384349E-01_dp
      fit_coef(3, 26) = 3.1884980883658924E-02_dp
      fit_coef(4, 26) = -2.3339754433360369E-02_dp
      fit_coef(5, 26) = -5.6037573210389264E-03_dp
      fit_coef(6, 26) = 1.0403698986274485E-02_dp
      fit_coef(7, 26) = 2.8865920265869838E-03_dp
      fit_coef(8, 26) = -3.2033768610669569E-03_dp
      fit_coef(9, 26) = -4.5416607006448602E-03_dp
      fit_coef(10, 26) = 7.3043246118188162E-03_dp
      fit_coef(11, 26) = -3.9689298970288006E-03_dp
      fit_coef(12, 26) = 8.1056551195733142E-04_dp
      fit_coef(0, 27) = 8.8000603036746139E-01_dp
      fit_coef(1, 27) = 4.1182532361787932E-01_dp
      fit_coef(2, 27) = 2.1497937217091584E-01_dp
      fit_coef(3, 27) = 2.1108047064593377E-02_dp
      fit_coef(4, 27) = -3.2191882260831459E-02_dp
      fit_coef(5, 27) = -2.5724974235411323E-03_dp
      fit_coef(6, 27) = 1.3955680448193437E-02_dp
      fit_coef(7, 27) = 1.8152463536498431E-03_dp
      fit_coef(8, 27) = -5.8233926497395817E-03_dp
      fit_coef(9, 27) = -2.9482077294913736E-03_dp
      fit_coef(10, 27) = 7.9261411558925292E-03_dp
      fit_coef(11, 27) = -4.8132854344172145E-03_dp
      fit_coef(12, 27) = 1.0440601147544259E-03_dp
      fit_coef(0, 28) = 1.4460627242259527E+00_dp
      fit_coef(1, 28) = 4.9271419953901929E-01_dp
      fit_coef(2, 28) = 2.3743318711644951E-01_dp
      fit_coef(3, 28) = 9.9076815053936729E-03_dp
      fit_coef(4, 28) = -3.9191794449869852E-02_dp
      fit_coef(5, 28) = 1.3099534113565326E-03_dp
      fit_coef(6, 28) = 1.6894034568497723E-02_dp
      fit_coef(7, 28) = 3.7043982265874470E-04_dp
      fit_coef(8, 28) = -8.3968062748629741E-03_dp
      fit_coef(9, 28) = -6.6028855186588131E-04_dp
      fit_coef(10, 28) = 7.8140284723942129E-03_dp
      fit_coef(11, 28) = -5.3528674153974880E-03_dp
      fit_coef(12, 28) = 1.2297837544111832E-03_dp
      fit_coef(0, 29) = 2.4319559958545311E+00_dp
      fit_coef(1, 29) = 5.6229813567117626E-01_dp
      fit_coef(2, 29) = 2.5870101698993758E-01_dp
      fit_coef(3, 29) = 2.8306941439785287E-03_dp
      fit_coef(4, 29) = -4.4047429705758694E-02_dp
      fit_coef(5, 29) = 4.3384122634558238E-03_dp
      fit_coef(6, 29) = 1.8978692315301807E-02_dp
      fit_coef(7, 29) = -7.9790741322960021E-04_dp
      fit_coef(8, 29) = -1.0330971785215522E-02_dp
      fit_coef(9, 29) = 1.2390333347493327E-03_dp
      fit_coef(10, 29) = 7.5464073935277154E-03_dp
      fit_coef(11, 29) = -5.6851817992854548E-03_dp
      fit_coef(12, 29) = 1.3586302192394280E-03_dp
      fit_coef(0, 30) = 4.5666853091127635E+00_dp
      fit_coef(1, 30) = 6.4009785574659228E-01_dp
      fit_coef(2, 30) = 2.9587958616654764E-01_dp
      fit_coef(3, 30) = 5.3258996536484551E-03_dp
      fit_coef(4, 30) = -4.8921450862390199E-02_dp
      fit_coef(5, 30) = 4.8519904409378041E-03_dp
      fit_coef(6, 30) = 2.1138966562008546E-02_dp
      fit_coef(7, 30) = -9.9509468782633297E-04_dp
      fit_coef(8, 30) = -1.1615586429846661E-02_dp
      fit_coef(9, 30) = 1.4899160811089988E-03_dp
      fit_coef(10, 30) = 8.3912654254730754E-03_dp
      fit_coef(11, 30) = -6.3525637249113569E-03_dp
      fit_coef(12, 30) = 1.5228141862253625E-03_dp

   END SUBROUTINE get_fit_coef_k15_mi

! **************************************************************************************************
!> \brief ...
!> \param k ...
!> \param fit_coef ...
!> \param L_b ...
!> \param U_b ...
! **************************************************************************************************
   SUBROUTINE get_fit_coef_k15_lo(k, fit_coef, L_b, U_b)
      INTEGER                                            :: k
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef
      REAL(KIND=dp)                                      :: L_b, U_b

      L_b = 9.3527642781376642E+01_dp
      U_b = 1.1750000000000000E+06_dp
      fit_coef(0, 1) = 3.2793279985494519E-06_dp
      fit_coef(1, 1) = 7.2761221861298733E-06_dp
      fit_coef(2, 1) = 6.6879597248076629E-04_dp
      fit_coef(3, 1) = -3.3377387242640721E-03_dp
      fit_coef(4, 1) = 3.6932190908574580E-02_dp
      fit_coef(5, 1) = -1.4859062400069553E-01_dp
      fit_coef(6, 1) = 4.4629627390711790E-01_dp
      fit_coef(7, 1) = -9.3682363031214910E-01_dp
      fit_coef(8, 1) = 1.3525709274164297E+00_dp
      fit_coef(9, 1) = -1.3157143856436124E+00_dp
      fit_coef(10, 1) = 8.2321722743838222E-01_dp
      fit_coef(11, 1) = -2.9905777894458391E-01_dp
      fit_coef(12, 1) = 4.7927187586119996E-02_dp
      fit_coef(0, 2) = 3.6131623637423831E-05_dp
      fit_coef(1, 2) = 4.7106962333648227E-05_dp
      fit_coef(2, 2) = 4.3052245932532632E-03_dp
      fit_coef(3, 2) = -2.3337397086091179E-02_dp
      fit_coef(4, 2) = 2.3805631448383252E-01_dp
      fit_coef(5, 2) = -9.8141547342868762E-01_dp
      fit_coef(6, 2) = 2.9514193901676489E+00_dp
      fit_coef(7, 2) = -6.1764264357918330E+00_dp
      fit_coef(8, 2) = 8.8882482060089814E+00_dp
      fit_coef(9, 2) = -8.6216924305612981E+00_dp
      fit_coef(10, 2) = 5.3820183799751744E+00_dp
      fit_coef(11, 2) = -1.9515484948325068E+00_dp
      fit_coef(12, 2) = 3.1228748667302658E-01_dp
      fit_coef(0, 3) = 1.9135313464109315E-04_dp
      fit_coef(1, 3) = 1.7460510061267730E-04_dp
      fit_coef(2, 3) = 1.5413483600038124E-02_dp
      fit_coef(3, 3) = -8.7803353742691062E-02_dp
      fit_coef(4, 3) = 8.0147578597001090E-01_dp
      fit_coef(5, 3) = -3.4506170283047766E+00_dp
      fit_coef(6, 3) = 1.0597559824370764E+01_dp
      fit_coef(7, 3) = -2.2472933987251839E+01_dp
      fit_coef(8, 3) = 3.2637358630788739E+01_dp
      fit_coef(9, 3) = -3.1863935196267601E+01_dp
      fit_coef(10, 3) = 1.9983399971045149E+01_dp
      fit_coef(11, 3) = -7.2707478368141611E+00_dp
      fit_coef(12, 3) = 1.1664152307259732E+00_dp
      fit_coef(0, 4) = 7.3883732820766948E-04_dp
      fit_coef(1, 4) = 5.0733833618914323E-04_dp
      fit_coef(2, 4) = 4.4370796269257232E-02_dp
      fit_coef(3, 4) = -2.6481196747393215E-01_dp
      fit_coef(4, 4) = 2.2595142950561815E+00_dp
      fit_coef(5, 4) = -1.0181630118704184E+01_dp
      fit_coef(6, 4) = 3.1679590172936724E+01_dp
      fit_coef(7, 4) = -6.7336030116199026E+01_dp
      fit_coef(8, 4) = 9.7677953648704644E+01_dp
      fit_coef(9, 4) = -9.5149794831930734E+01_dp
      fit_coef(10, 4) = 5.9530024414743906E+01_dp
      fit_coef(11, 4) = -2.1611138637657191E+01_dp
      fit_coef(12, 4) = 3.4601757819555536E+00_dp
      fit_coef(0, 5) = 2.3666313203208780E-03_dp
      fit_coef(1, 5) = 1.2786885886632092E-03_dp
      fit_coef(2, 5) = 1.1076208254179676E-01_dp
      fit_coef(3, 5) = -6.7839594619191412E-01_dp
      fit_coef(4, 5) = 5.5114285011328183E+00_dp
      fit_coef(5, 5) = -2.5395097951060229E+01_dp
      fit_coef(6, 5) = 7.9307253014224088E+01_dp
      fit_coef(7, 5) = -1.6859211837607822E+02_dp
      fit_coef(8, 5) = 2.4452261050037967E+02_dp
      fit_coef(9, 5) = -2.3820141119829208E+02_dp
      fit_coef(10, 5) = 1.4905893976672397E+02_dp
      fit_coef(11, 5) = -5.4127909713308576E+01_dp
      fit_coef(12, 5) = 8.6691001852426552E+00_dp
      fit_coef(0, 6) = 6.6782258468515501E-03_dp
      fit_coef(1, 6) = 2.9071241132673617E-03_dp
      fit_coef(2, 6) = 2.4975031064405542E-01_dp
      fit_coef(3, 6) = -1.5562125797528150E+00_dp
      fit_coef(4, 6) = 1.2185079476552678E+01_dp
      fit_coef(5, 6) = -5.6926325045922916E+01_dp
      fit_coef(6, 6) = 1.7855104053066674E+02_dp
      fit_coef(7, 6) = -3.8050660514368747E+02_dp
      fit_coef(8, 6) = 5.5288905790672220E+02_dp
      fit_coef(9, 6) = -5.3932598394454135E+02_dp
      fit_coef(10, 6) = 3.3782483380582988E+02_dp
      fit_coef(11, 6) = -1.2276125245641836E+02_dp
      fit_coef(12, 6) = 1.9671382604642698E+01_dp
      fit_coef(0, 7) = 1.7168691373333111E-02_dp
      fit_coef(1, 7) = 6.0978405270416770E-03_dp
      fit_coef(2, 7) = 5.2025306168908569E-01_dp
      fit_coef(3, 7) = -3.2821022955670651E+00_dp
      fit_coef(4, 7) = 2.4958349575372374E+01_dp
      fit_coef(5, 7) = -1.1765023272930763E+02_dp
      fit_coef(6, 7) = 3.7056008209422993E+02_dp
      fit_coef(7, 7) = -7.9184187030173371E+02_dp
      fit_coef(8, 7) = 1.1527200840436421E+03_dp
      fit_coef(9, 7) = -1.1258727391755413E+03_dp
      fit_coef(10, 7) = 7.0584825509402651E+02_dp
      fit_coef(11, 7) = -2.5665682670985007E+02_dp
      fit_coef(12, 7) = 4.1145879034693607E+01_dp
      fit_coef(0, 8) = 4.1066740219417780E-02_dp
      fit_coef(1, 8) = 1.1965752322210791E-02_dp
      fit_coef(2, 8) = 1.0148791234031664E+00_dp
      fit_coef(3, 8) = -6.4624725502258276E+00_dp
      fit_coef(4, 8) = 4.7976535067054421E+01_dp
      fit_coef(5, 8) = -2.2750659978981074E+02_dp
      fit_coef(6, 8) = 7.1929276146624795E+02_dp
      fit_coef(7, 8) = -1.5410667800156759E+03_dp
      fit_coef(8, 8) = 2.2474354642090384E+03_dp
      fit_coef(9, 8) = -2.1978239495400421E+03_dp
      fit_coef(10, 8) = 1.3791127707578614E+03_dp
      fit_coef(11, 8) = -5.0179288602294997E+02_dp
      fit_coef(12, 8) = 8.0484577079207980E+01_dp
      fit_coef(0, 9) = 9.2699687800883265E-02_dp
      fit_coef(1, 9) = 2.2155471885187787E-02_dp
      fit_coef(2, 9) = 1.8695902697850499E+00_dp
      fit_coef(3, 9) = -1.1991899383860448E+01_dp
      fit_coef(4, 9) = 8.7242978410239118E+01_dp
      fit_coef(5, 9) = -4.1536938867885857E+02_dp
      fit_coef(6, 9) = 1.3176623596141608E+03_dp
      fit_coef(7, 9) = -2.8299350356509945E+03_dp
      fit_coef(8, 9) = 4.1340902247522199E+03_dp
      fit_coef(9, 9) = -4.0476690586057271E+03_dp
      fit_coef(10, 9) = 2.5420687410683036E+03_dp
      fit_coef(11, 9) = -9.2553103476606759E+02_dp
      fit_coef(12, 9) = 1.4852229538367158E+02_dp
      fit_coef(0, 10) = 1.9947743279902694E-01_dp
      fit_coef(1, 10) = 3.8899312951373388E-02_dp
      fit_coef(2, 10) = 3.2680572849945402E+00_dp
      fit_coef(3, 10) = -2.1084081716614804E+01_dp
      fit_coef(4, 10) = 1.5075557198241123E+02_dp
      fit_coef(5, 10) = -7.1967411375827533E+02_dp
      fit_coef(6, 10) = 2.2897159192976455E+03_dp
      fit_coef(7, 10) = -4.9285687810081845E+03_dp
      fit_coef(8, 10) = 7.2112684192795978E+03_dp
      fit_coef(9, 10) = -7.0684961344198391E+03_dp
      fit_coef(10, 10) = 4.4429042040785398E+03_dp
      fit_coef(11, 10) = -1.6185943473816480E+03_dp
      fit_coef(12, 10) = 2.5986231052596014E+02_dp
      fit_coef(0, 11) = 4.1232198870536146E-01_dp
      fit_coef(1, 11) = 6.4893474864739994E-02_dp
      fit_coef(2, 11) = 5.4310004234894391E+00_dp
      fit_coef(3, 11) = -3.5204004815099751E+01_dp
      fit_coef(4, 11) = 2.4798012213498384E+02_dp
      fit_coef(5, 11) = -1.1858431117066064E+03_dp
      fit_coef(6, 11) = 3.7825064907174756E+03_dp
      fit_coef(7, 11) = -8.1581937192205814E+03_dp
      fit_coef(8, 11) = 1.1954131106568948E+04_dp
      fit_coef(9, 11) = -1.1729786611671518E+04_dp
      fit_coef(10, 11) = 7.3784592386581926E+03_dp
      fit_coef(11, 11) = -2.6896086736495204E+03_dp
      fit_coef(12, 11) = 4.3200462908615089E+02_dp
      fit_coef(0, 12) = 8.2375109973285454E-01_dp
      fit_coef(1, 12) = 1.0280784827845055E-01_dp
      fit_coef(2, 12) = 8.5755836263465604E+00_dp
      fit_coef(3, 12) = -5.5800839220514590E+01_dp
      fit_coef(4, 12) = 3.8804437171318665E+02_dp
      fit_coef(5, 12) = -1.8575547029234676E+03_dp
      fit_coef(6, 12) = 5.9379886844425691E+03_dp
      fit_coef(7, 12) = -1.2830180852831514E+04_dp
      fit_coef(8, 12) = 1.8824783357291468E+04_dp
      fit_coef(9, 12) = -1.8489266233281007E+04_dp
      fit_coef(10, 12) = 1.1638691611847828E+04_dp
      fit_coef(11, 12) = -4.2448366640433251E+03_dp
      fit_coef(12, 12) = 6.8208814863958401E+02_dp
      fit_coef(0, 13) = 1.6005401729604098E+00_dp
      fit_coef(1, 13) = 1.5427671373341209E-01_dp
      fit_coef(2, 13) = 1.2832947943920599E+01_dp
      fit_coef(3, 13) = -8.3759357589093241E+01_dp
      fit_coef(4, 13) = 5.7621942259264915E+02_dp
      fit_coef(5, 13) = -2.7598336160779349E+03_dp
      fit_coef(6, 13) = 8.8381623335643399E+03_dp
      fit_coef(7, 13) = -1.9126056849476197E+04_dp
      fit_coef(8, 13) = 2.8094638893959345E+04_dp
      fit_coef(9, 13) = -2.7617301955762628E+04_dp
      fit_coef(10, 13) = 1.7395609257161428E+04_dp
      fit_coef(11, 13) = -6.3475238166994641E+03_dp
      fit_coef(12, 13) = 1.0203400395817024E+03_dp
      fit_coef(0, 14) = 3.0540489317639166E+00_dp
      fit_coef(1, 14) = 2.1865700868129184E-01_dp
      fit_coef(2, 14) = 1.8148202718663860E+01_dp
      fit_coef(3, 14) = -1.1872710055025075E+02_dp
      fit_coef(4, 14) = 8.0984788084197976E+02_dp
      fit_coef(5, 14) = -3.8796556939052798E+03_dp
      fit_coef(6, 14) = 1.2441638354776058E+04_dp
      fit_coef(7, 14) = -2.6957427016371610E+04_dp
      fit_coef(8, 14) = 3.9635425931543767E+04_dp
      fit_coef(9, 14) = -3.8989060841935432E+04_dp
      fit_coef(10, 14) = 2.4571211983659101E+04_dp
      fit_coef(11, 14) = -8.9693909543608097E+03_dp
      fit_coef(12, 14) = 1.4422380562945843E+03_dp
      fit_coef(0, 15) = 5.8836496088258903E+00_dp
      fit_coef(1, 15) = 2.9333916159882772E-01_dp
      fit_coef(2, 15) = 2.4309148431266522E+01_dp
      fit_coef(3, 15) = -1.5928663857255199E+02_dp
      fit_coef(4, 15) = 1.0800144511418321E+03_dp
      fit_coef(5, 15) = -5.1742775717293125E+03_dp
      fit_coef(6, 15) = 1.6609418151702128E+04_dp
      fit_coef(7, 15) = -3.6019371990459913E+04_dp
      fit_coef(8, 15) = 5.2994637899369998E+04_dp
      fit_coef(9, 15) = -5.2156366955415026E+04_dp
      fit_coef(10, 15) = 3.2881619526078182E+04_dp
      fit_coef(11, 15) = -1.2006412191447664E+04_dp
      fit_coef(12, 15) = 1.9310054306610489E+03_dp
      fit_coef(0, 16) = 1.0490015508948160E-05_dp
      fit_coef(1, 16) = 1.9245918297480574E-05_dp
      fit_coef(2, 16) = 1.7514019628420299E-03_dp
      fit_coef(3, 16) = -8.7825475436631283E-03_dp
      fit_coef(4, 16) = 9.8132135498064849E-02_dp
      fit_coef(5, 16) = -4.0139357711585455E-01_dp
      fit_coef(6, 16) = 1.2151507236457464E+00_dp
      fit_coef(7, 16) = -2.5608071661985146E+00_dp
      fit_coef(8, 16) = 3.7049308475400640E+00_dp
      fit_coef(9, 16) = -3.6079085144794738E+00_dp
      fit_coef(10, 16) = 2.2586411450421502E+00_dp
      fit_coef(11, 16) = -8.2071784242647206E-01_dp
      fit_coef(12, 16) = 1.3153768317036563E-01_dp
      fit_coef(0, 17) = 6.9018219767429348E-05_dp
      fit_coef(1, 17) = 6.9648188439671224E-05_dp
      fit_coef(2, 17) = 6.2459333601416062E-03_dp
      fit_coef(3, 17) = -3.5441557914769072E-02_dp
      fit_coef(4, 17) = 3.2906922405173999E-01_dp
      fit_coef(5, 17) = -1.3703072930999050E+00_dp
      fit_coef(6, 17) = 4.1206399805320482E+00_dp
      fit_coef(7, 17) = -8.6245400433435595E+00_dp
      fit_coef(8, 17) = 1.2424942794117596E+01_dp
      fit_coef(9, 17) = -1.2071284705492808E+01_dp
      fit_coef(10, 17) = 7.5482859678406982E+00_dp
      fit_coef(11, 17) = -2.7416411252527824E+00_dp
      fit_coef(12, 17) = 4.3940704181003898E-01_dp
      fit_coef(0, 18) = 2.8267117174337328E-04_dp
      fit_coef(1, 18) = 2.0340536143980756E-04_dp
      fit_coef(2, 18) = 1.7678726727191417E-02_dp
      fit_coef(3, 18) = -1.0545298250093989E-01_dp
      fit_coef(4, 18) = 8.9650844672191476E-01_dp
      fit_coef(5, 18) = -4.0700645201639842E+00_dp
      fit_coef(6, 18) = 1.2778228701375353E+01_dp
      fit_coef(7, 18) = -2.7365134031898545E+01_dp
      fit_coef(8, 18) = 3.9913427510994921E+01_dp
      fit_coef(9, 18) = -3.9026449660683376E+01_dp
      fit_coef(10, 18) = 2.4477897493830259E+01_dp
      fit_coef(11, 18) = -8.9007828894477790E+00_dp
      fit_coef(12, 18) = 1.4266184467905938E+00_dp
      fit_coef(0, 19) = 9.1918738520005492E-04_dp
      fit_coef(1, 19) = 5.0032919937804867E-04_dp
      fit_coef(2, 19) = 4.3364999825556638E-02_dp
      fit_coef(3, 19) = -2.6781408638943244E-01_dp
      fit_coef(4, 19) = 2.1566125289358995E+00_dp
      fit_coef(5, 19) = -1.0022123878208944E+01_dp
      fit_coef(6, 19) = 3.1317396694353242E+01_dp
      fit_coef(7, 19) = -6.6447554440483501E+01_dp
      fit_coef(8, 19) = 9.6141728832817350E+01_dp
      fit_coef(9, 19) = -9.3449718908958914E+01_dp
      fit_coef(10, 19) = 5.8372166225769881E+01_dp
      fit_coef(11, 19) = -2.1166765338703282E+01_dp
      fit_coef(12, 19) = 3.3863816402153009E+00_dp
      fit_coef(0, 20) = 2.5882618681224243E-03_dp
      fit_coef(1, 20) = 1.1101963800374212E-03_dp
      fit_coef(2, 20) = 9.5153090359631570E-02_dp
      fit_coef(3, 20) = -5.9667092259469967E-01_dp
      fit_coef(4, 20) = 4.6158300825488254E+00_dp
      fit_coef(5, 20) = -2.1702960038752170E+01_dp
      fit_coef(6, 20) = 6.8104607605264903E+01_dp
      fit_coef(7, 20) = -1.4509864087700180E+02_dp
      fit_coef(8, 20) = 2.1082070081369434E+02_dp
      fit_coef(9, 20) = -2.0567820908774044E+02_dp
      fit_coef(10, 20) = 1.2886369297377382E+02_dp
      fit_coef(11, 20) = -4.6839217031726335E+01_dp
      fit_coef(12, 20) = 7.5072963372495414E+00_dp
      fit_coef(0, 21) = 6.5882379082633398E-03_dp
      fit_coef(1, 21) = 2.2621691339110277E-03_dp
      fit_coef(2, 21) = 1.9246964638130989E-01_dp
      fit_coef(3, 21) = -1.2217276753530955E+00_dp
      fit_coef(4, 21) = 9.1698242992620767E+00_dp
      fit_coef(5, 21) = -4.3465513365240064E+01_dp
      fit_coef(6, 21) = 1.3710631049681336E+02_dp
      fit_coef(7, 21) = -2.9321690733250756E+02_dp
      fit_coef(8, 21) = 4.2710092920695109E+02_dp
      fit_coef(9, 21) = -4.1732542801681706E+02_dp
      fit_coef(10, 21) = 2.6170639243954514E+02_dp
      fit_coef(11, 21) = -9.5176513369112968E+01_dp
      fit_coef(12, 21) = 1.5259801532746318E+01_dp
      fit_coef(0, 22) = 1.5543900572102644E-02_dp
      fit_coef(1, 22) = 4.3029374070070658E-03_dp
      fit_coef(2, 22) = 3.6379929804956312E-01_dp
      fit_coef(3, 22) = -2.3299690981198706E+00_dp
      fit_coef(4, 22) = 1.7060773313313661E+01_dp
      fit_coef(5, 22) = -8.1273230384261126E+01_dp
      fit_coef(6, 22) = 2.5746624440333989E+02_dp
      fit_coef(7, 22) = -5.5230172621897168E+02_dp
      fit_coef(8, 22) = 8.0612603717585921E+02_dp
      fit_coef(9, 22) = -7.8876694465068897E+02_dp
      fit_coef(10, 22) = 4.9512974619778498E+02_dp
      fit_coef(11, 22) = -1.8020197227270759E+02_dp
      fit_coef(12, 22) = 2.8909017650939049E+01_dp
      fit_coef(0, 23) = 3.4544493883799224E-02_dp
      fit_coef(1, 23) = 7.7057465719973899E-03_dp
      fit_coef(2, 23) = 6.4801783632712762E-01_dp
      fit_coef(3, 23) = -4.1791996544114722E+00_dp
      fit_coef(4, 23) = 2.9972887604114607E+01_dp
      fit_coef(5, 23) = -1.4323654947892643E+02_dp
      fit_coef(6, 23) = 4.5543253426208781E+02_dp
      fit_coef(7, 23) = -9.7965395620159916E+02_dp
      fit_coef(8, 23) = 1.4326296354673486E+03_dp
      fit_coef(9, 23) = -1.4036958339933367E+03_dp
      fit_coef(10, 23) = 8.8201906011468520E+02_dp
      fit_coef(11, 23) = -3.2125145658138058E+02_dp
      fit_coef(12, 23) = 5.1566695258306126E+01_dp
      fit_coef(0, 24) = 7.3122924052638857E-02_dp
      fit_coef(1, 24) = 1.3047345576674729E-02_dp
      fit_coef(2, 24) = 1.0921563698572587E+00_dp
      fit_coef(3, 24) = -7.0828159501016552E+00_dp
      fit_coef(4, 24) = 4.9900195235442546E+01_dp
      fit_coef(5, 24) = -2.3893035920690767E+02_dp
      fit_coef(6, 24) = 7.6209216068295359E+02_dp
      fit_coef(7, 24) = -1.6433368571765945E+03_dp
      fit_coef(8, 24) = 2.4074464958116155E+03_dp
      fit_coef(9, 24) = -2.3618333223611826E+03_dp
      fit_coef(10, 24) = 1.4854588983048657E+03_dp
      fit_coef(11, 24) = -5.4141801556510109E+02_dp
      fit_coef(12, 24) = 8.6954202781474834E+01_dp
      fit_coef(0, 25) = 1.4863184827325007E-01_dp
      fit_coef(1, 25) = 2.0897369993486958E-02_dp
      fit_coef(2, 25) = 1.7422013268701704E+00_dp
      fit_coef(3, 25) = -1.1349892975305522E+01_dp
      fit_coef(4, 25) = 7.8731198009107700E+01_dp
      fit_coef(5, 25) = -3.7738777230483873E+02_dp
      fit_coef(6, 25) = 1.2069682708983435E+03_dp
      fit_coef(7, 25) = -2.6084306788827607E+03_dp
      fit_coef(8, 25) = 3.8275166346901069E+03_dp
      fit_coef(9, 25) = -3.7594561324724668E+03_dp
      fit_coef(10, 25) = 2.3665593844071400E+03_dp
      fit_coef(11, 25) = -8.6313167371072871E+02_dp
      fit_coef(12, 25) = 1.3869371175524410E+02_dp
      fit_coef(0, 26) = 2.9197038671042685E-01_dp
      fit_coef(1, 26) = 3.1553064146021624E-02_dp
      fit_coef(2, 26) = 2.6212810102963688E+00_dp
      fit_coef(3, 26) = -1.7140795787714239E+01_dp
      fit_coef(4, 26) = 1.1729936023843838E+02_dp
      fit_coef(5, 26) = -5.6251944415470587E+02_dp
      fit_coef(6, 26) = 1.8032102744423644E+03_dp
      fit_coef(7, 26) = -3.9047643960741279E+03_dp
      fit_coef(8, 26) = 5.7382582448223930E+03_dp
      fit_coef(9, 26) = -5.6424075077151520E+03_dp
      fit_coef(10, 26) = 3.5547702376861953E+03_dp
      fit_coef(11, 26) = -1.2973000963117572E+03_dp
      fit_coef(12, 26) = 2.0855898009649175E+02_dp
      fit_coef(0, 27) = 5.5768130937544869E-01_dp
      fit_coef(1, 27) = 4.4572128359602628E-02_dp
      fit_coef(2, 27) = 3.6916856472540940E+00_dp
      fit_coef(3, 27) = -2.4213381713827740E+01_dp
      fit_coef(4, 27) = 1.6378750132966567E+02_dp
      fit_coef(5, 27) = -7.8545555955911163E+02_dp
      fit_coef(6, 27) = 2.5226651741402102E+03_dp
      fit_coef(7, 27) = -5.4722074342493534E+03_dp
      fit_coef(8, 27) = 8.0523685521856414E+03_dp
      fit_coef(9, 27) = -7.9256723613351096E+03_dp
      fit_coef(10, 27) = 4.9969395283607182E+03_dp
      fit_coef(11, 27) = -1.8246427440641146E+03_dp
      fit_coef(12, 27) = 2.9346534721488899E+02_dp
      fit_coef(0, 28) = 1.0450173925119002E+00_dp
      fit_coef(1, 28) = 5.8265207981760117E-02_dp
      fit_coef(2, 28) = 4.8143926353121813E+00_dp
      fit_coef(3, 28) = -3.1648523514665637E+01_dp
      fit_coef(4, 28) = 2.1213334922148832E+02_dp
      fit_coef(5, 28) = -1.0169903912741421E+03_dp
      fit_coef(6, 28) = 3.2709939445720015E+03_dp
      fit_coef(7, 28) = -7.1053218328636922E+03_dp
      fit_coef(8, 28) = 1.0466780738282137E+04_dp
      fit_coef(9, 28) = -1.0310448956212835E+04_dp
      fit_coef(10, 28) = 6.5044600754848825E+03_dp
      fit_coef(11, 28) = -2.3762315782073229E+03_dp
      fit_coef(12, 28) = 3.8231992824865711E+02_dp
      fit_coef(0, 29) = 1.9639167282286836E+00_dp
      fit_coef(1, 29) = 6.9961149923850899E-02_dp
      fit_coef(2, 29) = 5.7725439445353004E+00_dp
      fit_coef(3, 29) = -3.7996769718433761E+01_dp
      fit_coef(4, 29) = 2.5327571805969885E+02_dp
      fit_coef(5, 29) = -1.2137919470627476E+03_dp
      fit_coef(6, 29) = 3.9071985303035226E+03_dp
      fit_coef(7, 29) = -8.4944480629188183E+03_dp
      fit_coef(8, 29) = 1.2521478684961125E+04_dp
      fit_coef(9, 29) = -1.2340746128541903E+04_dp
      fit_coef(10, 29) = 7.7883138533160081E+03_dp
      fit_coef(11, 29) = -2.8461039477634558E+03_dp
      fit_coef(12, 29) = 4.5802666012412590E+02_dp
      fit_coef(0, 30) = 4.0342493016500347E+00_dp
      fit_coef(1, 30) = 7.9647259264362835E-02_dp
      fit_coef(2, 30) = 6.5710181771750698E+00_dp
      fit_coef(3, 30) = -4.3255572360835060E+01_dp
      fit_coef(4, 30) = 2.8820713216838652E+02_dp
      fit_coef(5, 30) = -1.3810306152021003E+03_dp
      fit_coef(6, 30) = 4.4457135475846326E+03_dp
      fit_coef(7, 30) = -9.6658445565099519E+03_dp
      fit_coef(8, 30) = 1.4249062004096166E+04_dp
      fit_coef(9, 30) = -1.4044077033410062E+04_dp
      fit_coef(10, 30) = 8.8636369152484167E+03_dp
      fit_coef(11, 30) = -3.2391609850813516E+03_dp
      fit_coef(12, 30) = 5.2129451486951780E+02_dp

   END SUBROUTINE get_fit_coef_k15_lo

   SUBROUTINE get_minimax_coeff_k35(k,E_range, aw)
      INTEGER, INTENT(IN)                                :: k
      REAL(KIND=dp), INTENT(IN)                          :: E_range
      REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw

      CHARACTER(LEN=*), PARAMETER :: routineN = 'get_minimax_coeff_k35', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, icoef, kkk
      REAL(KIND=dp)                                      :: L_b, R_max, R_min, R_transf, U_b, E_ratio
      REAL(KIND=dp), DIMENSION(0:12, 2*k)                :: fit_coef

      CALL timeset(routineN, handle)

      SELECT CASE (k)
      CASE (14)
        CALL get_coeff_14(k,E_range, aw) 
        E_ratio = 1.0_dp
      CASE (20)
        CALL get_coeff_20(k,E_range, aw, E_ratio)
     CASE (26)
        CALL get_coeff_26(k,E_range, aw, E_ratio)
     CASE (28)
        CALL get_coeff_28(k,E_range, aw, E_ratio)
      CASE (30)
        CALL get_coeff_30(k,E_range, aw, E_ratio)
     CASE (32)
        CALL get_coeff_32(k,E_range, aw, E_ratio)
      CASE (34)
        CALL get_coeff_34(k,E_range, aw, E_ratio)
      CASE (40)
        CALL get_coeff_40(k,E_range, aw, E_ratio)
      CASE DEFAULT
        CPABORT("The grid size you choose is not available.")
      END SELECT

      CALL rescale_grid(aw, E_range, E_ratio)

      CALL timestop(handle)

   END SUBROUTINE get_minimax_coeff_k35

   SUBROUTINE rescale_grid(aw, E_range, E_ratio)
      REAL(KIND=dp), DIMENSION(:), INTENT(INOUT)         :: aw
      REAL(KIND=dp), INTENT(IN)                          :: E_range
      REAL(KIND=dp), INTENT(INOUT)                          :: E_ratio

        IF(E_ratio > 1.5_dp) THEN
           E_ratio = E_ratio/1.5_dp
        ELSE
           E_ratio = 1.0_dp
        END IF
        aw(:) = aw(:)*E_ratio
   END SUBROUTINE 

SUBROUTINE get_coeff_14(k,E_range, aw)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw

IF(E_range < 50.0_dp) THEN
aw(:) = (/ &
0.007033988371446608_dp, &
0.037678504999091575_dp, &
0.09548956929385226_dp, &
0.18576410803886806_dp, &
0.31779312425467404_dp, &
0.5068349720755181_dp, &
0.7766836181145695_dp, &
1.1628502651609418_dp, &
1.7169189049668556_dp, &
2.513446536933386_dp, &
3.6622655064876106_dp, &
5.333971089569839_dp, &
7.826870009792442_dp, &
11.830960416296039_dp, &
0.018120413151890245_dp, &
0.04360100180643627_dp, &
0.07288193546283167_dp, &
0.10918256185633184_dp, &
0.15737550912443138_dp, &
0.22462209746145878_dp, &
0.32092428074278245_dp, &
0.45989714348645805_dp, &
0.6604758550584349_dp, &
0.9505077347494197_dp, &
1.374650674447938_dp, &
2.0155110635874958_dp, &
3.069017646840953_dp, &
5.294222991564442_dp /)
ELSE IF(E_range < 70.0_dp) THEN
aw(:) = (/ &
0.005391739468383946_dp, &
0.028957973986523645_dp, &
0.0737635957484961_dp, &
0.14466061865089944_dp, &
0.2503475221317881_dp, &
0.4053440365234253_dp, &
0.6324986286416501_dp, &
0.9661482804281842_dp, &
1.456643036277785_dp, &
2.1775896461863633_dp, &
3.2385114950733938_dp, &
4.810489359565789_dp, &
7.192640772635373_dp, &
11.073307414356625_dp, &
0.013898118700938116_dp, &
0.033619560778216415_dp, &
0.056777870172026554_dp, &
0.08643262370919962_dp, &
0.12731941382060033_dp, &
0.18643081014866653_dp, &
0.273528857166679_dp, &
0.4020467586083139_dp, &
0.5910136673964349_dp, &
0.868753030246825_dp, &
1.280697527302208_dp, &
1.9101755501348812_dp, &
2.952757221284849_dp, &
5.161070005995557_dp /)
ELSE IF(E_range < 100.0_dp) THEN
aw(:) = (/ &
0.004056225721001927_dp, &
0.021853980756421913_dp, &
0.056010262388666146_dp, &
0.11092332380989224_dp, &
0.19465691497094606_dp, &
0.32087042569825136_dp, &
0.5112373338421633_dp, &
0.7986247335644407_dp, &
1.231854342106747_dp, &
1.8833270763532868_dp, &
2.8620732519320047_dp, &
4.339611162384821_dp, &
6.616052064053293_dp, &
10.378662932191107_dp, &
0.010463064243724859_dp, &
0.025471259475826978_dp, &
0.04355170754788517_dp, &
0.06758573094614688_dp, &
0.10211827198226289_dp, &
0.15387011617076077_dp, &
0.23226499084209778_dp, &
0.35052472852518396_dp, &
0.5277975650994836_dp, &
0.7929347444960495_dp, &
1.1922412514014464_dp, &
1.8099471809409415_dp, &
2.8415517904786345_dp, &
5.0338612773645_dp /)
ELSE IF(E_range < 150.0_dp) THEN
aw(:) = (/ &
0.0029266875231949217_dp, &
0.015832900639212024_dp, &
0.04090592916621096_dp, &
0.0820604388934972_dp, &
0.1466458204874079_dp, &
0.2472812109466511_dp, &
0.4041743469930729_dp, &
0.6483762669700555_dp, &
1.0268300157107442_dp, &
1.6104050888338957_dp, &
2.507393818232554_dp, &
3.8896889677065363_dp, &
6.058615752782893_dp, &
9.700887782842155_dp, &
0.007556357553370067_dp, &
0.018547248084085394_dp, &
0.03222889351582765_dp, &
0.05127407861898456_dp, &
0.07996241074785554_dp, &
0.12462423757363507_dp, &
0.19425555448146978_dp, &
0.3018431126033537_dp, &
0.46666447671905603_dp, &
0.7181432527747482_dp, &
1.103582208191214_dp, &
1.7083522104651254_dp, &
2.728194733457147_dp, &
4.904359240038412_dp /)
ELSE IF(E_range < 300.0_dp) THEN
aw(:) = (/ &
0.0016663612973034454_dp, &
0.009091717259223682_dp, &
0.02389076529827524_dp, &
0.04924571718453551_dp, &
0.09133751379391733_dp, &
0.16096837788144863_dp, &
0.2757734566832361_dp, &
0.46366231422195536_dp, &
0.7682882010207237_dp, &
1.2576817182975237_dp, &
2.0385092477578888_dp, &
3.2829554315525318_dp, &
5.294434379381133_dp, &
8.759899195427739_dp, &
0.004310505477369055_dp, &
0.010762996322301278_dp, &
0.01934470025182372_dp, &
0.03236862863070855_dp, &
0.053580265879908114_dp, &
0.08856509136964401_dp, &
0.14560788021089902_dp, &
0.23730645268536366_dp, &
0.383059153361771_dp, &
0.6131011365593926_dp, &
0.97634830552257_dp, &
1.5602761817745894_dp, &
2.5616606875462824_dp, &
4.714407511566948_dp /)
ELSE IF(E_range < 600.0_dp) THEN
aw(:) = (/ &
0.0009445998990812008_dp, &
0.005210651155173852_dp, &
0.013999402414600478_dp, &
0.02988166051456169_dp, &
0.057975133883678906_dp, &
0.1073721376238129_dp, &
0.19332195057980348_dp, &
0.3408201894032817_dp, &
0.5903550237743922_dp, &
1.0070204521515524_dp, &
1.695519594299679_dp, &
2.827895454690515_dp, &
4.709483696284409_dp, &
8.028442009626652_dp, &
0.002449376620556013_dp, &
0.006252542853267622_dp, &
0.011735382985043822_dp, &
0.020859383986502725_dp, &
0.03680823107584084_dp, &
0.06446909805807077_dp, &
0.11150332879929417_dp, &
0.19009448335407098_dp, &
0.3195868998353066_dp, &
0.5307744054805075_dp, &
0.8739865496048677_dp, &
1.4388463236357882_dp, &
2.4237121361472624_dp, &
4.557297798241484_dp /)
ELSE IF(E_range < 1000.0_dp) THEN
aw(:) = (/ &
0.0006207562813744647_dp, &
0.003459377863029188_dp, &
0.009488853456821827_dp, &
0.02090056537030721_dp, &
0.042115248700416004_dp, &
0.08110531624723105_dp, &
0.1515715986179278_dp, &
0.2765845840370867_dp, &
0.4944573625309604_dp, &
0.868149511162554_dp, &
1.5007907852924771_dp, &
2.5640831549559073_dp, &
4.364606432116151_dp, &
7.591761535524144_dp, &
0.0016132137383178026_dp, &
0.004203321620266427_dp, &
0.008205175630212513_dp, &
0.0153343881811187_dp, &
0.02838882910760442_dp, &
0.0518157628140097_dp, &
0.09286688482943865_dp, &
0.16339325961277845_dp, &
0.28259950605373146_dp, &
0.48154088288418917_dp, &
0.8114387566425757_dp, &
1.3634496533377698_dp, &
2.3373122626452476_dp, &
4.458980283848654_dp /)
ELSE IF(E_range < 1500.0_dp) THEN
aw(:) = (/ &
0.0004447117240649359_dp, &
0.002502620491254362_dp, &
0.007001217127396729_dp, &
0.015869138034069044_dp, &
0.03303252051770302_dp, &
0.06567630102923491_dp, &
0.1264106898345732_dp, &
0.2369221605741954_dp, &
0.433913322188407_dp, &
0.7787042235860775_dp, &
1.3731435087964499_dp, &
2.388549935479365_dp, &
4.1323745368100315_dp, &
7.295122378518082_dp, &
0.0011581417012826308_dp, &
0.0030770724446050005_dp, &
0.006227720254628213_dp, &
0.012142862774813431_dp, &
0.02334653227145054_dp, &
0.043983678992645645_dp, &
0.08100358344341836_dp, &
0.14598199685699062_dp, &
0.2579678079108053_dp, &
0.4481469223672607_dp, &
0.7683567661774289_dp, &
1.3109138762728403_dp, &
2.2767229786363394_dp, &
4.390060321284838_dp /)
ELSE IF(E_range < 3000.0_dp) THEN
aw(:) = (/ &
0.00025163198847530084_dp, &
0.0014463969416761689_dp, &
0.004219287423148291_dp, &
0.010119568024871942_dp, &
0.022352106772967215_dp, &
0.046965180678614486_dp, &
0.09497614329353518_dp, &
0.18599499900409297_dp, &
0.35423025910205114_dp, &
0.6583665510745594_dp, &
1.198074653757926_dp, &
2.1438543879816323_dp, &
3.8044080833851703_dp, &
6.872224069166493_dp, &
0.0006582781934131662_dp, &
0.0018238959809057428_dp, &
0.003968925506923664_dp, &
0.008347152859937043_dp, &
0.017091740119353008_dp, &
0.03391226583525337_dp, &
0.06528124819859281_dp, &
0.12229912316054362_dp, &
0.22368974115005452_dp, &
0.4007368476435938_dp, &
0.7061526953087729_dp, &
1.2340830417614292_dp, &
2.1874721941502657_dp, &
4.288559195580037_dp /)
ELSE IF(E_range < 10000.0_dp) THEN
aw(:) = (/ &
9.432511324026308e-05_dp, &
0.0005743353815977165_dp, &
0.001856650422988761_dp, &
0.005008545726183479_dp, &
0.012329135220747828_dp, &
0.028441738345270376_dp, &
0.062298909504289575_dp, &
0.13069908990466905_dp, &
0.2643171695899576_dp, &
0.5179084117040162_dp, &
0.9876522941820988_dp, &
1.8423984016594872_dp, &
3.392384106263474_dp, &
6.333419958135755_dp, &
0.0002497977615891995_dp, &
0.000772139390410612_dp, &
0.001961758632072843_dp, &
0.004707823627981478_dp, &
0.010671857313366648_dp, &
0.022989131222192767_dp, &
0.0474281951951997_dp, &
0.09432099604798555_dp, &
0.18176506471568102_dp, &
0.3409650914972661_dp, &
0.625693511144038_dp, &
1.1327364993622373_dp, &
2.0684084501370754_dp, &
4.15314908898696_dp /)
ELSE
aw(:) = (/ &
3.930087896834751e-05_dp, &
0.0002623153029429105_dp, &
0.0009655315231040697_dp, &
0.0029321053602831514_dp, &
0.007940136661372303_dp, &
0.019772653898919544_dp, &
0.04610737554045433_dp, &
0.10192723077591104_dp, &
0.21552367898796046_dp, &
0.43887309414092707_dp, &
0.8655244728542174_dp, &
1.6628682898817153_dp, &
3.1419925614121733_dp, &
6.001262665573696_dp, &
0.00010626069042738675_dp, &
0.00038450843879575536_dp, &
0.0011446717296045412_dp, &
0.0030688085105842634_dp, &
0.007547735704169931_dp, &
0.01734252877364762_dp, &
0.037728050299063745_dp, &
0.07846081827343339_dp, &
0.15711016798804078_dp, &
0.3046803544836245_dp, &
0.5755275554309414_dp, &
1.0682374844273534_dp, &
1.9917202734688404_dp, &
4.065896088787167_dp /)
END IF
END SUBROUTINE

SUBROUTINE get_coeff_20(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio

E_ratio = 1.0_dp

IF(E_range < 160.0_dp) THEN
E_ratio = 160.0_dp/E_range
aw(:) = (/ &
0.001914630397440357_dp, &
0.010213576613401968_dp, &
0.025680891057421672_dp, &
0.0493449319236052_dp, &
0.08293685134693599_dp, &
0.12919624646173086_dp, &
0.19228285335550452_dp, &
0.2782760660904403_dp, &
0.39576304374993926_dp, &
0.5565788644848801_dp, &
0.7768133573403438_dp, &
1.078218225901317_dp, &
1.490179997259346_dp, &
2.05252868829164_dp, &
2.8197112496475043_dp, &
3.867494708835042_dp, &
5.305079876713698_dp, &
7.300875843119114_dp, &
10.15147891967995_dp, &
14.555017842792799_dp, &
0.004927656223630981_dp, &
0.011758515904440842_dp, &
0.019344088213568048_dp, &
0.028267449850217268_dp, &
0.03936606124214366_dp, &
0.05383767100019924_dp, &
0.07333913313314458_dp, &
0.10006644199150115_dp, &
0.1368635582237694_dp, &
0.18742663507099822_dp, &
0.25663534145442135_dp, &
0.35102647523788305_dp, &
0.4794682812719825_dp, &
0.6541976749238269_dp, &
0.892609049205798_dp, &
1.2207785955858315_dp, &
1.6815214839963029_dp, &
2.356430029250044_dp, &
3.4442813701617325_dp, &
5.727108804374391_dp /)
ELSE IF(E_range < 200.0_dp) THEN
aw(:) = (/ &
0.0015934118083825806_dp, &
0.008509063899457114_dp, &
0.021438168906707876_dp, &
0.04132186256649876_dp, &
0.06976220373713116_dp, &
0.10932128608557805_dp, &
0.16391921285449756_dp, &
0.23931259785276118_dp, &
0.3436664997737205_dp, &
0.48829754434795736_dp, &
0.6887047833087856_dp, &
0.9660158225052187_dp, &
1.349010682139282_dp, &
1.876993803313851_dp, &
2.6040407880564294_dp, &
3.605771589889266_dp, &
4.991497711230883_dp, &
6.929920412883439_dp, &
9.717725946130605_dp, &
14.051264969666578_dp, &
0.00410193748953601_dp, &
0.00980905245350431_dp, &
0.016202681787789122_dp, &
0.023826737638106887_dp, &
0.033474256616571156_dp, &
0.0462892956900987_dp, &
0.06385666391145611_dp, &
0.08827885704274316_dp, &
0.12229816206369343_dp, &
0.16952335915562816_dp, &
0.23477924911915304_dp, &
0.32459276577928703_dp, &
0.4478793216120171_dp, &
0.6169923625376843_dp, &
0.8495179832070735_dp, &
1.1717837086333052_dp, &
1.6268520987741988_dp, &
2.296386678135975_dp, &
3.3785363803003445_dp, &
5.651133301145076_dp /)
ELSE IF(E_range < 260.0_dp) THEN
aw(:) = (/ &
0.0012826603057563102_dp, &
0.0068587668174154475_dp, &
0.01732436797752842_dp, &
0.03352575629145031_dp, &
0.05692284719721156_dp, &
0.08987844528789371_dp, &
0.1360371747143393_dp, &
0.20077865488399435_dp, &
0.29177372565168536_dp, &
0.41973393936327597_dp, &
0.5994685778078079_dp, &
0.851370768653423_dp, &
1.2034954614185025_dp, &
1.6945031865392985_dp, &
2.377995017380597_dp, &
3.3293762143330468_dp, &
4.658057758626519_dp, &
6.533087585658005_dp, &
9.25133327731215_dp, &
13.507327774875803_dp, &
0.0033029757503794653_dp, &
0.007919729203352899_dp, &
0.013149299092943129_dp, &
0.0194919631120812_dp, &
0.027689636064854993_dp, &
0.03882150421363747_dp, &
0.054383526381991135_dp, &
0.07636569576619472_dp, &
0.10739236075305046_dp, &
0.15097187360710013_dp, &
0.2118622303874801_dp, &
0.2965700471632371_dp, &
0.4140528877951041_dp, &
0.5767873242476079_dp, &
0.8025758446055958_dp, &
1.1180439213635482_dp, &
1.5665659762545827_dp, &
2.2299380047180746_dp, &
3.3056758845372864_dp, &
5.567044533414704_dp /)
ELSE IF(E_range < 340.0_dp) THEN
aw(:) = (/ &
0.0010264504520625885_dp, &
0.00549683691935651_dp, &
0.013923537057194292_dp, &
0.027064335354502266_dp, &
0.046244583174484356_dp, &
0.07363428927862167_dp, &
0.11260531721148681_dp, &
0.1681595294038165_dp, &
0.24747280128089963_dp, &
0.36065226956640584_dp, &
0.5218158488053328_dp, &
0.7506122529239347_dp, &
1.0743481991029966_dp, &
1.5310024789311656_dp, &
2.17365271086166_dp, &
3.0774393380007683_dp, &
4.351838203868952_dp, &
6.166247581716586_dp, &
8.817790032549357_dp, &
12.999411677186895_dp, &
0.002644099460350311_dp, &
0.006358730129774125_dp, &
0.010617864369369328_dp, &
0.015879974975063757_dp, &
0.022836181425139004_dp, &
0.03249818247175067_dp, &
0.04626870726650271_dp, &
0.0660233342341586_dp, &
0.09427000407395004_dp, &
0.13441695883240462_dp, &
0.19115015382463443_dp, &
0.2709446768704422_dp, &
0.3827855855236461_dp, &
0.539258738081746_dp, &
0.7583757287107051_dp, &
1.0670660316353924_dp, &
1.5090418848411256_dp, &
2.1662834672628906_dp, &
3.2357675876251375_dp, &
5.486470388225014_dp /)
ELSE IF(E_range < 440.0_dp) THEN
aw(:) = (/ &
0.0008278623594831022_dp, &
0.004440143416805031_dp, &
0.01128003587129483_dp, &
0.022028032667051767_dp, &
0.03789020505509771_dp, &
0.06086193467530397_dp, &
0.09406309027051278_dp, &
0.14214406770495588_dp, &
0.21182122672254763_dp, &
0.31264097380596306_dp, &
0.4580758820114268_dp, &
0.6670709107703453_dp, &
0.9662125691628374_dp, &
1.3928075446835457_dp, &
1.999398370177953_dp, &
2.8608300883005477_dp, &
4.086604504026533_dp, &
5.846446628460555_dp, &
8.437777512375702_dp, &
12.552247765285534_dp, &
0.00213328557069995_dp, &
0.005146085588574016_dp, &
0.008644133670212282_dp, &
0.013048431012303657_dp, &
0.01900286013451288_dp, &
0.02745396370025683_dp, &
0.039714746998698626_dp, &
0.05755410449194266_dp, &
0.08337327676104966_dp, &
0.12048626756528526_dp, &
0.17350552959972318_dp, &
0.24886544660526308_dp, &
0.3555627201776054_dp, &
0.5062720274646171_dp, &
0.7191928570131636_dp, &
1.0215433649108963_dp, &
1.4573743510597394_dp, &
2.1088839879773635_dp, &
3.1726250587097207_dp, &
5.413783636472752_dp /)
ELSE IF(E_range < 560.0_dp) THEN
aw(:) = (/ &
0.0006765896832495788_dp, &
0.003634394085386253_dp, &
0.009260530719612768_dp, &
0.018169752658108284_dp, &
0.03146515341551808_dp, &
0.050988272535331905_dp, &
0.07963350697485291_dp, &
0.1217357830919052_dp, &
0.1835999837821604_dp, &
0.27426969542730256_dp, &
0.40663421160125035_dp, &
0.5989941130666432_dp, &
0.8772653269419931_dp, &
1.2781153183683012_dp, &
1.853561810334941_dp, &
2.678139838889795_dp, &
3.8613430246727076_dp, &
5.5731890537663755_dp, &
8.111415950629196_dp, &
12.166643844965888_dp, &
0.001744085534574605_dp, &
0.004220260226030547_dp, &
0.007131610577612271_dp, &
0.010866459413033155_dp, &
0.016026022819158336_dp, &
0.02349633931892759_dp, &
0.034508040063482036_dp, &
0.05073473657802367_dp, &
0.07448275058135573_dp, &
0.1089792057134355_dp, &
0.15876384778629773_dp, &
0.2302241599837806_dp, &
0.332355181520261_dp, &
0.47790039062791473_dp, &
0.6852225479188921_dp, &
0.9818045868913531_dp, &
1.4120229509638687_dp, &
2.0583115818361715_dp, &
3.1169029192999083_dp, &
5.349708551352753_dp /)
ELSE IF(E_range < 700.0_dp) THEN
aw(:) = (/ &
0.0005610728470394199_dp, &
0.0030184735435080496_dp, &
0.007713925764506502_dp, &
0.015206673875434671_dp, &
0.026511588660652154_dp, &
0.04333598414806455_dp, &
0.06837526629597918_dp, &
0.10568575085620181_dp, &
0.1612095192600486_dp, &
0.24354540087900603_dp, &
0.3650622527190533_dp, &
0.543478534168398_dp, &
0.8040951855844733_dp, &
1.1829825071783888_dp, &
1.731655055993821_dp, &
2.524334438621591_dp, &
3.6704791974285818_dp, &
5.3403641584517585_dp, &
7.8320463776255975_dp, &
11.835332307041483_dp, &
0.0014468097107413669_dp, &
0.0035116724039719513_dp, &
0.0059696932867474494_dp, &
0.009180942727492066_dp, &
0.013708539558394253_dp, &
0.020383477135754603_dp, &
0.030362622969538497_dp, &
0.04523605791649382_dp, &
0.06722634472750327_dp, &
0.09948099701646662_dp, &
0.14646924511428006_dp, &
0.21452819094851275_dp, &
0.3126412187281207_dp, &
0.4536035006902304_dp, &
0.6559176905693149_dp, &
0.9473057460219979_dp, &
1.372450650425951_dp, &
2.0140284663151338_dp, &
3.0680354176477675_dp, &
5.293568767070731_dp /)
ELSE IF(E_range < 1000.0_dp) THEN
aw(:) = (/ &
0.00041560421429592847_dp, &
0.0022418213242648387_dp, &
0.005758919791027741_dp, &
0.01144723621422771_dp, &
0.020193846705899098_dp, &
0.033507461983925436_dp, &
0.05378604262498878_dp, &
0.08466990073555877_dp, &
0.13155968308081434_dp, &
0.2023862592183233_dp, &
0.3087273646000756_dp, &
0.4674028243887369_dp, &
0.7027470769716094_dp, &
1.049871620190003_dp, &
1.559459158936237_dp, &
2.305183154542849_dp, &
3.3963931351407965_dp, &
5.003743319066008_dp, &
7.425838541826663_dp, &
11.351430214987415_dp, &
0.0010723394308518957_dp, &
0.0026167108295615284_dp, &
0.0044949873202556065_dp, &
0.007025855035654186_dp, &
0.010714434790090814_dp, &
0.016306786237111653_dp, &
0.024848663667522333_dp, &
0.037806540695016204_dp, &
0.05727652861002731_dp, &
0.08628011670186121_dp, &
0.12916822918318413_dp, &
0.1921852425532451_dp, &
0.2842775740585458_dp, &
0.41830025350191025_dp, &
0.6129567439235276_dp, &
0.8963364108719717_dp, &
1.3136174965892848_dp, &
1.947903440388086_dp, &
2.994921282155757_dp, &
5.209663723814196_dp /)
ELSE IF(E_range < 1400.0_dp) THEN
aw(:) = (/ &
0.00031286955750815576_dp, &
0.0016923580735978944_dp, &
0.004371282335548313_dp, &
0.008765574086694813_dp, &
0.01565550439603156_dp, &
0.026379948198065384_dp, &
0.04308062162687878_dp, &
0.06904028024723377_dp, &
0.10919325396802464_dp, &
0.17088922417992009_dp, &
0.2650068737779714_dp, &
0.407557514776118_dp, &
0.6219883180977661_dp, &
0.9425096109188739_dp, &
1.4189974311212061_dp, &
2.124562918652659_dp, &
3.1683968574584798_dp, &
4.7214722907969815_dp, &
7.082939256028212_dp, &
10.940793859502318_dp, &
0.0008077688195695463_dp, &
0.00198217907214279_dp, &
0.003442623190312048_dp, &
0.005472727115904357_dp, &
0.008526326287248634_dp, &
0.013274083368636464_dp, &
0.020666062246111463_dp, &
0.03206313216437523_dp, &
0.049449526388029234_dp, &
0.075729284506798_dp, &
0.115136915108001_dp, &
0.17381798603932938_dp, &
0.26066593218318485_dp, &
0.38856838846527497_dp, &
0.5763925947653621_dp, &
0.8525553749973281_dp, &
1.2627024319582374_dp, &
1.8903771533011222_dp, &
2.9311611560083013_dp, &
5.136578694592024_dp /)
ELSE IF(E_range < 2400.0_dp) THEN
aw(:) = (/ &
0.00019827924757951374_dp, &
0.0010780149952817663_dp, &
0.002812851626363984_dp, &
0.005733085307279989_dp, &
0.010472592090438678_dp, &
0.018132055487802043_dp, &
0.0304927641048669_dp, &
0.05033548822527872_dp, &
0.08193518503798143_dp, &
0.13180675469394101_dp, &
0.20980486496289558_dp, &
0.330730820064991_dp, &
0.5166727229258016_dp, &
0.800423257149481_dp, &
1.2305450135295029_dp, &
1.879178588117511_dp, &
2.8551601446227886_dp, &
4.329889768979011_dp, &
6.603416745038167_dp, &
10.362925760797669_dp, &
0.0005125033637486195_dp, &
0.0012706365288583487_dp, &
0.002252024006575626_dp, &
0.0036913680406142216_dp, &
0.00596768443335258_dp, &
0.009642980896647496_dp, &
0.01553381875613038_dp, &
0.02485193072868798_dp, &
0.03941530545876881_dp, &
0.06194402047568934_dp, &
0.09648086561048903_dp, &
0.1489966939762812_dp, &
0.2282704567062106_dp, &
0.3471993819576423_dp, &
0.5248628649710425_dp, &
0.790158669486661_dp, &
1.1894686829242134_dp, &
1.807093779555356_dp, &
2.838567904998341_dp, &
5.030578836489572_dp /)
ELSE IF(E_range < 3600.0_dp) THEN
aw(:) = (/ &
0.00014059271306170835_dp, &
0.0007678052276054606_dp, &
0.0020214454597474916_dp, &
0.004179433657126495_dp, &
0.007783077817148623_dp, &
0.013779744802789188_dp, &
0.02371898428729375_dp, &
0.04005840059696467_dp, &
0.0666432004521081_dp, &
0.10943384278942726_dp, &
0.1775904530928529_dp, &
0.2850755740119256_dp, &
0.45301322888085055_dp, &
0.7131639998184849_dp, &
1.113102510684593_dp, &
1.7242004030109677_dp, &
2.6549580323641053_dp, &
4.077028706944089_dp, &
6.291138446566697_dp, &
9.984122686799962_dp, &
0.0003637583285696214_dp, &
0.0009100112025278443_dp, &
0.0016417579827546728_dp, &
0.0027620931122371824_dp, &
0.004599979363190032_dp, &
0.007646660869866119_dp, &
0.01263359249237722_dp, &
0.02067423309446183_dp, &
0.03347107915887389_dp, &
0.053610731252129035_dp, &
0.08499108363616721_dp, &
0.13344321531224232_dp, &
0.20764158045301837_dp, &
0.320461116278926_dp, &
0.4911027073155423_dp, &
0.7487883585484156_dp, &
1.1404346721206706_dp, &
1.7509389336184635_dp, &
2.775923343651192_dp, &
4.958949968100653_dp /)
ELSE IF(E_range < 5500.0_dp) THEN
aw(:) = (/ &
9.811663396420181e-05_dp, &
0.0005387153388274525_dp, &
0.0014337324614504023_dp, &
0.003015475649069391_dp, &
0.00574242092592876_dp, &
0.010423553997436157_dp, &
0.018399712080256196_dp, &
0.03183594079136726_dp, &
0.054183034755727444_dp, &
0.09088394337446928_dp, &
0.15043915259366922_dp, &
0.2460004768427999_dp, &
0.397742789601456_dp, &
0.6363895630963665_dp, &
1.0084973731304356_dp, &
1.58461462515818_dp, &
2.472840066417251_dp, &
3.8450318812987394_dp, &
6.002608315811181_dp, &
9.632227959921071_dp, &
0.0002541593700342826_dp, &
0.0006427296045497553_dp, &
0.0011844269388975941_dp, &
0.0020534673481220505_dp, &
0.0035323770096513207_dp, &
0.0060486166381445595_dp, &
0.010257020310097073_dp, &
0.017179062167323594_dp, &
0.02840492188803033_dp, &
0.0463876574741554_dp, &
0.074876061766774_dp, &
0.1195516589811091_dp, &
0.18896785479097591_dp, &
0.295954218597094_dp, &
0.4598070254340165_dp, &
0.7100522049104695_dp, &
1.0941411622149884_dp, &
1.6976051976548783_dp, &
2.716249516175453_dp, &
4.890778387724059_dp /)
ELSE IF(E_range < 8000.0_dp) THEN
aw(:) = (/ &
7.13833187688556e-05_dp, &
0.00039407124822848477_dp, &
0.0010604003024174474_dp, &
0.002268881540747128_dp, &
0.0044151693299503475_dp, &
0.008202937097573668_dp, &
0.01481435405742254_dp, &
0.02619019968929153_dp, &
0.04547446807839964_dp, &
0.07770105168267799_dp, &
0.13084070019640795_dp, &
0.21738378803883784_dp, &
0.35671803090973053_dp, &
0.578691162487099_dp, &
0.9289811811665666_dp, &
1.4774029766398036_dp, &
2.3316667282273356_dp, &
3.6637643925932917_dp, &
5.77570509069075_dp, &
9.354119201211928_dp, &
0.00018513022698605404_dp, &
0.0004733144979171716_dp, &
0.0008910319001903643_dp, &
0.0015901132720062632_dp, &
0.002817026526408489_dp, &
0.004950987551991881_dp, &
0.00858807831282116_dp, &
0.014676442712350758_dp, &
0.024714034905529713_dp, &
0.04104181118323569_dp, &
0.06728048077253167_dp, &
0.10897899916886848_dp, &
0.17457702888159804_dp, &
0.2768488574187613_dp, &
0.4351510694991959_dp, &
0.6792490750144262_dp, &
1.05704310176894_dp, &
1.654625316923605_dp, &
2.66802454490393_dp, &
4.835727135810798_dp /)
ELSE IF(E_range < 11000.0_dp) THEN
aw(:) = (/ &
5.447865824224468e-05_dp, &
0.0003023252833436619_dp, &
0.0008222028799151143_dp, &
0.0017879791707883972_dp, &
0.0035487169071266567_dp, &
0.006729981428581205_dp, &
0.012396154693500248_dp, &
0.02231994729563704_dp, &
0.039412371274235955_dp, &
0.06839247718640222_dp, &
0.11681794886521146_dp, &
0.19665672113357646_dp, &
0.3266664792319806_dp, &
0.5359838367017119_dp, &
0.8695609961282702_dp, &
1.3965918547432172_dp, &
2.224437922826948_dp, &
3.525172811949831_dp, &
5.601287870515528_dp, &
9.13946477027197_dp, &
0.0001414493975419132_dp, &
0.0003654550117807148_dp, &
0.000702039683835934_dp, &
0.0012860992922918333_dp, &
0.0023370421382291296_dp, &
0.004198508502697536_dp, &
0.0074223074265970575_dp, &
0.012899566080348283_dp, &
0.02205503225498209_dp, &
0.037139217883663346_dp, &
0.061667645535398816_dp, &
0.10107779439722338_dp, &
0.16370967430591404_dp, &
0.2622817076716638_dp, &
0.4161860961861409_dp, &
0.6553710609430115_dp, &
1.0280988828074795_dp, &
1.620933455708893_dp, &
2.6301295998009273_dp, &
4.792491706479763_dp /)
ELSE IF(E_range < 16000.0_dp) THEN
aw(:) = (/ &
3.9650391457476895e-05_dp, &
0.00022158028472701614_dp, &
0.0006112114647923708_dp, &
0.0013575095488041788_dp, &
0.002761831542692714_dp, &
0.005369892359523151_dp, &
0.010125281664288239_dp, &
0.01862645222338442_dp, &
0.03353965840202935_dp, &
0.05924893350243236_dp, &
0.1028670071194057_dp, &
0.175792226768633_dp, &
0.296086876640148_dp, &
0.49209213781649114_dp, &
0.8079355076501964_dp, &
1.3120896126448494_dp, &
2.1114907233261895_dp, &
3.3782746327884623_dp, &
5.415473393257776_dp, &
8.909899582793463_dp, &
0.00010310476798567228_dp, &
0.0002701445181166159_dp, &
0.000532877592275924_dp, &
0.0010085116940090583_dp, &
0.0018886202572384375_dp, &
0.003480600088494636_dp, &
0.0062898935173928425_dp, &
0.011146278372867174_dp, &
0.01939438715492993_dp, &
0.033184306984395034_dp, &
0.05591288068846412_dp, &
0.09288916872140662_dp, &
0.15233437891683738_dp, &
0.24689315815619667_dp, &
0.395983526921868_dp, &
0.6297456476831799_dp, &
0.9968438011799718_dp, &
1.5843863193809051_dp, &
2.588926649886331_dp, &
4.745504704109157_dp /)
ELSE IF(E_range < 25000.0_dp) THEN
aw(:) = (/ &
2.717615579704421e-05_dp, &
0.0001533571888014568_dp, &
0.0004314089621086483_dp, &
0.0009855129983075699_dp, &
0.0020690362703024243_dp, &
0.004147562024030517_dp, &
0.008042659775306937_dp, &
0.015174233986015624_dp, &
0.02795389960551952_dp, &
0.05041209452241241_dp, &
0.0891854696033202_dp, &
0.15505475209023395_dp, &
0.2653178725200425_dp, &
0.4474288857984473_dp, &
0.7445797834635748_dp, &
1.22440583751609_dp, &
1.993324990447638_dp, &
3.223505217987872_dp, &
5.218581644674421_dp, &
8.66559898782423_dp, &
7.081544983542917e-05_dp, &
0.00018918815391733118_dp, &
0.0003867212572274142_dp, &
0.000762422950412275_dp, &
0.0014799057292581785_dp, &
0.0028101541541442733_dp, &
0.005210328062027634_dp, &
0.009444597928361062_dp, &
0.016770457544315606_dp, &
0.029227069062833674_dp, &
0.05007790744827613_dp, &
0.08448452825165322_dp, &
0.14052692488766483_dp, &
0.23075377939722663_dp, &
0.3745946021444043_dp, &
0.6023877662255268_dp, &
0.9632415613085664_dp, &
1.5448917744876522_dp, &
2.5442808578887193_dp, &
4.6946167261750755_dp /)
ELSE IF(E_range < 40000.0_dp) THEN
aw(:) = (/ &
1.82768685914661e-05_dp, &
0.00010441735132806156_dp, &
0.00030100097884435366_dp, &
0.0007108479100968664_dp, &
0.0015457139233036375_dp, &
0.003201819623271307_dp, &
0.006393975098919852_dp, &
0.012383252897810644_dp, &
0.023351012360050015_dp, &
0.043003080061775764_dp, &
0.07753277683957492_dp, &
0.13713767335568874_dp, &
0.23838383061179366_dp, &
0.4078631497260737_dp, &
0.6878427663389639_dp, &
1.145110973779935_dp, &
1.8855370076270503_dp, &
3.0812814312518766_dp, &
5.036563676320201_dp, &
8.438735784909705_dp, &
4.775105170315456e-05_dp, &
0.00013072535751322757_dp, &
0.0002788397223892091_dp, &
0.0005749660315378945_dp, &
0.0011586417040428197_dp, &
0.0022690318028462603_dp, &
0.004319418539129222_dp, &
0.00801292542398754_dp, &
0.014524632365889833_dp, &
0.02578712381287132_dp, &
0.04493348686003409_dp, &
0.07697794475687925_dp, &
0.12985482658883507_dp, &
0.21600603761696566_dp, &
0.35485484356242103_dp, &
0.5769158865490484_dp, &
0.9317239155196531_dp, &
1.507644325935015_dp, &
2.5020531430500816_dp, &
4.6465062068562695_dp /)
ELSE IF(E_range < 60000.0_dp) THEN
aw(:) = (/ &
1.2996456455083999e-05_dp, &
7.520394153832135e-05_dp, &
0.00022219672220177028_dp, &
0.0005415907631617633_dp, &
0.0012155029482221118_dp, &
0.0025907135024754754_dp, &
0.00530495258969124_dp, &
0.010502772892217318_dp, &
0.020194089168230928_dp, &
0.037839725116622466_dp, &
0.06929410944259942_dp, &
0.1243035205065752_dp, &
0.21886063137419387_dp, &
0.378873081895905_dp, &
0.6458630906008894_dp, &
1.0859234436620155_dp, &
1.8044559829817377_dp, &
2.973587525330566_dp, &
4.897999127084281_dp, &
8.265341669279621_dp, &
3.404721078241363e-05_dp, &
9.55719982888862e-05_dp, &
0.00021238064319698303_dp, &
0.00045565622189792205_dp, &
0.0009479324859284201_dp, &
0.0019053019198898493_dp, &
0.0037081546030584183_dp, &
0.00701301384980832_dp, &
0.012931170420230086_dp, &
0.023311564708109894_dp, &
0.041183395525769606_dp, &
0.07144135728423114_dp, &
0.12189849536468766_dp, &
0.2049026796474422_dp, &
0.3398600188228418_dp, &
0.5574133279003208_dp, &
0.9074318562550266_dp, &
1.4787943531508865_dp, &
2.4692596985566686_dp, &
4.609157984004122_dp /)
ELSE IF(E_range < 90000.0_dp) THEN
aw(:) = (/ &
9.25565591635203e-06_dp, &
5.4379793871950615e-05_dp, &
0.00016529041249751665_dp, &
0.0004168866047271283_dp, &
0.0009665590315332995_dp, &
0.002119701890324356_dp, &
0.004448656432436941_dp, &
0.008997704835679786_dp, &
0.017627313815903538_dp, &
0.033582228458602224_dp, &
0.06241471879404543_dp, &
0.11346444904349712_dp, &
0.2022020456835865_dp, &
0.3539053121009028_dp, &
0.6094023131307444_dp, &
1.034127349774388_dp, &
1.7330266988502332_dp, &
2.8781735622671407_dp, &
4.774671355320365_dp, &
8.110487276691934_dp, &
2.4325594418551197e-05_dp, &
7.032339781107698e-05_dp, &
0.00016342502104344062_dp, &
0.00036496196734052573_dp, &
0.0007833717870455653_dp, &
0.0016150107069731901_dp, &
0.0032114034906267487_dp, &
0.006187606435750351_dp, &
0.011597481073915308_dp, &
0.0212137722840156_dp, &
0.03796983999609342_dp, &
0.06664845750728544_dp, &
0.11494672304794436_dp, &
0.19511876400396244_dp, &
0.32654532756194626_dp, &
0.5399778215970892_dp, &
0.8855899018059838_dp, &
1.4527432413091068_dp, &
2.43957936533579_dp, &
4.575364183005572_dp /)
ELSE IF(E_range < 130000.0_dp) THEN
aw(:) = (/ &
6.814190272661025e-06_dp, &
4.069869051280147e-05_dp, &
0.00012737186292629312_dp, &
0.0003320226711903597_dp, &
0.0007932380758797745_dp, &
0.0017847550827699414_dp, &
0.0038282243222301777_dp, &
0.00788917626715274_dp, &
0.01570933707259002_dp, &
0.030359967326833645_dp, &
0.05714837487503194_dp, &
0.10508150808776481_dp, &
0.1891988830212347_dp, &
0.3342531973315357_dp, &
0.5804875943957438_dp, &
0.9927739293709223_dp, &
1.6756599162156112_dp, &
2.8011572912649165_dp, &
4.674719598942103_dp, &
7.984607594319158_dp, &
1.797150327322172e-05_dp, &
5.359985968339942e-05_dp, &
0.00013011002767413486_dp, &
0.00030130796361535977_dp, &
0.0006649389030619321_dp, &
0.001401893613810303_dp, &
0.0028406234000407605_dp, &
0.005562651151717404_dp, &
0.010574910172820922_dp, &
0.019587230337479018_dp, &
0.035452962885294766_dp, &
0.06286025084771178_dp, &
0.10940644613911357_dp, &
0.18726231582689418_dp, &
0.31578042935311507_dp, &
0.5257957373608585_dp, &
0.8677330061830318_dp, &
1.431364173453237_dp, &
2.4151720865879605_dp, &
4.547581134621985_dp /)
ELSE IF(E_range < 200000.0_dp) THEN
aw(:) = (/ &
4.7716265094763185e-06_dp, &
2.9166926555296936e-05_dp, &
9.488230396934703e-05_dp, &
0.0002575987504500076_dp, &
0.0006375738523465402_dp, &
0.0014774253841498193_dp, &
0.003248281876300821_dp, &
0.0068361862386404765_dp, &
0.013861680047468216_dp, &
0.027217194330710702_dp, &
0.05195517971142119_dp, &
0.09673336158813088_dp, &
0.17613481552787563_dp, &
0.31435129085132746_dp, &
0.5509948694809227_dp, &
0.9503229523022058_dp, &
1.616438266636084_dp, &
2.7212693591950243_dp, &
4.570641044122416_dp, &
7.853156552955678_dp, &
1.2647105161868332e-05_dp, &
3.9371013426940826e-05_dp, &
0.00010088612878844911_dp, &
0.00024366946359126982_dp, &
0.0005550207649994575_dp, &
0.0012002190364340662_dp, &
0.002484036402975316_dp, &
0.004953212877428526_dp, &
0.00956552131403576_dp, &
0.017964228953314672_dp, &
0.03291715355843139_dp, &
0.05901007747740615_dp, &
0.10373076556694004_dp, &
0.1791557412726208_dp, &
0.3046003739565861_dp, &
0.5109816376949545_dp, &
0.8489895178903741_dp, &
1.4088416891120237_dp, &
2.38940805816802_dp, &
4.518258620381161_dp /)
ELSE IF(E_range < 400000.0_dp) THEN
aw(:) = (/ &
2.7085780074763008e-06_dp, &
1.7382495556237058e-05_dp, &
6.080080265426087e-05_dp, &
0.00017679399321086206_dp, &
0.00046288083403377915_dp, &
0.001122499243799274_dp, &
0.002561980867183537_dp, &
0.00556376582571754_dp, &
0.01158821206139618_dp, &
0.02328836569293624_dp, &
0.04537152658422795_dp, &
0.0860172289138871_dp, &
0.15917675836779016_dp, &
0.28825678371879787_dp, &
0.5119753641421063_dp, &
0.8937060728187617_dp, &
1.536894937050839_dp, &
2.6133227988586185_dp, &
4.42932853536345_dp, &
7.674043939074506_dp, &
7.256724979426217e-06_dp, &
2.4610926999432795e-05_dp, &
6.912840577395855e-05_dp, &
0.00017828058314008386_dp, &
0.00042625566892881024_dp, &
0.0009579365113877111_dp, &
0.0020465973824414043_dp, &
0.004192094540853455_dp, &
0.008285096729880728_dp, &
0.01587687830483599_dp, &
0.0296155319303513_dp, &
0.05394149162759069_dp, &
0.09618399441437563_dp, &
0.16827869509337126_dp, &
0.28947646517971315_dp, &
0.49079626217953604_dp, &
0.8232936419170709_dp, &
1.3778229154713002_dp, &
2.353834965384448_dp, &
4.477779801378415_dp /)
ELSE IF(E_range < 1000000.0_dp) THEN
aw(:) = (/ &
1.3068203488599816e-06_dp, &
9.212301266702123e-06_dp, &
3.6077785810084084e-05_dp, &
0.00011500743375480786_dp, &
0.0003229429198398811_dp, &
0.0008270002068881072_dp, &
0.001971982044918491_dp, &
0.004439944614704386_dp, &
0.00953329585301145_dp, &
0.019665309964283544_dp, &
0.0391925858707856_dp, &
0.07580219423503365_dp, &
0.14278622266628385_dp, &
0.262721249810341_dp, &
0.47336628743918624_dp, &
0.8371297881550314_dp, &
1.456719526154987_dp, &
2.5037192622718836_dp, &
4.285003587485378_dp, &
7.490325205909173_dp, &
3.5816664156775507e-06_dp, &
1.4101265568215938e-05_dp, &
4.478948057844214e-05_dp, &
0.00012516101942682066_dp, &
0.0003171748836272924_dp, &
0.0007458653611380562_dp, &
0.0016532487268825003_dp, &
0.0034918788589282487_dp, &
0.007083652194632769_dp, &
0.01388412442627512_dp, &
0.02641488898485549_dp, &
0.04896020939381398_dp, &
0.0886753147317904_dp, &
0.1573355767805416_dp, &
0.27410773841195346_dp, &
0.47010188920080226_dp, &
0.7967524373831676_dp, &
1.3456030406971342_dp, &
2.3167695492236375_dp, &
4.435612422382358_dp /)
ELSE
aw(:) = (/ &
2.8057975042657943e-07_dp, &
2.9504841012167483e-06_dp, &
1.529449108035861e-05_dp, &
5.835738812143725e-05_dp, &
0.00018549636803432827_dp, &
0.0005205715421922644_dp, &
0.001332741709633465_dp, &
0.0031773466764213163_dp, &
0.007152812404460419_dp, &
0.015356440829864506_dp, &
0.03167411288231694_dp, &
0.06312050342894393_dp, &
0.1220727351170772_dp, &
0.2299360282547445_dp, &
0.4230917616548655_dp, &
0.7625318197450621_dp, &
1.3498427130340422_dp, &
2.356254436812211_dp, &
4.089380405073521_dp, &
7.239955448462523_dp, &
8.800290760673824e-07_dp, &
5.546603080912524e-06_dp, &
2.2360869953602353e-05_dp, &
7.205981239759e-05_dp, &
0.00020168204190107724_dp, &
0.0005110941710247_dp, &
0.0012017109415185414_dp, &
0.002663283670778342_dp, &
0.00562453491805382_dp, &
0.011408822070675338_dp, &
0.022359646995251253_dp, &
0.04253683853853289_dp, &
0.07883916269382399_dp, &
0.14279601051705815_dp, &
0.25342702929276556_dp, &
0.4419404134586229_dp, &
0.7602899768847363_dp, &
1.3010198298394264_dp, &
2.2652747116311405_dp, &
4.3770397670783465_dp /)
END IF
END SUBROUTINE

SUBROUTINE get_coeff_26(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio
E_ratio = 1.0_dp
IF(E_range < 799.0_dp) THEN
E_ratio = 799.0_dp/E_range
aw(:) = (/ &
0.00037867004555797054_dp, &
0.002019373143542281_dp, &
0.005074452374668418_dp, &
0.009741339051894624_dp, &
0.016351249099763064_dp, &
0.025426574118427834_dp, &
0.03775841037254848_dp, &
0.05450081932003549_dp, &
0.0772810957375601_dp, &
0.1083363446608876_dp, &
0.15069589349866513_dp, &
0.20843076717336742_dp, &
0.2869911650766578_dp, &
0.3936564974164155_dp, &
0.5381309695975518_dp, &
0.7333307052011886_dp, &
0.996427988952117_dp, &
1.3502505778649674_dp, &
1.8251937089224979_dp, &
2.4619239536701665_dp, &
3.3154239654767634_dp, &
4.461579856601731_dp, &
6.009262672635442_dp, &
8.12630074937428_dp, &
11.10929348730981_dp, &
15.660501844751538_dp, &
0.0009745069197489692_dp, &
0.002323920525926464_dp, &
0.003818499918990434_dp, &
0.005569525222864092_dp, &
0.007736056032946723_dp, &
0.01054486610207425_dp, &
0.014309337959120037_dp, &
0.01944447015531593_dp, &
0.026485944451580964_dp, &
0.03612533116765179_dp, &
0.04926670257709911_dp, &
0.06710398765236568_dp, &
0.09122022393749045_dp, &
0.12371484118454482_dp, &
0.16736960998240946_dp, &
0.2258689406847251_dp, &
0.30409900038431836_dp, &
0.40856816694415354_dp, &
0.548030961075663_dp, &
0.7344894715293502_dp, &
0.9849728438705977_dp, &
1.3251021306723128_dp, &
1.7972852200464182_dp, &
2.4831111234607426_dp, &
3.582848552118199_dp, &
5.8876320817174905_dp /)
ELSE IF(E_range < 995.0_dp) THEN
aw(:) = (/ &
0.00031416617596844365_dp, &
0.0016768070133771967_dp, &
0.004220387402167522_dp, &
0.008122019548232636_dp, &
0.013681535061623657_dp, &
0.021375962641606915_dp, &
0.03193242611639194_dp, &
0.04641604327774391_dp, &
0.0663338368572065_dp, &
0.09376694799993701_dp, &
0.13155072067060036_dp, &
0.18352265415559924_dp, &
0.2548586347835214_dp, &
0.3525227140309385_dp, &
0.4858648853075939_dp, &
0.667414665262307_dp, &
0.9139380450657599_dp, &
1.2478576632771592_dp, &
1.699195061182984_dp, &
2.3083135493129188_dp, &
3.13000622971369_dp, &
4.240119346476375_dp, &
5.747638276053859_dp, &
7.820581475731412_dp, &
10.755506760817696_dp, &
15.2531099085639_dp, &
0.0008086633399478686_dp, &
0.001931716707111731_dp, &
0.0031843549489835176_dp, &
0.0046679761695975755_dp, &
0.006529332077378777_dp, &
0.008979243990088065_dp, &
0.012309502604654843_dp, &
0.016906748893981475_dp, &
0.023272740862285972_dp, &
0.03206195264016037_dp, &
0.04413963056111703_dp, &
0.060659468387992695_dp, &
0.08316331669118938_dp, &
0.11371010407957528_dp, &
0.15504510373245767_dp, &
0.21082547573489732_dp, &
0.28592663287586667_dp, &
0.38687165073350616_dp, &
0.5224647206101295_dp, &
0.7047999539903935_dp, &
0.951040871895601_dp, &
1.286968657510818_dp, &
1.7551267580277088_dp, &
2.437087428172129_dp, &
3.5325565753903243_dp, &
5.829332987897013_dp /)
ELSE IF(E_range < 1293.0_dp) THEN
aw(:) = (/ &
0.000251216262630506_dp, &
0.0013422657589387847_dp, &
0.0033852831378795366_dp, &
0.006535679555855439_dp, &
0.01105946790435755_dp, &
0.017384205344456737_dp, &
0.02616631610050564_dp, &
0.03837171615696803_dp, &
0.05537294469791371_dp, &
0.07907682374899021_dp, &
0.11210158563747868_dp, &
0.15802208108125124_dp, &
0.22170328658555793_dp, &
0.3097484547351596_dp, &
0.43109811316492785_dp, &
0.5978298541062674_dp, &
0.826228906983008_dp, &
1.1382318360717252_dp, &
1.5634041077266412_dp, &
2.1417299191833186_dp, &
2.9277541939220253_dp, &
3.9972537939524453_dp, &
5.459336924183542_dp, &
7.482255655971381_dp, &
10.362568049094785_dp, &
14.79927594838438_dp, &
0.0006467895131307655_dp, &
0.0015483789662149572_dp, &
0.002562990408244027_dp, &
0.00378128864955077_dp, &
0.005336413195947385_dp, &
0.0074211270754708444_dp, &
0.010302498712647285_dp, &
0.014334977955021215_dp, &
0.01998237358756267_dp, &
0.027857715884000284_dp, &
0.038782119259036095_dp, &
0.05386231670249176_dp, &
0.07459076970770032_dp, &
0.10297648694866778_dp, &
0.1417181535853474_dp, &
0.19443583384059848_dp, &
0.26598602063677856_dp, &
0.36290208001606944_dp, &
0.4940399248668749_dp, &
0.671597172994337_dp, &
0.9128949281590978_dp, &
1.2439090969479831_dp, &
1.7073582868867894_dp, &
2.3848234540233704_dp, &
3.4754029453297375_dp, &
5.763153512206444_dp /)
ELSE IF(E_range < 1738.0_dp) THEN
aw(:) = (/ &
0.0001950376405774879_dp, &
0.0010434591619187578_dp, &
0.002638226324824606_dp, &
0.005113292156951265_dp, &
0.00870086598207974_dp, &
0.013778342596377403_dp, &
0.020929407226154117_dp, &
0.031017048599119232_dp, &
0.04527407877688829_dp, &
0.06542623021788685_dp, &
0.09386550782858209_dp, &
0.1338911807639374_dp, &
0.19003887760085253_dp, &
0.26852546690556994_dp, &
0.3778478282417773_dp, &
0.5295878007623821_dp, &
0.739496104873784_dp, &
1.0289606432203449_dp, &
1.427022553455354_dp, &
1.9732190338166054_dp, &
2.721788596346289_dp, &
3.7484019634398744_dp, &
5.162286381612833_dp, &
7.131961441982249_dp, &
9.954041146838616_dp, &
14.32582187195048_dp, &
0.0005022997712497699_dp, &
0.0012056339745159519_dp, &
0.0020057022650528562_dp, &
0.0029823455742529125_dp, &
0.004254672143777001_dp, &
0.0059963576674781395_dp, &
0.008448134528374854_dp, &
0.011930608474038368_dp, &
0.016868224194481772_dp, &
0.023830776674747995_dp, &
0.03359223660591877_dp, &
0.04720772204797841_dp, &
0.06611395904656978_dp, &
0.09226211932358148_dp, &
0.12829508657243463_dp, &
0.17778587596291995_dp, &
0.24556239855995798_dp, &
0.3381606409098896_dp, &
0.4644850568708377_dp, &
0.6368412465900805_dp, &
0.8727231844202037_dp, &
1.198329166055601_dp, &
1.6565902038491545_dp, &
2.329130996306286_dp, &
3.414441737018883_dp, &
5.692646250049292_dp /)
ELSE IF(E_range < 2238.0_dp) THEN
aw(:) = (/ &
0.0001570075292014777_dp, &
0.0008410019818081912_dp, &
0.002131222620982682_dp, &
0.004145569234358363_dp, &
0.007090677510918224_dp, &
0.011305441203226727_dp, &
0.01731697219552165_dp, &
0.02590771238094241_dp, &
0.03820108916428952_dp, &
0.05578097500511556_dp, &
0.08086136666211902_dp, &
0.11652297811189155_dp, &
0.16703767866605654_dp, &
0.23830959694579532_dp, &
0.33847251734047956_dp, &
0.47869779524697986_dp, &
0.6742880172162485_dp, &
0.9461646148413084_dp, &
1.3229154365464835_dp, &
1.8436824667954244_dp, &
2.5624225632785795_dp, &
3.5546935048363073_dp, &
4.929807210724086_dp, &
6.856512184015644_dp, &
9.631510529999963_dp, &
13.950799086937636_dp, &
0.0004044673162152666_dp, &
0.0009731515874605558_dp, &
0.0016264516115000137_dp, &
0.002435950210582391_dp, &
0.0035098078378292706_dp, &
0.005006471432500407_dp, &
0.0071456103103081205_dp, &
0.010221155980366915_dp, &
0.01462678603467543_dp, &
0.02089809034158199_dp, &
0.02977080559865047_dp, &
0.042257136734592034_dp, &
0.059746455952240894_dp, &
0.08413972088086426_dp, &
0.11803003373030124_dp, &
0.1649465354635964_dp, &
0.22968730801733273_dp, &
0.31878354856943114_dp, &
0.44117315496449167_dp, &
0.6092463295611383_dp, &
0.8406397450171417_dp, &
1.161742484167732_dp, &
1.6156774111835495_dp, &
2.2841323787971173_dp, &
3.36513821401701_dp, &
5.635684641787822_dp /)
ELSE IF(E_range < 3009.0_dp) THEN
aw(:) = (/ &
0.00012173433282636892_dp, &
0.000653044843078349_dp, &
0.0016597094198001826_dp, &
0.0032432122902587715_dp, &
0.005583731099072486_dp, &
0.00897974690437299_dp, &
0.01389829106398585_dp, &
0.021035914001762628_dp, &
0.03139940878203811_dp, &
0.046421111007702805_dp, &
0.06812376578216424_dp, &
0.09935124887637425_dp, &
0.1440867725430818_dp, &
0.20788864079466032_dp, &
0.29848482929563186_dp, &
0.42658280937440796_dp, &
0.6069727081348725_dp, &
0.8600354677167223_dp, &
1.213825493528422_dp, &
1.7070100456167858_dp, &
2.393195963336712_dp, &
3.3477843441213753_dp, &
4.680165709088339_dp, &
6.5593534066075065_dp, &
9.282192058087176_dp, &
13.54332387201712_dp, &
0.00031370730222591084_dp, &
0.0007570695990484558_dp, &
0.0012727251500342805_dp, &
0.0019236429401252328_dp, &
0.0028062999069873792_dp, &
0.004062578033501044_dp, &
0.00588929507373728_dp, &
0.00855186255429553_dp, &
0.012411138722817702_dp, &
0.017965594962442332_dp, &
0.025908461038187995_dp, &
0.03720320806718786_dp, &
0.05318444551258281_dp, &
0.0756939715032581_dp, &
0.10726485539108421_dp, &
0.1513713746489183_dp, &
0.21277116226058032_dp, &
0.297982246389678_dp, &
0.41597262011480124_dp, &
0.5792223633874778_dp, &
0.8055284041136077_dp, &
1.1215023110944804_dp, &
1.5705011278690364_dp, &
2.234313406989394_dp, &
3.310497107984999_dp, &
5.57261913216792_dp /)
ELSE IF(E_range < 4377.0_dp) THEN
aw(:) = (/ &
8.814790587980884e-05_dp, &
0.0004738555483452858_dp, &
0.001209162224838196_dp, &
0.0023779784718090265_dp, &
0.004131730232456246_dp, &
0.006724209366814066_dp, &
0.010555208090639503_dp, &
0.01622487926307075_dp, &
0.02460909648355756_dp, &
0.036969426022926015_dp, &
0.05511124918766332_dp, &
0.08160627715945637_dp, &
0.12010201314867164_dp, &
0.17574961207193285_dp, &
0.2557933011845474_dp, &
0.37038042948632577_dp, &
0.5336738618256365_dp, &
0.7653829566132899_dp, &
1.092887558492255_dp, &
1.5542408605150586_dp, &
2.2025803308040364_dp, &
3.113074354623969_dp, &
4.395181735531908_dp, &
6.218242520601409_dp, &
8.87933263258288_dp, &
13.071612343557268_dp, &
0.0002272629071152242_dp, &
0.0005507552663901378_dp, &
0.000933441836691016_dp, &
0.001428839311231769_dp, &
0.002120237386026644_dp, &
0.0031304749854564546_dp, &
0.0046303259129662575_dp, &
0.006853192947474817_dp, &
0.01012290436688526_dp, &
0.014894995829624644_dp, &
0.021812185159477044_dp, &
0.03177875816649789_dp, &
0.046061554071504564_dp, &
0.06642770441212868_dp, &
0.09533262717257988_dp, &
0.1361769876776064_dp, &
0.19366000664329205_dp, &
0.27427258699605633_dp, &
0.38700765699040135_dp, &
0.5444449877078907_dp, &
0.7645728576669742_dp, &
1.0742808530249828_dp, &
1.517233705376216_dp, &
2.175383868205045_dp, &
3.245782737215517_dp, &
5.498017132429984_dp /)
ELSE IF(E_range < 6256.0_dp) THEN
aw(:) = (/ &
6.476469600717215e-05_dp, &
0.00034891783555501116_dp, &
0.0008941569546994692_dp, &
0.001770490093242039_dp, &
0.003106199645499256_dp, &
0.0051184335340408015_dp, &
0.008151285145167414_dp, &
0.012724953439821175_dp, &
0.019606591631387022_dp, &
0.029915037575018_dp, &
0.04527205344500337_dp, &
0.06801655695473027_dp, &
0.10150518730808944_dp, &
0.1505317627030831_dp, &
0.22191039986921127_dp, &
0.32528375402693044_dp, &
0.4742415082162232_dp, &
0.6878698417118773_dp, &
0.9929114555738491_dp, &
1.426826349178114_dp, &
2.0422826887022167_dp, &
2.914194983288514_dp, &
4.152054985305785_dp, &
5.925500508772846_dp, &
8.531868435284434_dp, &
12.663118219514553_dp, &
0.00016705896835219716_dp, &
0.0004066434009563924_dp, &
0.0006951433499719869_dp, &
0.0010783855735630424_dp, &
0.0016285976932411544_dp, &
0.002452428674095891_dp, &
0.003698822692463904_dp, &
0.005574655476173014_dp, &
0.008372611367675033_dp, &
0.012511058447330109_dp, &
0.01858785984880576_dp, &
0.027453700962179424_dp, &
0.04031295449132496_dp, &
0.058862610090479196_dp, &
0.08548341120436852_dp, &
0.12350282964756092_dp, &
0.17755835622889385_dp, &
0.25410556179535176_dp, &
0.3621484874152968_dp, &
0.5143475290193675_dp, &
0.728860394738634_dp, &
1.0328353641593397_dp, &
1.470237966477004_dp, &
2.1232087233584265_dp, &
3.1884016755763915_dp, &
5.431943610280663_dp /)
ELSE IF(E_range < 9034.0_dp) THEN
aw(:) = (/ &
4.714267869850108e-05_dp, &
0.00025461261481751654_dp, &
0.0006556761719344113_dp, &
0.001308468824945838_dp, &
0.002321150197926946_dp, &
0.003878461874698709_dp, &
0.006274894154406174_dp, &
0.009959418801224819_dp, &
0.01560204953263909_dp, &
0.02419303582580116_dp, &
0.03718681238787033_dp, &
0.056707527706112246_dp, &
0.08584014778037409_dp, &
0.12904063424870618_dp, &
0.19271144671607915_dp, &
0.2860061442236756_dp, &
0.4219516424878554_dp, &
0.6190135338939492_dp, &
0.903289679923842_dp, &
1.3116275192031754_dp, &
1.8961959695325519_dp, &
2.7316203672239747_dp, &
3.927395419802043_dp, &
5.653448502577673_dp, &
8.207416056488858_dp, &
12.280210966173971_dp, &
0.00012167157119425256_dp, &
0.00029765224752192887_dp, &
0.0005138446741296063_dp, &
0.0008093181914490766_dp, &
0.0012462946084647903_dp, &
0.0019166924309280847_dp, &
0.002949916160920999_dp, &
0.004529176144148134_dp, &
0.006918766955709136_dp, &
0.010502242217013015_dp, &
0.01583456242948638_dp, &
0.02371430118668526_dp, &
0.03528415663667773_dp, &
0.05217068929971722_dp, &
0.07667811497799712_dp, &
0.11205676769912475_dp, &
0.16287592945291882_dp, &
0.23554669511627027_dp, &
0.3390729578857381_dp, &
0.4861840636834821_dp, &
0.6951983617349604_dp, &
0.9935218371108636_dp, &
1.4254333766864864_dp, &
2.073293734219426_dp, &
3.133425653975384_dp, &
5.368705541720747_dp /)
ELSE IF(E_range < 15564.0_dp) THEN
aw(:) = (/ &
2.9443353786715123e-05_dp, &
0.00015968865887835848_dp, &
0.0004146519996900374_dp, &
0.0008385653229066657_dp, &
0.0015154710444555315_dp, &
0.002590553782265674_dp, &
0.004297458120933096_dp, &
0.0069977275055487075_dp, &
0.011241365987483699_dp, &
0.017857778834567364_dp, &
0.02808897944065483_dp, &
0.043782322061949784_dp, &
0.06766742773244867_dp, &
0.1037518951041952_dp, &
0.1578839462681274_dp, &
0.2385489270378087_dp, &
0.3579930679384087_dp, &
0.5338068043298525_dp, &
0.791161589813274_dp, &
1.1660045885058554_dp, &
1.7097467191304874_dp, &
2.496543051760848_dp, &
3.635839789357852_dp, &
5.297957374408996_dp, &
7.7810148585358885_dp, &
11.774676503198675_dp, &
7.606243073747058e-05_dp, &
0.000187654111092642_dp, &
0.00032937906982910266_dp, &
0.0005320965908695892_dp, &
0.0008454741794288284_dp, &
0.0013429974722876687_dp, &
0.0021300795419892347_dp, &
0.003360638155325681_dp, &
0.005262693302642904_dp, &
0.008173980176265023_dp, &
0.0125917789247931_dp, &
0.01924336061080778_dp, &
0.029185562811316818_dp, &
0.043945018416338716_dp, &
0.0657148857554778_dp, &
0.09763017054646428_dp, &
0.1441532549628787_dp, &
0.2116173831358168_dp, &
0.3090075514896036_dp, &
0.4491309451933545_dp, &
0.6505180076674358_dp, &
0.9409371471889534_dp, &
1.3651313188527026_dp, &
2.0058258806612064_dp, &
3.0589785756980588_dp, &
5.283171742769873_dp /)
ELSE IF(E_range < 19500.0_dp) THEN
aw(:) = (/ &
2.4221939880648916e-05_dp, &
0.00013162518783252502_dp, &
0.0003431050538880391_dp, &
0.000698189132475488_dp, &
0.001272589233863707_dp, &
0.0021976433259721703_dp, &
0.0036856546266738404_dp, &
0.00606744928063478_dp, &
0.009850534125738348_dp, &
0.015806672884776662_dp, &
0.025100760660057786_dp, &
0.03947837634367904_dp, &
0.06153686837154947_dp, &
0.09511493193970912_dp, &
0.14584951983758487_dp, &
0.22196820786998617_dp, &
0.335412339591795_dp, &
0.5034260546260962_dp, &
0.7508088209241689_dp, &
1.1131398816244866_dp, &
1.6415127034352714_dp, &
2.409875693114849_dp, &
3.527638075049888_dp, &
5.165269756856301_dp, &
7.621098250903887_dp, &
11.584360794128232_dp, &
6.260081073117763e-05_dp, &
0.0001550481237336088_dp, &
0.00027425459381644776_dp, &
0.0004482075757040794_dp, &
0.0007220779147665958_dp, &
0.0011628086873198505_dp, &
0.0018674017508154028_dp, &
0.0029793398493427954_dp, &
0.004713349903984844_dp, &
0.00739000908391551_dp, &
0.011484659947886356_dp, &
0.017697084080236986_dp, &
0.0270505598375803_dp, &
0.04103207563980036_dp, &
0.06178995726179555_dp, &
0.09241161287658448_dp, &
0.13731382654384652_dp, &
0.20279428302236951_dp, &
0.2978245159544139_dp, &
0.43523615879690214_dp, &
0.6336386632046422_dp, &
0.9209431632237726_dp, &
1.3420830429471422_dp, &
1.9799449671624767_dp, &
3.0303736903126226_dp, &
5.250335959421828_dp /)
ELSE IF(E_range < 22300.0_dp) THEN
aw(:) = (/ &
2.1564828392874976e-05_dp, &
0.00011732900716735743_dp, &
0.00030658459101179467_dp, &
0.0006263105810174703_dp, &
0.0011476627771254888_dp, &
0.001994368166545555_dp, &
0.0033669949117497793_dp, &
0.005579455105811666_dp, &
0.009115754087433222_dp, &
0.014715595238257124_dp, &
0.023500724387693545_dp, &
0.037159447358704495_dp, &
0.058214268280899095_dp, &
0.09040780191303678_dp, &
0.1392561892112249_dp, &
0.21283882379303048_dp, &
0.32292082611394995_dp, &
0.48654501997143224_dp, &
0.728293233827829_dp, &
1.0835278065593024_dp, &
1.6031531143244464_dp, &
2.360992062122655_dp, &
3.466427749396279_dp, &
5.090015465029667_dp, &
7.53020860247837_dp, &
11.47601339726731_dp, &
5.574871333711635e-05_dp, &
0.00013841643320123116_dp, &
0.0002460246878681958_dp, &
0.00040498060578160683_dp, &
0.0006579585426571892_dp, &
0.0010682911223034433_dp, &
0.0017283542359863323_dp, &
0.002775834478515568_dp, &
0.0044179824445092115_dp, &
0.006965640555090498_dp, &
0.010881623599671406_dp, &
0.01684992291318925_dp, &
0.025874425224763156_dp, &
0.03941906165439105_dp, &
0.05960588810571241_dp, &
0.08949416080199679_dp, &
0.13347332726212882_dp, &
0.1978191873760846_dp, &
0.29149390319797547_dp, &
0.4273416663588288_dp, &
0.6240165631691874_dp, &
0.9095125059692878_dp, &
1.3288753869034111_dp, &
1.9650901753832044_dp, &
3.013944114046813_dp, &
5.231485480677545_dp /)
ELSE IF(E_range < 24783.0_dp) THEN
aw(:) = (/ &
1.968079412594506e-05_dp, &
0.00010718493641102105_dp, &
0.0002806353581711279_dp, &
0.0005751273705017238_dp, &
0.0010584284934163331_dp, &
0.001848589008502448_dp, &
0.003137423022955394_dp, &
0.005226209962089254_dp, &
0.008581348030605833_dp, &
0.013918427654363913_dp, &
0.022326614770472653_dp, &
0.03545081033521797_dp, &
0.0557565940780925_dp, &
0.0869132477657862_dp, &
0.1343443806186961_dp, &
0.20601549407747916_dp, &
0.31355578173538573_dp, &
0.4738522008187966_dp, &
0.7113173707719063_dp, &
1.061144148921342_dp, &
1.574088149869358_dp, &
2.323872276609698_dp, &
3.4198568581275923_dp, &
5.032662209058464_dp, &
7.460841068394957_dp, &
11.393228561803523_dp, &
5.0889409017427654e-05_dp, &
0.00012660472263169967_dp, &
0.0002259210980998603_dp, &
0.0003740655095418279_dp, &
0.0006118384940478349_dp, &
0.000999874816429794_dp, &
0.0016270976310772938_dp, &
0.0026268356436642305_dp, &
0.004200673380690335_dp, &
0.006652036004218306_dp, &
0.010434155159748769_dp, &
0.016218897955811638_dp, &
0.02499520138410832_dp, &
0.038209141684542064_dp, &
0.05796233141783732_dp, &
0.08729199004516089_dp, &
0.13056598592688773_dp, &
0.19404256395643912_dp, &
0.2866758412418721_dp, &
0.4213188646546868_dp, &
0.6166595702412992_dp, &
0.900755794469128_dp, &
1.3187413549017004_dp, &
1.9536794891939209_dp, &
3.00131680399958_dp, &
5.21699996325917_dp /)
ELSE IF(E_range < 41198.0_dp) THEN
aw(:) = (/ &
1.2674056285680165e-05_dp, &
6.938931077849002e-05_dp, &
0.0001836100520485683_dp, &
0.0003826779493233582_dp, &
0.0007202040757498056_dp, &
0.0012903838734080182_dp, &
0.002248236525495258_dp, &
0.00384173613615444_dp, &
0.006462366726733287_dp, &
0.010722120046078792_dp, &
0.017568831387303644_dp, &
0.028457372697558436_dp, &
0.04560190325370437_dp, &
0.07234504966297221_dp, &
0.11369479964375559_dp, &
0.1771007906182374_dp, &
0.2735712154870606_dp, &
0.4192744014302692_dp, &
0.637834936286395_dp, &
0.9636469742972286_dp, &
1.4467538592483378_dp, &
2.1603857309953622_dp, &
3.2137696680006567_dp, &
4.777814618678099_dp, &
7.151554564591586_dp, &
11.0231252332385_dp, &
3.2810000033385146e-05_dp, &
8.24959243779711e-05_dp, &
0.00015031831931530549_dp, &
0.0002565196061853303_dp, &
0.0004339110517853371_dp, &
0.0007317640998323022_dp, &
0.0012244457171152918_dp, &
0.0020265452915832964_dp, &
0.0033148331774103383_dp, &
0.005359851622721186_dp, &
0.008571929299584299_dp, &
0.013568145467107103_dp, &
0.021269239209990003_dp, &
0.033039016114448296_dp, &
0.05088380967787895_dp, &
0.07773666161595308_dp, &
0.11786137108586371_dp, &
0.177428416908905_dp, &
0.26534598308886587_dp, &
0.3944982542877601_dp, &
0.583721306443492_dp, &
0.8613658745137791_dp, &
1.2729807136647084_dp, &
1.9020150846924464_dp, &
2.9440730604569967_dp, &
5.151372469471086_dp /)
ELSE IF(E_range < 94407.0_dp) THEN
aw(:) = (/ &
6.186479912580181e-06_dp, &
3.4233291323971405e-05_dp, &
9.256119767321995e-05_dp, &
0.00019951435503347473_dp, &
0.000391782808208518_dp, &
0.0007348447007241284_dp, &
0.0013393430685668998_dp, &
0.0023881407715703996_dp, &
0.004179359736833124_dp, &
0.007193105406189777_dp, &
0.012193603475754204_dp, &
0.020384081500457202_dp, &
0.03363963393153156_dp, &
0.054854535910931194_dp, &
0.0884563818462964_dp, &
0.1411619416118249_dp, &
0.22308161106495417_dp, &
0.3493255541517698_dp, &
0.5423345187791336_dp, &
0.835274429372801_dp, &
1.277059229828611_dp, &
1.9400953490632433_dp, &
2.933324383841704_dp, &
4.428041894721328_dp, &
6.724054874045832_dp, &
10.508722000292966_dp, &
1.605270034684411e-05_dp, &
4.1236777606757607e-05_dp, &
7.834646084720584e-05_dp, &
0.00014151002746774427_dp, &
0.0002536786392750537_dp, &
0.00045044627173442046_dp, &
0.0007883280940822237_dp, &
0.0013578695679215443_dp, &
0.0023029240253660617_dp, &
0.0038494157006882636_dp, &
0.006348322282988488_dp, &
0.010339483088961496_dp, &
0.016645606665935727_dp, &
0.026509781185170438_dp, &
0.04179532069895724_dp, &
0.06527461060266015_dp, &
0.10104501326885583_dp, &
0.15512766861607885_dp, &
0.23633651295613511_dp, &
0.3575717157671407_dp, &
0.5378622500972828_dp, &
0.8059829678217727_dp, &
1.2081212436199635_dp, &
1.8283703742955213_dp, &
2.862256950464903_dp, &
5.057682060420358_dp /)
ELSE IF(E_range < 189080.0_dp) THEN
aw(:) = (/ &
3.398642214192787e-06_dp, &
1.9025312202899538e-05_dp, &
5.266161038810626e-05_dp, &
0.00011755445811565422_dp, &
0.00024054583577387794_dp, &
0.00047036774204517446_dp, &
0.0008915989199018158_dp, &
0.0016480620414409948_dp, &
0.00298048536052862_dp, &
0.005285920675047294_dp, &
0.009210353187127669_dp, &
0.01579153884075367_dp, &
0.026677173353753962_dp, &
0.044455103994883116_dp, &
0.07314891087057689_dp, &
0.11895589362798491_dp, &
0.19133835932619525_dp, &
0.3046280564326354_dp, &
0.48037697428284243_dp, &
0.7508072430736709_dp, &
1.1639388014245502_dp, &
1.7914938071058262_dp, &
2.7421307814095353_dp, &
4.18739539739007_dp, &
6.427708444690209_dp, &
10.150038983021496_dp, &
8.840904815568487e-06_dp, &
2.3243411245798182e-05_dp, &
4.61419948505059e-05_dp, &
8.798822647972772e-05_dp, &
0.0001658923006259312_dp, &
0.0003074688101295367_dp, &
0.0005583537180990002_dp, &
0.0009937301893068423_dp, &
0.001735787847791623_dp, &
0.002980424755862107_dp, &
0.005037848905051201_dp, &
0.008393759041169294_dp, &
0.013800737687712429_dp, &
0.022413678806554142_dp, &
0.03598900168633015_dp, &
0.057175813550115986_dp, &
0.08993936716796855_dp, &
0.1401758115607063_dp, &
0.21660925146824841_dp, &
0.33212798532910376_dp, &
0.5058820683958832_dp, &
0.7669510501463104_dp, &
1.1620122321700583_dp, &
1.775691260253936_dp, &
2.8035588004499257_dp, &
4.990540570288636_dp /)
ELSE IF(E_range < 457444.0_dp) THEN
aw(:) = (/ &
1.5924512725077461e-06_dp, &
9.091861526358144e-06_dp, &
2.617525769293245e-05_dp, &
6.171102340282755e-05_dp, &
0.00013395463159750724_dp, &
0.00027703236128610273_dp, &
0.0005524440573740519_dp, &
0.0010685815369178362_dp, &
0.002012774020332085_dp, &
0.0037030176407984823_dp, &
0.006670370391832498_dp, &
0.011788596543354984_dp, &
0.020475865129561125_dp, &
0.03500529115221168_dp, &
0.058978461638664194_dp, &
0.09804109351785203_dp, &
0.16095595007242294_dp, &
0.26120029165707_dp, &
0.41933283009911826_dp, &
0.6664993354060413_dp, &
1.0496740893941348_dp, &
1.6397469396988835_dp, &
2.54498795937397_dp, &
3.9371777982755236_dp, &
6.11745128618483_dp, &
9.772519404832513_dp, &
4.1599386041678434e-06_dp, &
1.1373735340879785e-05_dp, &
2.420740702971419e-05_dp, &
4.980537853880897e-05_dp, &
0.00010018420860309139_dp, &
0.00019590630897285626_dp, &
0.0003724648960236611_dp, &
0.0006901948951647627_dp, &
0.0012498546464018527_dp, &
0.0022170391172995803_dp, &
0.0038600176308915607_dp, &
0.006607746433407134_dp, &
0.011137935349822067_dp, &
0.018509521510997304_dp, &
0.030360279138384877_dp, &
0.04919938780336098_dp, &
0.07883793047658633_dp, &
0.12502011000770785_dp, &
0.19635093903737422_dp, &
0.3056814182538295_dp, &
0.47227226126957916_dp, &
0.7255282521113733_dp, &
1.112683210316189_dp, &
1.7190059451835937_dp, &
2.7402160914459284_dp, &
4.918151105513321_dp /)
ELSE IF(E_range < 2101965.0_dp) THEN
aw(:) = (/ &
4.376362889407548e-07_dp, &
2.643443650482503e-06_dp, &
8.433114286038146e-06_dp, &
2.2439275175400065e-05_dp, &
5.458751916823435e-05_dp, &
0.00012469615312218161_dp, &
0.00027090154346575427_dp, &
0.0005643774226213476_dp, &
0.001134437804378847_dp, &
0.0022105491534259005_dp, &
0.004191467146338867_dp, &
0.007757292968144634_dp, &
0.014048601297429358_dp, &
0.024949196870902246_dp, &
0.04352732410856656_dp, &
0.07471686858797275_dp, &
0.12635898905420878_dp, &
0.21078138430975438_dp, &
0.3471766909066865_dp, &
0.5651735745080563_dp, &
0.9102276762862731_dp, &
1.451965060215205_dp, &
2.2979916853626663_dp, &
3.6203271934806356_dp, &
5.721128905491445_dp, &
9.287036012890901_dp, &
1.15698350323011e-06_dp, &
3.5239168925237668e-06_dp, &
8.790907369106611e-06_dp, &
2.0802473841355257e-05_dp, &
4.666225571530836e-05_dp, &
9.967519249166123e-05_dp, &
0.00020419454443829576_dp, &
0.0004036329368680413_dp, &
0.0007736164897230748_dp, &
0.0014432538211834921_dp, &
0.002629031662183179_dp, &
0.004688130413199208_dp, &
0.008201289641287034_dp, &
0.014100202125235302_dp, &
0.023861417673799547_dp, &
0.03979881543108159_dp, &
0.06550130495364699_dp, &
0.10648414776164583_dp, &
0.17115716235174305_dp, &
0.27227855433947445_dp, &
0.42921718067635106_dp, &
0.6717958687439105_dp, &
1.048026560686365_dp, &
1.6441450980118213_dp, &
2.656245656029467_dp, &
4.82228579092364_dp /)
ELSE IF(E_range < 14140999.0_dp) THEN
aw(:) = (/ &
9.23626066730038e-08_dp, &
6.508986506667573e-07_dp, &
2.548256199033377e-06_dp, &
8.121244988671366e-06_dp, &
2.280050217208556e-05_dp, &
5.838011937544951e-05_dp, &
0.00013919237871694767_dp, &
0.0003133658603796936_dp, &
0.0006727986081208618_dp, &
0.0013877632092122146_dp, &
0.002765631333703304_dp, &
0.005348711634689046_dp, &
0.010074592422860395_dp, &
0.018535050517138794_dp, &
0.03338852029027121_dp, &
0.05900918410430352_dp, &
0.10249708343303293_dp, &
0.17523445903983162_dp, &
0.2952639242558341_dp, &
0.4909040641368631_dp, &
0.8062587841364091_dp, &
1.3097798505481377_dp, &
2.1083909166084847_dp, &
3.3742290423813217_dp, &
5.410341564313699_dp, &
8.903545764950158_dp, &
2.5312239739000636e-07_dp, &
9.960996067355009e-07_dp, &
3.1628525158730217e-06_dp, &
8.836550056217969e-06_dp, &
2.2389728532275477e-05_dp, &
5.264557648617341e-05_dp, &
0.0001166810769940904_dp, &
0.00024642731945934143_dp, &
0.0004998727804625503_dp, &
0.0009797084412907111_dp, &
0.0018638236468942935_dp, &
0.00345441392005227_dp, &
0.0062559619517267165_dp, &
0.011097497158462704_dp, &
0.019322067389699597_dp, &
0.03307729968505806_dp, &
0.05575685817980752_dp, &
0.09266621527033955_dp, &
0.15202316630925258_dp, &
0.24647013442575488_dp, &
0.39542567933234785_dp, &
0.6290352095921129_dp, &
0.9959743276132635_dp, &
1.583367038071677_dp, &
2.5877759334182198_dp, &
4.74419254971056_dp /)
ELSE
aw(:) = (/ &
1.6637944502107968e-08_dp, &
1.8303408054343632e-07_dp, &
9.687195306133241e-07_dp, &
3.7389122528334316e-06_dp, &
1.1973447021251308e-05_dp, &
3.378096920589386e-05_dp, &
8.68338189070593e-05_dp, &
0.00020767968945116668_dp, &
0.00046874696407475663_dp, &
0.001008553647726778_dp, &
0.0020841046186256115_dp, &
0.0041598936699362606_dp, &
0.008056312504669037_dp, &
0.015193061067365586_dp, &
0.02798246122774749_dp, &
0.050456461866690214_dp, &
0.08925401954572831_dp, &
0.15515900119886547_dp, &
0.2654733104360881_dp, &
0.44765566677475566_dp, &
0.7449030289121367_dp, &
1.2248551649261998_dp, &
1.9939328236196414_dp, &
3.2243037705041826_dp, &
5.219599665831332_dp, &
8.666863277748455_dp, &
5.318821639180232e-08_dp, &
3.4945656362875746e-07_dp, &
1.4305518660978924e-06_dp, &
4.650456118314243e-06_dp, &
1.30922605762957e-05_dp, &
3.33205562494026e-05_dp, &
7.860562991572967e-05_dp, &
0.00017467466389322248_dp, &
0.0003697061169016621_dp, &
0.0007513106176187417_dp, &
0.0014748157046539873_dp, &
0.0028095645215491876_dp, &
0.005213528655212034_dp, &
0.009451865748843455_dp, &
0.016782910646881463_dp, &
0.029246610819691533_dp, &
0.05010725392331038_dp, &
0.08452724443086411_dp, &
0.14058736988255555_dp, &
0.23083686098028847_dp, &
0.37470521981828925_dp, &
0.6025298048228135_dp, &
0.9634165096071392_dp, &
1.5450976370872966_dp, &
2.5445132697805617_dp, &
4.694880341120372_dp /)
END IF
END SUBROUTINE

!SUBROUTINE get_coeff_26(k,E_range, aw, E_ratio)
!INTEGER, INTENT(IN)                                :: k
!REAL(KIND=dp), INTENT(IN)                          :: E_range
!REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
!REAL(KIND=dp), INTENT(OUT)                          :: E_ratio
!E_ratio = 1.0_dp
!IF(E_range < 799.0_dp) THEN
!E_ratio = 799.0_dp/E_range 
!aw(:) = (/ &
!0.00037867004555797054_dp, &
!0.002019373143542281_dp, &
!0.005074452374668418_dp, &
!0.009741339051894624_dp, &
!0.016351249099763064_dp, &
!0.025426574118427834_dp, &
!0.03775841037254848_dp, &
!0.05450081932003549_dp, &
!0.0772810957375601_dp, &
!0.1083363446608876_dp, &
!0.15069589349866513_dp, &
!0.20843076717336742_dp, &
!0.2869911650766578_dp, &
!0.3936564974164155_dp, &
!0.5381309695975518_dp, &
!0.7333307052011886_dp, &
!0.996427988952117_dp, &
!1.3502505778649674_dp, &
!1.8251937089224979_dp, &
!2.4619239536701665_dp, &
!3.3154239654767634_dp, &
!4.461579856601731_dp, &
!6.009262672635442_dp, &
!8.12630074937428_dp, &
!11.10929348730981_dp, &
!15.660501844751538_dp, &
!0.0009745069197489692_dp, &
!0.002323920525926464_dp, &
!0.003818499918990434_dp, &
!0.005569525222864092_dp, &
!0.007736056032946723_dp, &
!0.01054486610207425_dp, &
!0.014309337959120037_dp, &
!0.01944447015531593_dp, &
!0.026485944451580964_dp, &
!0.03612533116765179_dp, &
!0.04926670257709911_dp, &
!0.06710398765236568_dp, &
!0.09122022393749045_dp, &
!0.12371484118454482_dp, &
!0.16736960998240946_dp, &
!0.2258689406847251_dp, &
!0.30409900038431836_dp, &
!0.40856816694415354_dp, &
!0.548030961075663_dp, &
!0.7344894715293502_dp, &
!0.9849728438705977_dp, &
!1.3251021306723128_dp, &
!1.7972852200464182_dp, &
!2.4831111234607426_dp, &
!3.582848552118199_dp, &
!5.8876320817174905_dp /)
!ELSE IF(E_range < 995.0_dp) THEN
!aw(:) = (/ &
!0.00031416617596844365_dp, &
!0.0016768070133771967_dp, &
!0.004220387402167522_dp, &
!0.008122019548232636_dp, &
!0.013681535061623657_dp, &
!0.021375962641606915_dp, &
!0.03193242611639194_dp, &
!0.04641604327774391_dp, &
!0.0663338368572065_dp, &
!0.09376694799993701_dp, &
!0.13155072067060036_dp, &
!0.18352265415559924_dp, &
!0.2548586347835214_dp, &
!0.3525227140309385_dp, &
!0.4858648853075939_dp, &
!0.667414665262307_dp, &
!0.9139380450657599_dp, &
!1.2478576632771592_dp, &
!1.699195061182984_dp, &
!2.3083135493129188_dp, &
!3.13000622971369_dp, &
!4.240119346476375_dp, &
!5.747638276053859_dp, &
!7.820581475731412_dp, &
!10.755506760817696_dp, &
!15.2531099085639_dp, &
!0.0008086633399478686_dp, &
!0.001931716707111731_dp, &
!0.0031843549489835176_dp, &
!0.0046679761695975755_dp, &
!0.006529332077378777_dp, &
!0.008979243990088065_dp, &
!0.012309502604654843_dp, &
!0.016906748893981475_dp, &
!0.023272740862285972_dp, &
!0.03206195264016037_dp, &
!0.04413963056111703_dp, &
!0.060659468387992695_dp, &
!0.08316331669118938_dp, &
!0.11371010407957528_dp, &
!0.15504510373245767_dp, &
!0.21082547573489732_dp, &
!0.28592663287586667_dp, &
!0.38687165073350616_dp, &
!0.5224647206101295_dp, &
!0.7047999539903935_dp, &
!0.951040871895601_dp, &
!1.286968657510818_dp, &
!1.7551267580277088_dp, &
!2.437087428172129_dp, &
!3.5325565753903243_dp, &
!5.829332987897013_dp /)
!ELSE IF(E_range < 1293.0_dp) THEN
!aw(:) = (/ &
!0.000251216262630506_dp, &
!0.0013422657589387847_dp, &
!0.0033852831378795366_dp, &
!0.006535679555855439_dp, &
!0.01105946790435755_dp, &
!0.017384205344456737_dp, &
!0.02616631610050564_dp, &
!0.03837171615696803_dp, &
!0.05537294469791371_dp, &
!0.07907682374899021_dp, &
!0.11210158563747868_dp, &
!0.15802208108125124_dp, &
!0.22170328658555793_dp, &
!0.3097484547351596_dp, &
!0.43109811316492785_dp, &
!0.5978298541062674_dp, &
!0.826228906983008_dp, &
!1.1382318360717252_dp, &
!1.5634041077266412_dp, &
!2.1417299191833186_dp, &
!2.9277541939220253_dp, &
!3.9972537939524453_dp, &
!5.459336924183542_dp, &
!7.482255655971381_dp, &
!10.362568049094785_dp, &
!14.79927594838438_dp, &
!0.0006467895131307655_dp, &
!0.0015483789662149572_dp, &
!0.002562990408244027_dp, &
!0.00378128864955077_dp, &
!0.005336413195947385_dp, &
!0.0074211270754708444_dp, &
!0.010302498712647285_dp, &
!0.014334977955021215_dp, &
!0.01998237358756267_dp, &
!0.027857715884000284_dp, &
!0.038782119259036095_dp, &
!0.05386231670249176_dp, &
!0.07459076970770032_dp, &
!0.10297648694866778_dp, &
!0.1417181535853474_dp, &
!0.19443583384059848_dp, &
!0.26598602063677856_dp, &
!0.36290208001606944_dp, &
!0.4940399248668749_dp, &
!0.671597172994337_dp, &
!0.9128949281590978_dp, &
!1.2439090969479831_dp, &
!1.7073582868867894_dp, &
!2.3848234540233704_dp, &
!3.4754029453297375_dp, &
!5.763153512206444_dp /)
!ELSE IF(E_range < 1738.0_dp) THEN
!aw(:) = (/ &
!0.00019503948178662084_dp, &
!0.0010434689928647827_dp, &
!0.002638251061502539_dp, &
!0.005113339659677273_dp, &
!0.008700945623371394_dp, &
!0.01377846615011094_dp, &
!0.020929590132331147_dp, &
!0.031017311417159306_dp, &
!0.045274448736969015_dp, &
!0.06542674293351501_dp, &
!0.09386620906729039_dp, &
!0.13389212805274583_dp, &
!0.19004014176220868_dp, &
!0.2685271340668203_dp, &
!0.377850003355214_dp, &
!0.5295906140381842_dp, &
!0.7394997208371181_dp, &
!1.0289652688271176_dp, &
!1.4270284407794631_dp, &
!1.9732264763083136_dp, &
!2.721797922074536_dp, &
!3.748413528798017_dp, &
!5.162300567395363_dp, &
!7.131978662209736_dp, &
!9.954061908511639_dp, &
!14.325847031936839_dp, &
!0.0005023045111882079_dp, &
!0.001205645302440015_dp, &
!0.002005720903775344_dp, &
!0.002982372714059457_dp, &
!0.004254709674118375_dp, &
!0.005996408570346338_dp, &
!0.008448203182614166_dp, &
!0.011930700751986554_dp, &
!0.016868347596357594_dp, &
!0.02383094048001708_dp, &
!0.033592451809102954_dp, &
!0.047208001153773035_dp, &
!0.06611431621378733_dp, &
!0.09226257124098161_dp, &
!0.12829565471537094_dp, &
!0.17778658987626775_dp, &
!0.24556329719382933_dp, &
!0.3381617692032787_dp, &
!0.46448645895791835_dp, &
!0.6368429605901209_dp, &
!0.8727252415918798_dp, &
!1.19833159205944_dp, &
!1.6565930235234085_dp, &
!2.329134258703709_dp, &
!3.41444560424951_dp, &
!5.692651404356115_dp /)
!ELSE IF(E_range < 2238.0_dp) THEN
!aw(:) = (/ &
!0.00015700744438413856_dp, &
!0.0008410015394614405_dp, &
!0.0021312215547787303_dp, &
!0.004145567315278482_dp, &
!0.00709067459392277_dp, &
!0.01130543737244799_dp, &
!0.01731696799068697_dp, &
!0.025907708951358625_dp, &
!0.038201088680563586_dp, &
!0.05578098176408108_dp, &
!0.08086138835374904_dp, &
!0.11652302642738899_dp, &
!0.16703776837300044_dp, &
!0.2383097430965355_dp, &
!0.33847273237964964_dp, &
!0.478698088475134_dp, &
!0.6742883954381591_dp, &
!0.9461650804919068_dp, &
!1.3229159835766504_dp, &
!1.843683077594644_dp, &
!2.562423206657593_dp, &
!3.5546941289197296_dp, &
!4.929807721254272_dp, &
!6.85651240675138_dp, &
!9.63151015771553_dp, &
!13.950797550420011_dp, &
!0.00040446709905684556_dp, &
!0.0009731510923965196_dp, &
!0.0016264508648019713_dp, &
!0.0024359492637851146_dp, &
!0.0035098068247113295_dp, &
!0.00500647069683762_dp, &
!0.0071456103963289525_dp, &
!0.010221157599416823_dp, &
!0.01462679065931212_dp, &
!0.020898100778432942_dp, &
!0.029770825719871098_dp, &
!0.04225717041142872_dp, &
!0.05974650511131695_dp, &
!0.08413978412410032_dp, &
!0.11803010770982635_dp, &
!0.16494661750291933_dp, &
!0.22968739528440968_dp, &
!0.31878363471753624_dp, &
!0.44117322961627403_dp, &
!0.6092463802533127_dp, &
!0.8406397563067061_dp, &
!1.1617424275144823_dp, &
!1.6156772275544151_dp, &
!2.284131965422472_dp, &
!3.365137399027895_dp, &
!5.635683006283261_dp /)
!ELSE IF(E_range < 3009.0_dp) THEN
!aw(:) = (/ &
!0.00012173433282636892_dp, &
!0.000653044843078349_dp, &
!0.0016597094198001826_dp, &
!0.0032432122902587715_dp, &
!0.005583731099072486_dp, &
!0.00897974690437299_dp, &
!0.01389829106398585_dp, &
!0.021035914001762628_dp, &
!0.03139940878203811_dp, &
!0.046421111007702805_dp, &
!0.06812376578216424_dp, &
!0.09935124887637425_dp, &
!0.1440867725430818_dp, &
!0.20788864079466032_dp, &
!0.29848482929563186_dp, &
!0.42658280937440796_dp, &
!0.6069727081348725_dp, &
!0.8600354677167223_dp, &
!1.213825493528422_dp, &
!1.7070100456167858_dp, &
!2.393195963336712_dp, &
!3.3477843441213753_dp, &
!4.680165709088339_dp, &
!6.5593534066075065_dp, &
!9.282192058087176_dp, &
!13.54332387201712_dp, &
!0.00031370730222591084_dp, &
!0.0007570695990484558_dp, &
!0.0012727251500342805_dp, &
!0.0019236429401252328_dp, &
!0.0028062999069873792_dp, &
!0.004062578033501044_dp, &
!0.00588929507373728_dp, &
!0.00855186255429553_dp, &
!0.012411138722817702_dp, &
!0.017965594962442332_dp, &
!0.025908461038187995_dp, &
!0.03720320806718786_dp, &
!0.05318444551258281_dp, &
!0.0756939715032581_dp, &
!0.10726485539108421_dp, &
!0.1513713746489183_dp, &
!0.21277116226058032_dp, &
!0.297982246389678_dp, &
!0.41597262011480124_dp, &
!0.5792223633874778_dp, &
!0.8055284041136077_dp, &
!1.1215023110944804_dp, &
!1.5705011278690364_dp, &
!2.234313406989394_dp, &
!3.310497107984999_dp, &
!5.57261913216792_dp /)
!ELSE IF(E_range < 4377.0_dp) THEN
!aw(:) = (/ &
!8.814798852846974e-05_dp, &
!0.00047385601650182306_dp, &
!0.0012091635256913527_dp, &
!0.00237798129797295_dp, &
!0.004131735600921136_dp, &
!0.006724218659489275_dp, &
!0.010555223117605267_dp, &
!0.01622490234306836_dp, &
!0.024609130469194388_dp, &
!0.036969473964326045_dp, &
!0.05511131378702244_dp, &
!0.08160636091339551_dp, &
!0.12010211906510675_dp, &
!0.17574974333525284_dp, &
!0.2557934591366403_dp, &
!0.37038061132573763_dp, &
!0.5336740586707772_dp, &
!0.7653831510177921_dp, &
!1.0928877203185414_dp, &
!1.554240937533324_dp, &
!2.2025802278286_dp, &
!3.113073901350633_dp, &
!4.395180645931341_dp, &
!6.218240326595125_dp, &
!8.879328547943995_dp, &
!13.071604900598574_dp, &
!0.00022726312287597392_dp, &
!0.0005507558437273343_dp, &
!0.0009334429662410583_dp, &
!0.0014288412873932998_dp, &
!0.0021202405544394317_dp, &
!0.003130479736808873_dp, &
!0.004630332714659572_dp, &
!0.006853202347928061_dp, &
!0.010122916817408381_dp, &
!0.014895011218089086_dp, &
!0.02181220303758935_dp, &
!0.03177877870332068_dp, &
!0.04606157792471934_dp, &
!0.06642773094696337_dp, &
!0.09533265331311527_dp, &
!0.1361770083100393_dp, &
!0.1936600146159572_dp, &
!0.2742725721016829_dp, &
!0.3870076032574855_dp, &
!0.5444448648761351_dp, &
!0.7645726080680213_dp, &
!1.074280383110475_dp, &
!1.517232873363435_dp, &
!2.175382440935495_dp, &
!3.2457802833401543_dp, &
!5.498012547834142_dp /)
!ELSE IF(E_range < 6256.0_dp) THEN
!aw(:) = (/ &
!6.476501325165667e-05_dp, &
!0.0003489195644131545_dp, &
!0.0008941614534819361_dp, &
!0.0017704990846086513_dp, &
!0.003106215270875755_dp, &
!0.0051184583066800245_dp, &
!0.008151321905396521_dp, &
!0.012725005409661648_dp, &
!0.019606662602541135_dp, &
!0.029915132446101094_dp, &
!0.04527217925861043_dp, &
!0.06801672344353663_dp, &
!0.10150540639451665_dp, &
!0.1505320490421629_dp, &
!0.22191077397826406_dp, &
!0.3252842451465471_dp, &
!0.47424215018241783_dp, &
!0.68787066046948_dp, &
!0.9929124585786613_dp, &
!1.4268275212139947_dp, &
!2.042283983468357_dp, &
!2.914196304359969_dp, &
!4.152056157791546_dp, &
!5.925501234806853_dp, &
!8.531868206909438_dp, &
!12.663116079707873_dp, &
!0.00016705978905450625_dp, &
!0.0004066454411086679_dp, &
!0.0006951469119322108_dp, &
!0.0010783910686938938_dp, &
!0.0016286055279901323_dp, &
!0.0024524391864277316_dp, &
!0.0036988362171391674_dp, &
!0.005574672460740364_dp, &
!0.008372632557555961_dp, &
!0.01251108542263018_dp, &
!0.01858789523471959_dp, &
!0.027453747272404343_dp, &
!0.04031301377684959_dp, &
!0.05886268627839266_dp, &
!0.08548351213716461_dp, &
!0.12350296368570872_dp, &
!0.17755852252016796_dp, &
!0.2541057458406055_dp, &
!0.3621486681799956_dp, &
!0.5143476813428599_dp, &
!0.7288604796469867_dp, &
!1.032835319018796_dp, &
!1.4702376945502111_dp, &
!2.1232080691504955_dp, &
!3.1884003549497018_dp, &
!5.431940899454501_dp /)
!ELSE IF(E_range < 9034.0_dp) THEN
!aw(:) = (/ &
!4.7143308106155734e-05_dp, &
!0.0002546160297198654_dp, &
!0.0006556849915321849_dp, &
!0.001308486317589335_dp, &
!0.002321180528643633_dp, &
!0.00387851039714272_dp, &
!0.006274968017379002_dp, &
!0.009959528054236609_dp, &
!0.015602208675286744_dp, &
!0.02419326568455119_dp, &
!0.03718714187701131_dp, &
!0.05670799513076175_dp, &
!0.08584080260381065_dp, &
!0.12904154057230738_dp, &
!0.19271268775071665_dp, &
!0.2860078270209378_dp, &
!0.4219539039351692_dp, &
!0.6190165464585299_dp, &
!0.9032936535864567_dp, &
!1.3116326969167833_dp, &
!1.8962026105735332_dp, &
!2.731628719534309_dp, &
!3.9274056926834273_dp, &
!5.6534608469569685_dp, &
!8.207430563278576_dp, &
!12.280227735780262_dp, &
!0.0001216731978032898_dp, &
!0.00029765625617544494_dp, &
!0.0005138515834050829_dp, &
!0.0008093287778000748_dp, &
!0.001246309890847173_dp, &
!0.0019167138219384776_dp, &
!0.002949945927403237_dp, &
!0.00452921789279372_dp, &
!0.006918826042313398_dp, &
!0.010502325922280422_dp, &
!0.015834679544717904_dp, &
!0.023714461765780157_dp, &
!0.03528437327807364_dp, &
!0.05217097884790169_dp, &
!0.07667849879558422_dp, &
!0.11205727237433083_dp, &
!0.16287658808361544_dp, &
!0.23554754506032946_dp, &
!0.33907403585363055_dp, &
!0.4861853967624974_dp, &
!0.6951999533883163_dp, &
!0.9935236614313603_dp, &
!1.425435383549766_dp, &
!2.0732958594988293_dp, &
!3.133427851009757_dp, &
!5.368707926259402_dp /)
!ELSE IF(E_range < 15564.0_dp) THEN
!aw(:) = (/ &
!2.9443346537201843e-05_dp, &
!0.00015968862646421845_dp, &
!0.0004146519539836927_dp, &
!0.0008385653413963636_dp, &
!0.0015154712602891832_dp, &
!0.002590554358727564_dp, &
!0.004297459311925436_dp, &
!0.00699772969683458_dp, &
!0.011241369599203199_dp, &
!0.01785778451619208_dp, &
!0.028088988720722612_dp, &
!0.04378233724788014_dp, &
!0.06766745014703779_dp, &
!0.10375192365887785_dp, &
!0.15788397699593362_dp, &
!0.23854895177740337_dp, &
!0.3579930733907525_dp, &
!0.5338067690466547_dp, &
!0.7911614759088416_dp, &
!1.1660043326643865_dp, &
!1.7097462188341677_dp, &
!2.4965421329192803_dp, &
!3.6358381495699823_dp, &
!5.2979544926024476_dp, &
!7.781009835617463_dp, &
!11.774667668631706_dp, &
!7.606241270032894e-05_dp, &
!0.00018765408353617493_dp, &
!0.0003293790826976408_dp, &
!0.0005320967163836746_dp, &
!0.0008454744514330148_dp, &
!0.0013429979364498713_dp, &
!0.0021300803335643497_dp, &
!0.0033606393645051213_dp, &
!0.005262694960984085_dp, &
!0.008173982814950133_dp, &
!0.012591783653297441_dp, &
!0.019243367519600263_dp, &
!0.029185569928311852_dp, &
!0.043945023133826665_dp, &
!0.06571488466295647_dp, &
!0.09763015883017251_dp, &
!0.14415322687574247_dp, &
!0.21161732699122274_dp, &
!0.30900744619426596_dp, &
!0.44913076043380107_dp, &
!0.6505176919859489_dp, &
!0.9409366048953239_dp, &
!1.3651303838199809_dp, &
!2.0058242712612966_dp, &
!3.0589757886307543_dp, &
!5.283166557998725_dp /)
!ELSE IF(E_range < 24783.0_dp) THEN
!aw(:) = (/ &
!1.9680677968787818e-05_dp, &
!0.00010718429961200576_dp, &
!0.00028063369699425877_dp, &
!0.0005751240880955878_dp, &
!0.0010584229459454301_dp, &
!0.0018485803926858477_dp, &
!0.003137409980916865_dp, &
!0.005226190408208622_dp, &
!0.008581319789353085_dp, &
!0.013918389625070615_dp, &
!0.022326567153525306_dp, &
!0.03545075249936955_dp, &
!0.055756520521443115_dp, &
!0.08691314520509733_dp, &
!0.1343442294932413_dp, &
!0.20601527180136564_dp, &
!0.31355545994115364_dp, &
!0.4738517368085432_dp, &
!0.7113167007261072_dp, &
!1.0611431840271215_dp, &
!1.5740867657116888_dp, &
!2.3238702896001606_dp, &
!3.4198539990669645_dp, &
!5.032658096852942_dp, &
!7.460835156838619_dp, &
!11.393219939193445_dp, &
!5.08891079758869e-05_dp, &
!0.00012660396743151422_dp, &
!0.00022591978741528656_dp, &
!0.00037406357239589875_dp, &
!0.000611835881523469_dp, &
!0.0009998711964009784_dp, &
!0.0016270922615530513_dp, &
!0.002626827962152557_dp, &
!0.004200663882037452_dp, &
!0.006652026187910315_dp, &
!0.010434145701018064_dp, &
!0.016218886149508362_dp, &
!0.024995180362107215_dp, &
!0.03820910360067658_dp, &
!0.057962271996807264_dp, &
!0.08729190632278588_dp, &
!0.13056586817338334_dp, &
!0.19404239367659404_dp, &
!0.28667559539845255_dp, &
!0.4213185151511697_dp, &
!0.6166590713784111_dp, &
!0.9007550732823644_dp, &
!1.3187403136722335_dp, &
!1.9536779982378873_dp, &
!3.001314644694727_dp, &
!5.2169965073036355_dp /)
!ELSE IF(E_range < 41198.0_dp) THEN
!aw(:) = (/ &
!1.2674109022592262e-05_dp, &
!6.938961424250064e-05_dp, &
!0.0001836109117095651_dp, &
!0.0003826798510826298_dp, &
!0.000720207804815144_dp, &
!0.0012903906777303414_dp, &
!0.002248248117076268_dp, &
!0.003841754605276488_dp, &
!0.006462394751163588_dp, &
!0.01072216143228886_dp, &
!0.017568891891104072_dp, &
!0.02845746098996065_dp, &
!0.04560203213902193_dp, &
!0.07234523770656233_dp, &
!0.11369507364927639_dp, &
!0.1771011885774494_dp, &
!0.2735717887355965_dp, &
!0.41927521216040264_dp, &
!0.6378360441380629_dp, &
!0.9636484216701742_dp, &
!1.4467556672020363_dp, &
!2.16038789497357_dp, &
!3.213772136975622_dp, &
!4.777817256471035_dp, &
!7.151557096116501_dp, &
!11.023127187346768_dp, &
!3.281013826130679e-05_dp, &
!8.249630450619683e-05_dp, &
!0.0001503190811628769_dp, &
!0.00025652097657410326_dp, &
!0.00043391341684158855_dp, &
!0.0007317679654007801_dp, &
!0.0012244514824103281_dp, &
!0.0020265533742082715_dp, &
!0.003314844387972562_dp, &
!0.005359867459588524_dp, &
!0.008571952178977296_dp, &
!0.013568178848330984_dp, &
!0.021269287959126144_dp, &
!0.03303908704226431_dp, &
!0.05088391255416086_dp, &
!0.07773680891852074_dp, &
!0.11786157641848992_dp, &
!0.17742868614688642_dp, &
!0.2653463049634412_dp, &
!0.39449860788459795_dp, &
!0.5837216700811525_dp, &
!0.8613662154615228_dp, &
!1.2729809688101992_dp, &
!1.9020151440554_dp, &
!2.9440727560738815_dp, &
!5.151371570579308_dp /)
!ELSE IF(E_range < 94407.0_dp) THEN
!aw(:) = (/ &
!6.1864873133800315e-06_dp, &
!3.42333321344293e-05_dp, &
!9.256130883872102e-05_dp, &
!0.0001995146059753597_dp, &
!0.00039178329741312076_dp, &
!0.0007348455222213412_dp, &
!0.0013393442963486831_dp, &
!0.0023881425887916493_dp, &
!0.004179362781511009_dp, &
!0.007193111088055994_dp, &
!0.012193614356274658_dp, &
!0.020384101548649825_dp, &
!0.0336396676142171_dp, &
!0.054854588656627214_dp, &
!0.08845646415130372_dp, &
!0.14116207201192815_dp, &
!0.2230818181620526_dp, &
!0.34932588584464863_dp, &
!0.5423350524867223_dp, &
!0.8352752656174427_dp, &
!1.277060465461711_dp, &
!1.940097049357285_dp, &
!2.9333265745101387_dp, &
!4.428044567563_dp, &
!6.724058005061396_dp, &
!10.508725651867227_dp, &
!1.6052719564954786e-05_dp, &
!4.1236826436924364e-05_dp, &
!7.834655867964726e-05_dp, &
!0.00014151021463892604_dp, &
!0.00025367892790421795_dp, &
!0.00045044664188715174_dp, &
!0.0007883285497964747_dp, &
!0.0013578703683878816_dp, &
!0.0023029258087636224_dp, &
!0.0038494193775215953_dp, &
!0.0063483292722449425_dp, &
!0.01033949449874093_dp, &
!0.016645622593311638_dp, &
!0.026509804263420932_dp, &
!0.041795358134028064_dp, &
!0.06527467086334215_dp, &
!0.10104510956132512_dp, &
!0.15512782675831138_dp, &
!0.2363367630969623_dp, &
!0.3575720701278201_dp, &
!0.5378626888119358_dp, &
!0.8059834516089599_dp, &
!1.2081217349019628_dp, &
!1.8283708436585544_dp, &
!2.862257407220438_dp, &
!5.0576827215027595_dp /)
!ELSE IF(E_range < 189080.0_dp) THEN
!aw(:) = (/ &
!3.3986206749267677e-06_dp, &
!1.902518391414751e-05_dp, &
!5.2661237083094025e-05_dp, &
!0.00011755363584667629_dp, &
!0.00024054425131242777_dp, &
!0.0004703649078539346_dp, &
!0.0008915941462073683_dp, &
!0.0016480542439914962_dp, &
!0.0029804726581487105_dp, &
!0.005285900202059813_dp, &
!0.009210321167109866_dp, &
!0.01579149035429721_dp, &
!0.026677101018615843_dp, &
!0.04445499620476586_dp, &
!0.07314875020986866_dp, &
!0.11895565544283601_dp, &
!0.1913380136358924_dp, &
!0.3046275682444795_dp, &
!0.48037629526556186_dp, &
!0.7508062947438406_dp, &
!1.1639374578031345_dp, &
!1.7914918915391436_dp, &
!2.7421280756252386_dp, &
!4.1873916307035435_dp, &
!6.427703228089019_dp, &
!10.150031607977615_dp, &
!8.840847793682273e-06_dp, &
!2.324324578616612e-05_dp, &
!4.6141659827954215e-05_dp, &
!8.79876463334303e-05_dp, &
!0.00016589132538121333_dp, &
!0.0003074672552501645_dp, &
!0.000558351334396754_dp, &
!0.000993726386236742_dp, &
!0.001735781672681653_dp, &
!0.0029804152406535532_dp, &
!0.0050378351578707235_dp, &
!0.00839373945692239_dp, &
!0.013800708853496916_dp, &
!0.02241363578850322_dp, &
!0.03598893766310146_dp, &
!0.05717572165982958_dp, &
!0.08993924338093535_dp, &
!0.14017564828256868_dp, &
!0.21660902810458357_dp, &
!0.3321276618732292_dp, &
!0.5058815925186743_dp, &
!0.7669503755350563_dp, &
!1.162011318543182_dp, &
!1.7756900338995414_dp, &
!2.8035570812969_dp, &
!4.990537806043327_dp /)
!ELSE IF(E_range < 457444.0_dp) THEN
!aw(:) = (/ &
!1.5924593122033248e-06_dp, &
!9.091912376803294e-06_dp, &
!2.617541550484418e-05_dp, &
!6.171139320571212e-05_dp, &
!0.00013395538711265_dp, &
!0.0002770337626188815_dp, &
!0.0005524463258968383_dp, &
!0.001068584922831763_dp, &
!0.0020127793539792362_dp, &
!0.003703027242240496_dp, &
!0.006670388431813212_dp, &
!0.011788628295122639_dp, &
!0.020475916475588652_dp, &
!0.03500536896338614_dp, &
!0.058978575947671516_dp, &
!0.09804126384191222_dp, &
!0.1609562140010296_dp, &
!0.26120070619246305_dp, &
!0.4193334613096446_dp, &
!0.6665002448586814_dp, &
!1.0496753145003448_dp, &
!1.6397484503431214_dp, &
!2.544989612656627_dp, &
!3.9371793139954767_dp, &
!6.117452207466385_dp, &
!9.772518958696445_dp, &
!4.159960277698401e-06_dp, &
!1.1373804071332653e-05_dp, &
!2.420755936335119e-05_dp, &
!4.980565992350797e-05_dp, &
!0.00010018471872308093_dp, &
!0.00019590707898534102_dp, &
!0.00037246585760723267_dp, &
!0.0006901962557824386_dp, &
!0.0012498574384655588_dp, &
!0.0022170451967030076_dp, &
!0.0038600286020127068_dp, &
!0.006607762974687318_dp, &
!0.0111379581375363_dp, &
!0.01850955210435411_dp, &
!0.03036032305211595_dp, &
!0.04919945906832491_dp, &
!0.07883804989723764_dp, &
!0.12502029334390763_dp, &
!0.1963511883011235_dp, &
!0.3056817220164481_dp, &
!0.472272577604597_dp, &
!0.7255284874947008_dp, &
!1.1126832375494256_dp, &
!1.7190056146581998_dp, &
!2.7402151858607713_dp, &
!4.918149142282079_dp /)
!ELSE IF(E_range < 2101965.0_dp) THEN
!aw(:) = (/ &
!4.3763782671313535e-07_dp, &
!2.6434563079962914e-06_dp, &
!8.433158677763288e-06_dp, &
!2.243939079033379e-05_dp, &
!5.458781239810723e-05_dp, &
!0.0001246968887131766_dp, &
!0.0002709032279176203_dp, &
!0.0005643809087532476_dp, &
!0.0011344444384498063_dp, &
!0.0022105609486985167_dp, &
!0.0041914876817813265_dp, &
!0.007757329010348183_dp, &
!0.014048664192222513_dp, &
!0.024949303411478253_dp, &
!0.04352749770037569_dp, &
!0.07471713938183712_dp, &
!0.1263593964480408_dp, &
!0.2107819821054517_dp, &
!0.34717755412498486_dp, &
!0.5651748081810047_dp, &
!0.9102294273524391_dp, &
!1.4519675215115317_dp, &
!2.2979950661737902_dp, &
!3.6203316604049878_dp, &
!5.721134535152919_dp, &
!9.287042786628001_dp, &
!1.1569880004787854e-06_dp, &
!3.523936820119075e-06_dp, &
!8.790953146508159e-06_dp, &
!2.0802580807693856e-05_dp, &
!4.6662530125936186e-05_dp, &
!9.967584450657681e-05_dp, &
!0.00020419584362552093_dp, &
!0.00040363532938237734_dp, &
!0.000773620486252419_dp, &
!0.001443260395829438_dp, &
!0.00262904310158851_dp, &
!0.0046881507710751745_dp, &
!0.008201323868377521_dp, &
!0.014100256343177825_dp, &
!0.023861498617344504_dp, &
!0.03979893039493727_dp, &
!0.06550146555914389_dp, &
!0.10648437147720813_dp, &
!0.1711574744018939_dp, &
!0.27227899024551494_dp, &
!0.4292177878516702_dp, &
!0.671796685405684_dp, &
!1.0480275756498045_dp, &
!1.6441462396403197_dp, &
!2.656246822368235_dp, &
!4.822286914504179_dp /)
!ELSE IF(E_range < 14140999.0_dp) THEN
!aw(:) = (/ &
!9.23624644128009e-08_dp, &
!6.508984293721556e-07_dp, &
!2.548260568559413e-06_dp, &
!8.121270144006962e-06_dp, &
!2.280058581644762e-05_dp, &
!5.8380296942493665e-05_dp, &
!0.00013919273921343152_dp, &
!0.00031336665128641237_dp, &
!0.0006728000565273929_dp, &
!0.0013877655298111305_dp, &
!0.002765635340538635_dp, &
!0.005348719455124167_dp, &
!0.010074608651547055_dp, &
!0.018535082073824607_dp, &
!0.033388575519434034_dp, &
!0.059009270960414834_dp, &
!0.10249720956254744_dp, &
!0.17523464318825271_dp, &
!0.295264210048237_dp, &
!0.49090452525406647_dp, &
!0.8062595222094527_dp, &
!1.3097809890691694_dp, &
!2.1083926132817887_dp, &
!3.374231527201458_dp, &
!5.410345151018989_dp, &
!8.903550827247484_dp, &
!2.531219703616332e-07_dp, &
!9.961008767194328e-07_dp, &
!3.162861663769258e-06_dp, &
!8.836587416639661e-06_dp, &
!2.238980678804299e-05_dp, &
!5.264569178724729e-05_dp, &
!0.0001166813619734366_dp, &
!0.0002464278889725518_dp, &
!0.0004998735078142292_dp, &
!0.0009797095669868766_dp, &
!0.0018638260888474189_dp, &
!0.0034544195528519897_dp, &
!0.006255973531337034_dp, &
!0.01109751648307314_dp, &
!0.01932209536383217_dp, &
!0.033077334643883766_dp, &
!0.05575690357295425_dp, &
!0.09266629027020219_dp, &
!0.15202329968174552_dp, &
!0.24647035652570645_dp, &
!0.3954260144539292_dp, &
!0.6290356804120132_dp, &
!0.9959749856781216_dp, &
!1.5833679713133475_dp, &
!2.5877772123240255_dp, &
!4.744194244583471_dp /)
!ELSE
!aw(:) = (/ &
!1.6637778692707924e-08_dp, &
!1.830325030227858e-07_dp, &
!9.68711373372054e-07_dp, &
!3.7388846034199107e-06_dp, &
!1.1973378330546111e-05_dp, &
!3.3780851112299e-05_dp, &
!8.683366737789759e-05_dp, &
!0.00020767934360549128_dp, &
!0.0004687458937392337_dp, &
!0.0010085507285261743_dp, &
!0.0020840977236299507_dp, &
!0.004159880282674519_dp, &
!0.008056290348835145_dp, &
!0.015193026788229584_dp, &
!0.027982410796583163_dp, &
!0.05045639323084222_dp, &
!0.08925393411518674_dp, &
!0.1551589049169124_dp, &
!0.2654732187906678_dp, &
!0.44765562892410937_dp, &
!0.7449031470687724_dp, &
!1.2248555688817835_dp, &
!1.9939336234567542_dp, &
!3.2243050783371476_dp, &
!5.219601772431826_dp, &
!8.66686700632047_dp, &
!5.3187750083822314e-08_dp, &
!3.49453598783478e-07_dp, &
!1.4305400504162423e-06_dp, &
!4.650427738791991e-06_dp, &
!1.3092207927518681e-05_dp, &
!3.332051922014777e-05_dp, &
!7.860557022994064e-05_dp, &
!0.00017467427142731318_dp, &
!0.0003697049726865861_dp, &
!0.000751307870394084_dp, &
!0.001474810428344357_dp, &
!0.002809556909505189_dp, &
!0.005213518516350395_dp, &
!0.009451851458463308_dp, &
!0.01678289302538877_dp, &
!0.029246592582298937_dp, &
!0.05010723935290573_dp, &
!0.08452723843244708_dp, &
!0.14058739074942753_dp, &
!0.23083695752353275_dp, &
!0.3747054404015692_dp, &
!0.6025301514677136_dp, &
!0.9634169520729101_dp, &
!1.5450982380645015_dp, &
!2.5445143484550763_dp, &
!4.694882729929951_dp /)
!END IF
!END SUBROUTINE

SUBROUTINE get_coeff_28(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio
E_ratio = 1.0_dp
IF(E_range < 1545.0_dp) THEN
E_ratio = 1545.0_dp/E_range
aw(:) = (/ &
0.00019912749139307142_dp, &
0.0010623606743332262_dp, &
0.002671743519332409_dp, &
0.00513532276440788_dp, &
0.008635195625198482_dp, &
0.013459847519256876_dp, &
0.020047630847673463_dp, &
0.029039389901087887_dp, &
0.04134031090156478_dp, &
0.05819756722342521_dp, &
0.0813050912035362_dp, &
0.11294734291964816_dp, &
0.156193889583069_dp, &
0.21515890144425773_dp, &
0.2953443526600425_dp, &
0.4040922882505582_dp, &
0.5511803231110346_dp, &
0.7496069714265129_dp, &
1.016632455219541_dp, &
1.3751727396754474_dp, &
1.8557042554346754_dp, &
2.4989586763196456_dp, &
3.3599578911008456_dp, &
4.514593597566258_dp, &
6.071708534726409_dp, &
8.199087984832916_dp, &
11.19334614657497_dp, &
15.75712024734422_dp, &
0.0005125044577801448_dp, &
0.0012232237231645984_dp, &
0.00201318827111783_dp, &
0.0029437832219448537_dp, &
0.004103311596300065_dp, &
0.005618171589079525_dp, &
0.007663188401742024_dp, &
0.010469929287815255_dp, &
0.01433805703894519_dp, &
0.01965640843594298_dp, &
0.026936172517646257_dp, &
0.036855644579373535_dp, &
0.05031748334107659_dp, &
0.06852207829804384_dp, &
0.09306268185731213_dp, &
0.12604985461633972_dp, &
0.1702755103243793_dp, &
0.22943160937864643_dp, &
0.30840755489079047_dp, &
0.4137086894086003_dp, &
0.5540782272946391_dp, &
0.741497238773558_dp, &
0.9929644420040254_dp, &
1.3340651862992305_dp, &
1.8071782497597113_dp, &
2.4938998954135005_dp, &
3.5946345041711467_dp, &
5.901304226915114_dp /)
ELSE IF(E_range < 2002.0_dp) THEN
aw(:) = (/ &
0.00015937273545708875_dp, &
0.0008511031358508971_dp, &
0.0021444475784000504_dp, &
0.004133806965295614_dp, &
0.006979927658415096_dp, &
0.010939953042854741_dp, &
0.016407375175241108_dp, &
0.02396000912675003_dp, &
0.034417229669099075_dp, &
0.04891406525328216_dp, &
0.06900316741446434_dp, &
0.09679560984224389_dp, &
0.13515199610246875_dp, &
0.18793838914519878_dp, &
0.2603666777263224_dp, &
0.3594458206099125_dp, &
0.494579482258333_dp, &
0.6783583159834518_dp, &
0.927614454599119_dp, &
1.2648378420472655_dp, &
1.7201130600602677_dp, &
2.3338552489887388_dp, &
3.1608903963862245_dp, &
4.2770721150582895_dp, &
5.7913663447139845_dp, &
7.871758093638662_dp, &
10.814809416909151_dp, &
15.321476179970999_dp, &
0.00041027809665924954_dp, &
0.0009811723452880425_dp, &
0.001620915236450962_dp, &
0.0023840944111134044_dp, &
0.0033503311465523967_dp, &
0.004634503191430588_dp, &
0.006395709602901231_dp, &
0.008844913900028784_dp, &
0.012257024616981905_dp, &
0.016993101307172367_dp, &
0.02353379728419698_dp, &
0.03252363542798268_dp, &
0.04482787026200198_dp, &
0.06160612640371776_dp, &
0.08440876245542936_dp, &
0.11530375228655552_dp, &
0.15704466612137769_dp, &
0.21329510792248377_dp, &
0.28893382637364007_dp, &
0.39048261973383075_dp, &
0.5267381633327601_dp, &
0.7097795124230211_dp, &
0.9567475153369529_dp, &
1.2933958092367415_dp, &
1.762243934642854_dp, &
2.4448656022179347_dp, &
3.5410603027618426_dp, &
5.839188641520334_dp /)
ELSE IF(E_range < 2600.0_dp) THEN
aw(:) = (/ &
0.00012724350537371705_dp, &
0.0006802420548887551_dp, &
0.0017174005833444652_dp, &
0.0033210405572595314_dp, &
0.005632832416553322_dp, &
0.008881578544535886_dp, &
0.013419698685145537_dp, &
0.019766846961404197_dp, &
0.028663060595610065_dp, &
0.041139669406250204_dp, &
0.058618369068556976_dp, &
0.08304861028830621_dp, &
0.11709467333563459_dp, &
0.16438742365500805_dp, &
0.2298611471707931_dp, &
0.32020292243210463_dp, &
0.44445136041510835_dp, &
0.6147946527884036_dp, &
0.8476375115705392_dp, &
1.1650387509190727_dp, &
1.5966797391113152_dp, &
2.1826420752053064_dp, &
2.9775354931417777_dp, &
4.057156308900191_dp, &
5.530582551166736_dp, &
7.5660060038749855_dp, &
10.45997985096501_dp, &
14.911921279047414_dp, &
0.00032764611790808316_dp, &
0.0007852281095230586_dp, &
0.001302499732766977_dp, &
0.001927928018050743_dp, &
0.0027331746862386277_dp, &
0.0038223319490506384_dp, &
0.005339653436386936_dp, &
0.0074768267093066575_dp, &
0.010485860861582824_dp, &
0.01470204224626301_dp, &
0.020577182689456824_dp, &
0.028723198501018832_dp, &
0.03996852241361688_dp, &
0.05543195597157087_dp, &
0.07662014166796133_dp, &
0.10555670439639143_dp, &
0.14495398753342661_dp, &
0.1984431170697882_dp, &
0.27088687577758624_dp, &
0.36881736713632424_dp, &
0.501078546468407_dp, &
0.6798426309902373_dp, &
0.9223911021644032_dp, &
1.2546500752337444_dp, &
1.719292335657813_dp, &
2.3978939243834243_dp, &
3.489702010303218_dp, &
5.779704886427998_dp /)
ELSE IF(E_range < 3300.0_dp) THEN
aw(:) = (/ &
0.0001035812061619439_dp, &
0.000554312958948722_dp, &
0.0014022178351092648_dp, &
0.0027199186238307522_dp, &
0.004633626332072742_dp, &
0.007348901862042884_dp, &
0.011184100569752935_dp, &
0.016610348266822107_dp, &
0.02430141997172798_dp, &
0.03520198814732897_dp, &
0.05062394356636958_dp, &
0.0723803681708885_dp, &
0.10296857818003698_dp, &
0.14581768397450023_dp, &
0.2056217187942231_dp, &
0.28878665461351954_dp, &
0.4040292970175386_dp, &
0.5631795391897031_dp, &
0.7822574690472879_dp, &
1.0829291748663241_dp, &
1.4945032140919154_dp, &
2.056746330708999_dp, &
2.8240560994099004_dp, &
3.87216062398999_dp, &
5.31022885971396_dp, &
7.306643557351196_dp, &
10.157983272563953_dp, &
14.562388075641433_dp, &
0.000266779348158189_dp, &
0.0006406783010995833_dp, &
0.0010669475657666657_dp, &
0.0015890577230093321_dp, &
0.0022720547662622424_dp, &
0.003210897076964011_dp, &
0.004537212534376952_dp, &
0.006426464787961385_dp, &
0.009111540381231737_dp, &
0.012906037854636046_dp, &
0.018237038812244945_dp, &
0.025687953276114096_dp, &
0.036054499502532374_dp, &
0.05041868971863853_dp, &
0.070247180509793_dp, &
0.0975222869562315_dp, &
0.134916933126155_dp, &
0.1860296678962944_dp, &
0.2557045471811927_dp, &
0.3504778528812439_dp, &
0.4792309189974202_dp, &
0.6542157222577841_dp, &
0.8928389893219679_dp, &
1.2211855078507017_dp, &
1.6820762830637361_dp, &
2.357109509558417_dp, &
3.4450757129111187_dp, &
5.728067480861656_dp /)
ELSE IF(E_range < 4000.0_dp) THEN
aw(:) = (/ &
8.771437194659891e-05_dp, &
0.00046981282949284237_dp, &
0.0011904567715385075_dp, &
0.0023152670914235655_dp, &
0.003959200577898668_dp, &
0.006310722793342447_dp, &
0.009662905002781183_dp, &
0.01445070779560698_dp, &
0.021298505075191456_dp, &
0.031086296591640365_dp, &
0.04504371996352949_dp, &
0.06488109904164477_dp, &
0.09296905461591826_dp, &
0.13258247030343728_dp, &
0.1882303453180191_dp, &
0.26610050569058297_dp, &
0.3746580755427752_dp, &
0.5254504247939676_dp, &
0.7341916622638109_dp, &
1.0222323068269783_dp, &
1.418577696146867_dp, &
1.9627342158223895_dp, &
2.7089187483169885_dp, &
3.732793667761796_dp, &
5.14359311285564_dp, &
7.109854244554064_dp, &
9.928196448146211_dp, &
14.295810062691444_dp, &
0.0002259584596775983_dp, &
0.0005436014065931008_dp, &
0.0009083522994476501_dp, &
0.001360022022949313_dp, &
0.0019587330256291995_dp, &
0.0027925453027995255_dp, &
0.003983542841235872_dp, &
0.005695043198214815_dp, &
0.008145671980235551_dp, &
0.011632693203529134_dp, &
0.01656425108725137_dp, &
0.023501595060412054_dp, &
0.03321470050609412_dp, &
0.04675631004359308_dp, &
0.06556089291033432_dp, &
0.09157703564672734_dp, &
0.1274448334217016_dp, &
0.17673475828550625_dp, &
0.2442731132825866_dp, &
0.33659611474367285_dp, &
0.4626113339076691_dp, &
0.6346313704421029_dp, &
0.8701615672163753_dp, &
1.1954150913540855_dp, &
1.6533375661453558_dp, &
2.325557793916785_dp, &
3.410528405591019_dp, &
5.688122811462863_dp /)
ELSE IF(E_range < 5000.0_dp) THEN
aw(:) = (/ &
7.231152188478278e-05_dp, &
0.00038772726865815517_dp, &
0.0009844837786035514_dp, &
0.0019209126836657545_dp, &
0.003300157154162172_dp, &
0.005292557663606624_dp, &
0.008164163748376724_dp, &
0.012311171926660728_dp, &
0.018304991861179568_dp, &
0.026956176212698496_dp, &
0.03940566382368068_dp, &
0.057252298754729467_dp, &
0.08272831849932695_dp, &
0.11893897414817305_dp, &
0.17018831821850053_dp, &
0.2424208413642749_dp, &
0.34381888188337045_dp, &
0.48560994949960345_dp, &
0.68315888981334_dp, &
0.9574526915358277_dp, &
1.3371435383436294_dp, &
1.8614300851160115_dp, &
2.5843101667934096_dp, &
3.5813587287033566_dp, &
4.961876490500659_dp, &
6.894579341534884_dp, &
9.676154843478928_dp, &
14.002777296643712_dp, &
0.00018632498066864197_dp, &
0.00044921855213653585_dp, &
0.0007537638296544074_dp, &
0.0011359077842782092_dp, &
0.0016504930186307747_dp, &
0.002378087896208448_dp, &
0.0034304274436559386_dp, &
0.004957778020056624_dp, &
0.007163453510160269_dp, &
0.010326973559848835_dp, &
0.0148355888498722_dp, &
0.02122577871325675_dp, &
0.030238442499465622_dp, &
0.04289294382786412_dp, &
0.06058666241299866_dp, &
0.08522882557656473_dp, &
0.11942055225842414_dp, &
0.16669801861501937_dp, &
0.23186426244014152_dp, &
0.32145183283272727_dp, &
0.44439450297432675_dp, &
0.6130709874752484_dp, &
0.8450980113616935_dp, &
1.1668375144951064_dp, &
1.621384532244084_dp, &
2.290416624110755_dp, &
3.3720270933942214_dp, &
5.643641293866133_dp /)
ELSE IF(E_range < 5800.0_dp) THEN
aw(:) = (/ &
6.358752692633583e-05_dp, &
0.0003412046362365837_dp, &
0.0008676052314903817_dp, &
0.001696724618707933_dp, &
0.0029245257586060763_dp, &
0.004710233349785592_dp, &
0.007303215812339811_dp, &
0.011075692416141687_dp, &
0.016566304596355927_dp, &
0.024542556688053108_dp, &
0.03609017098515814_dp, &
0.05273820097743335_dp, &
0.07663171929303365_dp, &
0.11076846568614528_dp, &
0.15932179622752599_dp, &
0.22808006029257497_dp, &
0.32504298834606_dp, &
0.4612301705022223_dp, &
0.6517778001414668_dp, &
0.9174329759430668_dp, &
1.2866126278634253_dp, &
1.7983078281539955_dp, &
2.506365504219328_dp, &
3.4862959694100204_dp, &
4.847438621720606_dp, &
6.7586271477598405_dp, &
9.516604972737115_dp, &
13.816917061150304_dp, &
0.00016387366227471576_dp, &
0.0003956833555057065_dp, &
0.0006658661861317086_dp, &
0.0010080085317879625_dp, &
0.0014736812958346364_dp, &
0.002138763598196897_dp, &
0.00310853235933617_dp, &
0.004525172218371909_dp, &
0.006582500342505047_dp, &
0.0095488999069641_dp, &
0.013798331312490912_dp, &
0.01985134913295637_dp, &
0.028429999006008653_dp, &
0.04053182527263095_dp, &
0.05752974961085405_dp, &
0.08130676932209915_dp, &
0.11443765763930915_dp, &
0.16043487908983078_dp, &
0.2240845367272352_dp, &
0.311914711856157_dp, &
0.4328741422164402_dp, &
0.5993829816920151_dp, &
0.8291300891292139_dp, &
1.1485759543667498_dp, &
1.6009173789914857_dp, &
2.2678714811822975_dp, &
3.3473105643753676_dp, &
5.615101807506323_dp /)
ELSE IF(E_range < 7000.0_dp) THEN
aw(:) = (/ &
5.40208922943882e-05_dp, &
0.0002901574218009231_dp, &
0.0007392141971595068_dp, &
0.0014500283280455789_dp, &
0.002510177021344463_dp, &
0.0040657991238261465_dp, &
0.006346521670349522_dp, &
0.009696151598289092_dp, &
0.01461448675071523_dp, &
0.021817888400210955_dp, &
0.03232620376087572_dp, &
0.047584791601727995_dp, &
0.06963362235719665_dp, &
0.10134008926088696_dp, &
0.14671822700211287_dp, &
0.21136500303928102_dp, &
0.30305507802937337_dp, &
0.43255028902361986_dp, &
0.6147016140691686_dp, &
0.8699548471471977_dp, &
1.2264290200274872_dp, &
1.7228495271474296_dp, &
2.412866754832892_dp, &
3.371902180409746_dp, &
4.7093379807234985_dp, &
6.594155669520039_dp, &
9.323180769772668_dp, &
13.591212022178699_dp, &
0.00013925034927611353_dp, &
0.00033689758020930645_dp, &
0.0005691288628998808_dp, &
0.000866760965686723_dp, &
0.0012774771286669078_dp, &
0.0018715412342993915_dp, &
0.0027465275546555167_dp, &
0.004035039597245034_dp, &
0.0059195955044671665_dp, &
0.0086551765726244_dp, &
0.012599568036048066_dp, &
0.018253775832994857_dp, &
0.02631653618635365_dp, &
0.03775822887923194_dp, &
0.053921093309633066_dp, &
0.07665492497469895_dp, &
0.10850075471633389_dp, &
0.1529401192419044_dp, &
0.2147361484742657_dp, &
0.30040908293128193_dp, &
0.41892390451696465_dp, &
0.5827503975591537_dp, &
0.8096664957210812_dp, &
1.1262568684972982_dp, &
1.5758496574356475_dp, &
2.240219834134439_dp, &
3.316979586478506_dp, &
5.5800996039424104_dp /)
ELSE IF(E_range < 8500.0_dp) THEN
aw(:) = (/ &
4.564555265477108e-05_dp, &
0.0002454346685212948_dp, &
0.0006265784407134182_dp, &
0.001233158622447077_dp, &
0.0021448652153409154_dp, &
0.0034954182642228226_dp, &
0.005495608826652996_dp, &
0.00846210660098653_dp, &
0.012857585326682552_dp, &
0.01934938119006627_dp, &
0.028893873468546265_dp, &
0.04285529874482882_dp, &
0.06317112153384097_dp, &
0.0925808252186819_dp, &
0.13494115158016368_dp, &
0.19565897567855495_dp, &
0.2822840227867385_dp, &
0.405318874954234_dp, &
0.5793256641233082_dp, &
0.8244427232040404_dp, &
1.168482383587411_dp, &
1.6498933305311403_dp, &
2.322117411253598_dp, &
3.2604766389691147_dp, &
4.574390024972901_dp, &
6.432989341834178_dp, &
9.133196407474784_dp, &
13.369095389855303_dp, &
0.00011768969816977658_dp, &
0.00028534927032496446_dp, &
0.0004840726346740154_dp, &
0.0007420579710874475_dp, &
0.0011032591807000774_dp, &
0.0016325138217318865_dp, &
0.002419994753070502_dp, &
0.0035891515886961726_dp, &
0.005311642196826666_dp, &
0.007829394245850389_dp, &
0.011484246513613412_dp, &
0.016757755663488007_dp, &
0.024325289635545783_dp, &
0.03512979233715042_dp, &
0.05048226189921966_dp, &
0.07219835845965666_dp, &
0.1027839824530896_dp, &
0.14568785531725492_dp, &
0.20564777829998312_dp, &
0.2891735001868529_dp, &
0.4052436947742419_dp, &
0.5663760127138803_dp, &
0.7904374777760682_dp, &
1.1041398360339376_dp, &
1.5509491079980848_dp, &
2.21270825293388_dp, &
3.2867828515308894_dp, &
5.5452724054483005_dp /)
ELSE IF(E_range < 11000.0_dp) THEN
aw(:) = (/ &
3.6488006929306264e-05_dp, &
0.00019649151810584802_dp, &
0.000503107877439482_dp, &
0.000994818980265892_dp, &
0.0017419330458111762_dp, &
0.00286324584306399_dp, &
0.004546790823266163_dp, &
0.0070764254736432365_dp, &
0.01086988103008825_dp, &
0.01653489379720316_dp, &
0.024950193189557313_dp, &
0.0373800572984858_dp, &
0.05563473102373159_dp, &
0.08229382069514388_dp, &
0.12101608504744633_dp, &
0.1769674560118736_dp, &
0.25741052665895353_dp, &
0.3725145015640989_dp, &
0.5364671734722614_dp, &
0.7690049629152851_dp, &
1.0975352324996732_dp, &
1.5601365551801154_dp, &
2.209966148011392_dp, &
3.1222026085468673_dp, &
4.406302512891578_dp, &
6.231592566777427_dp, &
8.895138025566414_dp, &
13.09015522563044_dp, &
9.411058161764564e-05_dp, &
0.00022887458151702356_dp, &
0.00039057680221556677_dp, &
0.0006042814819193009_dp, &
0.0009094035245976271_dp, &
0.0013641332935019103_dp, &
0.002049648033435726_dp, &
0.003078322757949694_dp, &
0.004608548064372825_dp, &
0.006866051954909783_dp, &
0.010172626387128262_dp, &
0.014985141010389004_dp, &
0.021949052283394388_dp, &
0.03197188778843743_dp, &
0.046323942696241904_dp, &
0.0667759221190767_dp, &
0.09578682159309504_dp, &
0.13676065620689198_dp, &
0.19439935510505243_dp, &
0.27519525733539296_dp, &
0.3881406112452913_dp, &
0.5458114344320582_dp, &
0.7661884133896568_dp, &
1.0761497928442625_dp, &
1.519347409242135_dp, &
2.177726198492106_dp, &
3.2483564776686458_dp, &
5.5009815780981555_dp /)
ELSE IF(E_range < 14000.0_dp) THEN
aw(:) = (/ &
2.9588446308604547e-05_dp, &
0.00015957764714927812_dp, &
0.00040980061859442377_dp, &
0.0008141581268268158_dp, &
0.0014351942867949477_dp, &
0.002379218143911262_dp, &
0.003815130228038144_dp, &
0.0059991869926145835_dp, &
0.009311238414156518_dp, &
0.014308514375164259_dp, &
0.021803455514145805_dp, &
0.03297432754720971_dp, &
0.04952104060181237_dp, &
0.07388347387186651_dp, &
0.1095460460641794_dp, &
0.16146092263899975_dp, &
0.23663401107846255_dp, &
0.3449341441338246_dp, &
0.5002090250697223_dp, &
0.7218265976855273_dp, &
1.036819021323394_dp, &
1.4829173538285452_dp, &
2.113006031585928_dp, &
3.002118887270464_dp, &
4.25973634734891_dp, &
6.055362139002693_dp, &
8.68620971747294_dp, &
12.84476212055696_dp, &
7.634110961434247e-05_dp, &
0.00018622528846790902_dp, &
0.00031969127466536623_dp, &
0.0004991918529001643_dp, &
0.0007602900708133781_dp, &
0.0011555100382159716_dp, &
0.0017584338105163752_dp, &
0.0026721061529931525_dp, &
0.004043591457129711_dp, &
0.00608455281208252_dp, &
0.009099150682157118_dp, &
0.013522328657324851_dp, &
0.01997274372862434_dp, &
0.029325909214924503_dp, &
0.04281497070968849_dp, &
0.06216914194983893_dp, &
0.08980351557167271_dp, &
0.12907941549114965_dp, &
0.184663246594884_dp, &
0.2630278258762335_dp, &
0.3731736346682672_dp, &
0.5277260380561264_dp, &
0.7447670823860634_dp, &
1.051327938626868_dp, &
1.4912362410625788_dp, &
2.1465432019415753_dp, &
3.214074319529515_dp, &
5.461496007193266_dp /)
ELSE IF(E_range < 18000.0_dp) THEN
aw(:) = (/ &
2.3779414127343956e-05_dp, &
0.00012846422005240743_dp, &
0.0003309928831073909_dp, &
0.0006610841016675125_dp, &
0.001174111517851935_dp, &
0.001964735642046167_dp, &
0.0031839437046973372_dp, &
0.005062138452978201_dp, &
0.007943587365676456_dp, &
0.012337793711020667_dp, &
0.01899408520862755_dp, &
0.02900819673889348_dp, &
0.0439733623874622_dp, &
0.06619336766731646_dp, &
0.09898159274083276_dp, &
0.1470789732012216_dp, &
0.2172359337388872_dp, &
0.31902010927146557_dp, &
0.465935489961852_dp, &
0.676974407399601_dp, &
0.9787827564627841_dp, &
1.4087293747222094_dp, &
2.0194090526711066_dp, &
2.8856949554637783_dp, &
4.117081205191889_dp, &
5.883249820883165_dp, &
8.481580484406694_dp, &
12.603865500485162_dp, &
6.137647450224398e-05_dp, &
0.00015022918522230813_dp, &
0.00025961689799477495_dp, &
0.00040956404897364297_dp, &
0.0006319903102840593_dp, &
0.0009740518614949558_dp, &
0.0015022032585384295_dp, &
0.002310722155775527_dp, &
0.0035358702066416604_dp, &
0.005375683414999454_dp, &
0.008117044767788793_dp, &
0.012173222749539778_dp, &
0.018136165890969762_dp, &
0.026849214184554587_dp, &
0.03950783752538176_dp, &
0.057798731689314514_dp, &
0.08409142368704033_dp, &
0.1217021378276602_dp, &
0.17525857508663245_dp, &
0.2512102337617962_dp, &
0.3585618870750593_dp, &
0.5099851018896473_dp, &
0.7236623896941952_dp, &
1.0267809590048702_dp, &
1.4633527890610423_dp, &
2.1155495825866195_dp, &
3.1799713652899073_dp, &
5.422242322222096_dp /)
ELSE IF(E_range < 22000.0_dp) THEN
aw(:) = (/ &
1.996975719114667e-05_dp, &
0.00010803766024790043_dp, &
0.00027914946672295566_dp, &
0.0005600685201304379_dp, &
0.0010010455833094946_dp, &
0.0016883469263468778_dp, &
0.0027600190175514337_dp, &
0.004427777354772284_dp, &
0.0070100858816421925_dp, &
0.010981630036803518_dp, &
0.017045370162921285_dp, &
0.026235975775831418_dp, &
0.04006720477930825_dp, &
0.060740772707250645_dp, &
0.09144094582526643_dp, &
0.13674817658247623_dp, &
0.20321752055734807_dp, &
0.3001847449852204_dp, &
0.44088738725627163_dp, &
0.6440243823419953_dp, &
0.9359372859467849_dp, &
1.3537064463445627_dp, &
1.9496922567626198_dp, &
2.798633041764183_dp, &
4.010025863301125_dp, &
5.753689968035089_dp, &
8.327144493740365_dp, &
12.421679016019592_dp, &
5.156001410295644e-05_dp, &
0.00012656590579908283_dp, &
0.00021996531665695752_dp, &
0.00035003655094498806_dp, &
0.0005460402965281552_dp, &
0.0008512186113124776_dp, &
0.0013268718440107125_dp, &
0.0020609116400843674_dp, &
0.0031816393309671174_dp, &
0.004876914587104125_dp, &
0.007420595371219653_dp, &
0.011209484450019506_dp, &
0.016815094297439877_dp, &
0.02505596575962014_dp, &
0.037098296020097635_dp, &
0.05459544297770877_dp, &
0.07988080804958117_dp, &
0.11623432159124078_dp, &
0.16825175472580994_dp, &
0.2423620424395829_dp, &
0.3475704605516449_dp, &
0.4965817576093905_dp, &
0.7076548141522077_dp, &
1.0080986813480735_dp, &
1.4420727665642223_dp, &
2.09185113214404_dp, &
3.1538737345296277_dp, &
5.392217904025483_dp /)
ELSE IF(E_range < 30000.0_dp) THEN
aw(:) = (/ &
1.5245472100141073e-05_dp, &
8.267520676446116e-05_dp, &
0.0002146253100574728_dp, &
0.00043387740900171205_dp, &
0.0007836944027554964_dp, &
0.0013387930663292789_dp, &
0.0022193817580203824_dp, &
0.0036113909961874636_dp, &
0.005797519906650395_dp, &
0.009203847862920226_dp, &
0.01446808989009802_dp, &
0.022538313347044602_dp, &
0.034814713569521594_dp, &
0.05335209638799308_dp, &
0.0811475703355085_dp, &
0.12254732266532327_dp, &
0.18381921201396537_dp, &
0.27395570440350625_dp, &
0.40579685475434923_dp, &
0.5976003453364507_dp, &
0.875245772613848_dp, &
1.2753705825844541_dp, &
1.8499693672606738_dp, &
2.6735624079289133_dp, &
3.855638524201632_dp, &
5.566220151857392_dp, &
8.103052167483417_dp, &
12.156727650422702_dp, &
3.9383299458817714e-05_dp, &
9.713918675705757e-05_dp, &
0.000170421510611186_dp, &
0.0002751089399893439_dp, &
0.0004367526322332042_dp, &
0.0006931539398019088_dp, &
0.0010985049239842675_dp, &
0.0017318671480729233_dp, &
0.0027102869498213354_dp, &
0.004207039027136905_dp, &
0.006477134639756053_dp, &
0.009893360814081634_dp, &
0.014997191697894255_dp, &
0.02257042180249529_dp, &
0.03373549826844379_dp, &
0.05009547634855471_dp, &
0.07392862462849442_dp, &
0.10845863055245018_dp, &
0.15823054590725322_dp, &
0.22963859190692054_dp, &
0.33168416368000975_dp, &
0.47711716376192037_dp, &
0.6843078702870563_dp, &
0.9807484854746213_dp, &
1.4108257998358733_dp, &
2.0569814817847933_dp, &
3.1154409728311463_dp, &
5.348031605731007_dp /)
ELSE IF(E_range < 40000.0_dp) THEN
aw(:) = (/ &
1.1868158588050793e-05_dp, &
6.451496656973e-05_dp, &
0.00016828367389327262_dp, &
0.0003428138578813575_dp, &
0.0006257728657443166_dp, &
0.001082546664145262_dp, &
0.0018189122009131555_dp, &
0.002999899699939122_dp, &
0.004879035508307028_dp, &
0.007842398392055541_dp, &
0.012473504017331851_dp, &
0.019647823262643315_dp, &
0.030669526100067415_dp, &
0.04746815571454119_dp, &
0.07287994053111492_dp, &
0.11104807754173227_dp, &
0.16798957023910677_dp, &
0.2523945702323838_dp, &
0.37675011948627407_dp, &
0.5589184164507482_dp, &
0.8243607282457994_dp, &
1.209308359311203_dp, &
1.7654144169432378_dp, &
2.566987933665634_dp, &
3.7234969430883793_dp, &
5.405141983044582_dp, &
7.9098851591902335_dp, &
11.927751732780631_dp, &
3.067518035880976e-05_dp, &
7.602736152767003e-05_dp, &
0.00013466099350051709_dp, &
0.0002205156203375216_dp, &
0.00035609796830194946_dp, &
0.0005747749473125263_dp, &
0.0009249835582363322_dp, &
0.0014785254567841132_dp, &
0.0023430318256698933_dp, &
0.003679400712048183_dp, &
0.005726484701954953_dp, &
0.008836295353538897_dp, &
0.013524092658912494_dp, &
0.02053931124436599_dp, &
0.030965498069383788_dp, &
0.04636050178159656_dp, &
0.06895241489824191_dp, &
0.1019128924749107_dp, &
0.14973885102211165_dp, &
0.2187896005944566_dp, &
0.31805835198625043_dp, &
0.4603305714250192_dp, &
0.6640726964679041_dp, &
0.9569407810619277_dp, &
1.3835309605116226_dp, &
2.026448697565121_dp, &
3.0817527233540427_dp, &
5.309325429687417_dp /)
ELSE IF(E_range < 55000.0_dp) THEN
aw(:) = (/ &
8.994970387322457e-06_dp, &
4.9038651219381535e-05_dp, &
0.00012866021417541092_dp, &
0.00026454506965524195_dp, &
0.0004890218504779385_dp, &
0.0008585083873106919_dp, &
0.0014649004957659442_dp, &
0.002453066714438582_dp, &
0.00404818477662515_dp, &
0.006597100721999991_dp, &
0.010629669867507806_dp, &
0.016948840913414677_dp, &
0.026762034057250313_dp, &
0.041871547214080776_dp, &
0.06494887934178976_dp, &
0.09992774926372437_dp, &
0.15256426013426896_dp, &
0.2312316469143259_dp, &
0.3480438390880272_dp, &
0.520441341810499_dp, &
0.7734349455842019_dp, &
1.1428133579373827_dp, &
1.6798513048817534_dp, &
2.4586158447712725_dp, &
3.5885384348523433_dp, &
5.240004673003719_dp, &
7.7112226946250795_dp, &
11.691667476075265_dp, &
2.3263971840934205e-05_dp, &
5.799733961007195e-05_dp, &
0.00010391874644151258_dp, &
0.0001730995125618299_dp, &
0.0002850761818642179_dp, &
0.00046892999174830916_dp, &
0.0007675526539416234_dp, &
0.0012456339379473573_dp, &
0.0020014080059075953_dp, &
0.0031832610399331806_dp, &
0.005013553854404307_dp, &
0.007822909693482165_dp, &
0.012099368861470802_dp, &
0.018558470273632704_dp, &
0.02824263531163439_dp, &
0.04266142900862469_dp, &
0.06398871559990858_dp, &
0.09533906369189366_dp, &
0.14115536539598936_dp, &
0.20775581572804178_dp, &
0.3041200099484807_dp, &
0.4430661693538553_dp, &
0.6431592914657807_dp, &
0.9322296232686789_dp, &
1.355102086814327_dp, &
1.9945706936958254_dp, &
3.046542150566686_dp, &
5.268893998445498_dp /)
ELSE IF(E_range < 75000.0_dp) THEN
aw(:) = (/ &
6.8673479673260504e-06_dp, &
3.755590753670896e-05_dp, &
9.915139981458229e-05_dp, &
0.00020590986120619216_dp, &
0.00038570411969655687_dp, &
0.0006874225835250557_dp, &
0.0011912995282688835_dp, &
0.0020251938479038355_dp, &
0.0033901629916074914_dp, &
0.005599342672133319_dp, &
0.009136035642934201_dp, &
0.0147397005849582_dp, &
0.023532325447503034_dp, &
0.03720291582111918_dp, &
0.05827512485526659_dp, &
0.09049318420974399_dp, &
0.13937535368635193_dp, &
0.21300368025384273_dp, &
0.3231464746438852_dp, &
0.4868502227776408_dp, &
0.7287007191729763_dp, &
1.0840642579491595_dp, &
1.6038486600019533_dp, &
2.3618790984089095_dp, &
3.467539083950752_dp, &
5.09138219825452_dp, &
7.53185924341038_dp, &
11.47798005469671_dp, &
1.7773461657358013e-05_dp, &
4.458780785431919e-05_dp, &
8.088403265892346e-05_dp, &
0.00013715743270713957_dp, &
0.00023041539243485957_dp, &
0.00038613012816423253_dp, &
0.0006425163485643913_dp, &
0.001058147634119284_dp, &
0.0017230321030210073_dp, &
0.002774473725192771_dp, &
0.004420090523231148_dp, &
0.006971221027705693_dp, &
0.010891151319944746_dp, &
0.01686432115872218_dp, &
0.025895072083461815_dp, &
0.03944781877116618_dp, &
0.059645137430023996_dp, &
0.0895468281752171_dp, &
0.13354285706817015_dp, &
0.1979094322641605_dp, &
0.291608890942118_dp, &
0.42748519623553954_dp, &
0.6241916121496675_dp, &
0.9097205132707596_dp, &
1.3291156770565788_dp, &
1.965360170512719_dp, &
3.0142420912237333_dp, &
5.231825921337635_dp /)
ELSE IF(E_range < 100000.0_dp) THEN
aw(:) = (/ &
5.34722652436438e-06_dp, &
2.9335345788005576e-05_dp, &
7.794420445994551e-05_dp, &
0.0001635102977617264_dp, &
0.0003103380666014091_dp, &
0.0005612614423847721_dp, &
0.0009871305012123896_dp, &
0.0017020454525922442_dp, &
0.002887383643586052_dp, &
0.004828517011721044_dp, &
0.007970062065688042_dp, &
0.012998268728170836_dp, &
0.02096295289681554_dp, &
0.03345668031415268_dp, &
0.052876323139595754_dp, &
0.08280245241436206_dp, &
0.1285464334202786_dp, &
0.19793519202552176_dp, &
0.30243198793781245_dp, &
0.4587328138519049_dp, &
0.6910421452320579_dp, &
1.0343433822030452_dp, &
1.5392071335504423_dp, &
2.279230497059965_dp, &
3.3637429883435237_dp, &
4.96344347806973_dp, &
7.377008968134632_dp, &
11.293074733218035_dp, &
1.3848866712786266e-05_dp, &
3.496420645730977e-05_dp, &
6.422528074399328e-05_dp, &
0.00011085156828305275_dp, &
0.00018979249819623976_dp, &
0.0003236153897281577_dp, &
0.0005467478277421638_dp, &
0.000912713310709869_dp, &
0.001504624320381209_dp, &
0.0024504009658675406_dp, &
0.003945073258731463_dp, &
0.006283381657271798_dp, &
0.009907137885517564_dp, &
0.015473588963650402_dp, &
0.023953503710281564_dp, &
0.03677113226590553_dp, &
0.056002962059688105_dp, &
0.08465896050001478_dp, &
0.1270800799647972_dp, &
0.18950233376937137_dp, &
0.2808690547904898_dp, &
0.4140431391215438_dp, &
0.6077531842367035_dp, &
0.8901351666158656_dp, &
1.3064316369806335_dp, &
1.9398043844798791_dp, &
2.9859550412644587_dp, &
5.199382405454911_dp /)
ELSE IF(E_range < 140000.0_dp) THEN
aw(:) = (/ &
3.991479954961762e-06_dp, &
2.198720501971829e-05_dp, &
5.890526173006801e-05_dp, &
0.00012518031616235489_dp, &
0.00024153524567753135_dp, &
0.00044471016635781095_dp, &
0.0007960973675010948_dp, &
0.0013958401833713643_dp, &
0.0024051630337238337_dp, &
0.004080730226383008_dp, &
0.0068267899198047955_dp, &
0.011273593700680132_dp, &
0.018394388130399692_dp, &
0.029678633788751282_dp, &
0.04738664952634921_dp, &
0.07492145227987683_dp, &
0.11736834467836472_dp, &
0.18227349630314565_dp, &
0.2807619982420454_dp, &
0.42913828874392707_dp, &
0.65117720616762_dp, &
0.9814263687158781_dp, &
1.4700673778736164_dp, &
2.190427750165236_dp, &
3.2517634440541547_dp, &
4.824930255638227_dp, &
7.208868255000471_dp, &
11.091834807498854_dp, &
1.0346862342302347e-05_dp, &
2.6338249637987535e-05_dp, &
4.916412854301245e-05_dp, &
8.674834277534896e-05_dp, &
0.00015194758780428323_dp, &
0.00026440532385478294_dp, &
0.00045469618333409404_dp, &
0.0007711009718467_dp, &
0.0012894728448137978_dp, &
0.002127757911813279_dp, &
0.0034674947327675617_dp, &
0.005585492420042323_dp, &
0.008900179331682422_dp, &
0.014038952538275588_dp, &
0.02193542166778625_dp, &
0.033968992230147974_dp, &
0.05216418521519709_dp, &
0.07947410573244414_dp, &
0.12018284604662782_dp, &
0.18047837042907597_dp, &
0.26927867836352665_dp, &
0.3994632895161757_dp, &
0.5898412436801937_dp, &
0.8687080586761889_dp, &
1.2815327235597582_dp, &
1.91168824489995_dp, &
2.954800144592791_dp, &
5.163666162590388_dp /)
ELSE IF(E_range < 200000.0_dp) THEN
aw(:) = (/ &
2.9287097241919463e-06_dp, &
1.6211294489452685e-05_dp, &
4.3860626639419625e-05_dp, &
9.463310124321836e-05_dp, &
0.00018604994609815954_dp, &
0.0003493943765511619_dp, &
0.0006375682448360294_dp, &
0.0011380795201948043_dp, &
0.001993709375951546_dp, &
0.00343456904920158_dp, &
0.005827197217933671_dp, &
0.009749051421957427_dp, &
0.016100560249175083_dp, &
0.026272317263678777_dp, &
0.042392659447799845_dp, &
0.06769166630053426_dp, &
0.10703277867836304_dp, &
0.16768452950177593_dp, &
0.26043499640668427_dp, &
0.4011952042193507_dp, &
0.6133044483273886_dp, &
0.9308639628203111_dp, &
1.4036513854446626_dp, &
2.1047075037398466_dp, &
3.143200044406983_dp, &
4.6901364769963045_dp, &
7.044731962208944_dp, &
10.894906305627384_dp, &
7.599945789593838e-06_dp, &
1.9535340209099297e-05_dp, &
3.7160656136204285e-05_dp, &
6.722566694134913e-05_dp, &
0.00012069652213103665_dp, &
0.00021459840053591846_dp, &
0.0003759981286012543_dp, &
0.0006483001768113125_dp, &
0.0011005125730724488_dp, &
0.0018410833738741884_dp, &
0.0030385893339244744_dp, &
0.004952462825414513_dp, &
0.007978278215654642_dp, &
0.012714011523554197_dp, &
0.02005632833363888_dp, &
0.03133963591517064_dp, &
0.04853577656951428_dp, &
0.07453950556949848_dp, &
0.11357558363229604_dp, &
0.17178036896999477_dp, &
0.2580420507206907_dp, &
0.3852520170265286_dp, &
0.5722961014253114_dp, &
0.847629061135607_dp, &
1.2569519792288701_dp, &
1.8838624046566568_dp, &
2.923931086168498_dp, &
5.1282957309537895_dp /)
ELSE IF(E_range < 280000.0_dp) THEN
aw(:) = (/ &
2.187887975241367e-06_dp, &
1.217302105778866e-05_dp, &
3.3280277384650595e-05_dp, &
7.294651585196924e-05_dp, &
0.0001461477506435841_dp, &
0.0002798259055981046_dp, &
0.0005201025046636251_dp, &
0.0009442991970046586_dp, &
0.0016801646722205163_dp, &
0.0029359404461932336_dp, &
0.005046805926292481_dp, &
0.008545920293678393_dp, &
0.014272096616783494_dp, &
0.023531600046055763_dp, &
0.038339327744621156_dp, &
0.061775614958965434_dp, &
0.09851040024837446_dp, &
0.15556831587766237_dp, &
0.243439200558_dp, &
0.3776832590482275_dp, &
0.5812481524315448_dp, &
0.8878298248851654_dp, &
1.3468344739190594_dp, &
2.031033514764075_dp, &
3.0495035108172392_dp, &
4.573381799063933_dp, &
6.902137854735053_dp, &
10.723425645457905_dp, &
5.683864101635673e-06_dp, &
1.4761618602067145e-05_dp, &
2.863922084439926e-05_dp, &
5.31200028912619e-05_dp, &
9.765726111159933e-05_dp, &
0.00017719091014345152_dp, &
0.00031593721698955905_dp, &
0.0005532594185059319_dp, &
0.000952423486460829_dp, &
0.001613835283450162_dp, &
0.002695001718544619_dp, &
0.004440391442639413_dp, &
0.007225736946277607_dp, &
0.01162326005144731_dp, &
0.018497020133015078_dp, &
0.02914136949940742_dp, &
0.04548081139069081_dp, &
0.07035707490413387_dp, &
0.10794019701910468_dp, &
0.1643177452299961_dp, &
0.24834763708314758_dp, &
0.37292767817567907_dp, &
0.5570086406904123_dp, &
0.829186245882319_dp, &
1.2353723738001468_dp, &
1.8593753882752646_dp, &
2.896735710900212_dp, &
5.0971507815821235_dp /)
ELSE IF(E_range < 400000.0_dp) THEN
aw(:) = (/ &
1.6070031867276192e-06_dp, &
8.996223208889756e-06_dp, &
2.4903330256736718e-05_dp, &
5.5597271413213874e-05_dp, &
0.00011378123225419507_dp, &
0.00022251943620373886_dp, &
0.0004218452833351526_dp, &
0.0007798397462569304_dp, &
0.001410464737003683_dp, &
0.0025017001910060442_dp, &
0.0043594012767892145_dp, &
0.007474944527557248_dp, &
0.012628553354571935_dp, &
0.021045689348535686_dp, &
0.03463175327773065_dp, &
0.05632151430959264_dp, &
0.09059550579868533_dp, &
0.1442379871996754_dp, &
0.2274428250359289_dp, &
0.35541953739182136_dp, &
0.5507212109659467_dp, &
0.8466317459799636_dp, &
1.2921752687688919_dp, &
1.9598406408284006_dp, &
2.9586009577769_dp, &
4.459717588193107_dp, &
6.762921787948558_dp, &
10.555632767123296_dp, &
4.180340751600858e-06_dp, &
1.099130280826382e-05_dp, &
2.1822771688601605e-05_dp, &
4.1621175974422805e-05_dp, &
7.848461568768336e-05_dp, &
0.00014548449813908576_dp, &
0.00026422505480514413_dp, &
0.0004703011644186479_dp, &
0.0008215666081596167_dp, &
0.0014107803843716492_dp, &
0.002384836263299515_dp, &
0.0039737397649737246_dp, &
0.006533894710819361_dp, &
0.010612242095881175_dp, &
0.017040581368346817_dp, &
0.027073311442312262_dp, &
0.04258730399078407_dp, &
0.06637038094501468_dp, &
0.10253618298821561_dp, &
0.15712087871604702_dp, &
0.23894860519216926_dp, &
0.36091957483420634_dp, &
0.5420459855749467_dp, &
0.8110632955302136_dp, &
1.2140975460551333_dp, &
1.8351776848641823_dp, &
2.8698308488017825_dp, &
5.066349620899207_dp /)
ELSE IF(E_range < 700000.0_dp) THEN
aw(:) = (/ &
9.917609137977847e-07_dp, &
5.61571565668021e-06_dp, &
1.5905578931798044e-05_dp, &
3.668045117366575e-05_dp, &
7.779911321491552e-05_dp, &
0.0001574703243782781_dp, &
0.0003080352227521647_dp, &
0.0005857339949829306_dp, &
0.0010866136806263423_dp, &
0.0019719871993643633_dp, &
0.0035086606810778888_dp, &
0.006131766091278258_dp, &
0.010541856280946737_dp, &
0.017853461412239566_dp, &
0.029820286970471933_dp, &
0.049173663410606715_dp, &
0.08012714696416229_dp, &
0.1291233103698595_dp, &
0.20593171918618178_dp, &
0.32525471926745003_dp, &
0.5090693312937355_dp, &
0.7900513389150534_dp, &
1.216653152188614_dp, &
1.8609312710223198_dp, &
2.831687087034342_dp, &
4.300351181611193_dp, &
6.567047622277454_dp, &
10.318912411219488_dp, &
2.5862090723245625e-06_dp, &
6.95626109831625e-06_dp, &
1.4392071614748994e-05_dp, &
2.874925884264212e-05_dp, &
5.642970465568392e-05_dp, &
0.00010814505912642593_dp, &
0.0002021008939835649_dp, &
0.0003688922773676731_dp, &
0.0006591003182679095_dp, &
0.0011550942822986209_dp, &
0.0019891964272602783_dp, &
0.0033713558748037233_dp, &
0.0056309115068987_dp, &
0.00927906547653282_dp, &
0.015101574633332329_dp, &
0.024295245067122773_dp, &
0.03866754865528627_dp, &
0.060926803424727255_dp, &
0.09510216196058996_dp, &
0.14715087305959354_dp, &
0.22584200604915872_dp, &
0.34407208794479094_dp, &
0.5209358458550991_dp, &
0.7853686557167923_dp, &
1.1838124863793347_dp, &
1.8006332661854303_dp, &
2.8313699702065134_dp, &
5.0223450191392205_dp /)
ELSE IF(E_range < 1200000.0_dp) THEN
aw(:) = (/ &
6.243870742444705e-07_dp, &
3.5837547513507616e-06_dp, &
1.0424530634626105e-05_dp, &
2.4910967538724005e-05_dp, &
5.48216536959138e-05_dp, &
0.0001148065303999889_dp, &
0.000231490488679646_dp, &
0.00045215860524583124_dp, &
0.0008590837338469464_dp, &
0.0015927682432411259_dp, &
0.0028891299956185297_dp, &
0.005138253059658019_dp, &
0.008976135083193377_dp, &
0.015426433468499803_dp, &
0.026117293389759495_dp, &
0.043609983170991296_dp, &
0.0718927383781267_dp, &
0.11711702885780055_dp, &
0.18868749571084845_dp, &
0.3008660121451099_dp, &
0.47512403493520045_dp, &
0.7435972663838382_dp, &
1.154222692676153_dp, &
1.7786576652841575_dp, &
2.7255320837393326_dp, &
4.166412460588324_dp, &
6.401776404663582_dp, &
10.11856522638141_dp, &
1.6329106731735064e-06_dp, &
4.511304110662316e-06_dp, &
9.770509823422929e-06_dp, &
2.0454024246766166e-05_dp, &
4.172679855059072e-05_dp, &
8.25404415798886e-05_dp, &
0.00015847711044459983_dp, &
0.00029619555479324744_dp, &
0.0005404696680087612_dp, &
0.000965266351597196_dp, &
0.0016909815488903273_dp, &
0.0029109600865339395_dp, &
0.004931885217935304_dp, &
0.008234716202420396_dp, &
0.01356583086594651_dp, &
0.022072242618732946_dp, &
0.035500720158786525_dp, &
0.05648910800327906_dp, &
0.08899040639475064_dp, &
0.13888891315246274_dp, &
0.21489980129365738_dp, &
0.32990928766337146_dp, &
0.5030773832314516_dp, &
0.7635106326257592_dp, &
1.157931116384902_dp, &
1.771014718506682_dp, &
2.798340277794168_dp, &
4.984574132350961_dp /)
ELSE IF(E_range < 2000000.0_dp) THEN
aw(:) = (/ &
4.0370316213559476e-07_dp, &
2.3541711422913717e-06_dp, &
7.057573612279361e-06_dp, &
1.751162431119145e-05_dp, &
3.998197703771905e-05_dp, &
8.651758104084795e-05_dp, &
0.00017950101241220228_dp, &
0.00035945947508640393_dp, &
0.0006981140824218166_dp, &
0.001319810424458466_dp, &
0.0024361923309759396_dp, &
0.004401544013890147_dp, &
0.007800025204094134_dp, &
0.013581611201940205_dp, &
0.023271704375636543_dp, &
0.03929119474223725_dp, &
0.06544075124204898_dp, &
0.10762751156023793_dp, &
0.17494730054790258_dp, &
0.2812860986936203_dp, &
0.44768004577537496_dp, &
0.705795087398322_dp, &
1.1031142167130359_dp, &
1.7109362858917831_dp, &
2.6377277178978358_dp, &
4.055162275207615_dp, &
6.264028054745519_dp, &
9.951137323627291_dp, &
1.0593262130160056e-06_dp, &
3.0185707613834373e-06_dp, &
6.8658831794879865e-06_dp, &
1.5046778554922312e-05_dp, &
3.182893187575826e-05_dp, &
6.484710011709137e-05_dp, &
0.00012766407757989652_dp, &
0.000243860437909589_dp, &
0.00045361234593234906_dp, &
0.0008241579727338823_dp, &
0.001466238321223863_dp, &
0.002559617356616681_dp, &
0.004392266266516603_dp, &
0.0074199303035661236_dp, &
0.012355848225194005_dp, &
0.02030472924789548_dp, &
0.0329612503015901_dp, &
0.05290210344877982_dp, &
0.08401329609895744_dp, &
0.13211363792116823_dp, &
0.20586780977743016_dp, &
0.3181478896142563_dp, &
0.48816468770340044_dp, &
0.7451687687530838_dp, &
1.1361253599700833_dp, &
1.7459877749472628_dp, &
2.7703909434425067_dp, &
4.952626774509824_dp /)
ELSE IF(E_range < 4500000.0_dp) THEN
aw(:) = (/ &
2.0330900665342417e-07_dp, &
1.2259721995909777e-06_dp, &
3.9000493864165235e-06_dp, &
1.0346548674455454e-05_dp, &
2.5104573254227977e-05_dp, &
5.722321695341548e-05_dp, &
0.0001240915794841191_dp, &
0.00025812480274594484_dp, &
0.0005181574776824829_dp, &
0.0010084979886612006_dp, &
0.0019102647242876139_dp, &
0.003532146072512572_dp, &
0.0063914997207670885_dp, &
0.011342304908756772_dp, &
0.019774721642387954_dp, &
0.03392304253932974_dp, &
0.05733626485492922_dp, &
0.0955906794116426_dp, &
0.15736013123361753_dp, &
0.2560121769552906_dp, &
0.4119766499624099_dp, &
0.6562580036536524_dp, &
1.0356912360288544_dp, &
1.6210528978513052_dp, &
2.5205569278227222_dp, &
3.906010837038453_dp, &
6.07864411392436_dp, &
9.725147433824423_dp, &
5.372960626792664e-07_dp, &
1.6313777234005112e-06_dp, &
4.053582188809558e-06_dp, &
9.562645330678188e-06_dp, &
2.140029237042944e-05_dp, &
4.5628018410594356e-05_dp, &
9.332815320834305e-05_dp, &
0.00018423715674012751_dp, &
0.0003527071740962291_dp, &
0.0006573395283833412_dp, &
0.0011963286460791483_dp, &
0.0021315852595548933_dp, &
0.003726211411404738_dp, &
0.006402085885163769_dp, &
0.010827471574921069_dp, &
0.018049062427338197_dp, &
0.029689317308315376_dp, &
0.048239120452325215_dp, &
0.0774890829270111_dp, &
0.12316278045588097_dp, &
0.1938483490800667_dp, &
0.30239008881049667_dp, &
0.4680611192958136_dp, &
0.7203071572907666_dp, &
1.1064348627930074_dp, &
1.7118002420359797_dp, &
2.732149945273555_dp, &
4.908937928661797_dp /)
ELSE IF(E_range < 10000000.0_dp) THEN
aw(:) = (/ &
1.0453117858944355e-07_dp, &
6.608617412970565e-07_dp, &
2.263239779628554e-06_dp, &
6.456580799232289e-06_dp, &
1.6648773232985477e-05_dp, &
3.988538010811069e-05_dp, &
9.013431830037881e-05_dp, &
0.00019412595183632066_dp, &
0.00040148594866505004_dp, &
0.000801953768960446_dp, &
0.0015541047607406542_dp, &
0.002932474750123789_dp, &
0.0054037297272779_dp, &
0.0097481213186169_dp, &
0.01725076571642057_dp, &
0.029999477630049398_dp, &
0.05134367891607686_dp, &
0.08659470336954253_dp, &
0.14408515309896788_dp, &
0.2367592019652705_dp, &
0.3845463783654337_dp, &
0.6178989765753549_dp, &
0.9831034515051852_dp, &
1.5504849076462728_dp, &
2.4280266299570483_dp, &
3.7876319869184405_dp, &
5.930901453777237_dp, &
9.544473553973464_dp, &
2.7911359387760054e-07_dp, &
9.225699169205882e-07_dp, &
2.525970869540034e-06_dp, &
6.400448408750558e-06_dp, &
1.5105857774384374e-05_dp, &
3.360200325419909e-05_dp, &
7.118848251951222e-05_dp, &
0.00014478894137667428_dp, &
0.0002844248664965114_dp, &
0.0005421792095475284_dp, &
0.001006646872918385_dp, &
0.0018259034247645597_dp, &
0.003243545770723624_dp, &
0.0056545862514208655_dp, &
0.009691209964763319_dp, &
0.016353052767086353_dp, &
0.027203315594152448_dp, &
0.044661579803193296_dp, &
0.07243798977464_dp, &
0.11617410461859103_dp, &
0.18438950337219645_dp, &
0.2898984327270837_dp, &
0.45201785956261764_dp, &
0.7003495143727266_dp, &
1.0824839695902242_dp, &
1.6841236033891487_dp, &
2.701136203287264_dp, &
4.87352184807004_dp /)
ELSE IF(E_range < 50000000.0_dp) THEN
aw(:) = (/ &
2.883608355013751e-08_dp, &
2.1533521290478e-07_dp, &
8.922309745530507e-07_dp, &
2.9625449220607536e-06_dp, &
8.567451943237967e-06_dp, &
2.243432347635065e-05_dp, &
5.4443793898855104e-05_dp, &
0.00012434900771925045_dp, &
0.0002702071667839422_dp, &
0.0005630773934606482_dp, &
0.001132097168889575_dp, &
0.0022064607569192898_dp, &
0.004184511435918261_dp, &
0.007745737975074526_dp, &
0.014029823575166972_dp, &
0.024919304936911416_dp, &
0.04348066531140458_dp, &
0.07464540686704813_dp, &
0.12625156513534663_dp, &
0.2106228965724005_dp, &
0.34694729093744464_dp, &
0.5648480693819434_dp, &
0.9097754411882123_dp, &
1.4513508415022998_dp, &
2.297177622667234_dp, &
3.619275978541861_dp, &
5.719806682040045_dp, &
9.28540894538464_dp, &
8.026935386394373e-08_dp, &
3.431660125547769e-07_dp, &
1.1507249938591508e-06_dp, &
3.3255150480524145e-06_dp, &
8.631544858985621e-06_dp, &
2.0673889190163757e-05_dp, &
4.650257994599218e-05_dp, &
9.941521329490123e-05_dp, &
0.0002037413884778986_dp, &
0.00040284592225487954_dp, &
0.0007722776557802875_dp, &
0.0014410293829550617_dp, &
0.0026254207291881958_dp, &
0.0046823958218659345_dp, &
0.008192369379080482_dp, &
0.014086602667824237_dp, &
0.023841088169455484_dp, &
0.039769020002557186_dp, &
0.06545851893082533_dp, &
0.1064240118543133_dp, &
0.1710745747479101_dp, &
0.2721680114354786_dp, &
0.42907347105223687_dp, &
0.671615155080402_dp, &
1.0478076752969256_dp, &
1.6438903234547084_dp, &
2.6559588426369642_dp, &
4.821957836337722_dp /)
ELSE
aw(:) = (/ &
7.2246193648183315e-09_dp, &
7.947805205434051e-08_dp, &
4.2064242178143404e-07_dp, &
1.6235298675242196e-06_dp, &
5.1991703168705275e-06_dp, &
1.4668546651111404e-05_dp, &
3.7705424076918065e-05_dp, &
9.01796708136201e-05_dp, &
0.0002035415234868656_dp, &
0.0004379390624442907_dp, &
0.0009049701100068852_dp, &
0.0018063291005268592_dp, &
0.0034982504743316482_dp, &
0.006597202794116675_dp, &
0.012150671665557649_dp, &
0.021909404675932484_dp, &
0.03875606413936719_dp, &
0.06737255773930659_dp, &
0.11526829896699824_dp, &
0.1943501837935314_dp, &
0.3233020082512726_dp, &
0.5311764928093146_dp, &
0.8628395371729328_dp, &
1.3874103242942566_dp, &
2.212207252365476_dp, &
3.5093118591573305_dp, &
5.581272308767191_dp, &
9.114780312141413_dp, &
2.3095697901478273e-08_dp, &
1.5174284298730647e-07_dp, &
6.211817019789001e-07_dp, &
2.0193446034719265e-06_dp, &
5.684986670359784e-06_dp, &
1.4468629362219993e-05_dp, &
3.413251437719442e-05_dp, &
7.584802072328596e-05_dp, &
0.00016053557840876553_dp, &
0.0003262378591667722_dp, &
0.0006404016188164239_dp, &
0.0012199823490146536_dp, &
0.0022638430447923645_dp, &
0.0041042325427612_dp, &
0.007287543243404729_dp, &
0.012699528275383098_dp, &
0.021757390677572788_dp, &
0.03670161910421734_dp, &
0.061035614898498174_dp, &
0.10018369236083766_dp, &
0.16247380492926894_dp, &
0.26061721476056016_dp, &
0.41400956019249524_dp, &
0.6526199011491635_dp, &
1.0247530046294373_dp, &
1.617029280212683_dp, &
2.625732716141634_dp, &
4.78747600908045_dp /)
END IF
END SUBROUTINE

SUBROUTINE get_coeff_30(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio

E_ratio = 1.0_dp

IF(E_range < 2906.0_dp) THEN
E_ratio = 2906.0_dp/E_range
aw(:) = (/ &
0.00010719782131833758_dp, &
0.0005720967766785331_dp, &
0.001439669174001423_dp, &
0.002769849121659816_dp, &
0.004664000049778321_dp, &
0.007283241780546586_dp, &
0.010873000028774003_dp, &
0.015792629918665993_dp, &
0.02255038608921401_dp, &
0.031847782299297675_dp, &
0.04463988143266245_dp, &
0.06221821662158963_dp, &
0.08632310206096853_dp, &
0.11929354905847961_dp, &
0.16426574047985984_dp, &
0.2254347024266385_dp, &
0.30839852293765396_dp, &
0.42061064113062563_dp, &
0.5719741696129382_dp, &
0.7756243966360713_dp, &
1.0489645297600516_dp, &
1.415051852321483_dp, &
1.9044914330966105_dp, &
2.558117126903484_dp, &
3.4310117702125913_dp, &
4.5990752550643395_dp, &
6.171106184489477_dp, &
8.314824292369494_dp, &
11.326871889957468_dp, &
15.91048796379593_dp, &
0.0002759212187476099_dp, &
0.0006589920695996667_dp, &
0.001085937188673809_dp, &
0.0015910122106733914_dp, &
0.0022237295869651035_dp, &
0.003055163774041633_dp, &
0.004183672861951183_dp, &
0.005739573533927359_dp, &
0.007891868776677179_dp, &
0.010860734164881057_dp, &
0.014936866838559482_dp, &
0.020507383044696917_dp, &
0.02808894961545452_dp, &
0.038370325842103176_dp, &
0.052267550952526234_dp, &
0.07099592893538823_dp, &
0.09616415398590554_dp, &
0.1298977096485147_dp, &
0.17500152116499199_dp, &
0.235176717264172_dp, &
0.3153154828400979_dp, &
0.4219165107665381_dp, &
0.5637038218742615_dp, &
0.7526245591453196_dp, &
1.005629375577748_dp, &
1.3482480324070423_dp, &
1.8228149068796895_dp, &
2.510939848106289_dp, &
3.6132439344998795_dp, &
5.9228966909439205_dp /)
ELSE IF(E_range < 3236.0_dp) THEN
aw(:) = (/ &
9.76676051737233e-05_dp, &
0.000521439839796676_dp, &
0.0013131673080453233_dp, &
0.0025293903535565823_dp, &
0.0042661250059111795_dp, &
0.0066765801036800295_dp, &
0.009994783165295324_dp, &
0.014564000812497294_dp, &
0.020870492553693578_dp, &
0.02958685422269625_dp, &
0.04163138568046643_dp, &
0.058249942771901375_dp, &
0.08112692839971898_dp, &
0.11253369957795559_dp, &
0.15552549205786806_dp, &
0.2142017013001525_dp, &
0.29404913356355256_dp, &
0.4023941070276559_dp, &
0.5489978387767454_dp, &
0.7468418586806641_dp, &
1.0131691905482183_dp, &
1.3708791149003934_dp, &
1.8504330035916199_dp, &
2.4925505474173923_dp, &
3.3522463174174093_dp, &
4.505410618668513_dp, &
6.060890691327169_dp, &
8.186478838726869_dp, &
11.178786177384513_dp, &
15.74038405881577_dp, &
0.00025141352623809566_dp, &
0.0006009323005964832_dp, &
0.0009917485369505094_dp, &
0.0014564055840763494_dp, &
0.002042200574311218_dp, &
0.0028172406584227633_dp, &
0.003875825130607345_dp, &
0.005342947480157929_dp, &
0.00738116098048549_dp, &
0.010203263291670922_dp, &
0.014091611662835478_dp, &
0.019423788992341998_dp, &
0.026705494627716198_dp, &
0.03661296179988794_dp, &
0.05004819563708869_dp, &
0.06821123613467164_dp, &
0.09269486715309114_dp, &
0.12560901090026905_dp, &
0.16974491105273157_dp, &
0.22879407342136257_dp, &
0.30764598759629597_dp, &
0.4128070021079015_dp, &
0.5530226357541534_dp, &
0.7402778334602197_dp, &
0.9915767219176012_dp, &
1.332510891997349_dp, &
1.8054641559432192_dp, &
2.49203146209568_dp, &
3.5925935804922893_dp, &
5.89893592679766_dp /)
ELSE IF(E_range < 3810.0_dp) THEN
aw(:) = (/ &
8.47819458078811e-05_dp, &
0.00045292247232938765_dp, &
0.0011419488590203569_dp, &
0.0022036034544386665_dp, &
0.0037263081029974285_dp, &
0.005851972143509195_dp, &
0.008798246612654406_dp, &
0.012885194918659285_dp, &
0.018567320347579703_dp, &
0.026475458982841895_dp, &
0.03747475660944132_dp, &
0.05274484615883288_dp, &
0.07388876569982615_dp, &
0.103079000554941_dp, &
0.14325196624460337_dp, &
0.19836605747691327_dp, &
0.27374326191660053_dp, &
0.3765207493043709_dp, &
0.5162475872690495_dp, &
0.7056742360429414_dp, &
0.9618016216492792_dp, &
1.3072886416997354_dp, &
1.772376230882795_dp, &
2.39760681930879_dp, &
3.237886311434963_dp, &
4.369085529708981_dp, &
5.900122672480505_dp, &
7.998902446784312_dp, &
10.962002095892634_dp, &
15.491026703179935_dp, &
0.00021827426323305446_dp, &
0.0005223668404015434_dp, &
0.0008641224881204865_dp, &
0.0012736436118541284_dp, &
0.001795042273807457_dp, &
0.0024921178119473955_dp, &
0.00345325315447348_dp, &
0.004795711504484412_dp, &
0.006672726568911133_dp, &
0.009286420295701826_dp, &
0.012906984705977642_dp, &
0.017897959998054307_dp, &
0.02474875419566456_dp, &
0.0341168643096291_dp, &
0.046883167956700916_dp, &
0.06422456826044022_dp, &
0.08770953467687731_dp, &
0.11942394754852313_dp, &
0.16213754631449626_dp, &
0.21952613567742027_dp, &
0.2964736644552558_dp, &
0.39949638395652903_dp, &
0.5373706722452656_dp, &
0.7221373998987525_dp, &
0.9708812677144617_dp, &
1.309288667276434_dp, &
1.7798217912856973_dp, &
2.464060270734846_dp, &
3.5620371023687465_dp, &
5.8635023001897295_dp /)
ELSE IF(E_range < 4405.0_dp) THEN
aw(:) = (/ &
7.475415175359787e-05_dp, &
0.0003995782940364109_dp, &
0.0010085401523130205_dp, &
0.0019494530926649103_dp, &
0.003304487334168418_dp, &
0.005206188391856581_dp, &
0.007858547072014874_dp, &
0.011562196417307553_dp, &
0.01674502130839572_dp, &
0.02400284151062447_dp, &
0.03415615402045192_dp, &
0.04832877588615427_dp, &
0.06805486346691048_dp, &
0.09542279350363236_dp, &
0.13326739942725138_dp, &
0.18542592486136436_dp, &
0.25707802147355546_dp, &
0.3551966628680496_dp, &
0.48914575629830964_dp, &
0.6714729217979144_dp, &
0.9189652081581203_dp, &
1.2540675821639358_dp, &
1.706822983665265_dp, &
2.3176124905704714_dp, &
3.1412403932005213_dp, &
4.253555048851896_dp, &
5.763534315696107_dp, &
7.839184077942836_dp, &
10.777063260156016_dp, &
15.27796168145106_dp, &
0.00019248226458298544_dp, &
0.0004611671410305988_dp, &
0.0007645470771240003_dp, &
0.0011307064674945877_dp, &
0.001601097912084883_dp, &
0.002235885143791768_dp, &
0.003118439955104845_dp, &
0.00435950977219228_dp, &
0.006104510985566756_dp, &
0.008546583020796826_dp, &
0.011945581030119579_dp, &
0.016652984461283625_dp, &
0.02314408570724871_dp, &
0.03206002850298206_dp, &
0.044263120398424796_dp, &
0.06090977430142275_dp, &
0.08354673398267808_dp, &
0.11423814381228614_dp, &
0.15573394602932067_dp, &
0.2116949419787466_dp, &
0.28699873776286255_dp, &
0.38816866941911854_dp, &
0.5240067210209787_dp, &
0.7066019007188943_dp, &
0.9531097377337013_dp, &
1.2893015561438181_dp, &
1.7577121859945424_dp, &
2.439914418694296_dp, &
3.5356481389173076_dp, &
5.832916244207001_dp /)
ELSE IF(E_range < 5400.0_dp) THEN
aw(:) = (/ &
6.263685476987877e-05_dp, &
0.0003350864107617077_dp, &
0.0008471021406263553_dp, &
0.0016414732496533968_dp, &
0.0027923267352349805_dp, &
0.004420069092229887_dp, &
0.006710856247430801_dp, &
0.009939856329159545_dp, &
0.014500042255919464_dp, &
0.02094127936473627_dp, &
0.030025328730095918_dp, &
0.042802276529506474_dp, &
0.06071482338334902_dp, &
0.08573906273600417_dp, &
0.1205734666631485_dp, &
0.168891746293846_dp, &
0.23568035070007362_dp, &
0.3276881055094081_dp, &
0.4540246498488577_dp, &
0.6269572963379791_dp, &
0.8629754834557104_dp, &
1.1842241235652444_dp, &
1.6204657246854692_dp, &
2.2118496417292866_dp, &
3.0130302596964915_dp, &
4.099816741357715_dp, &
5.581265514191353_dp, &
7.625526452879397_dp, &
10.52915117245689_dp, &
14.99185357007055_dp, &
0.00016131234676709978_dp, &
0.0003871323803883455_dp, &
0.0006438639425224474_dp, &
0.0009569827684469891_dp, &
0.0013644640382973136_dp, &
0.0019216628235682472_dp, &
0.002705305661659533_dp, &
0.003817552287654078_dp, &
0.005393559577915718_dp, &
0.007614615789728071_dp, &
0.0107267793784036_dp, &
0.01506524759897237_dp, &
0.02108608211163133_dp, &
0.029407976403737046_dp, &
0.040867556214958384_dp, &
0.05659266977268607_dp, &
0.07809948510740791_dp, &
0.10742117687653448_dp, &
0.1472789549319935_dp, &
0.20131104861420793_dp, &
0.27438405593923587_dp, &
0.3730286538370253_dp, &
0.50607994537212_dp, &
0.6856919919465865_dp, &
0.9291183488974969_dp, &
1.2622503775094018_dp, &
1.7277292939880777_dp, &
2.407128771305674_dp, &
3.4998023815666013_dp, &
5.791398263742947_dp /)
ELSE IF(E_range < 6800.0_dp) THEN
aw(:) = (/ &
5.1260796624319686e-05_dp, &
0.00027449998212479283_dp, &
0.0006952559402364613_dp, &
0.0013512581469068891_dp, &
0.0023084695955214183_dp, &
0.0036748571064951344_dp, &
0.005618146035815099_dp, &
0.008387084892047715_dp, &
0.012338401496919359_dp, &
0.01797419704821914_dp, &
0.025994966614220735_dp, &
0.03737346294738748_dp, &
0.05345583818085436_dp, &
0.0760988446106003_dp, &
0.10785503448538585_dp, &
0.15222193676863333_dp, &
0.2139764359272559_dp, &
0.29962253812775563_dp, &
0.4179901608236006_dp, &
0.5810358869798488_dp, &
0.8049164836351872_dp, &
1.11143825413834_dp, &
1.5300435290071273_dp, &
2.1006127473331047_dp, &
2.8776209315806276_dp, &
3.93682226153259_dp, &
5.387354342884511_dp, &
7.397531351751853_dp, &
10.263919356557402_dp, &
14.685102820756857_dp, &
0.00013204473011131683_dp, &
0.0003175248324260899_dp, &
0.0005301211290638112_dp, &
0.000792646090309341_dp, &
0.0011394726365449568_dp, &
0.0016209075495792747_dp, &
0.0023066917116162182_dp, &
0.0032900298122596973_dp, &
0.004695415042196518_dp, &
0.0066917090579473775_dp, &
0.009510290459434545_dp, &
0.013468802311717997_dp, &
0.019002357781084556_dp, &
0.02670497738030662_dp, &
0.037384817040713114_dp, &
0.05213776033026327_dp, &
0.072445378072951_dp, &
0.10030529745315488_dp, &
0.13840505606579703_dp, &
0.19035539860923475_dp, &
0.2610076849435031_dp, &
0.3568973694825793_dp, &
0.48689310015691706_dp, &
0.6632189222745253_dp, &
0.9032371034193576_dp, &
1.2329755271865048_dp, &
1.695201274879912_dp, &
2.3715025399546223_dp, &
3.4608285109987857_dp, &
5.746290251951392_dp /)
ELSE IF(E_range < 8400.0_dp) THEN
aw(:) = (/ &
4.264827171158114e-05_dp, &
0.00022859921129368268_dp, &
0.0005800651655637003_dp, &
0.0011306618266057897_dp, &
0.001939659276932296_dp, &
0.0031047255791507274_dp, &
0.004778205945319036_dp, &
0.007186731964418128_dp, &
0.010656673525706605_dp, &
0.01565005735652159_dp, &
0.02281575734702256_dp, &
0.03306097534037487_dp, &
0.04764947470167637_dp, &
0.0683354702649706_dp, &
0.09754528626464402_dp, &
0.13862301432535307_dp, &
0.19616179200890152_dp, &
0.27644949046799433_dp, &
0.3880673249123362_dp, &
0.5426935306072133_dp, &
0.7561844283103211_dp, &
1.0500376608384598_dp, &
1.4534003880410857_dp, &
2.005901250325418_dp, &
2.7618425195091647_dp, &
3.796918364528742_dp, &
5.220333260426152_dp, &
7.200552079960992_dp, &
10.034172866297915_dp, &
14.41882124837872_dp, &
0.0001098833579123983_dp, &
0.0002647437855875365_dp, &
0.0004436472669265954_dp, &
0.0006672101026873376_dp, &
0.0009667895306084336_dp, &
0.0013884135997529708_dp, &
0.0019959070037335724_dp, &
0.0028749526284055213_dp, &
0.004141093118158889_dp, &
0.005952643864305207_dp, &
0.008528340123693306_dp, &
0.01217052007837369_dp, &
0.017295858304956554_dp, &
0.024476476682491477_dp, &
0.03449505017711748_dp, &
0.04841859097574648_dp, &
0.06769708781585487_dp, &
0.09429528482789237_dp, &
0.13086897935334732_dp, &
0.1810021324793278_dp, &
0.24952974948475865_dp, &
0.3429885311920854_dp, &
0.4702742651357602_dp, &
0.6436716531183656_dp, &
0.880640225748484_dp, &
1.2073331130934455_dp, &
1.6666371474517014_dp, &
2.3401654531537157_dp, &
3.4265253480645885_dp, &
5.706615197039722_dp /)
ELSE IF(E_range < 10000.0_dp) THEN
aw(:) = (/ &
3.6638445936143444e-05_dp, &
0.00019654898584533858_dp, &
0.0004995365456010348_dp, &
0.0009761629469610179_dp, &
0.001680692380508944_dp, &
0.0027030234477272513_dp, &
0.004183824314963024_dp, &
0.006332898868769434_dp, &
0.009453508422633638_dp, &
0.013977114132271049_dp, &
0.020513053005428484_dp, &
0.02991804681839935_dp, &
0.043392038467264346_dp, &
0.06260934488416951_dp, &
0.089897358898568_dp, &
0.12847922117893293_dp, &
0.18280238941596594_dp, &
0.25898237331386087_dp, &
0.36540085871444006_dp, &
0.5135113514430796_dp, &
0.718925939367594_dp, &
1.0028894169830442_dp, &
1.3943048957706679_dp, &
1.9325902027814568_dp, &
2.6718997828791795_dp, &
3.6878719059499683_dp, &
5.089761008341865_dp, &
7.046156150294659_dp, &
9.853694451194254_dp, &
14.209263474504755_dp, &
9.441684068693651e-05_dp, &
0.00022786014262271372_dp, &
0.00038307351711976584_dp, &
0.0005790223101422384_dp, &
0.0008447661561374952_dp, &
0.0012230417856923205_dp, &
0.0017731334822002557_dp, &
0.002574990545123123_dp, &
0.003737321923840727_dp, &
0.005410302324390409_dp, &
0.007802786243657_dp, &
0.011205034742082835_dp, &
0.016019061597614896_dp, &
0.02279945700439318_dp, &
0.032308353668577205_dp, &
0.04558931015753348_dp, &
0.06406643328220389_dp, &
0.08967722561429306_dp, &
0.1250507942222136_dp, &
0.17374800817672645_dp, &
0.24058882427995695_dp, &
0.3321088613980677_dp, &
0.4572237934265981_dp, &
0.6282658627647734_dp, &
0.8627730246821379_dp, &
1.187001575527606_dp, &
1.6439399241381616_dp, &
2.3152296444802993_dp, &
3.3992160109205325_dp, &
5.6750508010285134_dp /)
ELSE IF(E_range < 12000.0_dp) THEN
aw(:) = (/ &
3.125352048295285e-05_dp, &
0.00016781351172929327_dp, &
0.00042725280780375096_dp, &
0.0008372368297287864_dp, &
0.0014472482996344947_dp, &
0.0023397061489668085_dp, &
0.00364397706459525_dp, &
0.005553554215846271_dp, &
0.008349288769114196_dp, &
0.012432933353972381_dp, &
0.018375235671074245_dp, &
0.0269833920064158_dp, &
0.03939439104594155_dp, &
0.057203310896803944_dp, &
0.08263890933961444_dp, &
0.118803108858221_dp, &
0.16999660677619183_dp, &
0.24216036024367976_dp, &
0.3434728927615001_dp, &
0.4851575786991668_dp, &
0.6825748514741242_dp, &
0.956707159367906_dp, &
1.3362022902301018_dp, &
1.8602550523674573_dp, &
2.582860494442485_dp, &
3.5795923551771076_dp, &
4.959752013207221_dp, &
6.892057403111867_dp, &
9.673196873736643_dp, &
13.999332658018544_dp, &
8.055654719138439e-05_dp, &
0.0001947658525515697_dp, &
0.00032859686532794565_dp, &
0.0004994306825817649_dp, &
0.0007340947583837473_dp, &
0.0010721035282296792_dp, &
0.001568310846736925_dp, &
0.002297102568888374_dp, &
0.0033605346972129828_dp, &
0.00490077013024238_dp, &
0.0071168327918336724_dp, &
0.010286868707377024_dp, &
0.01479808821232481_dp, &
0.021187271898748812_dp, &
0.03019555503700218_dp, &
0.042842362400043486_dp, &
0.06052496436664667_dp, &
0.08515234916902564_dp, &
0.11932532276255173_dp, &
0.16657974696296404_dp, &
0.23171843775477596_dp, &
0.3212739304330216_dp, &
0.4441803236632917_dp, &
0.6128171239515162_dp, &
0.8448023886795951_dp, &
1.1664998500858539_dp, &
1.6210063652060114_dp, &
2.2900001311037035_dp, &
3.371570229987856_dp, &
5.643112922874634_dp /)
ELSE IF(E_range < 15000.0_dp) THEN
aw(:) = (/ &
2.572460674430881e-05_dp, &
0.00013828809639808773_dp, &
0.0003528803185321503_dp, &
0.0006939963784805861_dp, &
0.0012058419667301263_dp, &
0.00196250811153749_dp, &
0.003080710120197529_dp, &
0.004735659001322314_dp, &
0.007183065998906602_dp, &
0.010791253848864667_dp, &
0.016087325873200328_dp, &
0.02382213056319325_dp, &
0.03506059129972535_dp, &
0.051306532736392435_dp, &
0.07467445268948096_dp, &
0.10812503757927873_dp, &
0.15578699671299615_dp, &
0.2233955307214805_dp, &
0.31888824542436384_dp, &
0.45321390622010965_dp, &
0.6414306299530386_dp, &
0.9042033223791973_dp, &
1.2698679478183343_dp, &
1.7773436933847355_dp, &
2.4804246798636003_dp, &
3.4545978455659587_dp, &
4.80921488681201_dp, &
6.713149815396915_dp, &
9.463167209838238_dp, &
13.75460415331797_dp, &
6.632324623333035e-05_dp, &
0.00016073115506706547_dp, &
0.0002724186286744679_dp, &
0.00041700914328903177_dp, &
0.0006188187442380069_dp, &
0.000913709771883424_dp, &
0.0013515392555132155_dp, &
0.002000450746020119_dp, &
0.0029549917404197673_dp, &
0.004348167160940399_dp, &
0.006367634264735562_dp, &
0.009277411142785488_dp, &
0.013447324345434384_dp, &
0.01939310256171694_dp, &
0.027830892296757826_dp, &
0.03975118351588367_dp, &
0.056518799068157294_dp, &
0.08000790205760616_dp, &
0.11278427119956595_dp, &
0.15835215912014322_dp, &
0.2214916578911036_dp, &
0.3087289982044895_dp, &
0.42901765162073763_dp, &
0.5947915793660636_dp, &
0.8237640916060611_dp, &
1.142429482918768_dp, &
1.5940199744065582_dp, &
2.2602675888442345_dp, &
3.3389719364130626_dp, &
5.605477212236766_dp /)
ELSE IF(E_range < 20000.0_dp) THEN
aw(:) = (/ &
2.0010353507234397e-05_dp, &
0.0001077433422458143_dp, &
0.00027580023006758454_dp, &
0.0005451258515021734_dp, &
0.0009539529098224693_dp, &
0.0015668441726082855_dp, &
0.0024859562396115408_dp, &
0.003865438814677545_dp, &
0.005931997279200694_dp, &
0.009015206545406032_dp, &
0.013591237471033118_dp, &
0.020344680291656994_dp, &
0.03025506904355767_dp, &
0.04471729454363523_dp, &
0.0657084625235162_dp, &
0.09601821543882148_dp, &
0.13956549956219141_dp, &
0.2018327748807436_dp, &
0.2904595479070126_dp, &
0.41605219017177664_dp, &
0.5932887745998796_dp, &
0.8424313672051914_dp, &
1.1914160744119529_dp, &
1.678803642976038_dp, &
2.3581208348346263_dp, &
3.3047307785320212_dp, &
4.6280385004171025_dp, &
6.497115386386612_dp, &
9.208843054933286_dp, &
13.457588010635305_dp, &
5.1609539489171625e-05_dp, &
0.00012547938250585046_dp, &
0.00021401951207906163_dp, &
0.0003308508348857692_dp, &
0.0004973777574979467_dp, &
0.0007451991801870302_dp, &
0.0011183825299403462_dp, &
0.0016778698024847744_dp, &
0.0025094583225573333_dp, &
0.0037353071207149595_dp, &
0.005529432603001048_dp, &
0.008138731296629963_dp, &
0.011911786975259289_dp, &
0.0173383968709882_dp, &
0.02510367274551181_dp, &
0.03616185208550007_dp, &
0.05183670042056887_dp, &
0.07395778747131394_dp, &
0.10504533244757883_dp, &
0.14856147943385148_dp, &
0.20925449434054086_dp, &
0.2936386766166653_dp, &
0.41068755123213285_dp, &
0.5728998438821491_dp, &
0.7981069410955339_dp, &
1.1129693992413612_dp, &
1.5608972768887142_dp, &
2.223705100788062_dp, &
3.298855550400102_dp, &
5.559194240832182_dp /)
ELSE IF(E_range < 28000.0_dp) THEN
aw(:) = (/ &
1.4913158562437485e-05_dp, &
8.046249863674396e-05_dp, &
0.00020679269836799646_dp, &
0.00041135606175085387_dp, &
0.000726429457279576_dp, &
0.001206954515218603_dp, &
0.0019403041016435166_dp, &
0.0030592205520702234_dp, &
0.004760844972277948_dp, &
0.007334997531745175_dp, &
0.011205132165430936_dp, &
0.016986616421301493_dp, &
0.025568952229820977_dp, &
0.03823114994605764_dp, &
0.056802896237992014_dp, &
0.08388874591064059_dp, &
0.12317874197383723_dp, &
0.17987720132408636_dp, &
0.26129273365588884_dp, &
0.3776482327688284_dp, &
0.5431920469707082_dp, &
0.7777259064649447_dp, &
1.1087233079248426_dp, &
1.5743236370665366_dp, &
2.2277313150621567_dp, &
3.144149304477041_dp, &
4.433028932012796_dp, &
6.263665410565493_dp, &
8.933099495687644_dp, &
13.134684015624826_dp, &
3.848089108885685e-05_dp, &
9.39454195069688e-05_dp, &
0.00016153027576924122_dp, &
0.0002528432683507074_dp, &
0.0003863006748200446_dp, &
0.000589103857593431_dp, &
0.0008994092337673156_dp, &
0.0013708168552981645_dp, &
0.0020800635852347526_dp, &
0.0031378613819817206_dp, &
0.004703620407911122_dp, &
0.007005705166384156_dp, &
0.010369489117192241_dp, &
0.015256172448531693_dp, &
0.02231631000951314_dp, &
0.03246335096484727_dp, &
0.04697433098325742_dp, &
0.06762737851039459_dp, &
0.09688924910721004_dp, &
0.13817140362709024_dp, &
0.19618184162864452_dp, &
0.2774160356134996_dp, &
0.39086437274546876_dp, &
0.5490937144681445_dp, &
0.7700664939132916_dp, &
1.0806338759241187_dp, &
1.5244170903885788_dp, &
2.183343539374224_dp, &
3.254529601017961_dp, &
5.508095437648589_dp /)
ELSE IF(E_range < 38000.0_dp) THEN
aw(:) = (/ &
1.1418872198921082e-05_dp, &
6.173413750253069e-05_dp, &
0.00015929220328335313_dp, &
0.00031889503854739523_dp, &
0.0005682342649483396_dp, &
0.0009547547779672074_dp, &
0.0015542608919948567_dp, &
0.0024827272879236525_dp, &
0.003914044995780404_dp, &
0.006106526724782397_dp, &
0.009441451620849054_dp, &
0.01447827421468276_dp, &
0.022033096235727034_dp, &
0.03328960553532699_dp, &
0.04995516756713029_dp, &
0.07447945346316431_dp, &
0.11035934551634245_dp, &
0.16256247169247554_dp, &
0.23811344985212937_dp, &
0.3469031346310937_dp, &
0.5028042883523934_dp, &
0.7252121362263632_dp, &
1.0411867545570728_dp, &
1.4884852069632906_dp, &
2.120012518561155_dp, &
3.0108136851550755_dp, &
4.2703677397197355_dp, &
6.068165404997884_dp, &
8.701408700808228_dp, &
12.8626331750244_dp, &
2.947787098787719e-05_dp, &
7.225954134028555e-05_dp, &
0.0001252402865134724_dp, &
0.00019846473664424758_dp, &
0.00030798161384563543_dp, &
0.0004775038459907116_dp, &
0.0007405481153276563_dp, &
0.0011449321771295313_dp, &
0.0017601237436588445_dp, &
0.0026874815648166833_dp, &
0.004074324357547685_dp, &
0.0061335185303593494_dp, &
0.009170846038462778_dp, &
0.013623145317805922_dp, &
0.020111251332522298_dp, &
0.029513178507526473_dp, &
0.043064913823154624_dp, &
0.062498815014491525_dp, &
0.09023329634576005_dp, &
0.12963292767495377_dp, &
0.18536685282870904_dp, &
0.26390946533406945_dp, &
0.3742607782294386_dp, &
0.5290426355523791_dp, &
0.7463296721416969_dp, &
1.0531417380773798_dp, &
1.493293292834725_dp, &
2.1488272848196677_dp, &
3.216586556758033_dp, &
5.464388958502095_dp /)
ELSE IF(E_range < 50000.0_dp) THEN
aw(:) = (/ &
8.982546693806235e-06_dp, &
4.865818678309805e-05_dp, &
0.00012604148335636865_dp, &
0.00025390861901678955_dp, &
0.0004563995854196066_dp, &
0.0007750952045840026_dp, &
0.0012767265489376275_dp, &
0.002064101925548432_dp, &
0.0032927700122513484_dp, &
0.005195994959981289_dp, &
0.008121235122746024_dp, &
0.012582721029292343_dp, &
0.019336695756878207_dp, &
0.029488494866509965_dp, &
0.044644170175709534_dp, &
0.06712415308115761_dp, &
0.10026296333569998_dp, &
0.14882782019535817_dp, &
0.2196010954609523_dp, &
0.3221882442730784_dp, &
0.47013658664261176_dp, &
0.6824859983381489_dp, &
0.9859314820504801_dp, &
1.417888163895422_dp, &
2.030988118342932_dp, &
2.900125529064684_dp, &
4.134793295219112_dp, &
5.904650857570968_dp, &
8.50705581410213_dp, &
12.633884453832914_dp, &
2.3198709666506916e-05_dp, &
5.709307522017352e-05_dp, &
9.972792415428492e-05_dp, &
0.0001599273928986242_dp, &
0.0002518605559821019_dp, &
0.00039647689416403176_dp, &
0.0006236485188629505_dp, &
0.0009766180462365724_dp, &
0.0015189982582493745_dp, &
0.0023445004429086033_dp, &
0.0035904524574293063_dp, &
0.005456804956339633_dp, &
0.008232881015796438_dp, &
0.012334887697009843_dp, &
0.01835828332102448_dp, &
0.027150564163334993_dp, &
0.03991202944665125_dp, &
0.05833481297012722_dp, &
0.08479428042701978_dp, &
0.12261246982859136_dp, &
0.17642213365813839_dp, &
0.2526759210340022_dp, &
0.3603782883628389_dp, &
0.5121951678672628_dp, &
0.7262964913958404_dp, &
1.0298496814196687_dp, &
1.4668430564566999_dp, &
2.119432331068027_dp, &
3.1842446525056274_dp, &
5.427158434113688_dp /)
ELSE IF(E_range < 64000.0_dp) THEN
aw(:) = (/ &
7.238349702296008e-06_dp, &
3.928458373506955e-05_dp, &
0.00010214563688496669_dp, &
0.00020702207021802802_dp, &
0.0003752580104799391_dp, &
0.0006437818300993232_dp, &
0.001072114382774453_dp, &
0.0017525850080720873_dp, &
0.002826070350614258_dp, &
0.004505651346498813_dp, &
0.007111320652225213_dp, &
0.011120302372900251_dp, &
0.017239502228764955_dp, &
0.026509227663294707_dp, &
0.04045088318380033_dp, &
0.06127621597555976_dp, &
0.09218233054747209_dp, &
0.13776574927878416_dp, &
0.2046011860485274_dp, &
0.3020478300836994_dp, &
0.443370229666864_dp, &
0.6472971651730427_dp, &
0.9402012658200933_dp, &
1.3591924682162404_dp, &
1.956655369037565_dp, &
2.8073424394429862_dp, &
4.020750844752612_dp, &
5.766686098521229_dp, &
8.342653009498491_dp, &
12.439991369032668_dp, &
1.8702020468622055e-05_dp, &
4.6203219782190167e-05_dp, &
8.131762180666392e-05_dp, &
0.00013190143534228278_dp, &
0.00021061299324796057_dp, &
0.0003361903026363539_dp, &
0.000535605997735673_dp, &
0.0008484276122519987_dp, &
0.0013334887847408908_dp, &
0.0020781819161196638_dp, &
0.003211507934451115_dp, &
0.004922574199422692_dp, &
0.007486789650255751_dp, &
0.011302797934300408_dp, &
0.016944298667165195_dp, &
0.02523242357196762_dp, &
0.03733639323117268_dp, &
0.05491299477721061_dp, &
0.08029934786534572_dp, &
0.11677913471088484_dp, &
0.16895145657274943_dp, &
0.24324744526263437_dp, &
0.3486724489732179_dp, &
0.49792796808610074_dp, &
0.7092652162431419_dp, &
1.009980880755789_dp, &
1.4442192775687972_dp, &
2.094243792411222_dp, &
3.1565101709731045_dp, &
5.395251764843703_dp /)
ELSE IF(E_range < 84000.0_dp) THEN
aw(:) = (/ &
5.706429844293361e-06_dp, &
3.1040364829701505e-05_dp, &
8.107328661660111e-05_dp, &
0.00016550326188419503_dp, &
0.0003029763316477536_dp, &
0.0005259009041296794_dp, &
0.0008867852843333584_dp, &
0.0014677450288040894_dp, &
0.002395273728414361_dp, &
0.003862517134876226_dp, &
0.006162129729070569_dp, &
0.00973422779264423_dp, &
0.015235895683403678_dp, &
0.023641332048275565_dp, &
0.03638532798309582_dp, &
0.055567718600602306_dp, &
0.08424323176407787_dp, &
0.12683043751734444_dp, &
0.18968622248619452_dp, &
0.2819098113026998_dp, &
0.4164652673968472_dp, &
0.6117484077902107_dp, &
0.893783996100621_dp, &
1.2993493017998998_dp, &
1.8805551654919748_dp, &
2.7119921346520504_dp, &
3.9031533013982855_dp, &
5.623997981838861_dp, &
8.172198452849088_dp, &
12.238558563265647_dp, &
1.4751349028333799e-05_dp, &
3.660909835718544e-05_dp, &
6.501211560655282e-05_dp, &
0.00010687536539993863_dp, &
0.00017337123549154558_dp, &
0.00028107588099189824_dp, &
0.000454137569267576_dp, &
0.0007285000756826098_dp, &
0.001158211238344495_dp, &
0.0018242672088928367_dp, &
0.0028471735084972696_dp, &
0.004404897269066725_dp, &
0.006758454939132606_dp, &
0.01028819929336629_dp, &
0.01554501319618277_dp, &
0.023322182186959424_dp, &
0.03475584758118297_dp, &
0.05146485198022031_dp, &
0.07574483876932264_dp, &
0.11083732616584702_dp, &
0.16130359908425193_dp, &
0.23354923720333218_dp, &
0.3365774511982118_dp, &
0.48312462993088157_dp, &
0.6915266416555069_dp, &
0.9892183849743513_dp, &
1.420514754404278_dp, &
2.0678032116600904_dp, &
3.127373080182709_dp, &
5.361746811812883_dp /)
ELSE IF(E_range < 110000.0_dp) THEN
aw(:) = (/ &
4.507925303711096e-06_dp, &
2.4580451494530444e-05_dp, &
6.451263281768541e-05_dp, &
0.00013272058871030806_dp, &
0.0002455197572529811_dp, &
0.0004313910776655633_dp, &
0.0007367445062259567_dp, &
0.001234786889942828_dp, &
0.002039383465959761_dp, &
0.00332603532508658_dp, &
0.005363004530984576_dp, &
0.008557044500582709_dp, &
0.01352012174944254_dp, &
0.021166153619467505_dp, &
0.03285042829603399_dp, &
0.05056939685207151_dp, &
0.07724543673788668_dp, &
0.11713067199112848_dp, &
0.17637697527453522_dp, &
0.263837336282913_dp, &
0.39218930488595755_dp, &
0.5795089364205666_dp, &
0.8514842120807159_dp, &
1.2445671334250805_dp, &
1.8105970655011243_dp, &
2.6239991397280553_dp, &
3.7942542214271437_dp, &
5.491467032685958_dp, &
8.013480504379215_dp, &
12.05062051565182_dp, &
1.165943020924494e-05_dp, &
2.907711562944157e-05_dp, &
5.2135381292650646e-05_dp, &
8.692931914146652e-05_dp, &
0.00014332433921579002_dp, &
0.00023600989072825049_dp, &
0.0003866717964299752_dp, &
0.0006280479343018203_dp, &
0.0010098870748762492_dp, &
0.001607379741603303_dp, &
0.00253325997236095_dp, &
0.003955231709080428_dp, &
0.006120960061358173_dp, &
0.009393712743103142_dp, &
0.014302912174317341_dp, &
0.02161546778139668_dp, &
0.03243595167457546_dp, &
0.04834668628642002_dp, &
0.07160299218511915_dp, &
0.10540485530435402_dp, &
0.15427553890661685_dp, &
0.22459367842315606_dp, &
0.3253574621758251_dp, &
0.46933370101682803_dp, &
0.6749372715340488_dp, &
0.9697356715189667_dp, &
1.3982111886845372_dp, &
2.0428791235793997_dp, &
3.099885407679866_dp, &
5.330156088439341_dp /)
ELSE IF(E_range < 160000.0_dp) THEN
aw(:) = (/ &
3.2492504706355543e-06_dp, &
1.7782935938647487e-05_dp, &
4.702094527107548e-05_dp, &
9.788765130060788e-05_dp, &
0.00018394800566011574_dp, &
0.0003290226690843259_dp, &
0.0005722787919864911_dp, &
0.0009762934532250238_dp, &
0.0016397346725346576_dp, &
0.002716663917457915_dp, &
0.004445425742834286_dp, &
0.007191491449190322_dp, &
0.011510534331810433_dp, &
0.01824066242763528_dp, &
0.02863642038498904_dp, &
0.044562287024877995_dp, &
0.06877046943314344_dp, &
0.10529755457345026_dp, &
0.16002804271088766_dp, &
0.24149146592791995_dp, &
0.361986165199179_dp, &
0.5391615439925698_dp, &
0.7982521139849742_dp, &
1.1752670593405194_dp, &
1.7216703901842698_dp, &
2.511651208522495_dp, &
3.654660979891775_dp, &
5.320994285506342_dp, &
7.808735565003434_dp, &
11.80762543615471_dp, &
8.410837332223653e-06_dp, &
2.113245841372446e-05_dp, &
3.845111058827689e-05_dp, &
6.548419810172488e-05_dp, &
0.00011052592304298797_dp, &
0.00018602054544749434_dp, &
0.0003107143750285354_dp, &
0.0005134448260073493_dp, &
0.0008386432327700857_dp, &
0.0013542399463937198_dp, &
0.002163165180703502_dp, &
0.0034200707057045_dp, &
0.0053555078104983965_dp, &
0.008310672199526598_dp, &
0.012787041135722151_dp, &
0.019516894799063707_dp, &
0.029562998187469335_dp, &
0.044458853860523784_dp, &
0.0664052872217014_dp, &
0.09854535532252624_dp, &
0.14534905940024856_dp, &
0.2131554574377401_dp, &
0.3109515151654752_dp, &
0.4515398978356489_dp, &
0.653437304234083_dp, &
0.9443877723109725_dp, &
1.3691021780052584_dp, &
2.010279336423766_dp, &
3.063897971847719_dp, &
5.288820184644417_dp /)
ELSE IF(E_range < 220000.0_dp) THEN
aw(:) = (/ &
2.460537860443357e-06_dp, &
1.3513517677623214e-05_dp, &
3.598515389781885e-05_dp, &
7.575261086818076e-05_dp, &
0.0001444216294004159_dp, &
0.00026247903310060725_dp, &
0.0004639039840208787_dp, &
0.0008036114726657992_dp, &
0.001369205915452374_dp, &
0.0022989755449385045_dp, &
0.0038090333768977967_dp, &
0.006233886266199457_dp, &
0.01008661073890614_dp, &
0.016147469119526157_dp, &
0.025593516684082565_dp, &
0.04018693519947543_dp, &
0.06254703089827915_dp, &
0.0965408099871576_dp, &
0.14784085162679408_dp, &
0.22471839016516207_dp, &
0.33916659817804756_dp, &
0.5084886871639984_dp, &
0.7575478060875892_dp, &
1.1219864437669036_dp, &
1.6529529420726783_dp, &
2.4244318079577103_dp, &
3.5458393070593472_dp, &
5.187620119758077_dp, &
7.648065546110419_dp, &
11.61648316806262_dp, &
6.374120402922641e-06_dp, &
1.6128252577743594e-05_dp, &
2.9754113670763418e-05_dp, &
5.166430893732005e-05_dp, &
8.901502625196717e-05_dp, &
0.00015264214630912148_dp, &
0.0002591687585053359_dp, &
0.00043455383137211705_dp, &
0.000719240735481554_dp, &
0.0011756556233951615_dp, &
0.0018992288700902785_dp, &
0.003034542038483531_dp, &
0.004798831948386736_dp, &
0.007515980097531093_dp, &
0.011665360331427038_dp, &
0.017951632069288924_dp, &
0.027403921176393425_dp, &
0.041516070343442624_dp, &
0.06244414300326091_dp, &
0.0932837802259788_dp, &
0.13845969920620585_dp, &
0.20427585444540086_dp, &
0.2997063234922189_dp, &
0.43757882986144675_dp, &
0.6364895482287001_dp, &
0.9243252918542686_dp, &
1.3459866886820373_dp, &
1.9843322156142236_dp, &
3.035224685946777_dp, &
5.255903429150008_dp /)
ELSE IF(E_range < 370000.0_dp) THEN
aw(:) = (/ &
1.5637740408395342e-06_dp, &
8.644581438556428e-06_dp, &
2.3326026659580105e-05_dp, &
5.012244143843333e-05_dp, &
9.804997776691128e-05_dp, &
0.000183175648888439_dp, &
0.0003325825937550191_dp, &
0.0005909069030179598_dp, &
0.0010307262850542825_dp, &
0.001768639145774748_dp, &
0.0029898323691103764_dp, &
0.004985271727949821_dp, &
0.008207526024267482_dp, &
0.013353895620644004_dp, &
0.0214892751837372_dp, &
0.034226463130117915_dp, &
0.053989018153179254_dp, &
0.08439206340451581_dp, &
0.1307907673641504_dp, &
0.20106621537104283_dp, &
0.3067466011130968_dp, &
0.46460274135109153_dp, &
0.698920779817791_dp, &
1.0447664262134877_dp, &
1.5527831264532448_dp, &
2.296617666855539_dp, &
3.3856119528633637_dp, &
4.990434105875191_dp, &
7.4097119078577_dp, &
11.33215729315795_dp, &
4.056801137252665e-06_dp, &
1.0400304864325584e-05_dp, &
1.968307964805301e-05_dp, &
3.537211676997081e-05_dp, &
6.309758174077909e-05_dp, &
0.00011156144898932593_dp, &
0.0001945228908869725_dp, &
0.0003339605532184735_dp, &
0.000564711446670201_dp, &
0.0009413729224733118_dp, &
0.0015485976328911168_dp, &
0.0025163527439127716_dp, &
0.004042366217657837_dp, &
0.006424908889055276_dp, &
0.010110351754549359_dp, &
0.01576172006986062_dp, &
0.02435693199744164_dp, &
0.03732881710193292_dp, &
0.05676374129839576_dp, &
0.08568239505343193_dp, &
0.12843634401072082_dp, &
0.19127041101280928_dp, &
0.28313227875730185_dp, &
0.41688110772833514_dp, &
0.6112296495300013_dp, &
0.8942832614684186_dp, &
1.3112417204467246_dp, &
1.9452277892500305_dp, &
2.9919599852080703_dp, &
5.206267501211972_dp /)
ELSE IF(E_range < 520000.0_dp) THEN
aw(:) = (/ &
1.162904999603526e-06_dp, &
6.460220717190712e-06_dp, &
1.7606607480068412e-05_dp, &
3.8410731277459235e-05_dp, &
7.652869912290393e-05_dp, &
0.00014570244755436322_dp, &
0.00026937164379552174_dp, &
0.0004866802067040396_dp, &
0.0008620687596415197_dp, &
0.001500239074857704_dp, &
0.0025692139442485516_dp, &
0.004335527408651135_dp, &
0.00721744537817494_dp, &
0.011864777494471953_dp, &
0.01927761941951377_dp, &
0.030981700048419165_dp, &
0.04928550606403292_dp, &
0.07765484694198317_dp, &
0.12125517110926898_dp, &
0.18773242857967326_dp, &
0.2883321960866146_dp, &
0.43949878526491704_dp, &
0.6651609536255874_dp, &
1.0000229662476519_dp, &
1.4944068851808134_dp, &
2.2217383426459816_dp, &
3.2913012165448343_dp, &
4.873895944019108_dp, &
7.268367413333011_dp, &
11.163103194460474_dp, &
3.0200745063435063e-06_dp, &
7.819224105572315e-06_dp, &
1.5081156733800399e-05_dp, &
2.776801766506454e-05_dp, &
5.0700455074630125e-05_dp, &
9.145615638049602e-05_dp, &
0.0001622504928743556_dp, &
0.00028286367265657205_dp, &
0.0004849887879932366_dp, &
0.0008187828750348566_dp, &
0.0013627221411797283_dp, &
0.002238314348285273_dp, &
0.003631885398877431_dp, &
0.005826597760169261_dp, &
0.009249172128608455_dp, &
0.014537611686674936_dp, &
0.0226387619956523_dp, &
0.03494804598497715_dp, &
0.053508604450163914_dp, &
0.0812939971568708_dp, &
0.12260888215529167_dp, &
0.18365876903663883_dp, &
0.27337120427563866_dp, &
0.4046203209300098_dp, &
0.5961868839329513_dp, &
0.8763095433545912_dp, &
1.2903758814089559_dp, &
1.9216820534456363_dp, &
2.9658783009061263_dp, &
5.176364379511127_dp /)
ELSE IF(E_range < 700000.0_dp) THEN
aw(:) = (/ &
8.982200709109713e-07_dp, &
5.013746520038623e-06_dp, &
1.3797556086375253e-05_dp, &
3.053886905121258e-05_dp, &
6.188411270912267e-05_dp, &
0.00011984809712211666_dp, &
0.00022515174308936954_dp, &
0.00041280335546669807_dp, &
0.0007410557365297971_dp, &
0.0013054805064640572_dp, &
0.0022608163607798406_dp, &
0.003854540538433318_dp, &
0.006477967304857807_dp, &
0.010743336152506139_dp, &
0.01759914400488243_dp, &
0.028501352515678965_dp, &
0.04566568530951012_dp, &
0.07243687914225057_dp, &
0.11382566141804518_dp, &
0.17728511649027562_dp, &
0.2738276317065499_dp, &
0.41962643246873177_dp, &
0.6383115216108352_dp, &
0.9642825901967439_dp, &
1.4475879782153047_dp, &
2.1614613507414213_dp, &
3.2151308325865324_dp, &
4.779503431415253_dp, &
7.153609650797204_dp, &
11.025589359258367_dp, &
2.33509910395076e-06_dp, &
6.10398857813977e-06_dp, &
1.1988268517456785e-05_dp, &
2.2570731247056643e-05_dp, &
4.206867163196105e-05_dp, &
7.722227756857685e-05_dp, &
0.00013907449507452585_dp, &
0.0002457091386543816_dp, &
0.00042636927263886435_dp, &
0.0007277240308901886_dp, &
0.0012233619517539264_dp, &
0.0020280492037779985_dp, &
0.0033189631613667883_dp, &
0.005367062740884365_dp, &
0.008583092022062843_dp, &
0.0135845766573915_dp, &
0.021292762750234086_dp, &
0.033072042665578866_dp, &
0.05092942142345538_dp, &
0.07779867772338114_dp, &
0.11794435500606822_dp, &
0.17753757154997848_dp, &
0.26548687080483935_dp, &
0.3946762763128839_dp, &
0.5839408961339644_dp, &
0.8616294721408878_dp, &
1.27328786626255_dp, &
1.902362545255473_dp, &
2.944458259522062_dp, &
5.1518134931437265_dp /)
ELSE IF(E_range < 1100000.0_dp) THEN
aw(:) = (/ &
6.070765776691511e-07_dp, &
3.416871225040894e-06_dp, &
9.561919657171797e-06_dp, &
2.1683051083242104e-05_dp, &
4.515689173365717e-05_dp, &
8.982491421280511e-05_dp, &
0.00017296367691311632_dp, &
0.00032428348083868487_dp, &
0.0005940202656656541_dp, &
0.001065793876730367_dp, &
0.0018767883279394732_dp, &
0.0032490764599530236_dp, &
0.0055377486897541_dp, &
0.009304158830434279_dp, &
0.015426407047504075_dp, &
0.02526460961553787_dp, &
0.0409062041812229_dp, &
0.06552740429734295_dp, &
0.10392219781712_dp, &
0.16327176797255086_dp, &
0.25425762049527484_dp, &
0.3926657116824371_dp, &
0.6016960744289409_dp, &
0.9153062522487351_dp, &
1.3831427017478757_dp, &
2.0781516959653166_dp, &
3.1094698633638336_dp, &
4.64815163738497_dp, &
6.99350202546865_dp, &
10.833342615226073_dp, &
1.5810344568110277e-06_dp, &
4.201983132608034e-06_dp, &
8.509301559244783e-06_dp, &
1.660217499192091e-05_dp, &
3.1937948250709083e-05_dp, &
6.0196622838881555e-05_dp, &
0.00011090164490251824_dp, &
0.00019990225033881748_dp, &
0.0003531805009952824_dp, &
0.000612721804071078_dp, &
0.0010454964371682956_dp, &
0.0017570690249359784_dp, &
0.0029120339233426416_dp, &
0.004764440773172712_dp, &
0.007702737556158165_dp, &
0.012315669576916223_dp, &
0.019488241944529968_dp, &
0.03054057337997689_dp, &
0.0474276677063864_dp, &
0.0730254826754342_dp, &
0.11153946595828376_dp, &
0.1690888848289202_dp, &
0.2545515426791882_dp, &
0.3808215765752087_dp, &
0.5668083660036856_dp, &
0.8410170454494971_dp, &
1.2492234533619844_dp, &
1.8750991975085831_dp, &
2.91420232217261_dp, &
5.117153113908438_dp /)
ELSE IF(E_range < 1800000.0_dp) THEN
aw(:) = (/ &
3.9676424024146875e-07_dp, &
2.2573188915559954e-06_dp, &
6.453978574008018e-06_dp, &
1.5076014172691726e-05_dp, &
3.241300804980766e-05_dp, &
6.644316731632574e-05_dp, &
0.00013145859105193587_dp, &
0.0002525113816701914_dp, &
0.0004726891479758136_dp, &
0.0008648183450896165_dp, &
0.0015500448519401517_dp, &
0.00272698303228361_dp, &
0.004716941205527944_dp, &
0.008033374074156422_dp, &
0.013487530814882825_dp, &
0.022347744374551653_dp, &
0.03657761361377878_dp, &
0.0591894063537043_dp, &
0.09476465831169727_dp, &
0.15021601736547668_dp, &
0.23589569241536856_dp, &
0.3672011309219512_dp, &
0.5668972518738795_dp, &
0.8684896297613345_dp, &
1.3212084913196627_dp, &
1.9976957850654782_dp, &
3.00698161245029_dp, &
4.520261622699826_dp, &
6.837125864760079_dp, &
10.645115539809895_dp, &
1.0356874216884052e-06_dp, &
2.812088667709278e-06_dp, &
5.914445821533368e-06_dp, &
1.2021336664547268e-05_dp, &
2.3940050145157427e-05_dp, &
4.6430499419310365e-05_dp, &
8.765815047744608e-05_dp, &
0.00016143932217372545_dp, &
0.0002907528424437285_dp, &
0.0005132235070598931_dp, &
0.0008895970920506905_dp, &
0.0015167009659352298_dp, &
0.0025470709958953243_dp, &
0.004218406965056695_dp, &
0.006897410357753999_dp, &
0.011144514451186464_dp, &
0.017808764773569234_dp, &
0.028165958189904933_dp, &
0.044118529654135305_dp, &
0.0684833014292294_dp, &
0.1054043456969387_dp, &
0.16094570015015752_dp, &
0.24395004441521287_dp, &
0.3673168074143_dp, &
0.550025709656669_dp, &
0.8207374761240367_dp, &
1.2254629611475691_dp, &
1.8481116761650778_dp, &
2.884215749911384_dp, &
5.082816358644647_dp /)
ELSE IF(E_range < 3300000.0_dp) THEN
aw(:) = (/ &
2.3571329342099874e-07_dp, &
1.3632343165414816e-06_dp, &
4.023756596461991e-06_dp, &
9.795394424152668e-06_dp, &
2.195714519062327e-05_dp, &
4.674841254419276e-05_dp, &
9.563512339593217e-05_dp, &
0.00018918230274547437_dp, &
0.00036348321259617266_dp, &
0.0006806544569921736_dp, &
0.001245724390748882_dp, &
0.002233454990133142_dp, &
0.00393042938989855_dp, &
0.0068003647335751725_dp, &
0.011584430850915652_dp, &
0.019453873577319414_dp, &
0.032240162866761804_dp, &
0.052779185694875504_dp, &
0.08542204924536101_dp, &
0.136787796738067_dp, &
0.2168656674798965_dp, &
0.3406212481786099_dp, &
0.5303313641534245_dp, &
0.8189892522205988_dp, &
1.2553469341161432_dp, &
1.9116889104572583_dp, &
2.896909089589307_dp, &
4.382351904816623_dp, &
6.667935664555373_dp, &
10.44093549828013_dp, &
6.17436261768095e-07_dp, &
1.7313881986083482e-06_dp, &
3.841234035608705e-06_dp, &
8.228727173271725e-06_dp, &
1.709820341911912e-05_dp, &
3.4332974259821285e-05_dp, &
6.676429369632049e-05_dp, &
0.00012617607786624823_dp, &
0.0002325052400773221_dp, &
0.0004189073145508487_dp, &
0.000739676578430109_dp, &
0.0012824873185346604_dp, &
0.0021871159459730906_dp, &
0.0036737980596718025_dp, &
0.006085791120718465_dp, &
0.009952725340223944_dp, &
0.01608416328295351_dp, &
0.02570677247015947_dp, &
0.04066411629576856_dp, &
0.06370600198386596_dp, &
0.09890586632566792_dp, &
0.15226260044097462_dp, &
0.23257478024275838_dp, &
0.35274192070832605_dp, &
0.5318168800175598_dp, &
0.7986316416317844_dp, &
1.1994633292908101_dp, &
1.8185003473260537_dp, &
2.851271019221645_dp, &
5.045111765749473_dp /)
ELSE IF(E_range < 6000000.0_dp) THEN
aw(:) = (/ &
1.4156108142769783e-07_dp, &
8.358408090279267e-07_dp, &
2.563425668333653e-06_dp, &
6.532496957198105e-06_dp, &
1.529085429435563e-05_dp, &
3.381004239065925e-05_dp, &
7.145759572627074e-05_dp, &
0.0001454048583206386_dp, &
0.0002863676783658706_dp, &
0.0005481076173216848_dp, &
0.0010229120533241523_dp, &
0.0018664583231752616_dp, &
0.0033372345185815654_dp, &
0.005858297289124_dp, &
0.010112965971934416_dp, &
0.017191603273076256_dp, &
0.028814641239522522_dp, &
0.04766849746482397_dp, &
0.07790741829272711_dp, &
0.12589760019317675_dp, &
0.2013131761655887_dp, &
0.3187416881227427_dp, &
0.5000288557106496_dp, &
0.7777107650960449_dp, &
1.200107192058866_dp, &
1.8391727503056268_dp, &
2.803666178964093_dp, &
4.265054523983134_dp, &
6.5235529232359735_dp, &
10.266242055151874_dp, &
3.724391604535185e-07_dp, &
1.0869250980985111e-06_dp, &
2.5605516562597e-06_dp, &
5.7850356623393795e-06_dp, &
1.2528619638131183e-05_dp, &
2.6016291223482876e-05_dp, &
5.204832564381964e-05_dp, &
0.00010081228246405614_dp, &
0.00018982349811898816_dp, &
0.0003486351977763963_dp, &
0.000626280938172421_dp, &
0.0011028905979224082_dp, &
0.0019076162117445236_dp, &
0.0032460080919710396_dp, &
0.005441424633595424_dp, &
0.00899711363781583_dp, &
0.01468851247521531_dp, &
0.023699427657021866_dp, &
0.03782154672437762_dp, &
0.059744945747629535_dp, &
0.09347917704011861_dp, &
0.1449628687820927_dp, &
0.2229515840638538_dp, &
0.34033980700933175_dp, &
0.5162398701965166_dp, &
0.7796319970356158_dp, &
1.1770306615230532_dp, &
1.7928810285625805_dp, &
2.8227297861916565_dp, &
5.012462806749877_dp /)
ELSE IF(E_range < 18000000.0_dp) THEN
aw(:) = (/ &
5.619543540742037e-08_dp, &
3.5072928313745456e-07_dp, &
1.1785442886447028e-06_dp, &
3.3029924055616623e-06_dp, &
8.394882189434066e-06_dp, &
1.987896075698813e-05_dp, &
4.449626404673206e-05_dp, &
9.506973724339386e-05_dp, &
0.00019528307697601834_dp, &
0.0003877720465039687_dp, &
0.0007475784942017128_dp, &
0.001404160723388905_dp, &
0.002576876821658531_dp, &
0.004631430044514_dp, &
0.008168570072664031_dp, &
0.014161927673231893_dp, &
0.024169988544510326_dp, &
0.04065896918886588_dp, &
0.06749026098429482_dp, &
0.1106503055340562_dp, &
0.17933541972048098_dp, &
0.28755422315274737_dp, &
0.45648527437792397_dp, &
0.7179486868235401_dp, &
1.119577032884951_dp, &
1.7327879813628193_dp, &
2.6661030665007965_dp, &
4.091161983842176_dp, &
6.3086506803163065_dp, &
10.005419833797736_dp, &
1.4962249041277701e-07_dp, &
4.835223883864694e-07_dp, &
1.2927723889950915e-06_dp, &
3.221113081192073e-06_dp, &
7.5085382788298944e-06_dp, &
1.653889597388599e-05_dp, &
3.475621269796358e-05_dp, &
7.020789074538884e-05_dp, &
0.000137107896070371_dp, &
0.0002600215543190807_dp, &
0.00048059369103468994_dp, &
0.0008682158907070847_dp, &
0.0015367276921454673_dp, &
0.0026702645844449866_dp, &
0.004562839898703722_dp, &
0.007678366608423676_dp, &
0.012740848501038493_dp, &
0.02086877544772841_dp, &
0.03377384609132226_dp, &
0.05405280484981297_dp, &
0.08561371875294983_dp, &
0.13429708936338725_dp, &
0.20878452860262342_dp, &
0.32195328331951945_dp, &
0.4929981089833493_dp, &
0.7511227960385115_dp, &
1.1432128213684527_dp, &
1.7541296769802968_dp, &
2.7794877584199225_dp, &
4.963023556769552_dp /)
ELSE IF(E_range < 50000000.0_dp) THEN
aw(:) = (/ &
2.4343774571442308e-08_dp, &
1.6513085830933388e-07_dp, &
6.196802532079424e-07_dp, &
1.9114039978346155e-06_dp, &
5.2375690334248195e-06_dp, &
1.316277905469794e-05_dp, &
3.092030271194152e-05_dp, &
6.876819455159673e-05_dp, &
0.0001461441409773963_dp, &
0.0002988264714800754_dp, &
0.0005910302376020617_dp, &
0.0011354771724870155_dp, &
0.0021261659208679307_dp, &
0.0038911087717360972_dp, &
0.0069760886290369745_dp, &
0.012276085772085119_dp, &
0.021239226973771258_dp, &
0.03618004618429005_dp, &
0.060756103616981624_dp, &
0.10068683961443668_dp, &
0.16482922554469034_dp, &
0.2667764798918916_dp, &
0.4272231990996556_dp, &
0.6774637303634756_dp, &
1.064618278549332_dp, &
1.6596947503579922_dp, &
2.5710208535654164_dp, &
3.9703478472317504_dp, &
6.158711044835771_dp, &
9.822845740299956_dp, &
6.607779435757453e-08_dp, &
2.4531316662915766e-07_dp, &
7.456766806787507e-07_dp, &
2.026328956129065e-06_dp, &
5.0327732385235214e-06_dp, &
1.1652295901214502e-05_dp, &
2.5505955301673674e-05_dp, &
5.331453501898778e-05_dp, &
0.0001072058234797249_dp, &
0.00020853777085767683_dp, &
0.0003941273627591354_dp, &
0.0007262462873395727_dp, &
0.0013084426047463612_dp, &
0.00231027419159065_dp, &
0.004005535279633006_dp, &
0.006830718203653046_dp, &
0.011473559821092195_dp, &
0.019005957730023747_dp, &
0.031081856143062504_dp, &
0.0502296954993904_dp, &
0.08028203266369409_dp, &
0.12700459445500437_dp, &
0.19901982937944282_dp, &
0.3091853239783209_dp, &
0.4767481905112323_dp, &
0.7310696923220488_dp, &
1.1193069602208368_dp, &
1.72663775937696_dp, &
2.74875508795335_dp, &
4.927904659791831_dp /)
ELSE
aw(:) = (/ &
3.3134286206436378e-09_dp, &
3.608476731974466e-08_dp, &
1.9011771017402648e-07_dp, &
7.319641454121685e-07_dp, &
2.3402664590560514e-06_dp, &
6.59508695658557e-06_dp, &
1.6937909149354525e-05_dp, &
4.048237319663682e-05_dp, &
9.132011445721631e-05_dp, &
0.00019639172259169963_dp, &
0.0004056673541505352_dp, &
0.0008094360939409007_dp, &
0.001567130287522824_dp, &
0.0029545924806747004_dp, &
0.005440438931022152_dp, &
0.009807773172795757_dp, &
0.017345795847210112_dp, &
0.030148043017114434_dp, &
0.0515717640546498_dp, &
0.08693872933672857_dp, &
0.14459504048117927_dp, &
0.23750169179013617_dp, &
0.3856081802743787_dp, &
0.6193889538319287_dp, &
0.9851526108818242_dp, &
1.5532426287786882_dp, &
2.431651849778814_dp, &
3.7922799603529507_dp, &
5.93671227399739_dp, &
9.551588089060862_dp, &
1.0546633313159667e-08_dp, &
6.866238130433312e-08_dp, &
2.8014791592507246e-07_dp, &
9.089927672975243e-07_dp, &
2.5558210419120903e-06_dp, &
6.498693707400478e-06_dp, &
1.531991208406891e-05_dp, &
3.4023779895267246e-05_dp, &
7.197856727334818e-05_dp, &
0.00014621506144801751_dp, &
0.0002869199340733901_dp, &
0.0005464264697957605_dp, &
0.0010137006051227095_dp, &
0.0018373548659257801_dp, &
0.0032617470785987366_dp, &
0.005682944230668296_dp, &
0.009734554219352109_dp, &
0.016418075032354977_dp, &
0.02729906435610205_dp, &
0.04479995778261982_dp, &
0.07263414004010217_dp, &
0.11644650227746163_dp, &
0.1847594577631817_dp, &
0.2903885825624092_dp, &
0.4526492240843357_dp, &
0.7011369514153092_dp, &
1.0834309171758731_dp, &
1.6852193328076974_dp, &
2.7023644845108916_dp, &
4.874923027763007_dp /)
END IF
END SUBROUTINE

SUBROUTINE get_coeff_32(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio
E_ratio = 1.0_dp
IF(E_range < 4862.0_dp) THEN
E_ratio = 4862.0_dp/E_range
aw(:) = (/ &
6.390102821995496e-05_dp, &
0.0003410048490902184_dp, &
0.000858015893339577_dp, &
0.0016504328288047757_dp, &
0.002778247721262644_dp, &
0.004336751679872466_dp, &
0.0064710209548362084_dp, &
0.009393407912277671_dp, &
0.01340416788754177_dp, &
0.018917565422190717_dp, &
0.026497301468281335_dp, &
0.036905213045864395_dp, &
0.0511672148624311_dp, &
0.07066128158960995_dp, &
0.09723384839331731_dp, &
0.13335312533591778_dp, &
0.18231048817210543_dp, &
0.24848450523394644_dp, &
0.3376865863344385_dp, &
0.4576131552357413_dp, &
0.6184374387625521_dp, &
0.8335859250026539_dp, &
1.1207633263472254_dp, &
1.5033220996389998_dp, &
2.0121331479703484_dp, &
2.6882384204905283_dp, &
3.5868412256539663_dp, &
4.783852020373796_dp, &
6.387973693919086_dp, &
8.566792744620582_dp, &
11.61702931918549_dp, &
16.243244163196902_dp, &
0.0001644750107629364_dp, &
0.0003927652580734362_dp, &
0.0006470526340880139_dp, &
0.0009476004809462573_dp, &
0.0013236672585754251_dp, &
0.00181722804655121_dp, &
0.0024863648293148356_dp, &
0.0034080264366486593_dp, &
0.004681956510169319_dp, &
0.006437986099490384_dp, &
0.00884737093390095_dp, &
0.01213798662642807_dp, &
0.01661375946425665_dp, &
0.022679584239104586_dp, &
0.030873595466529645_dp, &
0.0419091652309378_dp, &
0.05672961944641861_dp, &
0.07657952416894345_dp, &
0.10309761695759492_dp, &
0.13843826485639682_dp, &
0.18543118617833296_dp, &
0.24779409719067566_dp, &
0.3304222201755898_dp, &
0.43979746242641576_dp, &
0.5846010059990558_dp, &
0.7767070228661042_dp, &
1.0329650164221098_dp, &
1.3787901136203442_dp, &
1.8564287417380754_dp, &
2.547529279131862_dp, &
3.6531886321446416_dp, &
5.969270620351697_dp /)
ELSE IF(E_range < 5846.0_dp) THEN
aw(:) = (/ &
5.442349002883837e-05_dp, &
0.0002906163628933693_dp, &
0.0007321293214401846_dp, &
0.0014109784539713872_dp, &
0.0023816365307190094_dp, &
0.0037311772208214787_dp, &
0.005592773584012946_dp, &
0.0081619152463338_dp, &
0.011715729132325964_dp, &
0.016637929524033673_dp, &
0.023453105417011018_dp, &
0.03287403908552556_dp, &
0.045865900699847945_dp, &
0.06373214396298486_dp, &
0.08822857905561622_dp, &
0.12171425464017217_dp, &
0.16735049859917261_dp, &
0.22936294868166074_dp, &
0.3133859506800179_dp, &
0.4269147709487938_dp, &
0.5798994400836499_dp, &
0.7855261781163106_dp, &
1.0612512378052221_dp, &
1.4301841263166695_dp, &
1.9229772371853768_dp, &
2.580501782235826_dp, &
3.4578625966778955_dp, &
4.630962088019041_dp, &
6.208582258228505_dp, &
8.358418878620247_dp, &
11.37712612667173_dp, &
15.968170544511102_dp, &
0.0001401015051695221_dp, &
0.00033499665110023704_dp, &
0.0005532522530987489_dp, &
0.0008133560062402562_dp, &
0.0011422450286066607_dp, &
0.0015787618621387204_dp, &
0.002176697085238142_dp, &
0.0030073676766694266_dp, &
0.004163649639092577_dp, &
0.005767402761174793_dp, &
0.007980687395399007_dp, &
0.011020625281672214_dp, &
0.015178459352483943_dp, &
0.020844174229750198_dp, &
0.028538588286649656_dp, &
0.03895532923923562_dp, &
0.05301575138551127_dp, &
0.07194074865569469_dp, &
0.09734467595148226_dp, &
0.1313584288969531_dp, &
0.17679160298883081_dp, &
0.23734854692003116_dp, &
0.31792227850976007_dp, &
0.42500882037583165_dp, &
0.5673248919963556_dp, &
0.7568049338284129_dp, &
1.0103817704846647_dp, &
1.3535647228420928_dp, &
1.8286721197765015_dp, &
2.517319618502594_dp, &
3.6202102280059703_dp, &
5.930981680526662_dp /)
ELSE IF(E_range < 6665.0_dp) THEN
aw(:) = (/ &
4.8544808464658475e-05_dp, &
0.0002593485712181486_dp, &
0.0006539518705617748_dp, &
0.0012620998863735048_dp, &
0.002134648029009164_dp, &
0.003353250117842424_dp, &
0.005043178699940835_dp, &
0.007388686080392843_dp, &
0.010651467508498257_dp, &
0.01519484467183565_dp, &
0.02151725732341947_dp, &
0.03029858754646109_dp, &
0.04246309490817034_dp, &
0.05926381239480956_dp, &
0.08239494336727046_dp, &
0.11414097839482948_dp, &
0.15757400511758252_dp, &
0.21681422548727597_dp, &
0.29737332757633034_dp, &
0.4066065418252456_dp, &
0.5543077141349219_dp, &
0.7534939940706474_dp, &
1.0214457127670413_dp, &
1.3810991072747199_dp, &
1.8629492914867638_dp, &
2.5077428833752067_dp, &
3.370510682262872_dp, &
4.52714545796296_dp, &
6.086483167699423_dp, &
8.216299008759131_dp, &
11.213210923746692_dp, &
15.779945298966286_dp, &
0.00012498176375406726_dp, &
0.000299130758824733_dp, &
0.0004949254877793316_dp, &
0.0007296849915859161_dp, &
0.0010288042473983914_dp, &
0.0014290235014882077_dp, &
0.0019812415860893163_dp, &
0.0027529950493393966_dp, &
0.00383257137438519_dp, &
0.005336492795437645_dp, &
0.007420608589016774_dp, &
0.01029470391501449_dp, &
0.014241296894904196_dp, &
0.01964004703493962_dp, &
0.02699971016577413_dp, &
0.03700007872604122_dp, &
0.050547012803328586_dp, &
0.06884458868981981_dp, &
0.09348967784530542_dp, &
0.12659613160197555_dp, &
0.17095862941207848_dp, &
0.2302711267303353_dp, &
0.3094239152326776_dp, &
0.4149216994326361_dp, &
0.555505084520509_dp, &
0.7431502415960526_dp, &
0.9948487920453783_dp, &
1.3361777826685788_dp, &
1.809509286088653_dp, &
2.496441460806496_dp, &
3.597410863544938_dp, &
5.9045255378523525_dp /)
ELSE IF(E_range < 7800.0_dp) THEN
aw(:) = (/ &
4.2319791684398945e-05_dp, &
0.0002262256833340046_dp, &
0.0005710759339812151_dp, &
0.0011041000499916853_dp, &
0.0018721273080030145_dp, &
0.002950744254164121_dp, &
0.004456329948077625_dp, &
0.006560445775378507_dp, &
0.009507343812482433_dp, &
0.013637283047013359_dp, &
0.01941907839765781_dp, &
0.027495223047590025_dp, &
0.03874330338636378_dp, &
0.054358588747561776_dp, &
0.07596440384363516_dp, &
0.10575909660902735_dp, &
0.14671121294566578_dp, &
0.20281810337064493_dp, &
0.27944792395127677_dp, &
0.38379131316561327_dp, &
0.525457692252208_dp, &
0.7172635733637082_dp, &
0.9762793617085143_dp, &
1.3252331952635001_dp, &
1.7944297509586076_dp, &
2.4244629711653745_dp, &
3.270270626899914_dp, &
4.407729666691881_dp, &
5.945738068832148_dp, &
8.052167730865541_dp, &
11.023603896745234_dp, &
15.561924493694205_dp, &
0.00010896980958811763_dp, &
0.0002611185815588283_dp, &
0.00043301817888623847_dp, &
0.0006406820460768114_dp, &
0.0009077674409158128_dp, &
0.001268624601649587_dp, &
0.0017708556976942175_dp, &
0.002477704567241762_dp, &
0.0034722665030748113_dp, &
0.004865004260125733_dp, &
0.006804656478190395_dp, &
0.009492536281842186_dp, &
0.013201017163803366_dp, &
0.01829768554187978_dp, &
0.025277118873799437_dp, &
0.034802761118126926_dp, &
0.04776205831199323_dp, &
0.06533896907314592_dp, &
0.08910928388293891_dp, &
0.1211660799438405_dp, &
0.16428554188642505_dp, &
0.22214824519628396_dp, &
0.29963998828115096_dp, &
0.40327442471190905_dp, &
0.5418192779326029_dp, &
0.7272994796496782_dp, &
0.9767765646718063_dp, &
1.315909458684106_dp, &
1.7871373292268897_dp, &
2.4720433230436907_dp, &
3.5707586649853615_dp, &
5.873612522584144_dp /)
ELSE IF(E_range < 10044.0_dp) THEN
aw(:) = (/ &
3.393255390660133e-05_dp, &
0.00018157221465841256_dp, &
0.0004592307676320871_dp, &
0.0008905276367436023_dp, &
0.0015164746316516245_dp, &
0.0024038182425049908_dp, &
0.003655880996707775_dp, &
0.0054255179504103905_dp, &
0.007931245931576318_dp, &
0.0114792650369607_dp, &
0.01649451279684121_dp, &
0.023563843280162002_dp, &
0.033494992691347265_dp, &
0.047396265575400595_dp, &
0.06678363974811272_dp, &
0.093724229945987_dp, &
0.13102792038715239_dp, &
0.18250270781166214_dp, &
0.2532942033927338_dp, &
0.3503362853107127_dp, &
0.4829488331079465_dp, &
0.6636312054412428_dp, &
0.9091194626135006_dp, &
1.241807411536441_dp, &
1.6916904363164948_dp, &
2.299110603722669_dp, &
3.1188473050465815_dp, &
4.226742880544712_dp, &
5.73178901479017_dp, &
7.802015996740442_dp, &
10.73398018079764_dp, &
15.228283328315147_dp, &
8.739335036060486e-05_dp, &
0.00020983760313852035_dp, &
0.0003493229451058535_dp, &
0.0005199664146008647_dp, &
0.0007428694438573918_dp, &
0.001048824687194698_dp, &
0.0014805150854819704_dp, &
0.0020948189118065326_dp, &
0.0029671598393268173_dp, &
0.004198990073584458_dp, &
0.005928343272370908_dp, &
0.008343628060712927_dp, &
0.011701624144711492_dp, &
0.0163512213147089_dp, &
0.022764888824347185_dp, &
0.03158039794694182_dp, &
0.04365605404949107_dp, &
0.06014369592103612_dp, &
0.08258509281019695_dp, &
0.11303931699899293_dp, &
0.1542516135087157_dp, &
0.20987914795565896_dp, &
0.28479788380848736_dp, &
0.38553266378218126_dp, &
0.5208913754190811_dp, &
0.7029743313415381_dp, &
0.9489538722825864_dp, &
1.2846216880545815_dp, &
1.7525303625505235_dp, &
2.4342520367725067_dp, &
3.5294590063092537_dp, &
5.82574638448254_dp /)
ELSE IF(E_range < 14058.0_dp) THEN
aw(:) = (/ &
2.5286922170354894e-05_dp, &
0.00013550300515152938_dp, &
0.00034365254863389185_dp, &
0.0006692829830793605_dp, &
0.0011467768099562_dp, &
0.0018326842825994733_dp, &
0.0028151146148173586_dp, &
0.004225026003742279_dp, &
0.006250780019705354_dp, &
0.009158586718943748_dp, &
0.013321604127882514_dp, &
0.01926054419553751_dp, &
0.02769941545264008_dp, &
0.03964139074601272_dp, &
0.056471580120530336_dp, &
0.08009577889831763_dp, &
0.11312722824816804_dp, &
0.15913731138320336_dp, &
0.22299123810302657_dp, &
0.3112966160076166_dp, &
0.43300212495019025_dp, &
0.6001966661049284_dp, &
0.8291790761739616_dp, &
1.141900704443846_dp, &
1.5679415214227161_dp, &
2.147298128643314_dp, &
2.934523705189343_dp, &
4.00539737112386_dp, &
5.469022915149277_dp, &
7.493643904207673_dp, &
10.375817412384437_dp, &
14.814601510443076_dp, &
6.514769920993364e-05_dp, &
0.00015687397983623111_dp, &
0.00026260059811101074_dp, &
0.000394267458907849_dp, &
0.000569987485627106_dp, &
0.0008163276402528582_dp, &
0.0011701237352853642_dp, &
0.001680758396099781_dp, &
0.002414639429717159_dp, &
0.0034624986342827427_dp, &
0.004949404187288667_dp, &
0.00704789552613599_dp, &
0.00999535091125554_dp, &
0.014117167183015276_dp, &
0.019857768872086755_dp, &
0.02782204007098351_dp, &
0.03883055990410155_dp, &
0.05399308565310966_dp, &
0.07480617533828689_dp, &
0.1032828703661818_dp, &
0.14212537017885443_dp, &
0.19495650401901887_dp, &
0.2666345532300576_dp, &
0.3636933937165475_dp, &
0.4949878042793592_dp, &
0.6727122940047311_dp, &
0.91418280277028_dp, &
1.2453685853158643_dp, &
1.7089820529103965_dp, &
2.3866034784779506_dp, &
3.4773514427582293_dp, &
5.765409584335823_dp /)
ELSE IF(E_range < 19114.0_dp) THEN
aw(:) = (/ &
1.932135319864962e-05_dp, &
0.00010368197789193328_dp, &
0.0002636661772942731_dp, &
0.0005157191470705282_dp, &
0.0008891130218806238_dp, &
0.0014324323017507073_dp, &
0.002221781063255329_dp, &
0.003370774273254239_dp, &
0.005043864110016273_dp, &
0.00747544974094905_dp, &
0.010997228176632299_dp, &
0.01607649385487372_dp, &
0.023369007932130658_dp, &
0.033791453456649585_dp, &
0.04862029962984213_dp, &
0.06962623074409648_dp, &
0.09925634818539467_dp, &
0.14088037853305638_dp, &
0.199122452009721_dp, &
0.2803071363072363_dp, &
0.39305808457961094_dp, &
0.549101229063986_dp, &
0.7643445772055675_dp, &
1.0603390883163382_dp, &
1.466283147413628_dp, &
2.0218493641999524_dp, &
2.781370504056042_dp, &
3.8205517999540817_dp, &
5.248586694740235_dp, &
7.23391404025796_dp, &
10.073125286074314_dp, &
14.464007598891323_dp, &
4.9794334064441845e-05_dp, &
0.0001202442466871008_dp, &
0.00020239142207391234_dp, &
0.00030648388928353133_dp, &
0.00044826394218884584_dp, &
0.0006508962458350465_dp, &
0.0009465213181019036_dp, &
0.001378560941846276_dp, &
0.002006235556189542_dp, &
0.002911588665713435_dp, &
0.004208981250342487_dp, &
0.006057652087479506_dp, &
0.008678530354123182_dp, &
0.012376893766972177_dp, &
0.01757291406967402_dp, &
0.02484274557305271_dp, &
0.034973648086824426_dp, &
0.04903775694722466_dp, &
0.0684906283923225_dp, &
0.09530279631811456_dp, &
0.1321356654214213_dp, &
0.1825779770666456_dp, &
0.2514677535724369_dp, &
0.34534170861863805_dp, &
0.473091186246975_dp, &
0.6469906291688993_dp, &
0.8844828931716822_dp, &
1.2116994219481383_dp, &
1.6715060390575855_dp, &
2.3455108368533817_dp, &
3.432378601360354_dp, &
5.713384147526647_dp /)
ELSE IF(E_range < 25870.0_dp) THEN
aw(:) = (/ &
1.4818639796865528e-05_dp, &
7.963886363515554e-05_dp, &
0.00020311223580610354_dp, &
0.00039911404589973475_dp, &
0.0006926333595966217_dp, &
0.0011254943882583104_dp, &
0.0017635291623248677_dp, &
0.0027054780460561163_dp, &
0.004095260696633318_dp, &
0.006139843660554759_dp, &
0.009134914466544177_dp, &
0.013500979005727987_dp, &
0.019833495884515773_dp, &
0.02897206838249403_dp, &
0.042095533039198385_dp, &
0.06085215831868936_dp, &
0.08753730146299196_dp, &
0.1253350272549887_dp, &
0.17864571781396288_dp, &
0.2535290933953418_dp, &
0.35830209305587557_dp, &
0.5043450664128841_dp, &
0.7071902955364795_dp, &
0.9879995840794562_dp, &
1.3755955013642696_dp, &
1.9093257812455384_dp, &
2.6432950751333295_dp, &
3.6531218956627853_dp, &
5.048075954994877_dp, &
6.99678689342338_dp, &
9.795906582751762_dp, &
14.142088288126018_dp, &
3.820307739174523e-05_dp, &
9.253211423715926e-05_dp, &
0.00015666144789002137_dp, &
0.00023941035175647392_dp, &
0.0003544795957843895_dp, &
0.0005220706038556393_dp, &
0.0007702567427026573_dp, &
0.0011373348215882129_dp, &
0.001676302988186351_dp, &
0.002461551034244103_dp, &
0.003597851119996275_dp, &
0.005232383925906072_dp, &
0.007571020540645877_dp, &
0.01090045623553385_dp, &
0.01561826230147036_dp, &
0.022273572546892197_dp, &
0.03162200007875945_dp, &
0.04469955727258082_dp, &
0.06292193738822648_dp, &
0.08821770562181855_dp, &
0.12320712191071415_dp, &
0.17144328423198826_dp, &
0.23774090556227068_dp, &
0.3286348439935499_dp, &
0.45304683721268363_dp, &
0.6233242933376546_dp, &
0.8570306244015214_dp, &
1.1804560303362448_dp, &
1.6366228282940412_dp, &
2.3071833124934_dp, &
3.3903998854623985_dp, &
5.664862602289436_dp /)
ELSE IF(E_range < 35180.0_dp) THEN
aw(:) = (/ &
1.1315989220620849e-05_dp, &
6.091484824470965e-05_dp, &
0.0001558556348063919_dp, &
0.0003078203984115187_dp, &
0.0005380983662976303_dp, &
0.0008826023610964867_dp, &
0.0013981206622491853_dp, &
0.0021702768627684973_dp, &
0.0033248504851841115_dp, &
0.005044444967959304_dp, &
0.007592519217726472_dp, &
0.011347339662602948_dp, &
0.016849450868861544_dp, &
0.024867669342019266_dp, &
0.03649043474384647_dp, &
0.05325176815491175_dp, &
0.07730430050505666_dp, &
0.1116561161775443_dp, &
0.16049387125009024_dp, &
0.22962230812156062_dp, &
0.3270606812998151_dp, &
0.46385106272587334_dp, &
0.6551545548992531_dp, &
0.9217445361800276_dp, &
1.2920638357850092_dp, &
1.8051264992628424_dp, &
2.514796276521399_dp, &
3.496590833178514_dp, &
4.85984548686659_dp, &
6.773380959237685_dp, &
9.533934166411905_dp, &
13.837118014874001_dp, &
2.918394667912216e-05_dp, &
7.092102053120252e-05_dp, &
0.0001208490607819795_dp, &
0.00018654379226223144_dp, &
0.00027989427379719806_dp, &
0.0004184502969887145_dp, &
0.000626673050017551_dp, &
0.0009383293128013663_dp, &
0.0014008572107836192_dp, &
0.002081684144431634_dp, &
0.0030767293366374213_dp, &
0.004521918459768059_dp, &
0.006608939226821494_dp, &
0.009606834282132472_dp, &
0.013891522205312699_dp, &
0.0199860217642558_dp, &
0.028615075699006714_dp, &
0.04077910134384887_dp, &
0.057854053870927825_dp, &
0.08172606721705969_dp, &
0.11497300539733377_dp, &
0.16111009858027334_dp, &
0.22492545113901904_dp, &
0.3129477675464925_dp, &
0.43412426576959084_dp, &
0.6008706467915673_dp, &
0.8308678912066612_dp, &
1.1505656278854524_dp, &
1.6031493439397482_dp, &
2.270331555451258_dp, &
3.3500083706514476_dp, &
5.6182166490122984_dp /)
ELSE IF(E_range < 58986.0_dp) THEN
aw(:) = (/ &
7.188720818598015e-06_dp, &
3.881661788798954e-05_dp, &
9.991514616563082e-05_dp, &
0.00019924670807668527_dp, &
0.0003530921393375012_dp, &
0.0005892390933533132_dp, &
0.0009519594342766803_dp, &
0.0015086954155209336_dp, &
0.0023599840025008696_dp, &
0.003654238190455114_dp, &
0.005609186735480063_dp, &
0.008542461005020551_dp, &
0.01291487453690728_dp, &
0.01939133485309362_dp, &
0.028926175238232574_dp, &
0.04288217331115784_dp, &
0.06319585740381348_dp, &
0.09260618135064007_dp, &
0.13496966452445772_dp, &
0.19569319485572675_dp, &
0.2823266941469131_dp, &
0.40537308432315655_dp, &
0.5793949454394274_dp, &
0.8245311284762115_dp, &
1.1685945009236132_dp, &
1.6500342331912714_dp, &
2.322292524496194_dp, &
3.260691519334744_dp, &
4.574650070688025_dp, &
6.433299519631953_dp, &
9.133561259222846_dp, &
13.369520380969645_dp, &
1.8552570543955455e-05_dp, &
4.5365280395614445e-05_dp, &
7.824430449992235e-05_dp, &
0.000123064841905601_dp, &
0.00018917614532398957_dp, &
0.00029038908164917863_dp, &
0.0004461302902503803_dp, &
0.0006838555287975979_dp, &
0.001043105615047013_dp, &
0.0015811894241096348_dp, &
0.0023809364846653818_dp, &
0.003561415027894291_dp, &
0.005292829489029218_dp, &
0.007817192602841475_dp, &
0.011476902733337355_dp, &
0.01675409107074492_dp, &
0.024324590249295588_dp, &
0.035131696353015066_dp, &
0.050486674490227534_dp, &
0.07220539175013203_dp, &
0.10279391074071534_dp, &
0.14570107833907187_dp, &
0.20566478637597496_dp, &
0.28919483053487227_dp, &
0.40526987079618043_dp, &
0.5664074675124494_dp, &
0.7904744673132619_dp, &
1.1041823544223337_dp, &
1.5509968448017435_dp, &
2.2127606942845244_dp, &
3.286839828729209_dp, &
5.545336982828468_dp /)
ELSE IF(E_range < 85052.0_dp) THEN
aw(:) = (/ &
5.212842853721373e-06_dp, &
2.8216485728556984e-05_dp, &
7.298119875698969e-05_dp, &
0.00014666687860403872_dp, &
0.00026275225618308033_dp, &
0.0004444091505151155_dp, &
0.0007287615881022521_dp, &
0.0011728502574399147_dp, &
0.0018626816535886259_dp, &
0.0029267741621968922_dp, &
0.004555909601317652_dp, &
0.007031533755260647_dp, &
0.01076629475072859_dp, &
0.016361592669832718_dp, &
0.02468887356341755_dp, &
0.03700391916268461_dp, &
0.055106793600662854_dp, &
0.08156471559554047_dp, &
0.12002134058150862_dp, &
0.17562433753378992_dp, &
0.2556145691836236_dp, &
0.370135979817623_dp, &
0.5333479158700739_dp, &
0.7649561302861829_dp, &
1.092336920188996_dp, &
1.5535402649553307_dp, &
2.201701148128128_dp, &
3.111986649331865_dp, &
4.39385578004782_dp, &
6.216650142887412_dp, &
8.877446914810541_dp, &
13.069399675595056_dp, &
1.3460633699481115e-05_dp, &
3.307685267989567e-05_dp, &
5.7604518102217696e-05_dp, &
9.195641910083116e-05_dp, &
0.00014400755053456355_dp, &
0.00022539775962319305_dp, &
0.0003526662335951264_dp, &
0.0005496236389883362_dp, &
0.0008511354053799554_dp, &
0.0013083814172165737_dp, &
0.0019961328695681448_dp, &
0.0030229467394267007_dp, &
0.00454547537639009_dp, &
0.00678848631990991_dp, &
0.010072747132319853_dp, &
0.014853695284454384_dp, &
0.021774843424492554_dp, &
0.03174125164120884_dp, &
0.04602026166587432_dp, &
0.06637923482681696_dp, &
0.0952736156331887_dp, &
0.13610397708422503_dp, &
0.19356942137149777_dp, &
0.27416081565736494_dp, &
0.3868712620397728_dp, &
0.5442810461868627_dp, &
0.7643794009427055_dp, &
1.0740573000929907_dp, &
1.5169810386781961_dp, &
2.1751039860564134_dp, &
3.2454752973208514_dp, &
5.497663131109453_dp /)
ELSE IF(E_range < 126612.0_dp) THEN
aw(:) = (/ &
3.6756703015465298e-06_dp, &
1.9954750677579616e-05_dp, &
5.191549207164041e-05_dp, &
0.00010531730513004429_dp, &
0.00019114894568566566_dp, &
0.0003284345206580502_dp, &
0.0005478574351058358_dp, &
0.0008970642726998205_dp, &
0.001448846563954588_dp, &
0.0023134440610659145_dp, &
0.00365660754627558_dp, &
0.005725809627536767_dp, &
0.008888016831791598_dp, &
0.01368381089343313_dp, &
0.020904518142872384_dp, &
0.031701561841387796_dp, &
0.04774073523681854_dp, &
0.07141882485947987_dp, &
0.10616643258939448_dp, &
0.15686955518574938_dp, &
0.23045435097627895_dp, &
0.33669592809026266_dp, &
0.48933536087816165_dp, &
0.7076244370403264_dp, &
1.0184762856769092_dp, &
1.4595108245228903_dp, &
2.0835241856266555_dp, &
2.965501956804563_dp, &
4.214929989807013_dp, &
6.001367447071038_dp, &
8.622077664029446_dp, &
12.76932426608392_dp, &
9.497600314818484e-06_dp, &
2.347760247868802e-05_dp, &
4.136863276834202e-05_dp, &
6.721947207800727e-05_dp, &
0.00010755652914793868_dp, &
0.00017204310109274583_dp, &
0.00027461020521950513_dp, &
0.00043573410701580053_dp, &
0.0006859083122114215_dp, &
0.0010704829571126738_dp, &
0.0016564711186942442_dp, &
0.0025422076964986724_dp, &
0.003871039755848531_dp, &
0.005850643295798488_dp, &
0.008780144703142966_dp, &
0.013088020403816946_dp, &
0.01938482007375741_dp, &
0.02853620387815694_dp, &
0.04176374049584569_dp, &
0.06078358076425594_dp, &
0.087996854859105_dp, &
0.12675113895160917_dp, &
0.18170115882627896_dp, &
0.2593129035941683_dp, &
0.3685886307380454_dp, &
0.5221684374385586_dp, &
0.7381657526850111_dp, &
1.0436599990100175_dp, &
1.4825352832265704_dp, &
2.136878793265541_dp, &
3.2034437811357703_dp, &
5.449257489548676_dp /)
ELSE IF(E_range < 247709.0_dp) THEN
aw(:) = (/ &
2.0394837579119305e-06_dp, &
1.1136428689231387e-05_dp, &
2.9311026954581983e-05_dp, &
6.057395938117077e-05_dp, &
0.00011273361636596131_dp, &
0.00019945262227642557_dp, &
0.00034307020413915784_dp, &
0.0005789986244949678_dp, &
0.0009626136389075716_dp, &
0.0015797047702976012_dp, &
0.0025620328180390273_dp, &
0.00411027194192014_dp, &
0.006527602589322277_dp, &
0.010268581467814463_dp, &
0.016009798955308914_dp, &
0.02475143770962476_dp, &
0.03796243042890685_dp, &
0.057786831913653906_dp, &
0.0873357482301367_dp, &
0.1310983577206706_dp, &
0.19551814727621097_dp, &
0.2897978980472306_dp, &
0.4270216153720778_dp, &
0.6257183147230295_dp, &
0.9120522680792509_dp, &
1.3229345707400562_dp, &
1.9105866829861018_dp, &
2.7496652071687744_dp, &
3.9496659363104385_dp, &
5.680486965184195_dp, &
8.239732353995374_dp, &
12.318415719148392_dp, &
5.276628830314947e-06_dp, &
1.3196668481729192e-05_dp, &
2.379458469172375e-05_dp, &
3.999888381841533e-05_dp, &
6.655277116013134e-05_dp, &
0.00011053727387478837_dp, &
0.00018248747181848964_dp, &
0.00029842809341232193_dp, &
0.0004828462223872788_dp, &
0.0007729214376067976_dp, &
0.001224631278196218_dp, &
0.001921577765363735_dp, &
0.0029876866735399437_dp, &
0.0046053732161687655_dp, &
0.007041380807008615_dp, &
0.010683336283217822_dp, &
0.016091204665789092_dp, &
0.02406937359680715_dp, &
0.03576720002987578_dp, &
0.05281872257231177_dp, &
0.07753624384014632_dp, &
0.11317828438072734_dp, &
0.16432147566237754_dp, &
0.2373819281103962_dp, &
0.34136408980953675_dp, &
0.48899079151871744_dp, &
0.6985643332266451_dp, &
0.9974643386084848_dp, &
1.4299370040199422_dp, &
2.0783190576822954_dp, &
3.138964327630411_dp, &
5.375073752754771_dp /)
ELSE IF(E_range < 452410.0_dp) THEN
aw(:) = (/ &
1.2029417634752254e-06_dp, &
6.609895643735541e-06_dp, &
1.761886326860991e-05_dp, &
3.714701435066809e-05_dp, &
7.09606454750247e-05_dp, &
0.00012924586800944736_dp, &
0.0002289177620104591_dp, &
0.00039735575327091586_dp, &
0.0006783049971089168_dp, &
0.001140918627613451_dp, &
0.0018934100817969098_dp, &
0.0031034810098270562_dp, &
0.0050286516433828725_dp, &
0.008060966671683137_dp, &
0.012792436500279173_dp, &
0.020110206423758527_dp, &
0.03133410211483852_dp, &
0.04841425570428367_dp, &
0.07421348496567018_dp, &
0.11290867950409522_dp, &
0.1705586323733128_dp, &
0.2559040285419376_dp, &
0.38149111837813177_dp, &
0.5652486706311622_dp, &
0.8327086526705059_dp, &
1.2201713489652732_dp, &
1.7793482729282073_dp, &
2.5845851027859807_dp, &
3.745354303892851_dp, &
5.4318267387256425_dp, &
7.941927006849686_dp, &
11.965772395891321_dp, &
3.1166015808204985e-06_dp, &
7.893584647537054e-06_dp, &
1.4590421523143668e-05_dp, &
2.540166835137321e-05_dp, &
4.38866273658644e-05_dp, &
7.544235506820306e-05_dp, &
0.00012836853819339647_dp, &
0.0002156510284757784_dp, &
0.00035755049600582473_dp, &
0.0005853791673295657_dp, &
0.0009470589489445969_dp, &
0.0015152726995104847_dp, &
0.0023993396733646157_dp, &
0.0037624022005128713_dp, &
0.005846143988201738_dp, &
0.009006130762738474_dp, &
0.01376205648174844_dp, &
0.02086880821287917_dp, &
0.03141648831640661_dp, &
0.04697057987136505_dp, &
0.06976767926454802_dp, &
0.10298830876584288_dp, &
0.1511376630263013_dp, &
0.2205811788480397_dp, &
0.320313759532162_dp, &
0.4631152216015924_dp, &
0.6674360776688134_dp, &
0.9609048344557064_dp, &
1.3880820200304151_dp, &
2.031544618845653_dp, &
3.087377771089118_dp, &
5.31578682998986_dp /)
ELSE IF(E_range < 1104308.0_dp) THEN
aw(:) = (/ &
5.514768391533348e-07_dp, &
3.0661463050989084e-06_dp, &
8.370607042832422e-06_dp, &
1.830793322995512e-05_dp, &
3.65866489044243e-05_dp, &
6.987140555359488e-05_dp, &
0.00012955245746241426_dp, &
0.00023469156977383316_dp, &
0.00041673066301572737_dp, &
0.0007268469188482661_dp, &
0.0012473021552990137_dp, &
0.002108784237240361_dp, &
0.0035166643543472435_dp, &
0.005790416868049286_dp, &
0.009422327737378779_dp, &
0.015164272167993003_dp, &
0.02415507435487594_dp, &
0.038106181361843135_dp, &
0.05957064947730835_dp, &
0.09233052534676135_dp, &
0.14195168557080626_dp, &
0.2165746550862792_dp, &
0.3280373760158957_dp, &
0.493465990037226_dp, &
0.7375325370826453_dp, &
1.0956894156446793_dp, &
1.6189196862917987_dp, &
2.381098691427177_dp, &
3.4916207819996172_dp, &
5.1210060943387665_dp, &
7.5676555925110085_dp, &
11.520669806145102_dp, &
1.4324492833187423e-06_dp, &
3.714945365379106e-06_dp, &
7.187993199450151e-06_dp, &
1.3287667789384625e-05_dp, &
2.435221365340668e-05_dp, &
4.406786544877652e-05_dp, &
7.839477316844559e-05_dp, &
0.00013700432533564813_dp, &
0.00023541888373746023_dp, &
0.0003982417143724247_dp, &
0.0006640223423546365_dp, &
0.0010925276402432504_dp, &
0.0017755203327472222_dp, &
0.0028526152739405167_dp, &
0.0045344377578051915_dp, &
0.007136218994535957_dp, &
0.011126225466246156_dp, &
0.01719515628222213_dp, &
0.026355033830470446_dp, &
0.04007940751399853_dp, &
0.06050127009572156_dp, &
0.09069160330701329_dp, &
0.13505126572297815_dp, &
0.19986523979796597_dp, &
0.2940997147198093_dp, &
0.4305938358212778_dp, &
0.6279833290723216_dp, &
0.9142278655051918_dp, &
1.3343266294937697_dp, &
1.971223529289826_dp, &
3.0207289388880563_dp, &
5.239269717116346_dp /)
ELSE IF(E_range < 2582180.0_dp) THEN
aw(:) = (/ &
2.6346628583772034e-07_dp, &
1.4870682611802608e-06_dp, &
4.18499777462391e-06_dp, &
9.565754079998203e-06_dp, &
2.0095438915823317e-05_dp, &
4.030712638205273e-05_dp, &
7.820310873891258e-05_dp, &
0.0001476189063429688_dp, &
0.00027206125557350966_dp, &
0.0004908276742423054_dp, &
0.0008686422030538531_dp, &
0.0015106576888148007_dp, &
0.0025855700579017434_dp, &
0.0043608859416256075_dp, &
0.007256247781115456_dp, &
0.011923381221441286_dp, &
0.01936500732097008_dp, &
0.031110395247887928_dp, &
0.049472724953768164_dp, &
0.07792391423249093_dp, &
0.12163719777397036_dp, &
0.18826820696011917_dp, &
0.2890741916661454_dp, &
0.44051299024821017_dp, &
0.6665282145336267_dp, &
1.0018392157983365_dp, &
1.4967815518731842_dp, &
2.2247902473537535_dp, &
3.295151802440413_dp, &
4.87866127574686_dp, &
7.274154425483326_dp, &
11.170031678415262_dp, &
6.865691087424058e-07_dp, &
1.8349687179642652e-06_dp, &
3.75363535847515e-06_dp, &
7.406302503889748e-06_dp, &
1.4385998863803715e-05_dp, &
2.7332856632656454e-05_dp, &
5.070293872960844e-05_dp, &
9.194731882902252e-05_dp, &
0.00016333066755487376_dp, &
0.0002847458523672339_dp, &
0.0004880328440654819_dp, &
0.0008235400684252141_dp, &
0.001370000313997599_dp, &
0.002249268185463258_dp, &
0.003648135807576948_dp, &
0.005850384133560756_dp, &
0.009283540163912576_dp, &
0.01458663604211278_dp, &
0.022707798579110823_dp, &
0.03504400012572891_dp, &
0.05364018073213217_dp, &
0.08147187176803067_dp, &
0.12284570471579664_dp, &
0.18396885882280342_dp, &
0.27376977104743744_dp, &
0.4051220237008429_dp, &
0.5968036286827566_dp, &
0.8770477402799991_dp, &
1.2912341021707199_dp, &
1.9226515113345246_dp, &
2.9669527288987725_dp, &
5.177596037459169_dp /)
ELSE IF(E_range < 10786426.0_dp) THEN
aw(:) = (/ &
7.696764691509552e-08_dp, &
4.5183075501644766e-07_dp, &
1.3712733197245156e-06_dp, &
3.452398063116915e-06_dp, &
7.99119238288884e-06_dp, &
1.7499890158507848e-05_dp, &
3.668225172908821e-05_dp, &
7.411367023252019e-05_dp, &
0.0001450609740104884_dp, &
0.0002761333471939282_dp, &
0.0005128366382617873_dp, &
0.0009316733536716902_dp, &
0.0016592776527650684_dp, &
0.0029023227583525545_dp, &
0.004993762465329047_dp, &
0.008463615785356367_dp, &
0.014146318917564842_dp, &
0.02334212828021178_dp, &
0.03805782856888779_dp, &
0.0613630042242151_dp, &
0.0979136503362365_dp, &
0.15471676332581316_dp, &
0.24224053787692895_dp, &
0.3760196167953944_dp, &
0.5789730031642062_dp, &
0.8847668378302425_dp, &
1.3427798213464026_dp, &
2.0257631969638736_dp, &
3.042786371879784_dp, &
4.564995894383141_dp, &
6.891880055654338_dp, &
10.711074682469496_dp, &
2.0224960490304485e-07_dp, &
5.837519191389681e-07_dp, &
1.3534097826701593e-06_dp, &
3.0163324999016133e-06_dp, &
6.464297707355912e-06_dp, &
1.331006273228438e-05_dp, &
2.6438493784462622e-05_dp, &
5.089377871159909e-05_dp, &
9.53135908336066e-05_dp, &
0.00017421956824809404_dp, &
0.0003116271902197154_dp, &
0.0005466663994455823_dp, &
0.0009422361350401303_dp, &
0.0015981946213501235_dp, &
0.0026712728244169207_dp, &
0.004404875333973468_dp, &
0.00717331310959028_dp, &
0.011546952692401411_dp, &
0.018387490277830258_dp, &
0.028986360940377882_dp, &
0.04526460807725196_dp, &
0.07006006488175816_dp, &
0.10753871373789183_dp, &
0.16378446149863707_dp, &
0.24765288245177752_dp, &
0.3720420890669611_dp, &
0.5559074354965244_dp, &
0.8278548563568215_dp, &
1.2338117262467065_dp, &
1.8576021658437778_dp, &
2.8947650724141107_dp, &
5.094894347251689_dp /)
ELSE IF(E_range < 72565710.0_dp) THEN
aw(:) = (/ &
1.5580960122828055e-08_dp, &
1.0192299261551187e-07_dp, &
3.657511654714275e-07_dp, &
1.0874230387694592e-06_dp, &
2.8971440093678986e-06_dp, &
7.1231112287328e-06_dp, &
1.6440376252868036e-05_dp, &
3.6035911594207885e-05_dp, &
7.564822282610485e-05_dp, &
0.0001530604093750471_dp, &
0.0002999674552596099_dp, &
0.0005716626169231291_dp, &
0.0010627819854068627_dp, &
0.0019325381417571955_dp, &
0.0034446462658281212_dp, &
0.006029758821968773_dp, &
0.010382043759726605_dp, &
0.017607085158632928_dp, &
0.02944628039743525_dp, &
0.04861435988540585_dp, &
0.07930297685700823_dp, &
0.12792652074627453_dp, &
0.20421936908257549_dp, &
0.3228415686086713_dp, &
0.5057217947095013_dp, &
0.785484422991801_dp, &
1.2105330871831885_dp, &
1.852886865611021_dp, &
2.8213315422092657_dp, &
4.2873111269979765_dp, &
6.550982976987269_dp, &
10.299461972660755_dp, &
4.192785461816036e-08_dp, &
1.4682586910958697e-07_dp, &
4.248767433871547e-07_dp, &
1.117972412187448e-06_dp, &
2.71251708562224e-06_dp, &
6.166664608205133e-06_dp, &
1.3299698724334126e-05_dp, &
2.7458360437668858e-05_dp, &
5.463565964132874e-05_dp, &
0.00010531494807998692_dp, &
0.0001974602213268062_dp, &
0.0003612942851684024_dp, &
0.0006468325897487988_dp, &
0.001135618706607295_dp, &
0.0019588032650003236_dp, &
0.003324710089253406_dp, &
0.005560469162003084_dp, &
0.009174344598623402_dp, &
0.014948285309521006_dp, &
0.024074306780617383_dp, &
0.03835406948576763_dp, &
0.060489181968959096_dp, &
0.0945015966761899_dp, &
0.1463417320393485_dp, &
0.22477371348415012_dp, &
0.34269334477792374_dp, &
0.5192018690921564_dp, &
0.7832512032052583_dp, &
1.1813099794431254_dp, &
1.7977731364021527_dp, &
2.828182155100705_dp, &
5.018697893750251_dp /)
ELSE
aw(:) = (/ &
1.4796558719757977e-09_dp, &
1.6277762282636202e-08_dp, &
8.615102100471783e-08_dp, &
3.3251166768869524e-07_dp, &
1.064831343891703e-06_dp, &
3.0042406722730696e-06_dp, &
7.722404719947578e-06_dp, &
1.8469628120696557e-05_dp, &
4.168720676495336e-05_dp, &
8.969407716668371e-05_dp, &
0.00018534670944540425_dp, &
0.00036995434452620024_dp, &
0.0007164776276992899_dp, &
0.0013511756947926401_dp, &
0.0024885840769342527_dp, &
0.004487272639985271_dp, &
0.007937634518341114_dp, &
0.013798542180092185_dp, &
0.023607852053926744_dp, &
0.03980353161167122_dp, &
0.06620914631345257_dp, &
0.10876176019481658_dp, &
0.17659514989891748_dp, &
0.28364164912494744_dp, &
0.45099127302133024_dp, &
0.7103683461225575_dp, &
1.109312567732835_dp, &
1.7191679096665389_dp, &
2.648421848786194_dp, &
4.06873544538514_dp, &
6.280858212318488_dp, &
9.971616727521422_dp, &
4.7301698039073534e-09_dp, &
3.1078236876760897e-08_dp, &
1.2722293347158736e-07_dp, &
4.135773629354279e-07_dp, &
1.1643319139692173e-06_dp, &
2.9633030798047164e-06_dp, &
6.990659700961397e-06_dp, &
1.5534394085699292e-05_dp, &
3.28792052583233e-05_dp, &
6.681670495514343e-05_dp, &
0.00013116076656301764_dp, &
0.0002498651223560153_dp, &
0.0004636589258144925_dp, &
0.0008405899279454757_dp, &
0.0014925644239063773_dp, &
0.0026009920182698686_dp, &
0.004456120959091793_dp, &
0.00751677674389606_dp, &
0.012500261250050913_dp, &
0.020516489982978243_dp, &
0.03326657617193578_dp, &
0.053334802223457274_dp, &
0.08461553683522265_dp, &
0.13293582624456382_dp, &
0.20696680596446562_dp, &
0.31958257452387945_dp, &
0.4899879264178536_dp, &
0.7474157763110646_dp, &
1.1388011986618074_dp, &
1.7490627041810483_dp, &
2.7738272660229994_dp, &
4.9565545529977175_dp /)
END IF
END SUBROUTINE


SUBROUTINE get_coeff_34(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio
E_ratio = 1.0_dp
IF(E_range < 9649.0_dp) THEN
E_ratio = 9649.0_dp/E_range
aw(:) = (/ &
3.284029233271286e-05_dp, &
0.0001753443602861097_dp, &
0.0004416382833363001_dp, &
0.0008508530310865113_dp, &
0.001435500605855746_dp, &
0.0022474948219129297_dp, &
0.003366173371324329_dp, &
0.004907942966760921_dp, &
0.007037748206337417_dp, &
0.009983831117905245_dp, &
0.014057973412231587_dp, &
0.019683414844193538_dp, &
0.02743271515183254_dp, &
0.03807838736551588_dp, &
0.052660091008446035_dp, &
0.07257342975355677_dp, &
0.09968697353631621_dp, &
0.13649612296545818_dp, &
0.18632500099011615_dp, &
0.2535908886542759_dp, &
0.3441500984275648_dp, &
0.46575004766668_dp, &
0.6286204407593304_dp, &
0.8462483877785438_dp, &
1.1364010501626234_dp, &
1.5224916494827863_dp, &
2.035445381255887_dp, &
2.716346186249304_dp, &
3.620421040738206_dp, &
4.823581438788581_dp, &
6.434509907627798_dp, &
8.620766164346367_dp, &
11.679090118639257_dp, &
16.31432816585662_dp, &
8.4538052355896e-05_dp, &
0.00020209318539166847_dp, &
0.00033361605486797254_dp, &
0.0004901319572083779_dp, &
0.0006876803644732437_dp, &
0.0009493703556928221_dp, &
0.0013072005549165328_dp, &
0.0018035890085942577_dp, &
0.0024937295684409103_dp, &
0.0034499397324379713_dp, &
0.00476826036404788_dp, &
0.0065772141488150675_dp, &
0.009049030640506878_dp, &
0.012414126010019873_dp, &
0.016979951815415794_dp, &
0.02315561592334051_dp, &
0.031484041042937_dp, &
0.04268391545844468_dp, &
0.05770434764295233_dp, &
0.07779602216248598_dp, &
0.10460388481675117_dp, &
0.14028819232415718_dp, &
0.18768361957484037_dp, &
0.2505110517811862_dp, &
0.33366599536844854_dp, &
0.4436264870131247_dp, &
0.5890644515470643_dp, &
0.7818386670916874_dp, &
1.0387776399865274_dp, &
1.3852730956246349_dp, &
1.8635541689521873_dp, &
2.5552790973254864_dp, &
3.661647560244881_dp, &
5.979097365202448_dp /)
ELSE IF(E_range < 15161.0_dp) THEN
aw(:) = (/ &
2.2094687262226044e-05_dp, &
0.00011816317101145329_dp, &
0.00029854362986100777_dp, &
0.0005779739014899124_dp, &
0.0009819112670486493_dp, &
0.0015516034173501654_dp, &
0.0023507187818878183_dp, &
0.0034732685500535195_dp, &
0.00505334279364789_dp, &
0.007278229326144087_dp, &
0.010406819589614869_dp, &
0.014795161155330522_dp, &
0.020931286923915934_dp, &
0.029482154546951286_dp, &
0.04135653654321724_dp, &
0.05778898092867048_dp, &
0.08045158704880996_dp, &
0.11160243319073283_dp, &
0.15428220578421106_dp, &
0.2125741178044686_dp, &
0.2919468573365305_dp, &
0.3997065271731107_dp, &
0.545592087818525_dp, &
0.7425611332686428_dp, &
1.0078318394632224_dp, &
1.3642789958888877_dp, &
1.8423416665523202_dp, &
2.48272174251732_dp, &
3.34042303033154_dp, &
4.49133411535853_dp, &
6.044309582675748_dp, &
8.167153038233048_dp, &
11.156471603344267_dp, &
15.714736908294974_dp, &
5.689777248462453e-05_dp, &
0.0001364647878286462_dp, &
0.0002266959149215116_dp, &
0.0003363227347021138_dp, &
0.0004783149170860289_dp, &
0.0006715498663212837_dp, &
0.0009421913049362618_dp, &
0.0013250426248297155_dp, &
0.0018660179676916586_dp, &
0.00262648268480602_dp, &
0.003689465255982615_dp, &
0.005167783411766887_dp, &
0.007214591621086576_dp, &
0.010037222243044243_dp, &
0.013915457617770334_dp, &
0.01922566776214666_dp, &
0.02647264919678102_dp, &
0.03633153858613665_dp, &
0.049702887216291285_dp, &
0.06778493328431712_dp, &
0.09216841458681874_dp, &
0.12496114099205212_dp, &
0.1689524368948172_dp, &
0.22783243637559455_dp, &
0.3064902609254688_dp, &
0.411433435566918_dp, &
0.5514107558972641_dp, &
0.7384129467740507_dp, &
0.9894523508058952_dp, &
1.3301301240547123_dp, &
1.8028378573543902_dp, &
2.489168614889156_dp, &
3.5894672490279578_dp, &
5.895310588871672_dp /)
ELSE IF(E_range < 29986.0_dp) THEN
aw(:) = (/ &
1.2132019239369146e-05_dp, &
6.506554392462203e-05_dp, &
0.00016528151320033429_dp, &
0.00032271987439700723_dp, &
0.0005549904907647636_dp, &
0.0008912186380916402_dp, &
0.0013769333304078336_dp, &
0.00207998381188429_dp, &
0.0030983248733590404_dp, &
0.004571097441508312_dp, &
0.006694456763746173_dp, &
0.009743704668641925_dp, &
0.0141037713916676_dp, &
0.02031087384545962_dp, &
0.029109193553804015_dp, &
0.04152772200656733_dp, &
0.05898411878582171_dp, &
0.08342465247911345_dp, &
0.11751221427114543_dp, &
0.1648782324800013_dp, &
0.23045938623741916_dp, &
0.3209467823295079_dp, &
0.44538447707509593_dp, &
0.6159672665711226_dp, &
0.8491072823533072_dp, &
1.1668710896432313_dp, &
1.5989474740844776_dp, &
2.1854242653839444_dp, &
2.9809152725739025_dp, &
4.061217864274067_dp, &
5.535407811609172_dp, &
7.571672639611632_dp, &
10.466564989305379_dp, &
14.91952959384994_dp, &
3.126217597574477e-05_dp, &
7.540587737410632e-05_dp, &
0.00012663926218204873_dp, &
0.00019110769108438634_dp, &
0.00027820581560489743_dp, &
0.000401742249839652_dp, &
0.0005808693019056328_dp, &
0.0008413761481210849_dp, &
0.0012182358494751703_dp, &
0.001759637722089599_dp, &
0.00253246317609668_dp, &
0.003629502407336348_dp, &
0.00517906799440672_dp, &
0.0073579025027302025_dp, &
0.01040852585005309_dp, &
0.014662503493069159_dp, &
0.020571571811010672_dp, &
0.028749161300949972_dp, &
0.04002565337439286_dp, &
0.05552176300115481_dp, &
0.07674587007131725_dp, &
0.10572312939601615_dp, &
0.14516718239506418_dp, &
0.19871016099008434_dp, &
0.2712154487579484_dp, &
0.3692151806287191_dp, &
0.5015525805409441_dp, &
0.6803982020135383_dp, &
0.9230308920142081_dp, &
1.2553734442165063_dp, &
1.7200955901449364_dp, &
2.3987730554757754_dp, &
3.4906629653715_dp, &
5.780815874478559_dp /)
ELSE IF(E_range < 49196.0_dp) THEN
aw(:) = (/ &
7.845531705064765e-06_dp, &
4.217579566703759e-05_dp, &
0.0001076255406955941_dp, &
0.0002116704074840506_dp, &
0.0003678026352786526_dp, &
0.000598637239972256_dp, &
0.000939800791185946_dp, &
0.0014447771966513766_dp, &
0.002191627535252642_dp, &
0.0032927972200630018_dp, &
0.004909216880554122_dp, &
0.007270144108734435_dp, &
0.010700753084513917_dp, &
0.015660261119126325_dp, &
0.02279439226673758_dp, &
0.033007300151380946_dp, &
0.04755981744820865_dp, &
0.06820320985515761_dp, &
0.09736066855556012_dp, &
0.13837281589809153_dp, &
0.19582885861759528_dp, &
0.2760121868265214_dp, &
0.3874989458973457_dp, &
0.5419617440352286_dp, &
0.7552508475941064_dp, &
1.048857678610048_dp, &
1.4519234318259462_dp, &
2.004071631268318_dp, &
2.759600998924181_dp, &
3.7942044073242056_dp, &
5.217087596451319_dp, &
7.196718451225344_dp, &
10.029695823161129_dp, &
14.413626849760721_dp, &
2.022742080929783e-05_dp, &
4.902127023830298e-05_dp, &
8.308838409631613e-05_dp, &
0.00012719755797047105_dp, &
0.00018877141608928316_dp, &
0.00027875729779708194_dp, &
0.00041237481957666014_dp, &
0.0006104274314639044_dp, &
0.0009017824485847016_dp, &
0.0013270526990110587_dp, &
0.0019435429478184875_dp, &
0.0028318806888766317_dp, &
0.004105011117474572_dp, &
0.0059204430431561185_dp, &
0.008496895591229079_dp, &
0.012136858469163465_dp, &
0.017257067102487582_dp, &
0.02442954195917895_dp, &
0.03443669419909579_dp, &
0.04834513188831433_dp, &
0.06760433082481605_dp, &
0.09417845588835595_dp, &
0.13072271761835513_dp, &
0.18082056576074482_dp, &
0.2493066841224739_dp, &
0.3427177893721846_dp, &
0.46995019394607707_dp, &
0.6432897947695325_dp, &
0.8801980582788952_dp, &
1.2068306280775711_dp, &
1.6660767737909483_dp, &
2.3395502248176463_dp, &
3.4258517245052684_dp, &
5.7058363778559205_dp /)
ELSE IF(E_range < 109833.0_dp) THEN
aw(:) = (/ &
3.865603619912879e-06_dp, &
2.087430645188353e-05_dp, &
5.3738065099243575e-05_dp, &
0.00010718454630525007_dp, &
0.0001900012319697737_dp, &
0.00031718955789873134_dp, &
0.0005126537329759226_dp, &
0.0008128168395638415_dp, &
0.0012719920833382853_dp, &
0.001970384125983329_dp, &
0.0030256963224515113_dp, &
0.00460969356562721_dp, &
0.006971643847773047_dp, &
0.010471317877795917_dp, &
0.01562522911317097_dp, &
0.023171141484273704_dp, &
0.034157679556524934_dp, &
0.05006830174859498_dp, &
0.0729921440414609_dp, &
0.1058585802344639_dp, &
0.1527581466742932_dp, &
0.2193802673426605_dp, &
0.31360877991591823_dp, &
0.4463309353744743_dp, &
0.6325368438825444_dp, &
0.8928196332263253_dp, &
1.2554443726415525_dp, &
1.759267106284206_dp, &
2.458035348630705_dp, &
3.427215114290294_dp, &
4.77616847618911_dp, &
6.673804689443379_dp, &
9.416907604833742_dp, &
13.70063546514659_dp, &
9.976454765718753e-06_dp, &
2.4397956666150202e-05_dp, &
4.2091633541468636e-05_dp, &
6.622948123983615e-05_dp, &
0.00010186030376456873_dp, &
0.00015644233845744574_dp, &
0.0002404690130604266_dp, &
0.0003687782406776715_dp, &
0.0005627493820533361_dp, &
0.0008533803933831028_dp, &
0.0012854870718438482_dp, &
0.0019235107198793406_dp, &
0.002859592750775737_dp, &
0.004224781787747757_dp, &
0.0062045295634641385_dp, &
0.009060029072853486_dp, &
0.0131574822853581_dp, &
0.019008095244250185_dp, &
0.02732254096178089_dp, &
0.03908488812876366_dp, &
0.055652686745107795_dp, &
0.07889222315477536_dp, &
0.11136127365662019_dp, &
0.15655676737113142_dp, &
0.21925337945032786_dp, &
0.30597553969409247_dp, &
0.42568066454670456_dp, &
0.5908146138910935_dp, &
0.8191119697759165_dp, &
1.137096600179962_dp, &
1.5880319126786595_dp, &
2.2536634913094686_dp, &
3.331728541732483_dp, &
5.597118014202835_dp /)
ELSE IF(E_range < 276208.0_dp) THEN
aw(:) = (/ &
1.7147362511403065e-06_dp, &
9.319802731855972e-06_dp, &
2.4302616807482706e-05_dp, &
4.948261860530359e-05_dp, &
9.026404229610094e-05_dp, &
0.0001560250715046401_dp, &
0.00026193159643869907_dp, &
0.0004316317370034808_dp, &
0.0007014310228413562_dp, &
0.0011265811323026616_dp, &
0.001790538099368303_dp, &
0.002818442710487939_dp, &
0.004396613985562227_dp, &
0.006800572618267542_dp, &
0.01043510374307695_dp, &
0.015891226467704685_dp, &
0.02402679235817231_dp, &
0.03607995947798214_dp, &
0.05382821091103021_dp, &
0.07981021388476216_dp, &
0.11763406509273454_dp, &
0.17240391938168603_dp, &
0.25130849721857534_dp, &
0.36443086121519086_dp, &
0.525861591666699_dp, &
0.7552321345423689_dp, &
1.0798433566303969_dp, &
1.5376756421659512_dp, &
2.1818096583688225_dp, &
3.087384050821325_dp, &
4.363863206811109_dp, &
6.180624421988004_dp, &
8.83477361053772_dp, &
13.019312520088521_dp, &
4.43185996007463e-06_dp, &
1.0980758773489005e-05_dp, &
1.943725634861249e-05_dp, &
3.1800270744211854e-05_dp, &
5.129667367655266e-05_dp, &
8.27095880626353e-05_dp, &
0.0001329766243411858_dp, &
0.00021237174266944332_dp, &
0.0003362840713927107_dp, &
0.0005277050042553842_dp, &
0.0008207467241987074_dp, &
0.001265656823748127_dp, &
0.001935949524221289_dp, &
0.002938498181312597_dp, &
0.004427745609890722_dp, &
0.00662561356009038_dp, &
0.009849268137973679_dp, &
0.01454967083492188_dp, &
0.021364880820575197_dp, &
0.031193465982641635_dp, &
0.04529525929852885_dp, &
0.06542926430703719_dp, &
0.09404211859223464_dp, &
0.13452588222462222_dp, &
0.19157264817324873_dp, &
0.27166955823017896_dp, &
0.3838116245690371_dp, &
0.5405893273969905_dp, &
0.7600124959108089_dp, &
1.0690028998164918_dp, &
1.5112619233396074_dp, &
2.168763561959426_dp, &
3.238506003394989_dp, &
5.489633460260332_dp /)
ELSE IF(E_range < 852991.0_dp) THEN
aw(:) = (/ &
6.35749973501056e-07_dp, &
3.492169089246581e-06_dp, &
9.302356973224121e-06_dp, &
1.9592566972831844e-05_dp, &
3.7377729183938204e-05_dp, &
6.798117415794394e-05_dp, &
0.00012023568563760358_dp, &
0.0002084236297872551_dp, &
0.00035534106272681936_dp, &
0.0005969900768905776_dp, &
0.0009896587375544644_dp, &
0.001620507511517086_dp, &
0.0026232770839453665_dp, &
0.004201426984363522_dp, &
0.006661986770200809_dp, &
0.010464757550535931_dp, &
0.01629338409394917_dp, &
0.025157418466823473_dp, &
0.038538075135755495_dp, &
0.05859528378548497_dp, &
0.088460352930584_dp, &
0.13264771717146925_dp, &
0.19763178422753508_dp, &
0.2926522443111407_dp, &
0.4308357751366649_dp, &
0.6307586891808089_dp, &
0.9186346830441839_dp, &
1.3314221217344053_dp, &
1.9213813610912545_dp, &
2.76319200354934_dp, &
3.9663504041897766_dp, &
5.700732751780247_dp, &
8.263919165830904_dp, &
12.346999008174697_dp, &
1.646994148539095e-06_dp, &
4.168704436960474e-06_dp, &
7.695536879782202e-06_dp, &
1.3374173961844378e-05_dp, &
2.3064311988245763e-05_dp, &
3.958321376425686e-05_dp, &
6.725628425988863e-05_dp, &
0.00011284251346096221_dp, &
0.0001868773135187132_dp, &
0.0003056293856662244_dp, &
0.0004939782736171523_dp, &
0.0007896314087630558_dp, &
0.0012492643080895956_dp, &
0.001957402304692391_dp, &
0.0030391875844591482_dp, &
0.004678623430285177_dp, &
0.007144500300310407_dp, &
0.010827043576140467_dp, &
0.016289458703073287_dp, &
0.024340088611569204_dp, &
0.03613298947254673_dp, &
0.05330759009262987_dp, &
0.07818208570144934_dp, &
0.11402099109647866_dp, &
0.16540631928634394_dp, &
0.23875782212611196_dp, &
0.3430802599354601_dp, &
0.49109150628018977_dp, &
0.7010818269759505_dp, &
1.0004112023843146_dp, &
1.4333015739126787_dp, &
2.0820719815510036_dp, &
3.1430997613534135_dp, &
5.379828456565977_dp /)
ELSE
aw(:) = (/ &
1.117489994479298e-09_dp, &
1.0437767934573854e-08_dp, &
5.084528505414382e-08_dp, &
1.872213850467172e-07_dp, &
5.813350591395648e-07_dp, &
1.6041946888512065e-06_dp, &
4.054564459005386e-06_dp, &
9.568204612679587e-06_dp, &
2.1360561665920905e-05_dp, &
4.553900601252074e-05_dp, &
9.336828090817695e-05_dp, &
0.00018510270750313838_dp, &
0.0003563527131882449_dp, &
0.0006684922084643366_dp, &
0.0012254246986170441_dp, &
0.002200236453707954_dp, &
0.0038770557829877106_dp, &
0.006716063190006149_dp, &
0.011453417554768287_dp, &
0.019253386304452922_dp, &
0.0319378960087355_dp, &
0.052330033425615166_dp, &
0.0847641054822857_dp, &
0.13583764990140365_dp, &
0.21551319196591645_dp, &
0.3387243798072503_dp, &
0.527711765290698_dp, &
0.8154303015314811_dp, &
1.2505960123470905_dp, &
1.9054661547111273_dp, &
2.8889237961629872_dp, &
4.37232400643935_dp, &
6.65560984092506_dp, &
10.426038695834771_dp, &
3.3457340053937633e-09_dp, &
1.873386519501546e-08_dp, &
7.20396391424945e-08_dp, &
2.259124367219339e-07_dp, &
6.20700398368113e-07_dp, &
1.5516519034576341e-06_dp, &
3.6098773524770243e-06_dp, &
7.932321947059267e-06_dp, &
1.6633942623129364e-05_dp, &
3.353876533808818e-05_dp, &
6.539287755778778e-05_dp, &
0.00012384295496778607_dp, &
0.00022861478182158466_dp, &
0.0004125495803926801_dp, &
0.0007294839100515848_dp, &
0.0012664381661580966_dp, &
0.002162272032202245_dp, &
0.003635959152010226_dp, &
0.006029053426126615_dp, &
0.009868936147464018_dp, &
0.015962272528537477_dp, &
0.025532104472748035_dp, &
0.040417626126145245_dp, &
0.06336363493503006_dp, &
0.09843823938908386_dp, &
0.1516353587304579_dp, &
0.23175010955983177_dp, &
0.35168177852781385_dp, &
0.5304884039781542_dp, &
0.7970145382697631_dp, &
1.197557199061726_dp, &
1.8163260289376923_dp, &
2.848850144582029_dp, &
5.042342061181277_dp /)
END IF
END SUBROUTINE

SUBROUTINE get_coeff_40(k,E_range, aw, E_ratio)
INTEGER, INTENT(IN)                                :: k
REAL(KIND=dp), INTENT(IN)                          :: E_range
REAL(KIND=dp), DIMENSION(2*k), INTENT(OUT)         :: aw
REAL(KIND=dp), INTENT(OUT)                          :: E_ratio

E_ratio = 1.0_dp

IF(E_range < 57856.0_dp) THEN
E_ratio = 57856.0_dp/E_range
aw(:) = (/ &
5.625164703631762e-06_dp, &
3.0057314707929964e-05_dp, &
7.581456715955565e-05_dp, &
0.00014639299371835805_dp, &
0.0002477800293143379_dp, &
0.0003896039386773891_dp, &
0.0005866518866042211_dp, &
0.0008606874916545427_dp, &
0.0012426428371535472_dp, &
0.0017755036713625403_dp, &
0.002518314919509984_dp, &
0.003551724759050896_dp, &
0.004985520405439102_dp, &
0.006968741207146314_dp, &
0.009703159185474989_dp, &
0.013461179169270084_dp, &
0.018609540533758993_dp, &
0.025640622632516107_dp, &
0.03521369441558543_dp, &
0.048209139752427646_dp, &
0.06579957734629878_dp, &
0.08954293300354034_dp, &
0.12150398353288118_dp, &
0.1644127688907646_dp, &
0.22187068662660148_dp, &
0.29861821954053525_dp, &
0.40088237881728134_dp, &
0.536827520003883_dp, &
0.717141013310299_dp, &
0.9557968845548823_dp, &
1.2710591975583654_dp, &
1.6868195160872768_dp, &
2.2344246534909913_dp, &
2.9552779389928863_dp, &
3.9047778130295248_dp, &
5.1588289603459465_dp, &
6.825944473288935_dp, &
9.073479881127453_dp, &
12.198374786218169_dp, &
16.907889055332152_dp, &
1.4482915271757386e-05_dp, &
3.467519102649355e-05_dp, &
5.740919321142222e-05_dp, &
8.472639971765309e-05_dp, &
0.00011962680279675794_dp, &
0.000166455373763227_dp, &
0.00023123210315609694_dp, &
0.00032194744189025393_dp, &
0.00044906174589933476_dp, &
0.0006264093770111557_dp, &
0.0008725302537805241_dp, &
0.0012124206033923789_dp, &
0.0016797883840064898_dp, &
0.002319986809094753_dp, &
0.0031938591674912755_dp, &
0.004382787784909753_dp, &
0.005995317816762346_dp, &
0.008175830015706832_dp, &
0.011115870779351075_dp, &
0.015068920499380183_dp, &
0.02036960303919392_dp, &
0.027458623009814043_dp, &
0.036915082633833536_dp, &
0.049498301382351756_dp, &
0.06620187932042156_dp, &
0.08832358165406289_dp, &
0.11755580175795301_dp, &
0.15610312371837842_dp, &
0.20683635840053374_dp, &
0.27349745289637084_dp, &
0.36097925980218354_dp, &
0.4757237047702205_dp, &
0.6263240728164432_dp, &
0.8245131573104876_dp, &
1.0869519279807407_dp, &
1.4388502802776517_dp, &
1.9223110911786352_dp, &
2.6190934773994248_dp, &
3.7312648401424546_dp, &
6.060020766373581_dp /)
ELSE IF(E_range < 79311.0_dp) THEN
aw(:) = (/ &
4.250536980156458e-06_dp, &
2.273654241848127e-05_dp, &
5.746653135010074e-05_dp, &
0.00011132024540199449_dp, &
0.00018928142122576544_dp, &
0.0002994382569938377_dp, &
0.000454288627622533_dp, &
0.000672299503582154_dp, &
0.000979832246004811_dp, &
0.0014137495720868076_dp, &
0.002025080331686625_dp, &
0.002884110197625469_dp, &
0.004087325104251008_dp, &
0.005766778037762526_dp, &
0.008102652786456355_dp, &
0.011340055969581531_dp, &
0.015811395964487288_dp, &
0.02196612810440917_dp, &
0.030410189283770143_dp, &
0.041958146879091635_dp, &
0.057701991689176224_dp, &
0.07910167157621399_dp, &
0.10810396603988766_dp, &
0.14729823857007074_dp, &
0.20012010495068866_dp, &
0.2711173056417601_dp, &
0.3662963437398923_dp, &
0.49357419553441345_dp, &
0.6633674109186503_dp, &
0.8893627154039357_dp, &
1.1895319275941671_dp, &
1.587486363643954_dp, &
2.1143270663453757_dp, &
2.811272663465346_dp, &
3.733624369137816_dp, &
4.957294124671807_dp, &
6.590896278962817_dp, &
8.801904457294249_dp, &
11.88713123840321_dp, &
16.55238483684307_dp, &
1.094638989809278e-05_dp, &
2.6264541167916206e-05_dp, &
4.366423701624123e-05_dp, &
6.485707900640055e-05_dp, &
9.239132118095284e-05_dp, &
0.00012997958448661634_dp, &
0.00018276773839222705_dp, &
0.00025760504750009657_dp, &
0.00036354398051613267_dp, &
0.0005127116159650827_dp, &
0.0007215497957439241_dp, &
0.0010124366335163807_dp, &
0.001415795043352708_dp, &
0.0019728647860616224_dp, &
0.002739367020466087_dp, &
0.003790350458028016_dp, &
0.005226588633300203_dp, &
0.007183005233622441_dp, &
0.009839743808395621_dp, &
0.01343667650026383_dp, &
0.01829237662590928_dp, &
0.024828876515768098_dp, &
0.0336039121463881_dp, &
0.045352847970262476_dp, &
0.06104312125664014_dp, &
0.08194491101177495_dp, &
0.10972294540871524_dp, &
0.14655615112533493_dp, &
0.19529470168413016_dp, &
0.2596689868212105_dp, &
0.3445744387802778_dp, &
0.4564753463949847_dp, &
0.6040123174743215_dp, &
0.798993008358502_dp, &
1.0581770447413485_dp, &
1.40688033743651_dp, &
1.887277553918721_dp, &
2.5810634779601025_dp, &
3.689783760886334_dp, &
6.01179163953066_dp /)
ELSE IF(E_range < 115761.0_dp) THEN
aw(:) = (/ &
3.0368483523366833e-06_dp, &
1.6266974194706644e-05_dp, &
4.12242033669947e-05_dp, &
8.019166810373411e-05_dp, &
0.0001371714648132165_dp, &
0.000218728028885675_dp, &
0.00033507234961707125_dp, &
0.0005013623526533384_dp, &
0.000739363225101565_dp, &
0.0010797635446663115_dp, &
0.0015654649769365116_dp, &
0.0022561684161992524_dp, &
0.0032346603912939064_dp, &
0.004615352375884511_dp, &
0.006555823716919149_dp, &
0.009272370969787912_dp, &
0.013060890951248518_dp, &
0.018324846273714744_dp, &
0.025612610229014288_dp, &
0.03566720010832204_dp, &
0.04949233241879972_dp, &
0.0684399313263429_dp, &
0.094325771015324_dp, &
0.1295819368053425_dp, &
0.17745738870935_dp, &
0.24228129488458383_dp, &
0.32980824883053056_dp, &
0.4476704404060215_dp, &
0.6059700961848395_dp, &
0.8180575196842615_dp, &
1.1015588777601752_dp, &
1.4797500617867587_dp, &
1.983433359686095_dp, &
2.6535973719065566_dp, &
3.5454149361674943_dp, &
4.734794435293425_dp, &
6.330463734446495_dp, &
8.500042905168106_dp, &
11.540228533157341_dp, &
16.15522809507984_dp, &
7.82325854928367e-06_dp, &
1.8823444129083424e-05_dp, &
3.14620780748313e-05_dp, &
4.7125707594704645e-05_dp, &
6.79097312124436e-05_dp, &
9.688391441936373e-05_dp, &
0.0001383041584752257_dp, &
0.0001978655315728085_dp, &
0.0002831968772318314_dp, &
0.00040467853946619345_dp, &
0.000576570999057209_dp, &
0.0008184930792544974_dp, &
0.0011573694269334603_dp, &
0.0016300217332690654_dp, &
0.0022866264378665017_dp, &
0.003195323516104103_dp, &
0.00444834488212991_dp, &
0.006170141680702589_dp, &
0.008528133604565441_dp, &
0.011746889493182871_dp, &
0.01612678889492019_dp, &
0.022068523492347825_dp, &
0.030105195225203187_dp, &
0.04094428573089771_dp, &
0.05552245131035983_dp, &
0.07507700169834372_dp, &
0.10123916928428976_dp, &
0.1361560989005508_dp, &
0.18265135286102935_dp, &
0.24443864210663185_dp, &
0.3264127277785362_dp, &
0.4350602162090508_dp, &
0.5790737277050019_dp, &
0.7703465996268644_dp, &
1.0257546136519846_dp, &
1.3707423427424361_dp, &
1.8475782498699762_dp, &
2.5378991422868293_dp, &
3.6426749233799782_dp, &
5.957056861015189_dp /)
ELSE IF(E_range < 1000000.0_dp) THEN
aw(:) = (/ &
4.4572138185844315e-07_dp, &
2.4133766193055494e-06_dp, &
6.245927429436962e-06_dp, &
1.2564391587696646e-05_dp, &
2.253962411588127e-05_dp, &
3.818611134369075e-05_dp, &
6.273353100197593e-05_dp, &
0.00010114965717706532_dp, &
0.00016093543664301592_dp, &
0.0002533158064145081_dp, &
0.0003949772257419352_dp, &
0.000610568462213574_dp, &
0.0009362738069429562_dp, &
0.0014248908689701706_dp, &
0.0021530104668083796_dp, &
0.0032311199139085753_dp, &
0.004817754268156738_dp, &
0.00713922937782429_dp, &
0.010517040789356705_dp, &
0.015405750432827419_dp, &
0.0224451686927125_dp, &
0.03253195347335625_dp, &
0.04691749468897535_dp, &
0.06734126646359322_dp, &
0.09621189177533797_dp, &
0.13685221917698803_dp, &
0.19383009136320542_dp, &
0.273403673446551_dp, &
0.38411996989021985_dp, &
0.5376188379394832_dp, &
0.7497150344788303_dp, &
1.0418633156365826_dp, &
1.4431695917845953_dp, &
1.9932270233242353_dp, &
2.7463131720403404_dp, &
3.778113074292076_dp, &
5.197839840908659_dp, &
7.173979132328904_dp, &
10.003134479263514_dp, &
14.382802880300657_dp, &
1.1510235486502786e-06_dp, &
2.8301629259030296e-06_dp, &
4.934830570408326e-06_dp, &
7.89229075610394e-06_dp, &
1.238793054403981e-05_dp, &
1.9434921643028237e-05_dp, &
3.047490419602064e-05_dp, &
4.75880028665218e-05_dp, &
7.382592651342797e-05_dp, &
0.00011367493434526954_dp, &
0.0001736974909398471_dp, &
0.0002634324672190349_dp, &
0.0003966600831073136_dp, &
0.000593173425371456_dp, &
0.0008812480653638947_dp, &
0.001301069342806644_dp, &
0.0019094686473009853_dp, &
0.0027864416556577327_dp, &
0.004044083812124213_dp, &
0.005838793353305664_dp, &
0.008387876312088368_dp, &
0.0119920644645412_dp, &
0.01706595086938794_dp, &
0.024178996150197186_dp, &
0.03411061750535862_dp, &
0.047924012774613556_dp, &
0.06706490596801937_dp, &
0.09349352702708642_dp, &
0.12986124903897087_dp, &
0.1797482296944362_dp, &
0.2479870531016808_dp, &
0.34111438047292614_dp, &
0.4680295147606695_dp, &
0.641025338185714_dp, &
0.8775747301933432_dp, &
1.2038482306754261_dp, &
1.6627496029710054_dp, &
2.335896196568727_dp, &
3.4218496429649594_dp, &
5.701207727799133_dp /)
ELSE IF(E_range < 10000000.0_dp) THEN
aw(:) = (/ &
5.784905699357422e-08_dp, &
3.197847227059498e-07_dp, &
8.628534718508499e-07_dp, &
1.8539643437757063e-06_dp, &
3.626459808749055e-06_dp, &
6.774345618312585e-06_dp, &
1.2298850971534939e-05_dp, &
2.18500237667714e-05_dp, &
3.811067585632014e-05_dp, &
6.539061815055223e-05_dp, &
0.00011053458940373901_dp, &
0.00018429657222883774_dp, &
0.000303402625870388_dp, &
0.0004936226832917519_dp, &
0.0007943104359834751_dp, &
0.0012650654868938422_dp, &
0.0019954446725526237_dp, &
0.0031190286332746543_dp, &
0.004833674445421699_dp, &
0.007430508568658357_dp, &
0.011335207791487448_dp, &
0.01716647253714365_dp, &
0.02581844261866524_dp, &
0.038576311913552645_dp, &
0.05727778948423146_dp, &
0.08453762482429264_dp, &
0.1240585731633163_dp, &
0.18106049481027023_dp, &
0.2628705900885491_dp, &
0.37973340793466137_dp, &
0.5459216887422058_dp, &
0.7812634311720118_dp, &
1.1132586814602379_dp, &
1.5800711773747709_dp, &
2.2349241732379_dp, &
3.153030259662271_dp, &
4.4438384827116755_dp, &
6.276631233785271_dp, &
8.9484392918494_dp, &
13.152670324501127_dp, &
1.5007354710371465e-07_dp, &
3.847240168397885e-07_dp, &
7.280525401291225e-07_dp, &
1.3082355183250643e-06_dp, &
2.333420104749775e-06_dp, &
4.125299172286422e-06_dp, &
7.1924772969395586e-06_dp, &
1.2347340911211979e-05_dp, &
2.0877488569818263e-05_dp, &
3.480079813183638e-05_dp, &
5.7245841622678835e-05_dp, &
9.301567501841688e-05_dp, &
0.0001494172836813467_dp, &
0.0002374729092417609_dp, &
0.0003736770608393821_dp, &
0.000582529019240851_dp, &
0.0009001612240169092_dp, &
0.0013795079661119084_dp, &
0.002097628205440898_dp, &
0.0031660274306166227_dp, &
0.00474513584875455_dp, &
0.0070645196874516295_dp, &
0.010450968811914_dp, &
0.015367367350655144_dp, &
0.022466270926174387_dp, &
0.03266347470659499_dp, &
0.04723869852544491_dp, &
0.06797303336386028_dp, &
0.0973363343361155_dp, &
0.1387430350515734_dp, &
0.19690354767411897_dp, &
0.2783145616789539_dp, &
0.3919656593459528_dp, &
0.5504199646285972_dp, &
0.7716325358298327_dp, &
1.0824436400415847_dp, &
1.5264622057751107_dp, &
2.1856086044013003_dp, &
3.2570177843215506_dp, &
5.510961573974098_dp /)
ELSE IF(E_range < 100000000.0_dp) THEN
aw(:) = (/ &
7.70640846068886e-09_dp, &
4.448258917490421e-08_dp, &
1.3080762276377026e-07_dp, &
3.169523848586334e-07_dp, &
7.072220328449342e-07_dp, &
1.4995800728437142e-06_dp, &
3.056823036796932e-06_dp, &
6.028078462050908e-06_dp, &
1.1550200941161362e-05_dp, &
2.1576018875684373e-05_dp, &
3.9401792476965267e-05_dp, &
7.050386217198551e-05_dp, &
0.00012384994957978493_dp, &
0.00021393184360306275_dp, &
0.0003638844302990059_dp, &
0.0006102265589083663_dp, &
0.0010100037090592896_dp, &
0.0016514606584494686_dp, &
0.0026698656201440464_dp, &
0.004270802074829893_dp, &
0.00676421807964012_dp, &
0.010613879617228771_dp, &
0.01650875608923568_dp, &
0.025465464296727323_dp, &
0.03897447084393801_dp, &
0.0592076527270105_dp, &
0.08931150904420052_dp, &
0.13381945117520944_dp, &
0.19922910504512237_dp, &
0.2948078527876033_dp, &
0.4337143546881619_dp, &
0.6345603342446002_dp, &
0.9235964505628683_dp, &
1.337816370435542_dp, &
1.9295094395846712_dp, &
2.7733723263513426_dp, &
3.9789016642109_dp, &
5.715957168737114_dp, &
8.282101115138765_dp, &
12.368480022044524_dp, &
2.0178107659639145e-08_dp, &
5.6367263755163653e-08_dp, &
1.2429756911327585e-07_dp, &
2.6475339343624617e-07_dp, &
5.476340262894731e-07_dp, &
1.0955916127624847e-06_dp, &
2.1238343189909884e-06_dp, &
4.002849084627152e-06_dp, &
7.358331190312341e-06_dp, &
1.3229121263101538e-05_dp, &
2.331386732889364e-05_dp, &
4.035186295308255e-05_dp, &
6.870487234055741e-05_dp, &
0.00011523753234545123_dp, &
0.00019063770116232648_dp, &
0.00031138004039752104_dp, &
0.0005026233933172349_dp, &
0.0008024546937357704_dp, &
0.0012680625134385962_dp, &
0.0019846594459929524_dp, &
0.0030782984853653895_dp, &
0.0047341751428670595_dp, &
0.007222618433054622_dp, &
0.01093580663098367_dp, &
0.01643937738198581_dp, &
0.02454463892885994_dp, &
0.036409168383726855_dp, &
0.053676432340477115_dp, &
0.07866903483313294_dp, &
0.11465596124675029_dp, &
0.16622322356602454_dp, &
0.23979326244249682_dp, &
0.3443710265724193_dp, &
0.4926706308373399_dp, &
0.702973296327583_dp, &
1.0026242874881994_dp, &
1.435827423004336_dp, &
2.084888606410416_dp, &
3.146202983961085_dp, &
5.3833964099143214_dp /)
ELSE
aw(:) = (/ &
9.244859262370341e-11_dp, &
9.688645482714197e-10_dp, &
5.0141911163399264e-09_dp, &
1.9114844822899284e-08_dp, &
6.072338332197972e-08_dp, &
1.7034082386726808e-07_dp, &
4.3595933052332385e-07_dp, &
1.039095738862811e-06_dp, &
2.33872594629414e-06_dp, &
5.020178455196909e-06_dp, &
1.035308758243267e-05_dp, &
2.0629103356489906e-05_dp, &
3.989105138561648e-05_dp, &
7.512779107540237e-05_dp, &
0.00013820343998499602_dp, &
0.0002489306968715815_dp, &
0.00043990662666871475_dp, &
0.0007640355919502909_dp, &
0.0013061109710126168_dp, &
0.0022004740232120595_dp, &
0.00365769567358094_dp, &
0.006004552297192447_dp, &
0.009743449616281806_dp, &
0.015640101358779388_dp, &
0.024851991158373003_dp, &
0.03911535108729149_dp, &
0.06101562961258857_dp, &
0.09437645297561886_dp, &
0.14481597085262468_dp, &
0.22053880202868112_dp, &
0.33345906861219865_dp, &
0.5007898737276563_dp, &
0.7472971932166872_dp, &
1.1085269729134606_dp, &
1.6355437589901618_dp, &
2.4022768908057746_dp, &
3.518131795990682_dp, &
5.1535915819478095_dp, &
7.607002861713278_dp, &
11.56756637445306_dp, &
2.895556289542478e-10_dp, &
1.8191615349764746e-09_dp, &
7.325069811525452e-09_dp, &
2.3589548740506355e-08_dp, &
6.599224781728753e-08_dp, &
1.6717805065771188e-07_dp, &
3.929748485323067e-07_dp, &
8.707470599599378e-07_dp, &
1.838598331323936e-06_dp, &
3.728874808825026e-06_dp, &
7.307135065161298e-06_dp, &
1.3899407967673643e-05_dp, &
2.575809372083427e-05_dp, &
4.664295094550909e-05_dp, &
8.273201398017582e-05_dp, &
0.00014403314948488176_dp, &
0.00024654777250497477_dp, &
0.0004155549321896223_dp, &
0.00069055129044841_dp, &
0.001132618058962_dp, &
0.0018353181252988874_dp, &
0.002940697347526154_dp, &
0.004662616644432808_dp, &
0.007320542840879426_dp, &
0.01138818282123991_dp, &
0.017563076730560113_dp, &
0.026865638325958154_dp, &
0.040779390272580325_dp, &
0.061448678838256214_dp, &
0.09195663405552101_dp, &
0.1367158806495642_dp, &
0.20202078358459333_dp, &
0.296841502199106_dp, &
0.43401169638772713_dp, &
0.6321477781269876_dp, &
0.9191736530880799_dp, &
1.3400399323943677_dp, &
1.9776481665483123_dp, &
3.0278338226087165_dp, &
5.247421174586073_dp /)
END IF
END SUBROUTINE

END MODULE minimax_exp_k15
