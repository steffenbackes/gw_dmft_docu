module constants             
   implicit none             
   integer,    parameter :: kr   = selected_real_kind(8)
   integer,    parameter :: ki   = selected_int_kind(8)
   real(kr),   parameter :: pi   = 4.0_kr*atan(1.0_kr)
   real(kr),   parameter :: fpi  = 4.0_kr*pi
   real(kr),   parameter :: zero = 0.0_kr
   real(kr),   parameter :: one  = 1.0_kr
   complex(kr),parameter :: ci   = (0.0_kr,1.0_kr)
 
   integer,parameter  :: norb   = 3
   integer,parameter  :: nprodstates = 6
   integer,parameter  :: nr     = 12
   integer,parameter  :: ntheta = 12
   integer,parameter  :: nphi   = 12
   real(kr),parameter :: rmax   = 1.5_kr
   real(kr),parameter :: dr     = rmax/nr
   real(kr),parameter :: dtheta = pi/ntheta
   real(kr),parameter :: dphi   = 2*pi/nphi
   real(kr),parameter :: hartree = 27.21138602 !eV

   real(kr),parameter :: beta   = 10.0_kr
   integer,parameter  :: nw     = 5
   real(kr),parameter :: wmax   = 5.0_kr
   real(kr),parameter :: dw     = wmax/nw
   real(kr),parameter :: Z = 26.0_kr
   real(kr),parameter :: LagPref = sqrt( 1.0/720.0 ) * ( 2*Z/3.0_kr )**3.5_kr 

      interface
         function func(r,t,p)
            implicit none             
            integer,    parameter :: kr   = selected_real_kind(8)
            real(kr)            :: func
            real(kr),intent(in) :: r,t,p
        end function func
      end interface
      type Contains_f_ptr
         procedure (func), pointer, nopass :: my_f_ptr
      end type Contains_f_ptr


     type (Contains_f_ptr), dimension (norb) :: spstates
     type (Contains_f_ptr), dimension (nprodstates) :: prodstates

   contains

      function wn(n)
         real(kr)               :: wn
         integer(ki),intent(in) :: n
        wn = (2*n+1)*pi/beta
      end function wn

   
     ! atomic single-particle t2g basis functions

      function psi1(r,t,p)
         real(kr)            :: psi1
         real(kr),intent(in) :: r,t,p
         psi1 = sin(2*p)*sin(t)**2  * 0.25*sqrt(15.0/pi)  * LagPref*exp(-r*Z/3.0)*r**2  ! 3dxy
      end function psi1

      function psi2(r,t,p)
         real(kr)             :: psi2
         real(kr),intent(in) :: r,t,p
         psi2 = cos(p)*sin(t)*cos(t) * sqrt(15/(4.0*pi)) * LagPref*exp(-r*Z/3.0)*r**2 ! 3dxz
      end function psi2

      function psi3(r,t,p)
         real(kr)            :: psi3
         real(kr),intent(in) :: r,t,p
         psi3 = sin(p)*sin(t)*cos(t) * sqrt(15/(4.0*pi)) * LagPref*exp(-r*Z/3.0)*r**2  ! 3dyz
      end function psi3

  ! simple product of all single-particle wave functions without duplicates 3*2 = 6 functions
      function Bf1(r,t,p)
         real(kr)            :: Bf1
         real(kr),intent(in) :: r,t,p
         Bf1 = psi1(r,t,p)*psi1(r,t,p)
      end function Bf1

      function Bf2(r,t,p)
         real(kr)            :: Bf2
         real(kr),intent(in) :: r,t,p
         Bf2 = psi1(r,t,p)*psi2(r,t,p)
      end function Bf2

      function Bf3(r,t,p)
         real(kr)            :: Bf3
         real(kr),intent(in) :: r,t,p
         Bf3 = psi1(r,t,p)*psi3(r,t,p)
      end function Bf3

      function Bf4(r,t,p)
         real(kr)            :: Bf4
         real(kr),intent(in) :: r,t,p
         Bf4 = psi2(r,t,p)*psi2(r,t,p)
      end function Bf4

      function Bf5(r,t,p)
         real(kr)            :: Bf5
         real(kr),intent(in) :: r,t,p
         Bf5 = psi2(r,t,p)*psi3(r,t,p)
      end function Bf5

      function Bf6(r,t,p)
         real(kr)            :: Bf6
         real(kr),intent(in) :: r,t,p
         Bf6 = psi3(r,t,p)*psi3(r,t,p)
      end function Bf6

      subroutine overlap_mat( basisarray, nstates, omat )
         integer(ki),intent(in) :: nstates
         real(kr),intent(in)    :: basisarray(nstates,nr,ntheta,nphi)
         real(kr),intent(inout) :: omat(nstates,nstates)
         integer  :: i,j,k,m1,m2
         real(kr) :: r,t, sint

         omat = (0.0_kr)
         do m1=1,nstates
            do m2=1,nstates
               do i=1,nr
                  r = (i-1)*dr
                  do j=1,ntheta
                     t = (j-1)*dtheta
                     sint = sin(t)
                     do k=1,nphi

                        omat(m1,m2) = omat(m1,m2) + basisarray(m1,i,j,k)*basisarray(m2,i,j,k)*sint*r**2

                     enddo
                  enddo
               enddo
            enddo
         enddo
         omat = omat*dr*dphi*dtheta
      end subroutine overlap_mat

      subroutine diagonalize(A, nstates, evals, evecs)
         real(kr),intent(in)    :: A(:,:)
         integer(ki),intent(in) :: nstates
         real(kr),intent(inout) :: evals(:)
         real(kr),intent(inout) :: evecs(:,:)

         double precision,allocatable :: matrix(:,:) 
         integer                      :: INFO_lpck, m1,m2
         double precision,allocatable :: W_lpck(:), WORK_lpck(:)
         external DSYEV

         allocate( matrix(nstates,nstates) )
         allocate( W_lpck(nstates) )
         allocate( WORK_lpck(3*nstates-1) )

         ! copy data
         do m1=1,nstates
            do m2=1,nstates
               matrix(m1,m2) = A(m1,m2)
            enddo
         enddo

         call DSYEV( 'V', 'U', nstates, matrix, nstates, W_lpck, WORK_lpck, 3*nstates-1, INFO_lpck )

         if (INFO_lpck /= 0) then
            write(*,'(A,I2)') 'ERROR: Diagonalization returned: ',INFO_lpck
            stop
         endif

         ! copy data
         do m1=1,nstates
            evals(m1) = W_lpck(m1)
            do m2=1,nstates
               evecs(m1,m2) = matrix(m1,m2)
            enddo
         enddo

         deallocate( matrix )
         deallocate( W_lpck )
         deallocate( WORK_lpck )

      end subroutine diagonalize

      subroutine invert(A, nstates)
         real(kr),intent(inout)    :: A(:,:)
         integer(ki),intent(in) :: nstates

         real(8),allocatable :: matrix(:,:)  ! matrix for DGETRI
         integer             :: INFO_lpck, m1,m2, ipiv(nstates)
         real(8)             :: WORK_lpck(nstates)
         external DGETRI
         external DGETRF

         allocate( matrix(nstates,nstates) )

         ! copy data
         do m1=1,nstates
            do m2=1,nstates
               matrix(m1,m2) = A(m1,m2)
            enddo
         enddo

         call DGETRF( nstates, nstates, matrix, nstates, ipiv, INFO_lpck )
         if (INFO_lpck /= 0) then
            write(*,'(A,I2)') 'ERROR: DGETRF returned: ',INFO_lpck
            stop
         endif

         call DGETRI( nstates, matrix, nstates, ipiv, WORK_lpck, nstates, INFO_lpck )
         if (INFO_lpck /= 0) then
            write(*,'(A,I2)') 'ERROR: DGETRI returned: ',INFO_lpck
            stop
         endif

         ! copy data
         do m1=1,nstates
            do m2=1,nstates
               A(m1,m2) = matrix(m1,m2)
            enddo
         enddo

         deallocate( matrix )

      end subroutine invert

      function Utensor_sp(m1,m2,m3,m4, basisarray )
         real(kr)               :: Utensor_sp
         integer(ki),intent(in) :: m1,m2,m3,m4
         real(kr),intent(in)    :: basisarray(norb,nr,ntheta,nphi)
         integer  :: i1,j1,k1,i2,j2,k2
         real(kr) :: r1,r2,t1,t2,p1,p2, sint1,sint2,cost1,cost2,cosp1p2, dist

         Utensor_sp = 0.0_kr

         do j1=1,ntheta
            t1 = (j1-1)*dtheta - dtheta/10
            sint1 = sin(t1)
            cost1 = cos(t1)

            do j2=1,ntheta
               t2 = (j2-1)*dtheta + dtheta/10
               sint2 = sin(t2)
               cost2 = cos(t2)

               do k1=1,nphi
                  p1 = (k1-1)*dphi - dphi/10
                  do k2=1,nphi
                     p2 = (k2-1)*dphi + dphi/10
                     cosp1p2 = cos(p1-p2)

                     do i1=1,nr
                        r1 = (i1-1)*dr - dr/10
                        do i2=1,nr
                           r2 = (i2-1)*dr + dr/10

                           dist = sqrt(r1*r1+r2*r2-2*r1*r2*( sint1*sint2*cosp1p2 + cost1*cost2 ) )
                     
                           Utensor_sp = Utensor_sp                                     &
                                &  + basisarray(m1,i1,j1,k1)*basisarray(m2,i2,j2,k2)   &
                                &   *basisarray(m3,i1,j1,k1)*basisarray(m4,i2,j2,k2)   &
                                &   *sint1*sint2*r1*r1*r2*r2                           &
                                &   /dist                                              

                        enddo
                     enddo
                  enddo      
               enddo
            enddo
         enddo
         Utensor_sp = Utensor_sp*(dr*dphi*dtheta)**2 * hartree
      end function Utensor_sp

      function Utensor_prod(m1,m2, basisarray )
         real(kr)               :: Utensor_prod
         integer(ki),intent(in) :: m1,m2
         real(kr),intent(in)    :: basisarray(nprodstates,nr,ntheta,nphi)
         integer  :: i1,j1,k1,i2,j2,k2
         real(kr) :: r1,r2,t1,t2,p1,p2, sint1,sint2,cost1,cost2,cosp1p2, dist

         Utensor_prod = 0.0_kr

         do j1=1,ntheta
            t1 = (j1-1)*dtheta - dtheta/10
            sint1 = sin(t1)
            cost1 = cos(t1)

            do j2=1,ntheta
               t2 = (j2-1)*dtheta + dtheta/10
               sint2 = sin(t2)
               cost2 = cos(t2)

               do k1=1,nphi
                  p1 = (k1-1)*dphi - dphi/10
                  do k2=1,nphi
                     p2 = (k2-1)*dphi + dphi/10
                     cosp1p2 = cos(p1-p2)

                     do i1=1,nr
                        r1 = (i1-1)*dr - dr/10
                        do i2=1,nr
                           r2 = (i2-1)*dr + dr/10

                           dist = sqrt(r1*r1+r2*r2-2*r1*r2*( sint1*sint2*cosp1p2 + cost1*cost2 ) )
                     
                           Utensor_prod = Utensor_prod           &
                                &  + basisarray(m1,i1,j1,k1)     &
                                &   *basisarray(m2,i2,j2,k2)     &
                                &   *sint1*sint2*r1*r1*r2*r2     &
                                &   /dist                                              

                        enddo
                     enddo
                  enddo      
               enddo
            enddo
         enddo
         Utensor_prod = Utensor_prod*(dr*dphi*dtheta)**2 * hartree
      end function Utensor_prod

      subroutine get_polarization(pol, basisarray )
         complex(kr),intent(inout) :: pol(nr,ntheta,nphi, nr,ntheta,nphi, nw) 
         real(kr),intent(in)       :: basisarray(norb,nr,ntheta,nphi)
         integer     :: i1,j1,k1,i2,j2,k2, m1,m2,n
         real(kr)    :: eps(norb)
         complex(kr) :: prefac

         eps(1) = -0.01_kr
         eps(2) =  0.0_kr
         eps(3) =  0.01_kr

         pol = (0.0_kr)

         do m1=1,norb
            do m2=1,norb
               do n=1,nw
                  prefac = ( 1.0_kr/(1.0_kr+exp(beta*eps(m1))) - 1.0_kr/(1.0_kr+exp(beta*eps(m2)))  ) &
                      &    /( (n-1)*dw + eps(m1) - eps(m2) + CMPLX(0.0_kr,0.01_kr,kind=kr) )

                  do i1=1,nr
                     do j1=1,ntheta
                        do k1=1,nphi
                           do i2=1,nr
                              do j2=1,ntheta
                                 do k2=1,nphi

                                    pol(i1,j1,k1, i2,j2,k2, n) = pol(i1,j1,k1, i2,j2,k2, n) &
                                          & + prefac*basisarray(m1,i1,j1,k1)                &
                                          &         *basisarray(m2,i1,j1,k1)                &
                                          &         *basisarray(m1,i2,j2,k2)                &
                                          &         *basisarray(m2,i2,j2,k2)              
                                 enddo
                              enddo
                           enddo

                        enddo
                     enddo
                  enddo      

               enddo
            enddo
         enddo
      end subroutine get_polarization

      recursive real(kr) function det(n,A) result(res)
         integer(ki),intent(in) :: n
         real(kr),intent(in)    :: A(n,n)

         integer(ki)          :: i, i1, m1,m2
         real(kr)             :: d
         real(kr),allocatable :: B(:,:)

         d=0.0_kr
         if (n==1) then
            d = A(1,1)
         else
            allocate( B(n-1,n-1) )
            do i=1,n

               ! create submatrix !!!
               i1 = 1
               do m1=1,n
                  if (m1/=i) then
                     do m2=1,n-1
                        B(i1,m2) = A(m1,m2+1)
                     enddo
                     i1 = i1+1
                  endif
               enddo             
               !!!!!!!!!!!!!!!!!!!!!!     
               d = d + (-1)**(i-1) * A(i,1) * det(n-1,B)
            end do
            deallocate( B )
         endif
         res = d
         return
      end function det

end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
   use constants
   implicit none
   integer                 :: i,j,k,l,m1,m2,m3,m4
   real(kr)                :: tmp,omatrix_sp(norb,norb),omatrix_prod(nprodstates,nprodstates)
   real(kr)                :: evals(nprodstates), evecs(nprodstates,nprodstates)
   real(kr)                :: Umat(nprodstates,nprodstates), Dmat(nprodstates,nprodstates), r,t,p,sint
   real(kr),allocatable    :: spbasis(:,:,:,: ), prodbasis(:,:,:,: ), proj_spprod(:,:,:)
   real(kr)                :: Vcoul_prod(nprodstates,nprodstates), Vcoul_sp(norb**2,norb**2)
   complex(kr),allocatable :: pol_rspace(:,:,:, :,:,:, : )

   allocate( spbasis(norb,nr,ntheta,nphi) )
   allocate( prodbasis(nprodstates,nr,ntheta,nphi) )
   allocate( proj_spprod(norb,norb,nprodstates) )
   allocate( pol_rspace(nr,ntheta,nphi, nr,ntheta,nphi, nw) )
   write(*,'(A,F7.3,A)') 'Allocated polarization in real space:', sizeof(pol_rspace)/1024.0_kr**3,' Gbyte'
   write(*,'(A,F7.3,A)') 'Allocated product basis:',(sizeof(spbasis)+sizeof(prodbasis))/1024.0_kr**3,' Gbyte'

   ! First assign function pointers
   spstates(1)%my_f_ptr => psi1
   spstates(2)%my_f_ptr => psi2
   spstates(3)%my_f_ptr => psi3
   prodstates(1)%my_f_ptr => Bf1
   prodstates(2)%my_f_ptr => Bf2
   prodstates(3)%my_f_ptr => Bf3
   prodstates(4)%my_f_ptr => Bf4
   prodstates(5)%my_f_ptr => Bf5
   prodstates(6)%my_f_ptr => Bf6

   write(*,'(A)') 'Fill single-particle basis array with values...'
   spbasis = (0.0_kr)
   do m1=1,norb
      do i=1,nr
         r = (i-1)*dr
         do j=1,ntheta
            t = (j-1)*dtheta
            do k=1,nphi
               p = (k-1)*dphi
               spbasis(m1,i,j,k) = spstates(m1)%my_f_ptr(r,t,p)
            enddo
         enddo
      enddo
   enddo

   write(*,'(A)') 'Fill the initial product basis array with values...'
   prodbasis = (0.0_kr)
   do m1=1,nprodstates
      do i=1,nr
         r = (i-1)*dr
         do j=1,ntheta
            t = (j-1)*dtheta
            do k=1,nphi
               p = (k-1)*dphi
               prodbasis(m1,i,j,k) = prodstates(m1)%my_f_ptr(r,t,p)
            enddo
         enddo
      enddo
   enddo

   write(*,'(A)') 'Single-particle basis  overlap matrix...'
   call overlap_mat( spbasis, norb, omatrix_sp )
   do m1=1,norb
      do m2=1,norb
         write(*,'(F7.4,2X)',advance='no') omatrix_sp(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Product-basis overlap matrix...'
   call overlap_mat( prodbasis, nprodstates,omatrix_prod )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F7.3,2X)',advance='no') omatrix_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Eigenvalues of the product-states overlap matrix...'
   call diagonalize(omatrix_prod, nprodstates, evals, evecs)
   Dmat = (0.0_kr)
   do m1=1,nprodstates
      write(*,'(ES10.3,2X)',advance='no') evals(m1)
      Dmat(m1,m1) = 1.0_kr/sqrt(evals(m1))
   enddo
   write(*,'(A)') ''

   ! choose the sign for the Evectors so that maximum entry is positive
   do m1=1,nprodstates
      tmp = 0.0_kr
      do m2=1,nprodstates
         if (abs(tmp)<abs(evecs(m2,m1))) tmp=evecs(m2,m1)
      enddo
      do m2=1,nprodstates
         evecs(m2,m1) = evecs(m2,m1) * ( tmp/abs(tmp) )
      enddo
   enddo

   write(*,'(A)') 'Eigenvectors of the product-states overlap matrix...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F7.3,2X)',advance='no') evecs(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   ! get the trafo
   write(*,'(A)') 'Basis-Trafo matrix for the product-basis...'
   Umat = matmul( evecs, Dmat )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F7.3,2X)',advance='no') Umat(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Fill the product basis array with values...'
   prodbasis = (0.0_kr)
   do m1=1,nprodstates
      do i=1,nr
         r = (i-1)*dr
         do j=1,ntheta
            t = (j-1)*dtheta
            do k=1,nphi
               p = (k-1)*dphi
               do m2=1,nprodstates
                  prodbasis(m1,i,j,k) = prodbasis(m1,i,j,k) &
                         & + prodstates(m2)%my_f_ptr(r,t,p)*Umat(m2,m1)
               enddo
            enddo
         enddo
      enddo
   enddo

   write(*,'(A)') 'Product-basis overlap matrix after orthonormalization...'
   call overlap_mat( prodbasis, nprodstates, omatrix_prod )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F7.3,2X)',advance='no') omatrix_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Generate projectors between single-particle and product basis...'
   proj_spprod = (0.0_kr)
   do i=1,nr
      r = (i-1)*dr
      do j=1,ntheta
         t = (j-1)*dtheta
         sint = sin(t)
         do k=1,nphi

            do m1=1,norb
               do m2=1,norb
                  do l=1,nprodstates
                     proj_spprod(m1,m2,l) = proj_spprod(m1,m2,l)    &
                                       &  + spbasis(m1,i,j,k)       & 
                                       &  * spbasis(m2,i,j,k)       & 
                                       &  * prodbasis(l,i,j,k) *sint*r**2
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   proj_spprod = proj_spprod*dr*dphi*dtheta

   do m1=1,norb
   do m2=1,norb
      do l=1,nprodstates
         write(*,'(F7.3,2X)',advance='no') proj_spprod(m1,m2,l)
      enddo
      write(*,'(A)') ''
   enddo
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'P^dagger P product (nprodstates x nprodstates)...'
   do m1=1,nprodstates
      do m2=1,nprodstates

         tmp = 0.0_kr
         do i=1,norb
            do j=1,norb
               tmp = tmp + proj_spprod(i,j,m1)*proj_spprod(i,j,m2)
            enddo
         enddo
         write(*,'(F7.3,2X)',advance='no') tmp

      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'P P^dagger product (norb**2 x norb**2)...'
   do i=1,norb
      do j=1,norb
         do k=1,norb
            do l=1,norb

               tmp = 0.0_kr
               do m1=1,nprodstates
                  tmp = tmp + proj_spprod(i,j,m1)*proj_spprod(k,l,m1)
               enddo
               write(*,'(F7.3,2X)',advance='no') tmp
            enddo
         enddo
         write(*,'(A)') ''
      enddo
   enddo
   write(*,'(A)') ''
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   write(*,'(A)') 'Generate polarization in real space...'
!   call get_polarization(pol_rspace, spbasis )
!write(*,*) pol_rspace(1,1,1, 1,1,1, :)
!stop 0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! single-particle convergence check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,'(I3,2X,F7.3,2X,I3,2X,I3,2X,F7.3,2X,F7.3,2X,F7.3,2X)') nr, rmax, ntheta,nphi, & 
!        & Utensor_sp(1,1,1,1, spbasis), Utensor_sp(1,2,1,2, spbasis), Utensor_sp(1,2,2,1, spbasis)
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! product basis convergence check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(*,'(I3,2X,F7.3,2X,I3,2X,I3,2X,F7.3,2X,F7.3,2X,F7.3,2X)') nr, rmax, ntheta,nphi, & 
!        & Utensor_prod(1,4, prodbasis ), Utensor_prod(4,1, prodbasis ), Utensor_prod(6,6, prodbasis )
!stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,'(A)') 'Calculate single-particle U_abab matrix...'
   do m1=1,norb
      do m2=1,norb
         write(*,'(F7.3,2X)',advance='no') Utensor_sp(m1,m2,m1,m2, spbasis )
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Calculate single-particle U_abba matrix...'
   do m1=1,norb
      do m2=1,norb
         write(*,'(F7.3,2X)',advance='no') Utensor_sp(m1,m2,m2,m1, spbasis )
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Calculate single-particle U_aabb matrix...'
   do m1=1,norb
      do m2=1,norb
         write(*,'(F7.3,2X)',advance='no') Utensor_sp(m1,m1,m2,m2, spbasis )
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Full Coulomb interaction matrix in product basis...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         Vcoul_prod(m1,m2) = Utensor_prod(m1,m2, prodbasis )
         write(*,'(F7.3,2X)',advance='no') Vcoul_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') 'Symmetrize Coulomb interaction matrix in product basis...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         tmp = 0.5*( Vcoul_prod(m1,m2)+Vcoul_prod(m2,m1) ) ! symmetrize to get rid of numerical errors
         Vcoul_prod(m1,m2) = tmp
         Vcoul_prod(m2,m1) = tmp
         write(*,'(F7.3,2X)',advance='no') Vcoul_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A,ES10.3)') 'Determinant:',det(nprodstates,Vcoul_prod)
   write(*,'(A)') ''

   write(*,'(A)') 'Full Coulomb matrix in index-combination basis...'
   do m1=1,norb
   do m2=1,norb
      do m3=1,norb
      do m4=1,norb
         Vcoul_sp((m1-1)*norb+m2, (m3-1)*norb+m4) = Utensor_sp(m1,m2,m3,m4, spbasis )
         write(*,'(F7.3,2X)',advance='no') Vcoul_sp((m1-1)*norb+m2, (m3-1)*norb+m4)
      enddo
      enddo
      write(*,'(A)') ''
   enddo
   enddo
   write(*,'(A,ES10.3)') 'Determinant:',det(norb**2,Vcoul_sp)
   write(*,'(A)') ''

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,'(A)') 'Full Coulomb matrix in index-combination generated from Product representation...'
   do m1=1,norb
   do m2=1,norb
      do m3=1,norb
      do m4=1,norb
         tmp = 0.0_kr
         do i=1,nprodstates
            do j=1,nprodstates
               tmp = tmp + proj_spprod(m1,m3,i) * Vcoul_prod(i,j) * proj_spprod(m4,m2,j)
            enddo
         enddo
         write(*,'(F7.3,2X)',advance='no') tmp
         Vcoul_sp((m1-1)*norb+m2, (m3-1)*norb+m4) = tmp
      enddo
      enddo
      write(*,'(A)') ''
   enddo
   enddo
   write(*,'(A)') 'Copy this to the index-combination tensor since it is better...' 
   write(*,'(A)') ''
 
   write(*,'(A)') 'Full Coulomb interaction matrix in product basis generated from 4-index tensor...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         tmp = 0.0_kr
     
         do i=1,norb
            do j=1,norb
               do k=1,norb
                  do l=1,norb
                     tmp = tmp + proj_spprod(i,k,m1) * Vcoul_sp((i-1)*norb+j, (k-1)*norb+l) * proj_spprod(l,j,m2)
                  enddo
               enddo
            enddo
         enddo

         write(*,'(F7.3,2X)',advance='no') tmp
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,'(A)') 'Eigenvalues of the Coulomb matrix in the product-basis...'
   call diagonalize(Vcoul_prod, nprodstates, evals, evecs)
   do m1=1,nprodstates
      write(*,'(ES10.3,2X)',advance='no') evals(m1)
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Eigenvalues of the Coulomb matrix in the index-combination basis...'
   call diagonalize(Vcoul_sp, norb**2, evals, evecs)
   do m1=1,norb**2
      write(*,'(ES10.3,2X)',advance='no') evals(m1)
   enddo
   write(*,'(A)') ''

   write(*,'(A)') ''

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   write(*,'(A)') 'Inverse Coulomb interaction matrix in product basis...'
   call invert(Vcoul_prod, nprodstates) 
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F7.3,2X)',advance='no') Vcoul_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Inverse Coulomb interaction matrix in index-combination basis...'
   call invert(Vcoul_sp, norb**2) 
   do m1=1,norb**2
      do m2=1,norb**2
         write(*,'(F7.3,2X)',advance='no') Vcoul_sp(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Inverse Coulomb matrix in index-combination generated from Product representation...'
   do m1=1,norb
   do m2=1,norb
      do m3=1,norb
      do m4=1,norb
         tmp = 0.0_kr
         do i=1,nprodstates
            do j=1,nprodstates
               tmp = tmp + proj_spprod(m1,m3,i) * Vcoul_prod(i,j) * proj_spprod(m4,m2,j)
            enddo
         enddo
         write(*,'(F7.3,2X)',advance='no') tmp
      enddo
      enddo
      write(*,'(A)') ''
   enddo
   enddo
   write(*,'(A)') ''

   deallocate( prodbasis )
   deallocate( spbasis )
   deallocate( proj_spprod )
   deallocate( pol_rspace )

end program main
