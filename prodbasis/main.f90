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
   integer,parameter  :: nr     = 30
   integer,parameter  :: ntheta = 25
   integer,parameter  :: nphi   = 25
   real(kr),parameter :: rmax   = 1.5_kr
   real(kr),parameter :: dr     = rmax/nr
   real(kr),parameter :: dtheta = pi/ntheta
   real(kr),parameter :: dphi   = 2*pi/nphi
   real(kr),parameter :: hartree = 27.21138602 !eV

   real(kr),parameter :: beta   = 40.0_kr
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
         real(kr),intent(in)    :: basisarray(nstates,nr,ntheta,nphi)
         integer(ki),intent(in) :: nstates
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

         real,allocatable :: matrix(:,:)  ! matrix for SSYEV
         integer            :: INFO_lpck, m1,m2
         real,allocatable :: W_lpck(:), WORK_lpck(:)
         external SSYEV

         allocate( matrix(nstates,nstates) )
         allocate( W_lpck(nstates) )
         allocate( WORK_lpck(3*nstates-1) )

         ! copy data
         do m1=1,nstates
            do m2=1,nstates
               matrix(m1,m2) = A(m1,m2)
            enddo
         enddo

            ! SSYEV( JOBZ, UPLO, N,        A,     LDA,     W,       WORK,      LWORK,      INFO )
         call SSYEV( 'V', 'U', nstates, matrix, nstates, W_lpck, WORK_lpck, 3*nstates-1, INFO_lpck )

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

end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
   use constants
   implicit none
   integer                 :: iounit=10, i,j,k,l,m1,m2,m3,m4
   real(kr)                :: tmp,tmp1,omatrix_sp(norb,norb),omatrix_prod(nprodstates,nprodstates)
   real(kr)                :: evals(nprodstates), evecs(nprodstates,nprodstates)
   real(kr)                :: Umat(nprodstates,nprodstates), Dmat(nprodstates,nprodstates), r,t,p
   real(kr),allocatable    :: spbasis(:,:,:,: ), prodbasis(:,:,:,: )

   allocate( spbasis(norb,nr,ntheta,nphi) )
   allocate( prodbasis(nprodstates,nr,ntheta,nphi) )
   write(*,'(A,F8.5,A)') 'Allocated product basis:',(sizeof(spbasis)+sizeof(prodbasis))/1024.0**3,' Gbyte'

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
         write(*,'(F10.5,3X)',advance='no') omatrix_sp(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Product-basis overlap matrix...'
   call overlap_mat( prodbasis, nprodstates,omatrix_prod )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F9.5,3X)',advance='no') omatrix_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Eigenvalues of the product-states overlap matrix...'
   call diagonalize(omatrix_prod, nprodstates, evals, evecs)
   Dmat = (0.0_kr)
   do m1=1,nprodstates
      write(*,'(ES12.5,3X)',advance='no') evals(m1)
      Dmat(m1,m1) = 1.0_kr/sqrt(evals(m1))
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Eigenvectors of the product-states overlap matrix...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F8.5,3X)',advance='no') evecs(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   ! get the trafo
   write(*,'(A)') 'Basis-Trafo matrix for the product-basis...'
   Umat = matmul( evecs, Dmat )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F9.5,3X)',advance='no') Umat(m1,m2)
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
         write(*,'(F8.5,3X)',advance='no') omatrix_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,'(A)') 'Calculate single-particle U_abab matrix...'
   do m1=1,norb
      do m2=1,norb
         write(*,'(F8.5,3X)',advance='no') Utensor_sp(m1,m2,m1,m2, spbasis )
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Calculate single-particle U_abba matrix...'
   do m1=1,norb
      do m2=1,norb
         write(*,'(F8.5,3X)',advance='no') Utensor_sp(m1,m2,m2,m1, spbasis )
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Calculate single-particle U_aabb matrix...'
   do m1=1,norb
      do m2=1,norb
         write(*,'(F8.5,3X)',advance='no') Utensor_sp(m1,m1,m2,m2, spbasis )
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Full Coulomb interaction matrix in product basis...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         !write(*,'(F8.5,3X)',advance='no') Utensor_prod(m1,m2, prodbasis )
    !!!!!!!!!!!!!! 
         omatrix_prod(m1,m2) = Utensor_prod(m1,m2, prodbasis )
         write(*,'(F8.5,3X)',advance='no') omatrix_prod(m1,m2)
   !!!!!!!!!!!!!!!!
      enddo
      write(*,'(A)') ''
   enddo

   write(*,'(A)') 'Diagonalize Coulomb interaction in product basis...'
   call diagonalize(omatrix_prod, nprodstates, evals, evecs)
   do m1=1,nprodstates
      write(*,'(F8.5,3X)') evals(m1)
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   deallocate( prodbasis )
   deallocate( spbasis )

end program main
