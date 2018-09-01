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
   integer,parameter  :: nr     = 25
   integer,parameter  :: ntheta = 20
   integer,parameter  :: nphi   = 20
   real(kr),parameter :: rmax   = 40.0_kr
   real(kr),parameter :: dr     = rmax/nr
   real(kr),parameter :: dtheta = pi/ntheta
   real(kr),parameter :: dphi   = 2*pi/nphi

   real(kr),parameter :: beta   = 40.0_kr
   real(kr),parameter :: LagPref = sqrt((2.0/3)**3 /( 6*120.0 ))*10

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
         psi1 = sin(2*p)*sin(t)**2   / ( 0.25*pi*sqrt(6.0) ) * LagPref*exp(-r/3.0)*(r*2/3.0)**2  ! 3dxy
      end function psi1

      function psi2(r,t,p)
         real(kr)             :: psi2
         real(kr),intent(in) :: r,t,p
         psi2 = cos(p)*sin(t)*cos(t) / ( 0.25*pi*sqrt(2.0) ) * LagPref*exp(-r/3.0)*(r*2/3.0)**2  ! 3dyz
      end function psi2

      function psi3(r,t,p)
         real(kr)            :: psi3
         real(kr),intent(in) :: r,t,p
         psi3 = sin(p)*sin(t)*cos(t) / ( 0.25*pi*sqrt(2.0) ) * LagPref*exp(-r/3.0)*(r*2/3.0)**2  ! 3dxz
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

      function overlap(f1,f2)
         real(kr)          :: overlap
         real(kr),external :: f1,f2
         integer  :: i,j,k
         real(kr) :: r,t,p

         overlap = 0.0_kr
         do i=1,nr
            r = i*dr
            do j=1,ntheta
               t = j*dtheta
               do k=1,nphi
                  p = k*dphi

                  overlap = overlap + f1(r,t,p)*f2(r,t,p)

               enddo
            enddo
         enddo
         overlap = overlap*dr*dphi*dtheta
      end function overlap

      subroutine overlap_mat( basisarray, omat )
         real(kr),intent(in)    :: basisarray(nprodstates,nr,ntheta,nphi)
         real(kr),intent(inout) :: omat(nprodstates,nprodstates)
         integer  :: i,j,k,m1,m2

         omat = 0.0_kr
         do m1=1,nprodstates
            do m2=1,nprodstates
               do i=1,nr
                  do j=1,ntheta
                     do k=1,nphi

                        omat(m1,m2) = omat(m1,m2) + basisarray(m1,i,j,k)*basisarray(m2,i,j,k)

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

      function Utensor(m1,m2,m3,m4, basisarray )
         real(kr)               :: Utensor
         integer(ki),intent(in) :: m1,m2,m3,m4
         real(kr),intent(in)    :: basisarray(nprodstates,nr,ntheta,nphi)
         integer  :: i1,j1,k1,i2,j2,k2
         real(kr) :: r1,r2,t1,t2, dist

         Utensor = 0.0_kr
         do i1=1,nr
            r1 = i1*dr
            do j1=1,ntheta
               t1 = j1*dtheta
               do k1=1,nphi

                  do i2=1,nr
                     r2 = i2*dr + dr/2
                     do j2=1,ntheta
                        t2 = j2*dtheta + dtheta/2
                        do k2=1,nphi
                     
                           dist = sqrt(r1**2+r2**2-2*r1*r2*cos(t1-t2))
                           Utensor = Utensor                                           &
                                &  + basisarray(m1,i1,j1,k1)*basisarray(m2,i2,j2,k2)   &
                                &   *( 1.0_kr/dist )                                   &
                                &   *basisarray(m3,i1,j1,k1)*basisarray(m4,i2,j2,k2)  

                        enddo
                     enddo
                  enddo      
               enddo
            enddo
         enddo
         Utensor = Utensor*(dr*dphi*dtheta)**2
      end function Utensor

end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
   use constants
   implicit none
   integer                 :: iounit=10, i,j,k,l,m1,m2
   real(kr)                :: tmp,tmp1,omatrix_sp(norb,norb),omatrix_prod(nprodstates,nprodstates)
   real(kr)                :: evals(nprodstates), evecs(nprodstates,nprodstates)
   real(kr)                :: Umat(nprodstates,nprodstates), Dmat(nprodstates,nprodstates), r,t,p
   real(kr),allocatable    :: prodbasis(:,:,:,: )

   allocate( prodbasis(nprodstates,nr,ntheta,nphi) )
   write(*,'(A,F8.5,A)') 'Allocated product basis:',sizeof(prodbasis)/1024.0**3,' Gbyte'

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

   write(*,'(A)') 'Calculate single-particle basis overlap...'
   do m1=1,norb
      do m2=1,norb
         omatrix_sp(m1,m2) = overlap( spstates(m1)%my_f_ptr, spstates(m2)%my_f_ptr )
         write(*,'(F8.5,3X)',advance='no') omatrix_sp(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   write(*,'(A)') 'Fill the initial product basis array with values...'
   prodbasis = (0.0_kr)
   do m1=1,nprodstates
      do i=1,nr
         r = i*dr
         do j=1,ntheta
            t = j*dtheta
            do k=1,nphi
               p = k*dphi
               prodbasis(m1,i,j,k) = prodstates(m1)%my_f_ptr(r,t,p)
            enddo
         enddo
      enddo
   enddo

   write(*,'(A)') 'Product-basis overlap matrix...'
   call overlap_mat( prodbasis, omatrix_prod )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F8.5,3X)',advance='no') omatrix_prod(m1,m2)
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
         r = i*dr
         do j=1,ntheta
            t = j*dtheta
            do k=1,nphi
               p = k*dphi
               do m2=1,nprodstates
                  prodbasis(m1,i,j,k) = prodbasis(m1,i,j,k) &
                         & + prodstates(m2)%my_f_ptr(r,t,p)*Umat(m2,m1)
               enddo
            enddo
         enddo
      enddo
   enddo

   write(*,'(A)') 'Product-basis overlap matrix after orthonormalization...'
   call overlap_mat( prodbasis, omatrix_prod )
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F8.5,3X)',advance='no') omatrix_prod(m1,m2)
      enddo
      write(*,'(A)') ''
   enddo
   write(*,'(A)') ''

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   write(*,'(A)') 'Calculate density-density terms of the Coulomb interaction...'
   do m1=1,nprodstates
      do m2=1,nprodstates
         write(*,'(F8.5,3X)',advance='no') Utensor(m1,m2,m1,m2, prodbasis )
      enddo
      write(*,'(A)') ''
   enddo

   deallocate( prodbasis )

end program main
