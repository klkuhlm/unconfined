Program Driver

  ! only needed for intel compiler
#ifdef INTEL
  use ifport, only : dbesj0, dbesj1     
#endif

  ! type definitions 
  use types
  
  ! constants and coefficients
  use constants, only : DP, PI, EP
 
  ! function to be evaluated in Laplace/Hankel space
  use lap_hank_soln, only : unconfined_wellbore_slug

  ! inverse Laplace transform routine (currently only de Hoog, et al)
  use inverse_Laplace_Transform, only : dehoog_invlap, dehoog_pvalues 

  ! some non-built-in functions
  use utilities, only : logspace

  implicit none

  type(invLaplace) :: lap
  type(invHankel) :: hank
  type(GaussLobatto) :: gl
  type(TanhSinh) :: ts
  type(well) :: w
  type(formation) :: f
  type(solution) :: s

  real(DP) :: x, dx

  integer :: i, numt, terms, minlogsplit, maxlogsplit, np

  
  ! vectors of results and intermediate steps

  real(DP), allocatable  :: ts(:), td(:), dt(:) 
  complex(EP), allocatable :: fa(:,:)
  complex(EP), allocatable :: finint(:), infint(:), totlap(:), p(:)
  real(EP) :: totint
  real(DP) :: tee

  !! k parameter in tanh-sinh quadrature (order is 2^k - 1)
  integer ::  m, nn, jj, j, numtfile, numtcomp

  ! read in data from file, do minor error checking
  ! and allocate some solution vectors
  call read_input(w,f,s,lap,hank,gl,ts)
  
  call write_header(w,f,s,lap,hank,gl,ts,20)
     
  ! loop over all requested times
  do i = 1, numt
     if (.not. quiet) then
        write(*,'(I4,A,ES11.4)') i,' td ',td(i)
     end if
     tsval = td(i)

     tskk(1:nst) = [(tsk-m, m=nst-1,0,-1)]
     tsNN(1:nst) = 2**tskk - 1
     tshh(1:nst) = 4.0_EP/(2**tskk)

     ! compute abcissa and for highest-order case (others are subsets)
     call tanh_sinh_setup(tskk(nst),j0zero(splitv(i))/rDw,tsw(1:tsNN(nst)),tsa(1:tsNN(nst)))

     ! compute solution at densest set of abcissa
     !$OMP PARALLEL DO PRIVATE(nn) SHARED(fa,tsa,tsNN,nst)
     do nn=1,tsNN(nst)
        fa(nn,1:np) = unconfined_wellbore_slug(tsa(nn))
     end do
     !$OMP END PARALLEL DO

     tmp(nst,1:np) = (j0zero(splitv(i))/rDw)/2.0_DP*sum(spread(tsw(:),2,np)*fa(:,:),dim=1)

     do j=nst-1,1,-1
        !  only need to re-compute weights for each subsequent step
        call tanh_sinh_setup(tskk(j),j0zero(splitv(i))/rDw,tsw(1:tsNN(j)),tsa(1:tsNN(j)))

        ! compute index vector, to slice up solution 
        ! for nst'th turn count regular integers
        ! for (nst-1)'th turn count even integers
        ! for (nst-2)'th turn count every 4th integer, etc...
        ii(1:tsNN(j)) = [( m* 2**(nst-j), m=1,tsNN(j) )]

        ! estimate integral with subset of function evaluations and appropriate weights
        tmp(j,1:np) = (j0zero(splitv(i))/rDw)/2.0_DP* &
             & sum(spread(tsw(1:tsNN(j)),2,np)*fa(ii(1:tsNN(j)),1:np),dim=1)
     end do

     if (nst > 1) then
        !$OMP PARALLEL DO PRIVATE(jj,PolErr) SHARED(tshh,tmp,finint)
        do jj=1,np
           call polint(tshh(1:nst),tmp(1:nst,jj),0.0_EP,finint(jj),PolErr)
        end do
        !$OMP END PARALLEL DO
     else
        finint(:) = tmp(1,1:np)
     end if


     !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! "infinite" portion of Hankel  integral for each time level
     ! integrate between zeros of Bessel function, extrapolate 
     ! area from series of areas of alternating sign

     if(.not. allocated(GLx)) then
        allocate(Glx(GLorder-2),GLw(GLorder-2))
        call gauss_lobatto_setup(GLorder,GLx,GLw)  ! get weights and abcissa
     end if

     tsval = td(i)

     infint(1:np) = inf_integral(unconfined_wellbore_slug,splitv(i),naccel,&
          &glorder,j0zero,GLx,GLw)

     ! perform numerical inverse Laplace transform and sum infinite and finite
     ! portions of Hankel integral
     tee = td(i)*2.0
     p = dehoog_pvalues(tee,lap)
     totlap(1:np) = finint(1:np) + infint(1:np) ! omega
     totlap(1:np) = beta(1) + beta(2)*p + gamma*totlap(1:np)/2.0 ! f
     totlap(1:np) = totlap(:)/(1.0 + p(:)*totlap(:))  ! Psi
     totint = dehoog_invlap(td(i),tee,totlap(1:np),lap) 

     ! write results to file (as we go)
     if (dimless) then
        write (20,'(3(1x,ES24.15E3))') td(i), totint
     else
        write (20,'(3(1x,ES24.15E3))') ts(i), totint
     end if
     
  end do
  deallocate(finint,infint,ts,td,splitv,dt)

contains

  !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
  !!   end of program flow
  !! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !! ###################################################
  function inf_integral(integrand,split,num,ord,zeros,GLx,GLw) result(integral)
    use constants, only : DP,EP
    use shared_data, only : rdw, lap
    implicit none

    interface 
       function integrand(a) result(H)
         use constants, only : EP
         use shared_data, only : lap
         real(EP), intent(in) :: a
         complex(EP), dimension(2*lap%M+1) :: H
       end function integrand
    end interface

    integer, intent(in) :: split
    integer, intent(in) :: num,ord
    real(DP), dimension(:), intent(in) :: zeros

    real(EP), dimension(ord-2) :: GLy
    complex(EP), dimension(ord-2,2*lap%M+1) :: GLz
    real(EP), dimension(:), intent(in) :: GLx,GLw
    complex(EP), dimension(num,2*lap%M+1) :: area
    real(EP) :: lob,hib,width
    complex(EP), dimension(2*lap%M+1) :: integral
    integer :: i, j, k, np

    np = 2*lap%M+1

    do j = split+1, split+num
       lob = real(zeros(j-1)/rdw,EP) ! lower bound
       hib = real(zeros(j)/rdw,EP)   ! upper bound
       width = real(hib,EP) - real(lob,EP)
       
       ! transform GL abcissa to global coordinates
       GLy(:) = (width*GLx(:) + (hib + lob))/2.0
          
       !$OMP PARALLEL DO PRIVATE(k) SHARED(GLz,GLy,ord,np)
       do k=1,ord-2
          GLz(k,1:np) = integrand( GLy(k) )
       end do
       !$OMP END PARALLEL DO
          
       area(j-split,1:np) = width/2.0* sum(GLz(1:ord-2,:)*spread(GLw(1:ord-2),2,np),dim=1)
    end do
 
    !$OMP PARALLEL DO PRIVATE(i) SHARED(area,integral,num)
    do i=1,np
       ! accelerate each series independently
       integral(i) = wynn_epsilon(area(1:num,i))
    end do
    !$OMP END PARALLEL DO
    
  end function inf_integral

  !! ###################################################
  subroutine tanh_sinh_setup(k,s,w,a)
    use constants, only : PIOV2EP
    implicit none
    
    integer, intent(in) :: k
    real(DP), intent(in) :: s
    real(EP), intent(out), dimension(2**k-1) :: w, a

    integer :: N,r,i
    real(EP) :: h
    real(EP), allocatable :: u(:,:)

    !! compute weights 
    N = 2**k-1
    r = (N-1)/2
    h = 4.0_EP/2**k
    allocate(u(2,N))
       
    u(1,1:N) = PIOV2EP*cosh(h*[(i, i=-r,r)])
    u(2,1:N) = PIOV2EP*sinh(h*[(i, i=-r,r)])
    
    a(1:N) = tanh(u(2,1:N))
    w(1:N) = u(1,1:N)/cosh(u(2,:))**2
    w(1:N) = 2.0_EP*w(1:N)/sum(w(1:N))

    ! map the -1<=x<=1 interval onto 0<=a<=s
    a = (a + 1.0_EP)*s/2.0_EP

  end subroutine tanh_sinh_setup


  !! ###################################################
  subroutine gauss_lobatto_setup(ord,abcissa,weight)
    use constants, only : PIEP, EP
    implicit none
    integer, intent(in) :: ord
    real(EP), intent(out), dimension(ord-2) :: weight, abcissa
    real(EP), dimension(ord,ord) :: P
    real(EP), dimension(ord) :: x,xold,w
    integer :: i, N, N1, k

    ! leave out the endpoints (abcissa = +1 & -1), since 
    ! they will be the zeros of the Bessel functions
    ! (therefore, e.g., 5th order integration only uses 3 points)
    ! code copied from Matlab routine by Greg von Winckel, at
    ! http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=4775

    N = ord-1
    N1 = N+1
    ! first guess
    x = cos(PIEP*[(i, i=0,N)]/N)
    ! initialize Vandermonde matrix
    P = 0.0
    xold = 2.0

    iter: do 
       if (maxval(abs(x-xold)) > spacing(1.0_EP)) then
          xold = x
          P(:,1) = 1.0
          P(:,2) = x
          
          do k = 2,N
             P(:,k+1) = ((2*k-1)*x*P(:,k) - (k-1)*P(:,k-1))/k
          end do
          
          x = xold - (x*P(:,N1)-P(:,N))/(N1*P(:,N1))
       else
          exit iter
       end if
    end do iter
    
    w = 2.0/(N*N1*P(:,N1)**2)

    ! leave off endpoints
    abcissa = x(2:ord-1)
    weight = w(2:ord-1)

  end subroutine gauss_lobatto_setup


  !! ###################################################
  !! wynn-epsilon acceleration of partial sums, given a series
  !! all intermediate sums / differences are done in extended precision
  function wynn_epsilon(series) result(accsum)
    use constants, only : EONE, EP, DP

    implicit none
    integer, parameter :: MINTERMS = 4
    complex(EP), dimension(1:), intent(in) :: series
    complex(DP) :: accsum
    complex(EP) :: denom
    integer :: np, i, j, m
    complex(EP), dimension(1:size(series),-1:size(series)-1) :: eps

    np = size(series)

    ! first column is intiallized to zero
    eps(:,-1) = 0.0

    ! build up partial sums, but check for problems
    check: do i=1,np
       if((abs(series(i)) > huge(0.0D0)).or.(series(i) /= series(i))) then
          ! +/- Infinity or NaN ? truncate series to avoid corrupting accelerated sum
          np = i-1
          if(np < MINTERMS) then
             write(*,'(A)',advance='no') 'not enough Wynn-Epsilon series to accelerate '
             accsum = -999999.9  ! make it clear answer is bogus
             goto 777
          else
             write(*,'(A,I3,A)',advance='no') 'Wynn-Epsilon series&
                  &, truncated to ',np,' terms. '
             exit check
          end if
       else
          ! term is good, continue
          eps(i,0) = sum(series(1:i))
       end if
    end do check
    
    ! build up epsilon table (each column has one less entry)
    do j=0,np-2 
       do m = 1,np-(j+1)
          denom = eps(m+1,j) - eps(m,j)
          if(abs(denom) > tiny(EONE)) then ! check for div by zero
             eps(m,j+1) = eps(m+1,j-1) + EONE/denom
          else
             accsum = real(eps(m+1,j),DP)
             write(*,'(A,I2,1X,I2,A)',advance='no') 'epsilon cancel ',m,j,' '
             goto 777
          end if
       end do
    end do

    ! if made all the way through table use "corner value" of triangle as answer
    if(mod(np,2) == 0) then
       accsum = real(eps(2,np-2),DP)  ! even number of terms - corner is acclerated value
    else
       accsum = real(eps(2,np-3),DP)  ! odd numbers, use one column in from corner
    end if
777 continue

  end function wynn_epsilon

  !! ###################################################
  ! polynomial extrapolation modified from numerical recipes f90 (section 3.1)
  SUBROUTINE polint(xa,ya,x,y,dy)
    ! xa and ya are given x and y locations to fit an nth degree polynomial
    ! through.  results is a value y at given location x, with error estimate dy

    use constants, only : DP, EP
    IMPLICIT NONE
    REAL(EP), DIMENSION(:), INTENT(IN) :: xa
    complex(EP), dimension(:), intent(in) :: ya
    REAL(EP), INTENT(IN) :: x
    complex(EP), INTENT(OUT) :: y,dy
    INTEGER :: m,n,ns
    COMPLEX(EP), DIMENSION(size(xa)) :: c,d,den
    REAL(EP), DIMENSION(size(xa)) :: ho
    integer, dimension(1) :: tmp

    n=size(xa)
    c=ya
    d=ya
    ho=xa-x
    tmp = minloc(abs(x-xa))
    ns = tmp(1)
    y=ya(ns)
    ns=ns-1
    do m=1,n-1
       den(1:n-m)=ho(1:n-m)-ho(1+m:n)
       if (any(abs(den(1:n-m)) == 0.0)) then
          print *, 'polint: calculation failure'
          stop
       end if
       den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
       d(1:n-m)=ho(1+m:n)*den(1:n-m)
       c(1:n-m)=ho(1:n-m)*den(1:n-m)
       if (2*ns < n-m) then
          dy=c(ns+1)
       else
          dy=d(ns)
          ns=ns-1
       end if
       y=y+dy
    end do
  END SUBROUTINE polint

end program Driver
