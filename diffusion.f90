program diffusion 
  implicit none
  ! define parameters and variables
  integer, parameter :: N=100, Nmax=1000
  real*4,  parameter :: D=2.,t0=2.
  real*4   delta_t, delta_x, A, alpha, t, acc, L, eps
  real*4,  dimension(N) :: x, rho_new, rho, analytical, u, ucheck, ones, p, q, r, s, pivot
  real*4,  dimension(N) :: acc_array
  real*4,  dimension(Nmax) :: integral, tlist
  real*4,  dimension(Nmax,N) :: res, errorres
  real*4,  dimension(N,N) :: M
  integer  i, j, k, counter, rc, stepplot
  integer  BC

  ! read in the boundary conditions
  write(6,*) 'Enter the boundary conditions for the  one-dimensional diffusion equation: Dirichlet = 1, Neumann = 2'
  read(5,*) BC

  L = 10.
  delta_x = L/(N-1) 
  rho = 0.
  rho_new = 0.    
  delta_t = (delta_x*delta_x)/(2*D) 
  eps = 1.0e-6  
  stepplot = 10

  ! set initial values
  do i=1,N
    x(i) = (i-1)*delta_x-L/2  
  enddo

  do i=1,N
    rho(i) = initial(x(i))
  enddo

  res(1,:) = rho
  counter = 1
  t = t0 + delta_t
  tlist(1) = 2.
  integral(1) = simpson(rho,-5.,5.,N-1)

 ! boundary conditions
  if (BC.eq.1) then   		! Dirichlet BC
  do while(t .le. 4.0)
	! define alpha and A = (1 + 2*alpha)
	alpha = delta_t*D/(2*delta_x*delta_x)	  
	A = (1 + 2*alpha) 	  
	
	! determin u the right-hand side
	do i=1,N
	if (i == 1) then 
		u(i) = density(x(1),t)
	else if (i == N) then 
		u(i) = density(x(N),t)	
	else
		u(i) = alpha*rho(i-1)+(1-2*alpha)*rho(i)+alpha*rho(i+1)
	endif 
	enddo

	! determin the matrix for each time step in matrix and vector form
	! matrix form 
	M = 0.
	M = reshape(M,[N,N])
	do i = 1,N
	M(i,i) = A
	if(i .lt. N) then 
		M(i,i+1) = -alpha
		M(i+1,i) = -alpha
	endif
	enddo

	M(1,1)  = 1.
	M(N,N)  = 1.
	M(1,2)  = 0.
	M(N,N-1)= 0.	 

	! vector form  
	ones =  1.
	p    = -alpha*ones
	q    =  A*ones
	r    = -alpha*ones
	p(1) =  0.
	r(N) =  0.
	q(1) =  1.
	q(N) =  1.
	r(1) =  0. 
	p(N) =  0.

	! solve the system by the thomas algorithm and by the lapack function for testing reasons
	rho_new = thomasVec(p,q,r,u)
	call sgesv(N, 1, M, N, pivot, u, N, rc)

	! run test between thomas algorithm and lapack function (eps = 1.0e-6 )
	do i=1,N
		if(abs(rho_new(i)-u(i)) .gt. eps) then
			print*, 'Rho check encountered an error!'
			stop
		endif  
	enddo 

	! calculate the analytical solution  
	do i=1,N
	analytical(i) = density(x(i),t)
	enddo   

	! calculate the mass by the simpson rule and store the values for each time step
	integral(counter+1) = simpson(rho_new,-5.0,5.0,N-1)
	! calulate the relative error for each time step and store the values
	errorres(counter,:) = abs(((rho_new-analytical)/analytical)*100.0)

	! calculate the accuracy
	do k=2,N-1    
	  acc_array(k) = abs(rho_new(k)-rho(k))/minval([rho(k),rho_new(k)])
	enddo
	acc = maxval(acc_array)

	! store the rho values for the new time step
	rho = rho_new

	! determine the accuracy configuration
	if (acc.lt. 0.15) then
	 	t = t + delta_t
		counter = counter + 1	
		res(counter,:) = rho_new	
		tlist(counter) = t
	else 
	  	delta_t = delta_t/2
	 	t = t + delta_t
	endif  
  enddo
  else if (BC.eq.2) then  	! Neumann BC
  do while(t .le. 4.0)
	! define alpha and A = (1 + 2*alpha)
	alpha = delta_t*D/(2*delta_x*delta_x)	  
	A = (1 + 2*alpha) 

	! determin u the right-hand side
	do i=1,N
		if (i == 1) then 
			u(i) = (1-2*alpha)*rho(i)+2*alpha*rho(i+1)
		else if (i == N) then 
			u(i) = 2*alpha*rho(i-1)+(1-2*alpha)*rho(i)
		else
			u(i) = alpha*rho(i-1)+(1-2*alpha)*rho(i)+alpha*rho(i+1)
		endif 
    enddo

	! determin the matrix for each time step in matrix and vector form
	! matrix form 
	M = 0.
	M = reshape(M,[N,N])
	do i = 1,(N)
		M(i,i) = A
		if(i .lt. N) then 
			M(i,i+1) = -alpha
			M(i+1,i) = -alpha
		endif 
	enddo

	M(1,2)  = -2.*alpha
	M(N,N-1)= -2.*alpha	 
	
	! vector form 
	ones =   1.
	p	 =  -alpha*ones
	q 	 =   A*ones
 	r 	 =  -alpha*ones
	p(1) =   0.
	r(N) =   0.
	r(1) =  -2.*alpha 
	p(N) =  -2.*alpha

	! solve the system by the thomas algorithm and by the lapack function for testing reasons
	rho_new = thomasVec(p,q,r,u)
	call sgesv(N, 1, M, N, pivot, u, N, rc)

	! run test between thomas algorithm and lapack function (eps = 1.0e-6)
	do i=1,N
		if(abs(rho_new(i)-u(i)) .gt. eps) then
			print*, 'Rho check encountered an error!'
			stop
		endif  
	enddo   

	! calculate the analytical solution  
	do i=1,N
		analytical(i) = density(x(i),t)
	enddo   

	! calculate the mass by the simpson rule and store the values for each time step
	integral(counter+1) = simpson(rho_new,-5.0,5.0,N-1)
	! calulate the relative error for each time step and store the values
	errorres(counter,:) = abs(((rho_new-analytical)/analytical)*100.0)

	! calculate the accuracy
	do k=2,N-1    
  	    acc_array(k) = abs(rho_new(k)-rho(k))/minval([rho(k),rho_new(k)])
	enddo
	acc = maxval(acc_array)

	! store the rho values for the new time step
	rho = rho_new

	! determine the accuracy configuration
	if (acc.lt. 0.15) then
	 	t = t + delta_t
		counter = counter + 1	
		res(counter,:) = rho_new	
		tlist(counter) = t
	else 
	  	delta_t = delta_t/2
	 	t = t + delta_t
	endif  
	
  enddo
  else
  	print*,'no valid boundary condition'
  endif 

  ! print steps, rho(x,t=4.0), analytical solution (t=4.0), relative errors(t=4.0) and total mass
  print*,'steps:'
  print*, counter
  print*,'numerical solution for t=4.0:'
  print*, res(counter,:) 
  print*,'analytical solution for t=4.0:'
  print*, analytical	
  print*,'relative errors for t=4.0:'
  print*, errorres(counter-1,:)
  print*,'total mass as a function of x (simpsons rule):'
  print*, integral(1:counter)
  !print*, tlist(1:counter)

! print out the plotting variables
  open(11,file='graphplot.txt',status='unknown')
  open(12,file='integral.txt',status='unknown')
  open(13,file='error.txt',status='unknown')
  k = 0.
  do j=1,counter-1
	if ((mod(j,stepplot) == 0) .or. (t .lt. 4.0)) then		
		do i=1,N
			write(11,*) x(i), tlist(j), res(j,i) 
			write(13,*) x(i), tlist(j), errorres(j,i)
		enddo	
	endif
	write(12,*) tlist(j), integral(j)
  enddo  

  open(14,file='graphplotana.txt',status='unknown')
  open(15,file='graphplotana_N.txt',status='unknown')
  do j=1,N
		write(14,*) x(j), analytical(j)
		write(15,*) x(j), res(counter,j)
  enddo

contains
! analytical solution
function density(x,t) result(dens)
	implicit none
	real*4, parameter :: t0 = 2.0, D = 2.0
	real*4 t, x, dens
	dens = sqrt(t0/t)*exp(-(x*x)/(4*D*t))
end function density

! initial values 
function initial(x) result(dens)
	implicit none
	real*4, parameter :: t0 = 2.0, D = 2.0
	real*4 x, dens
	dens = exp(-(x*x)/(4*D*t0))
end function initial

! thomas algorithm 
function thomasVec(p,q,r,s) result(res)
	implicit none
	integer i, k
	real*4, dimension(N) :: p, q, r, s, res, inits
    integer, parameter :: N = size(s)
    
	inits = s
	do k = 2, N
		a    = p(k)/q(k-1)
		q(k) = q(k) - a*r(k-1)
		inits(k) = inits(k) - a*inits(k-1)
	enddo
	res(N) = inits(N)/q(N)
	do k = N-1,1,-1		
		res(k) = (inits(k)-r(k)*res(k+1))/q(k)
	enddo 
end function thomasVec

! simpson's rule that combines the 1/3 rule and the 3/8 rule
function simpson(f,a,b,N) result(integral)
implicit none
integer N, i
real*4  a, b, integral, s
real*4  h, x
real*4, dimension(N) :: f

s = 0.0
h = (b-a)/N
do i=3, N-4, 2
   s   = s + 2.0*f(i) + 4.0*f(i+1)
end do
integral = (s + f(1) + 4.0*f(2))*h/3.0 + (3*h/8)*(f(N-3)+ 3*f(N-2) + 3*f(N-1) + f(N))
return
end function simpson
end program diffusion

