!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!
! MAE 6263 Computational Fluid Dynamics - Project 3 !
! Sk. Mashfiqur Rahman !
!-----------------------------------------------------------------------------!
!-----------------------------------------------------------------------------!

program channel
implicit none
integer :: nx,ny,i,j,opt,nitr,ns
real*8 :: Re,pe,Tmax,nu,dy,dx,dt,rho,tolt,m
real*8 :: t1,t2,lx,ly,tol,gamma,tol1,en
real*8,allocatable ::u(:,:),v(:,:),p(:,:),hx(:,:),hy(:,:)

common/Reynolds/Re
common/contol/tol
common/density/rho

open(9,file='project3input.txt')
read(9,*)lx		    !length in x -direction
read(9,*)ly	    	!length in y -direction
read(9,*)nx		    !resolution in x direction
read(9,*)ny	    	!resolution in y direction
read(9,*)ns	    	!number of snapshots
read(9,*)gamma		!value of gamma to measure time-step size
read(9,*)rho	    !density
read(9,*)nu	    	!kinematic viscosity
read(9,*)pe	    	!Peclet number
read(9,*)Tmax		!final time
read(9,*)tol        !tolerance for iterative solvers
read(9,*)tolt       !tolerance for time integration
read(9,*)Re         !Reynolds number
read(9,*)opt    	![0]PPE without momentum interpolation, [1]PPE with momentum interpolation
close(9)     

!Calculating grid spacing (spatial)
dx = lx/dfloat(nx)
dy = ly/dfloat(ny)
print*, dx
print*, dy

!Time step
!dt = ((dx*dx)*gamma) /nu
dt = Tmax/dfloat(ns)

!Initial Condition & Boundary Condition Setup
allocate(u(0:nx,0:ny))
allocate(v(0:nx,0:ny))
allocate(p(0:nx,0:ny))
allocate(hx(0:nx,0:ny))
allocate(hy(0:nx,0:ny))

do i = 0,nx
  do j = 0,ny
    u(i,j) = 0
    v(i,j) = 0
    p(i,j) = 0
  end do
end do

do j = 0,ny
    p(0,j) = 100.16d0
    p(nx,j) = 100.0d0    
    
end do

nitr = 0
!IC File Ouput
open(20,file="InitialCondition.plt")
write(20,*)'Title="IC Data set"'
write(20,*)'variables ="x","y","u","v","p"'
close(20)

open(20,file="InitialCondition.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.3)',advance="no")dfloat(i)/dfloat(ny),dfloat(j)/dfloat(ny),u(i,j),v(i,j),p(i,j)
      write(20,*) ''
    end do
  end do
close(20)

!Output file setup
open(20,file="ContourPlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="x","y","u","v","p"'
close(20)

open(20,file="LinePlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="x","p"'
close(20)

open(20,file="UprofileLinePlots.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="u","y"'
close(20)

open(20,file="energy.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="e","t"'
close(20)

open(20,file="mass.plt")
write(20,*)'Title="Transient data set"'
write(20,*)'variables ="m","t"'
close(20)

call cpu_time(t1)

!-----------------------------------------------------------------------------!
!Time integration
!-----------------------------------------------------------------------------!

tol1 = 1.0d0  
do while (tol1.gt.tolt)
  tol1 = 0.0d0
  nitr = nitr + 1
call RHS(u,v,hx,hy,nx,ny,dx,dy)

if (opt==0) then
call compact(p,hx,hy,nx,ny,dx,dy,tol)
else 
call momentum(p,u,v,dt,nx,ny,dx,dy,tol,nitr)
end if

call velocity(u,v,p,hx,hy,nx,ny,dx,dy,dt,tol1)

!Output .plt file for Tecplot  
if (mod(nitr,500)==0) then
 
open(20,file="ContourPlots.plt",position="append")
write(20,"(a,i8,a,i8,a)")'Zone I = ',nx+1,',J=',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',nitr
  do j = 0,ny
    do i = 0,nx
      write (20, '(1600F14.6,1600F14.6,1600F14.6,1600F14.6)',advance="no")dfloat(i)/dfloat(nx/2),&
      &dfloat(j)/dfloat(ny),u(i,j),v(i,j),p(i,j)
      write(20,*) ''
    end do
  end do
close(20)

open(20,file="LinePlots.plt",position="append")
write(20,"(a,i8,a)")'Zone I = ',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',nitr
    do j = 0,ny
      write (20, '(1600F14.6,1600F14.6)',advance="no")dfloat(j)/dfloat(ny),p(nx/2,j)
      write(20,*) ''
    end do
close(20)

open(20,file="UprofileLinePlots.plt",position="append")
write(20,"(a,i8,a)")'Zone I = ',ny+1,',F=POINT'
write(20,"(a,i8)")'StrandID=0,SolutionTime=',nitr
    do j = 0,ny
      write (20, '(1600F14.6,1600F14.6)',advance="no")u(nx/2,j),dfloat(j)/dfloat(ny)
      write(20,*) ''
    end do
close(20)

en = 0.0d0
m = 0.0d0

open(20,file="energy.plt",position="append")
    do i = 0,nx
      do j = 0,ny-1
        en = en + 0.5d0*rho*u(i,j)*u(i,j)
        end do
        end do
      write (20,*) dfloat(nitr),en
      write(20,*) ''
close(20)

open(20,file="mass.plt",position="append")
      do j = 1,ny-1
        m = m + 0.5d0*(u(nx/2+1,j)+u(nx/2+1,j+1))
        end do
      write (20,*) dfloat(nitr),m*dy
      write(20,*) ''
close(20)
end if

end do

call cpu_time(t2)

open(4,file='cpu.txt')
write(4,*)"cpu time (sec)=",(t2-t1)
close(4)

end


!-----------------------------------------------------------------------------!
!Convection-diffusion
!-----------------------------------------------------------------------------!
subroutine RHS(u,v,hx,hy,nx,ny,dx,dy)
implicit none

integer :: nx,ny,i,j
real*8::dx,dy,Re
real*8,dimension(0:nx,0:ny):: u,v,hx,hy
real*8,dimension(-1:nx+1,-1:ny+1)::u2,v2
real*8,dimension(0:nx,0:ny)::udx,udy,uddx,uddy,vdx,vdy,vddx,vddy

common/Reynolds/Re

do i = 0,nx
  do j = 0,ny
    u2(i,j) = u(i,j)
    v2(i,j) = v(i,j)
  end do
end do

do j = 0,ny
  u2(-1,j) = u(nx-1,j)
  u2(nx+1,j) = u(1,j)
  v2(nx+1,j) = v(1,j)
  v2(-1,j) = v(nx-1,j)
end do
 
do i = 0,nx
  	do j = 1,ny-1
    if(u(i,j)>=0) then
  		udx(i,j) = (u2(i,j)-u2(i-1,j))/dx
        vdx(i,j) = (v2(i,j)-v2(i-1,j))/dx
    else if(u(i,j)<0) then
        udx(i,j) = (u2(i+1,j)-u2(i,j))/dx
        vdx(i,j) = (v2(i+1,j)-v2(i,j))/dx
    end if   
	end do
    udx(i,0) = udx(i,1)
    vdx(i,0) = vdx(i,1)
    udx(i,ny) = udx(i,ny-1)
    vdx(i,ny) = vdx(i,ny-1)
end do

do i = 0,nx
  	do j = 1,ny-1
    if(v(i,j)>=0) then
  		udy(i,j) = (u2(i,j)-u2(i,j-1))/dy
        vdy(i,j) = (v2(i,j)-v2(i,j-1))/dy
    else if(v(i,j)<0) then
      udy(i,j) = (u2(i,j+1)-u2(i,j))/dy
      vdy(i,j) = (v2(i,j+1)-v2(i,j))/dy
    end if
	end do
    udy(i,0) = udy(i,1)
    vdy(i,0) = vdy(i,1)
    udy(i,ny) = udy(i,ny-1)
    vdy(i,ny) = vdy(i,ny-1)
end do

do i = 0,nx
  	do j = 1,ny-1
  		uddx(i,j) = (u2(i+1,j)+u2(i-1,j)-2.*u2(i,j))/(dx**2)
        vddx(i,j) = (v2(i+1,j)+v2(i-1,j)-2.*v2(i,j))/(dx**2)
	end do
    uddx(i,0) = uddx(i,1)
    vddx(i,0) = vddx(i,1)
    uddx(i,ny) = uddx(i,ny-1)
    vddx(i,ny) = vddx(i,ny-1)
end do

do i = 0,nx
  	do j = 1,ny-1
  		uddy(i,j) = (u2(i,j+1)+u2(i,j-1)-2d0*u2(i,j))/(dy**2)
        vddy(i,j) = (v2(i,j+1)+v2(i,j-1)-2d0*v2(i,j))/(dy**2)
	end do
    uddy(i,0) = uddy(i,1)
    vddy(i,0) = vddy(i,1)
    uddy(i,ny) = uddy(i,ny-1)
    vddy(i,ny) = vddy(i,ny-1)
end do

do i = 0,nx
	do j = 0,ny
    	hx(i,j) = -u2(i,j)*(udx(i,j))-v2(i,j)*(udy(i,j))+1/Re*(uddx(i,j)+uddy(i,j))
    	hy(i,j) = -u2(i,j)*(vdx(i,j))-v2(i,j)*(vdy(i,j))+1/Re*(vddx(i,j)+vddy(i,j))
        
	end do
end do

return
end

!-----------------------------------------------------------------------------!
!Pressure
!-----------------------------------------------------------------------------!
subroutine compact(p,hx,hy,nx,ny,dx,dy,tol)
implicit none

integer :: nx,ny,i,j
real*8::dy,dx,a,b,tol,err
real*8,dimension(0:nx,0:ny):: p,hx,hy,hxdx,hydy,hy2,e
real*8,dimension(0:nx,-1:ny+1)::p2,p1
real*8,dimension(-1:nx+1,0:ny)::hx2

err=1.0d0
  
p1(0:nx,0:ny) = p

do i = 0,nx
  do j = 0,ny
    hx2(i,j) = hx(i,j)
    hy2(i,j) = hy(i,j)   
  end do
end do

do j = 0,ny
  hx2(-1,j) = hx(nx-1,j)
  hx2(nx+1,j) = hx(1,j)
end do
   
do i = 0,nx
  do j = 1,ny-1
    hxdx(i,j) = (hx2(i+1,j)-hx2(i-1,j))/2.0d0
    hydy(i,j) = (hy2(i,j+1)-hy2(i,j-1))/2.0d0
  end do
    hxdx(i,0)=hxdx(i,1)
    hxdx(i,ny)=hxdx(i,ny-1)
    hydy(i,0)=hydy(i,1)
    hydy(i,ny)=hydy(i,ny-1)
end do  

do i=0,nx
  p1(i,-1) = p1(i,1)
  p1(i,ny+1) = p1(i,ny-1)
end do  
  
do while(err.gt.tol)
err=0.0d0
p2 = p1
do i = 1,nx-1
  do j = 0,ny
      a = p1(i+1,j)+p2(i-1,j)+p2(i,j-1)+p1(i,j+1)
      b = hxdx(i,j) + hydy(i,j)
      p2(i,j) = a/4-(b*dx)/4
  end do 
    p2(i,-1) = p2(i,1)
	p2(i,ny+1) = p2(i,ny-1)
end do	
           
    !compute L1 norm
    do j=0,ny
    do i=0,nx
    e(i,j) = dabs(p2(i,j)-p(i,j))
    err =  e(i,j) + err
    end do
    end do 	  


p1 = p2
do j = 0,ny
  do i = 0,nx
        p(i,j) = p1(i,j)        
  end do
end do
end do

return
end

!-----------------------------------------------------------------------------!
!velocity update
!-----------------------------------------------------------------------------!
subroutine velocity(u,v,p,hx,hy,nx,ny,dx,dy,dt,tol1)
implicit none

common/density/rho

real*8::dx,dy,dt,rho,tol1
integer::nx,ny,i,j
real*8,dimension(0:nx,0:ny)::u,v,p,hx,hy,p2,s
real*8,dimension(0:nx,0:ny)::gradpx,gradpy,uo,vo

do i = 0,nx
  do j = 0,ny
    p2(i,j) = p(i,j)
  end do
end do

do i=0,nx
  do j=0,ny
    gradpx(i,j)=0.0
    gradpy(i,j)=0.0
  end do
end do 

do j = 0,ny
  do i = 0,nx
        uo(i,j) = u(i,j)
        vo(i,j) = v(i,j)        
end do
end do 

do i = 1,nx-1
  do j = 1,ny-1
    gradpx(i,j) = (p2(i+1,j)-p2(i-1,j))/(2.*dx)
  end do
end do

do i = 1,nx-1
  	do j = 1,ny-1
    	gradpy(i,j) = (p2(i,j+1)-p2(i,j-1))/(2.*dy)
	end do
end do

do j = 0,ny
   	gradpx(0,j) = gradpx(1,j)
   	gradpx(nx,j) = gradpx(nx-1,j)         
end do

do i = 0,nx
  do j = 1,ny-1
	u(i,j) = u(i,j) + dt*(hx(i,j)-gradpx(i,j))
  	v(i,j) = v(i,j) + dt*(hy(i,j)-gradpy(i,j))
  end do
end do

    !compute L1 norm
    do j=0,ny
    do i=0,nx
    s(i,j) = dabs(u(i,j)-uo(i,j))+dabs(v(i,j)-vo(i,j))
    tol1 =  s(i,j) + tol1
    end do
    end do
    print*,tol1 	  

return
end
