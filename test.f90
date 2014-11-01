	program test
	use gnufor2
	implicit none
!***********************************************************************************
	integer, parameter	:: N1=50, N2=100, N3=200, N4=500, Nx=80, Ny=40, Nm=800
	real(kind=8)		:: x1(N1), x2(N2), x3(N3), x4(N4)
	real(kind=8)		:: f1(N1), f2(N2), f3(N3), f4(N4)
	real(kind=8)		:: u1(1000), u2(1000), u3(1000), t(1000)
	real(kind=8)		:: xx(Nx), yy(Ny), zz(Nx,Ny), xyz(3,Nx,Ny), u(Nx), v(Nx)
	real(kind=4)		:: z_fractal(Nm,Nm), x_fractal(Nm), y_fractal(Nm), &
		& xstart,ystart,xend,yend,zmax,escape,color_scale,z_power
	integer			:: rgb(3,Nm,Nm)
	complex(kind=8)		:: z, a
	integer			:: i,j,k, N0
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
! generate data for 2D plots
	do  i=1,N1
		x1(i)=5.0*i/N1
	end do
	do  i=1,N2
		x2(i)=6.0*i/N2
	end do
	do  i=1,N3
		x3(i)=7.0*i/N3
	end do
	do  i=1,N4
		x4(i)=8.0*i/N4
	end do
	f1=sin(2*x1)
	f2=4*x2/(1+x2**2)
	f3=4-x3
	f4=exp(x4/3)-5
!***********************************************************************************
! generate data for 3D plot of a curve
	do i=1,1000
		t(i)=80.0*i/1000
	end do
	u1=4*cos(t)-cos(5*t)
	u2=4*sin(t)-sin(5*t)
	u3=10*sin(t/7)
!***********************************************************************************
! generate data for surface plots
	do i=1,Nx
		xx(i)=6.0*(dble(i)/Nx-0.5)
	end do
	do i=1,Ny
		yy(i)=4.0*(dble(i)/Ny-0.5)
	end do
	do i=1,Nx
	do j=1,Ny
		zz(i,j)=1-exp(-(xx(i)-1)**2-yy(j)**2)
	end do
	end do
!***********************************************************************************
! generate data for surface plots
	do i=1,Nx
		u(i)=1.5*3.15*(dble(i)/Nx)
	end do
	do i=1,Ny
		v(i)=1.6*3.15*(dble(i)/Ny)+0.7
	end do
	do i=1,Nx
	do j=1,Ny
		xyz(1,i,j)=(3+cos(v(j)))*cos(u(i))
	 	xyz(2,i,j)=(3+cos(v(j)))*sin(u(i))
		xyz(3,i,j)=sin(v(j))
		xyz(1,i,j)=xyz(1,i,j)-0.5*xyz(2,i,j)-0.5*xyz(3,i,j)
		xyz(3,i,j)=xyz(2,i,j)-5*xyz(3,i,j)
	end do
	end do
!***********************************************************************************
! create date for fractal plot
	xstart=-2.0
	xend=1.0
	ystart=-1.5
	yend=1.5
	escape=0.85
	z_power=2.0
	color_scale=0.2
	zmax=exp(log(10.0)*escape)
	N0=150
	z_fractal=1.0
	do i=1,Nm
		x_fractal(i)=(xstart+(xend-xstart)*(real(i-1)/(Nm-1)))
	do j=1,Nm
		y_fractal(j)=(ystart+(yend-ystart)*(real(j-1)/(Nm-1)))
		a=x_fractal(i)+y_fractal(j)*(0.0,1.0)
		z=a
		k=1
		do while ((abs(z)<zmax).and.(k<N0))
			k=k+1
			z=z*z+a
		end do
		if (k<N0) then
			z_fractal(i,j)=k+(log(abs(z))/(escape*log(10.0))-1.0)/(1.0-z_power)
			z_fractal(i,j)=0.5*(1.0+sin(z_fractal(i,j)*color_scale))
		end if
		rgb(1,i,j)=floor(z_fractal(i,j)*255)
		rgb(2,i,j)=floor(255*(atan(aimag(z)/real(z))/1.58+1)/2.0)
		rgb(3,i,j)=floor(255*(1+cos(abs(z)/5.0))/2.0)
	end do
	end do
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 1: simple 2D graph'
	print *,'call plot(x1,f1)'
	call plot(x1,f1)
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 2: style parameter - plotting with lines and points'
	print *,'call plot(x1,f1,''10-'')'
	call plot(x1,f1,'10-')
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 3: style parameter - plotting with points only'
	print *,'call plot(x1,f1,'' 5.'')'
	call plot(x1,f1,' 5.')
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 4: pause 5 seconds, then clear the graph'
	print *,'call plot(x1,f1,''15-'',pause=5.0,persist=''no'')'
	call plot(x1,f1,'15-',pause=5.0,persist='no')
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 5: pause until you hit enter, do not clear the graph'
	print *,'call plot(x1,f1,''15-'',pause=-1.0)'
	call plot(x1,f1,'15-',pause=-1.0)
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 6: specifying color'
	print *,'call plot(x1,f1,''21-'',pause=-1.0,persist=''no'',color1=''red'')'
	call plot(x1,f1,'21-',pause=-1.0,persist='no',color1='red')
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 7: specifying terminal'
	print *,'call plot(x1,f1,''22-'',pause=5.0,persist=''no'',terminal=''x11'')'
	call plot(x1,f1,'22-',pause=-1.0,persist='no',terminal='x11')
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 8: saving to jpeg file'
	print *,'the output file is in the same directory as test.out' 
	print *,'call plot(x1,f1,''22-'',terminal=''jpeg'')'
	call plot(x1,f1,'22-',terminal='jpeg')
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 9: saving to ps file, specifying the output file name'
	print *,'the output file is in the same directory as test.out '
	print *,'call plot(x1,f1,''22-'',terminal=''ps'',filename=''picture'')'
	call plot(x1,f1,'22-',terminal='ps',filename='picture')
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 10: specifying the linewidth'
	print *,'call plot(x1,f1,'' 0-'',pause=-2.0,persist=''no'',linewidth=3.5)'
	call plot(x1,f1,' 0-',pause=-2.0,persist='no',linewidth=3.5)
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 11: specifying the specific command/data file name'
	print *,'call plot(x1,f1,'' 8-'',pause=-2.0,persist=''no'',input=''my_filename'')'
	call plot(x1,f1,' 8-',pause=-2.0,persist='no',input='my_filename')
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 12: 2D graphs in polar coordinates'
	print *,'call plot(x1,f1,'' 8-'',pause=-2.0,persist=''no'',polar=''yes'')'
	call plot(x1,f1,' 8-',pause=-2.0,persist='no',polar='yes')
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 13: several graphs in the same coordinate system'
	print *,'call plot(x1,f1,x2,f2,x3,f3,pause=-2.0,persist=''no'')'
	call plot(x1,f1,x2,f2,x3,f3,pause=-2.0,persist='no')
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 14: several graphs in the same coordinate system'
	print *,'call plot(x1,f1,x2,f2,x3,f3,x4,f4,&
		&'' 5.22- 0- 0-'',color2=''dark-yellow'',color1=''#40e0d0'',pause=-1.0,persist=''no'')'
	call plot(x1,f1,x2,f2,x3,f3,x4,f4,' 5.22- 0- 0-',color2='dark-yellow',color1='#40e0d0',pause=-1.0,persist='no')
! ***********************************************************************************
	print *,'********************************************************************'
	print *,'Example 15: several graphs in the same polar coordinate system'
	print *,'call plot(x1,f1,x2,f2,x3,f3,x4,f4,pause=-1.0,persist=''no'',polar=''yes'',linewidth=2.0)'
	call plot(x1,f1,x2,f2,x3,f3,x4,f4,pause=-1.0,persist='no',polar='yes',linewidth=2.0)
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 16: histogram with 10 bins'
	print *,'call hist(f1,10,color=''#779944'',pause=-1.0,persist=''no'')'
	call  hist(f1,10,color='#779944',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 17: 3D plot of a curve'
	print *,'call plot3D(u1,u2,u3,pause=-1.0,persist=''no'',linewidth=2.0)'
	call  plot3D(u1,u2,u3,pause=-1.0,persist='no',linewidth=2.0)
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 18: surface plot'
	print *,'call surf(xx,yy,zz)'
	call  surf(xx,yy,zz)
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 19: surface plot without specying x and y grid'
	print *,'call surf(zz)'
	call  surf(zz)
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 20: parametric surface plot'
	print *,'call surf(xyz)'
	call  surf(xyz)
	print *,'press ENTER to go to the next example'
	read  *
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 21: parametric surface plot'
	print *,'call surf(xyz,pm3d=''pm3d'',pause=-1.0,persist=''no'')'
	call  surf(xyz,pm3d='pm3d',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 22: parametric surface plot, different palette'
	print *,'call surf(xyz,pm3d=''pm3d'',palette=''RGB'',pause=-1.0,persist=''no'')'
	call  surf(xyz,pm3d='pm3d',palette='RGB',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 23: parametric surface plot, different palette'
	print *,'call surf(xyz,pm3d=''pm3d'',palette=''gray negative'',pause=-1.0,persist=''no'')'
	call  surf(xyz,pm3d='pm3d',palette='gray negative',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'*******************************************************************'
	print *,'Example 24: parametric surface plot, different palette'
	print *,'call surf(xyz,pm3d=''pm3d'',palette=''rgbformulae 23,28,3'',pause=-1.0,persist=''no'')'
	call  surf(xyz,pm3d='pm3d',palette='rgbformulae 23,28,3',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 25: surface plot'
	print *,'call surf(xx,yy,zz,pm3d=''pm3d at b'',palette=''HSV'',pause=-1.0,persist=''no'')'
	call  surf(xx,yy,zz,pm3d='pm3d at b',palette='HSV',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 26: surface plot'
	print *,'call surf(xx,yy,zz,pm3d=''pm3d implicit map'',palette=''rgbformulae 31,-11,32'',pause=-1.0,persist=''no'')'
	call  surf(xx,yy,zz,pm3d='pm3d implicit map',palette='rgbformulae 31,-11,32',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 27: surface plot with contour on the surface'
	print *,'call surf(xx,yy,zz,contour=''surface'',palette=''YIQ'',pause=-1.0,persist=''no'')'
	call  surf(xx,yy,zz,contour='surface',palette='YIQ',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 28: surface plot with contour on the XY plane'
	print *,'call surf(xx,yy,zz,contour=''xy'',pause=-1.0,persist=''no'')'
	call  surf(xx,yy,zz,contour='xy',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 29: surface plot with contour on the surface and the XY plane'
	print *,'call surf(xx,yy,zz,contour=''both'',pause=-1.0,persist=''no'')'
	call  surf(xx,yy,zz,contour='both',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 30: image plot'
	print *,'call image(z_fractal,palette=''gray'',pause=-1.0,persist=''no'')'
	call  image(z_fractal,palette='gray',pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 31: image plot with specified x and y axes'
	print *,'call image(x_fractal,y_fractal,z_fractal,pause=-1.0,persist=''no'')'
	call  image(x_fractal,y_fractal,z_fractal,pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 32: RGB image plot'
	print *,'call  image(rgb,pause=-1.0,persist=''no'')'
	call  image(rgb,pause=-1.0,persist='no')
! ***********************************************************************************
 	print *,'********************************************************************'
	print *,'Example 33: RGB image plot with specified x and y axes'
	print *,'call  image(x_fractal,y_fractal,rgb,pause=-1.0,persist=''no'')'
	call  image(x_fractal,y_fractal,rgb,pause=-1.0,persist='no')
!***********************************************************************************
!***********************************************************************************
!***********************************************************************************
	end program test