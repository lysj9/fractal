	include "for_plot.f90"
	use for_plot
	implicit none
	integer,parameter::n=100000
	real::x(n)
	integer::i,ios,m=0

	do i=1,n
		read (*,*,iostat=ios) x(i)
		if ( ios/=0 ) exit
		m=m+1
	end do

	call histbin(x(1:m),50)
	end
