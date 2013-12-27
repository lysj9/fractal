Module for_plot
    use ISO_FORTRAN_ENV
    integer,parameter,private::PC=kind(1.0)
    contains
    subroutine pl2d(x,y,bin,xlow,xhigh,cr,file_unit)
    implicit none
    real(kind=PC),dimension(:),intent(in)::x,y
    integer,intent(in)::bin
! x-range, optional. if not set, use the range: [min(x),max(x)]
    real(kind=PC),optional,intent(in)::xlow,xhigh
! critical number, optional. only the bin with number above will be chosen.
! if not set, use 0
    integer,optional,intent(in)::cr
    integer::cr_num
! output file unit. if not set, use the standard output unit OUTPUT_UNIT
    integer,optional,intent(in)::file_unit
    integer::out_unit
    real(kind=PC),dimension(bin)::xbin,ybin,ybin_err
    integer,dimension(bin)::num_bin
    real::xmin,xmax,xdel
    integer::i,j
    integer::xydim

    cr_num=0
    if ( present(cr) ) cr_num=cr

    out_unit=OUTPUT_UNIT
    if ( present(file_unit) ) out_unit=file_unit

    xydim=ubound(x,1)-lbound(x,1)+1
    if ( xydim/=ubound(y,1)-lbound(y,1)+1 ) then
        write (ERROR_UNIT,*) "x, y dimension unmacth!"
        stop
    end if

    xmin=minval(x)
    xmax=maxval(x)
    write (ERROR_UNIT,*) "minimun and maximum of x: ",xmin,xmax
    if ( present(xlow) ) then
!        if ( xmin<xlow ) xmin=xlow
        xmin=xlow
    end if
    if ( present(xhigh) ) then
!        if ( xmax>xhigh ) xmax=xhigh
        xmax=xhigh
    end if
    xdel=(xmax-xmin)/bin
    write (ERROR_UNIT,*) "x-range: ",xmin,xmax
    do j=1,bin
        xbin(j) = xmin+j*xdel
    end do

    num_bin=0
    ybin=0
    ybin_err=0
    do i=1,xydim
        do j=1,bin
            if ( x(i)<=xbin(j) .and. x(i)>=xmin ) then
                num_bin(j) = num_bin(j)+1
                ybin(j) = ybin(j)+y(i)
                ybin_err(j) = ybin_err(j)+y(i)*y(i)
                exit
            end if
        end do
    end do
    where ( num_bin/=0 )
!        ybin_err = sqrt( (ybin_err-ybin(:)**2/num_bin) )/num_bin
        ybin_err = sqrt( (ybin_err-ybin(:)**2/num_bin) /num_bin )
    end where
    where ( num_bin/=0 )
        ybin = ybin/num_bin
    end where
    do j=1,bin
        if ( num_bin(j)>=cr_num ) then
            write (out_unit,*) xbin(j),num_bin(j),ybin(j),ybin_err(j)
        end if
    end do

    end subroutine pl2d

    subroutine pl3d(x,y,z,bin1,bin2,xlow,xhigh,ylow,yhigh,cr,file_unit)
    implicit none
    real(kind=PC),dimension(:),intent(in)::x,y,z
    integer,intent(in)::bin1,bin2
! x-range, optional. if not set, use the range: [min(x),max(x)]
    real(kind=PC),optional,intent(in)::xlow,xhigh
! y-range, optional. if not set, use the range: [min(y),max(y)]
    real(kind=PC),optional,intent(in)::ylow,yhigh
! critical number, optional. only the bin with number above will be chosen.
! if not set, use 0
    integer,optional,intent(in)::cr
    integer::cr_num
! output file unit. if not set, use the standard output unit OUTPUT_UNIT
    integer,optional,intent(in)::file_unit
    integer::out_unit
    real(kind=PC),dimension(bin1)::xbin
    real(kind=PC),dimension(bin2)::ybin
    real(kind=PC),dimension(bin1,bin2)::zbin,zbin_err
    integer,dimension(bin1,bin2)::num_bin
    real::xmin,xmax,xdel,ymin,ymax,ydel
    integer::i,jx,jy
    integer::xdim,ydim,zdim
    cr_num=0
    if ( present(cr) ) cr_num=cr
    out_unit=OUTPUT_UNIT
    if ( present(file_unit) ) out_unit=file_unit
    xdim=ubound(x,1)-lbound(x,1)+1
    ydim=ubound(y,1)-lbound(y,1)+1
    zdim=ubound(z,1)-lbound(z,1)+1
    if ( xdim/=zdim .OR. ydim/=zdim ) then
        write (ERROR_UNIT,*) "x, y, z dimension unmacth!"
        stop
    end if

    xmin=minval(x)
    xmax=maxval(x)
    write (ERROR_UNIT,*) "minimun and maximum of x: ",xmin,xmax
    if ( present(xlow) ) then
!        if ( xmin<xlow ) xmin=xlow
        xmin=xlow
    end if
    if ( present(xhigh) ) then
!        if ( xmax>xhigh ) xmax=xhigh
        xmax=xhigh
    end if
    xdel=(xmax-xmin)/bin1
    write (ERROR_UNIT,*) "x-range: ",xmin,xmax
    do jx=1,bin1
        xbin(jx) = xmin+jx*xdel
    end do

    ymin=minval(y)
    ymax=maxval(y)
    write (ERROR_UNIT,*) "minimun and maximum of y: ",ymin,ymax
    if ( present(ylow) ) then
!        if ( ymin<ylow ) ymin=ylow
        ymin=ylow
    end if
    if ( present(yhigh) ) then
!        if ( ymax>yhigh ) ymax=yhigh
        ymax=yhigh
    end if
    ydel=(ymax-ymin)/bin2
    write (ERROR_UNIT,*) "y-range: ",ymin,ymax
    do jy=1,bin2
        ybin(jy) = ymin+jy*ydel
    end do

    num_bin=0
    zbin=0
    zbin_err=0
    do i=1,zdim
        x_mesh:do jx=1,bin1
            if ( x(i)<=xbin(jx) .and. x(i)>=xmin ) then
            y_mesh:do jy=1,bin2
                if ( y(i)<=ybin(jy) .and. y(i)>=ymin ) then
                    num_bin(jx,jy) = num_bin(jx,jy)+1
                    zbin(jx,jy) = zbin(jx,jy)+z(i)
                    zbin_err(jx,jy) = zbin_err(jx,jy)+z(i)*z(i)
                    exit
                end if
            end do y_mesh
            exit
            end if
        end do x_mesh
    end do

    where ( num_bin/=0 )
!        zbin_err = sqrt( (zbin_err-zbin*zbin/num_bin) )/num_bin
        zbin_err = sqrt( (zbin_err-zbin*zbin/num_bin) /num_bin )
    end where
    where ( num_bin/=0 )
        zbin = zbin/num_bin
    end where
    do jy=1,bin2
        do jx=1,bin1
            if ( num_bin(jx,jy)>=cr_num ) then
                write (out_unit,*) &
                &xbin(jx),ybin(jy),num_bin(jx,jy),&
                &zbin(jx,jy),zbin_err(jx,jy)
            end if
        end do
        write (out_unit,*)
    end do

    end subroutine pl3d

    subroutine histbin(array,bin,alow,ahigh,file_unit)
    implicit none
    real(kind=PC),dimension(:),intent(in)::array
    integer::n
    integer,intent(in)::bin
    real(kind=PC),optional,intent(in)::alow,ahigh
    integer,optional,intent(in)::file_unit
    integer::out_unit
    integer,dimension(bin)::num_bin
    real(kind=PC),dimension(bin)::arr_bin
    real(kind=PC)::amax,amin,adel
    integer::i,j

    n=ubound(array,1)-lbound(array,1)+1
    out_unit=OUTPUT_UNIT
    if ( present(file_unit) ) out_unit=file_unit

    amin=minval(array)
    amax=maxval(array)
    write (ERROR_UNIT,*) "array min and max: ",amin,amax
    if ( present(alow) ) then
!        if ( amin<alow ) amin=alow
        amin=alow
    end if
    if ( present(ahigh) ) then
!        if ( amax>ahigh ) amax=ahigh
        amax=ahigh
    end if
    adel=(amax-amin)/bin
    write (ERROR_UNIT,*) "hist range: ",amin,amax
    do j=1,bin
        arr_bin(j) = amin+j*adel
    end do

    num_bin=0
    do i=1,n
        do j=1,bin
            if ( array(i)<=arr_bin(j) .and. array(i)>=amin ) then
                num_bin(j) = num_bin(j)+1
                exit
            end if
        end do
    end do
!    write (ERROR_UNIT,'(2I8)') "plot ",sum(num_bin)," from ",n," points"
    write (ERROR_UNIT,*) sum(num_bin),n
    do j=1,bin
        write (out_unit,*) arr_bin(j),num_bin(j)
    end do
    end subroutine histbin

end Module for_plot
