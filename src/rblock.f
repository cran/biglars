        subroutine rblock( nr, r, lr, nrowx, ncolx, x, lx, v, w,
     *                     yr, yx)
c number of cols of X used should either equal or be one less than
c that of R (nr); y should be the appropriate block of y

        implicit none

        integer            nr, lr, nrowx, ncolx, lx

c	double precision   r(lr,nr), x(lx,nr), v(nrowx), w(nr)
        double precision   r(lr,*), x(lx,*)
        double precision   v(*), w(*)
        double precision   yr(*), yx(*)

        integer            i, j, ncols

        double precision   v1, temp, vmag, vnorm, rtemp

        double precision   zero, one, two
        parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

        double precision   dnrm2, ddot
        external           dnrm2, ddot

        if (ncolx .eq. nr) then
c do not add an intercept

          ncols = nr

          do j = 1, nr

            call dcopy( nrowx, x(1,j), 1, v, 1)

            vnorm = dnrm2( nrowx, v, 1)

            v1    = r(j,j)

c           temp = sqrt(v1*v1 + vnorm*vnorm)

            vmag = max(abs(v1),vnorm)

            if (vmag .eq. zero) goto 111

            if (vmag .le. one) then
              temp = sqrt(v1*v1 + vnorm*vnorm)
            else
              temp = sqrt((v1/vmag)**2 + (vnorm/vmag)**2)*vmag
            end if
 
            v1    = r(j,j) + sign(one,r(j,j))*temp

            vmag = max(abs(v1),vnorm)

            if (vmag .le. one) then
              temp = sqrt(v1*v1 + vnorm*vnorm)
            else
              temp = sqrt((v1/vmag)**2 + (vnorm/vmag)**2)*vmag
            end if

c use the unit vector for numerical stability
            v1 = v1/temp
            call dscal( nrowx, one/temp, v, 1)

            call dgemv( 'T', nrowx, ncols, one, x(1,j), lx,
     *                   v, 1, zero, w, 1)

            call daxpy( ncols, v1, r(j,j), lr, w, 1)

            call daxpy( ncols, (-two*v1), w, 1, r(j,j), lr)

            call dger( nrowx, ncols, (-two), v, 1, w, 1, 
     *                 x(1,j), lx)

            temp  = -two*(v1*yr(j)+ddot(nrowx,v,1,yx,1))

            yr(j) = yr(j) + temp*v1

            call daxpy( nrowx, temp, v, 1, yx, 1)

111         continue

            ncols = ncols - 1

          end do
   
        else if (ncolx+1 .eq. nr) then
c add an intercept

          call dcopy( nrowx, one, 0, v, 1)

          vnorm = sqrt(dble(nrowx))

          v1    = r(1,1)

c         temp = sqrt(v1*v1 + vnorm*vnorm)

          vmag = max(abs(v1),vnorm)

          if (vmag .le. one) then
            temp = sqrt(v1*v1 + vnorm*vnorm)
          else
            temp = sqrt((v1/vmag)**2 + (vnorm/vmag)**2)*vmag
          end if
 
          v1    = r(1,1) + sign(one,r(1,1))*temp

          vmag = max(abs(v1),vnorm)

          if (vmag .le. one) then
            temp = sqrt(v1*v1 + vnorm*vnorm)
          else
            temp = sqrt((v1/vmag)**2 + (vnorm/vmag)**2)*vmag
          end if
      
          v1    = v1/temp
          call dscal( nrowx, one/temp, v, 1)

          call dgemv( 'T', nrowx, ncolx, one, x(1,1), lx,
     *                 v, 1, zero, w(2), 1)

          w(1) = ddot( nrowx, one, 0, v, 1)

          call daxpy( nr, v1, r(1,1), lr, w, 1)

          call daxpy( nr, (-two*v1), w, 1, r(1,1), lr)

          call dger( nrowx, ncolx, (-two), v, 1, w(2), 1, 
     *               x(1,1), lx)

          temp  = -two*(v1*yr(1)+ddot(nrowx,v,1,yx,1))

          yr(1) = yr(1) + temp*v1

          call daxpy( nrowx, temp, v, 1, yx, 1)

          i     = 1

          ncols = ncolx

          do j = 2, nr

            call dcopy( nrowx, x(1,i), 1, v, 1)

            vnorm = dnrm2( nrowx, v, 1)

            v1    = r(j,j)

            vmag = max(abs(v1),vnorm)

C Zero-length vector:  nothing left to do.

            if (vmag .eq. zero) goto 222

            if (vmag .le. one) then
              temp = sqrt(v1*v1 + vnorm*vnorm)
            else
              temp = sqrt((v1/vmag)**2 + (vnorm/vmag)**2)*vmag
            end if
 
            v1    = r(j,j) + sign(one,r(j,j))*temp

            vmag = max(abs(v1),vnorm)

            if (vmag .le. one) then
              temp = sqrt(v1*v1 + vnorm*vnorm)
            else
              temp = sqrt((v1/vmag)**2 + (vnorm/vmag)**2)*vmag
            end if

c Zero-length vector:  nothing left to do.

            if (temp .eq. zero) goto 222

            rtemp = one / temp
 
c use the unit vector for numerical stability
            v1 = v1 * rtemp
            call dscal( nrowx, rtemp, v, 1)

            call dgemv( 'T', nrowx, ncols, one, x(1,i), lx,
     *                   v, 1, zero, w, 1)
 
            call daxpy( ncols, v1, r(j,j), lr, w, 1)

            call daxpy( ncols, (-two*v1), w, 1, r(j,j), lr)

            call dger( nrowx, ncols, (-two), v, 1, w, 1, 
     *                 x(1,i), lx)

            temp  = -two*(v1*yr(j)+ddot(nrowx,v,1,yx,1))

            yr(j) = yr(j) + temp*v1

            call daxpy( nrowx, temp, v, 1, yx, 1)

222         continue

            ncols = ncols - 1

            i     = j

          end do

        else
  
           call intpr("ncol(X) incompatible with ncol(R)", 33, j, 0)
 
        end if

        return
        end
