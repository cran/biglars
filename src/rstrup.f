        subroutine rstrup( nrowr, ncolr, r, lr, loc, v, w, y)
c loc are the locations of the columns in the original R (in proper order)
c R is reduced to the columns in loc (ncolr is length of loc)

        implicit none
        integer            nrowr, ncolr, lr, loc(ncolr)

        double precision   r(lr,ncolr), y(nrowr)
        double precision   v(*), w(*)

        integer            j, lj, nr, nc

        double precision   v1, temp, vnorm

        double precision   zero, one, two
        parameter         (zero = 0.d0, one = 1.d0, two = 2.d0)

        double precision   dnrm2, ddot
        external           dnrm2, ddot

          nc = ncolr

          do j = 1, ncolr

            lj = loc(j)

            if (lj .eq. j) goto 111

            nr = lj - j + 1

            call dcopy( nr, r(j,j), 1, v, 1)
 
            vnorm = dnrm2(nr, v, 1)
 
            if (vnorm .eq. zero) goto 111

            v1    = r(j,j) + sign(one,r(j,j))*vnorm

            v(1)  = v1

            vnorm = dnrm2(nr, v, 1)

            call dscal( nr, one/vnorm, v, 1)

            call dgemv( 'T', nr, nc, one, r(j,j), lr,
     *                   v, 1, zero, w, 1)
 
            call dger( nr, nc, (-two), v, 1, w, 1, 
     *                 r(j,j), lr)

            temp  = - two * ddot( nr, v, 1, y(j), 1)
            call daxpy( nr, temp, v, 1, y(j), 1)

111         continue
         
            nc = nc - 1

          end do

        return
        end

