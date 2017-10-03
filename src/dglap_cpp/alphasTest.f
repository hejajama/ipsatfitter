      PROGRAM alphasTest
      INTEGER i, order
      DOUBLE PRECISION IX, F, ALPHAS
      order = 0
      call INITALPHAS(order,1.D0,1.D0,0.5D0,1.4D0,4.75D0,1.D10)
      print *, 'order=', order
      do i=1,200,2
            IX = i
            F = ALPHAS(IX)
            print *, IX, F
      enddo
      end
