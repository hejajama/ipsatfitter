      PROGRAM testEvolFortran
      INTEGER i, j, order,  coupling
      DOUBLE PRECISION x, Q2, gluon, mc, mb, mu0, asmur_, Ag, lambdag, As, lambdas

      order = 0
      call INITALPHAS(order,1.D0,1.D0,0.5D0,1.4D0,4.75D0,1.D10) ! initializes alpha_s
                     !  0   FR2  MUR  ASMUR,MC,   MB,    MT
      call INIT  ! initializes evolution

      mc = 1.4D0
      mb = 4.75D0
      mu0 = 1.D0
      asmur = 0.5D0
      Ag =  2.29831
      lambdag = 0.09708
      As = 0
      lambdas = 0
      coupling = 0

      do i=1,1000,10
         do j=1, 1000, 10
            x = i/2000.
            Q2 = j
            call LO_evol(x, Q2, gluon, coupling, mc, mb, mu0, asmur,
     &              Ag, lambdag, As, lambdas)
            print *, 'x=', x, ', Q2=', Q2, ', value=', gluon
         enddo
      enddo
      end

