SUBROUTINE crwdriftn2ll(tau2y, tau2x, Qmat, Tmat, x, y, loctype, &
                          stay, ay, ax, Py, Px, lonadj, N, lly, llx)
  INTEGER N
  INTEGER loctype(N), stay(N)
  DOUBLE PRECISION lly, llx, tau2y(N), tau2x(N), Qmat(N,5), Tmat(N,4)
  DOUBLE PRECISION x(N), y(N), lonadj(N)
  DOUBLE PRECISION ay(3,1), ax(3,1), Py(3,3), Px(3,3)
  DOUBLE PRECISION vy, vx, Fy, Fx
  DOUBLE PRECISION Z(1,3)
  DOUBLE PRECISION Ky(3,1), Kx(3,1)
  DOUBLE PRECISION Qy(3,3), Qx(3,3), T(3,3), Lx(3,3), Ly(3,3)

  !! INITIAL VALUES !!
  Z = RESHAPE((/1.0, 0.0, 0.0/),(/1, 3/))
  T = 0.0
  T(1,1) = 1.0
  Qy = 0.0
  Qx = 0.0

 !! BEGIN FILTER LOOP !!
  DO i=1,N
   !! GENERATE Q and T MATRICES !!
    IF(stay(i)==1) THEN
      Qy = 0.0
      Qx = 0.0
      T(1,2) = 0.0
      T(2,2) = 0.0
      T(1,3) = 0.0
      T(3,3) = 0.0
    ELSE
      Qy(1,1) = Qmat(i,1)
      Qy(2,1) = Qmat(i,2)
      Qy(3,1) = Qmat(i,3)
      Qy(1,2) = Qy(2,1)
      Qy(2,2) = Qmat(i,4)
      Qy(1,3) = Qy(3,1)
      Qy(3,3) = Qmat(i,5)
      Qx = Qy/(lonadj(i)*lonadj(i))
      T(1,2) = Tmat(i,1)
      T(2,2) = Tmat(i,2)
      T(1,3) = Tmat(i,3)
      T(3,3) = Tmat(i,4)
    END IF

   !! GENERAL KF LOOP !!
    Fy = Py(1,1) + tau2y(i)
    Fx = Px(1,1) + tau2x(i)/(lonadj(i)*lonadj(i))
    IF(loctype(i)==1 .OR. Fy==0.0) THEN
      ay = MATMUL(T,ay)
      Py = MATMUL(MATMUL(T,Py),TRANSPOSE(T)) + Qy
    ELSE
      vy = y(i)-ay(1,1)
      lly = lly - (log(Fy) + vy*vy/Fy)/2
      Ky = MATMUL(MATMUL(T,Py),TRANSPOSE(Z))/Fy
      Ly = T - MATMUL(Ky,Z)
      ay = MATMUL(T,ay) + Ky*vy
      Py = MATMUL(MATMUL(T,Py),TRANSPOSE(Ly)) + Qy

    END IF
    IF(loctype(i)==1 .OR. Fx==0.0) THEN      
      ax = MATMUL(T,ax)
      Px = MATMUL(MATMUL(T,Px),TRANSPOSE(T)) + Qx
    ELSE
      vx = x(i)-ax(1,1)
      llx = llx - (log(Fx) + vx*vx/Fx)/2
      Kx = MATMUL(MATMUL(T,Px),TRANSPOSE(Z))/Fx
      Lx = T - MATMUL(Kx,Z)      
      ax = MATMUL(T,ax) + Kx*vx      
      Px = MATMUL(MATMUL(T,Px),TRANSPOSE(Lx)) + Qx
    END IF
  END DO
END SUBROUTINE crwdriftn2ll
