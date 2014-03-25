SUBROUTINE crwdrift_simulate( tau2y, tau2x, Qmat, Tmat, x, y, loctype, &
                             stay, a1y, a1x, P1y, P1x, lonadj, N,  lly, llx, &
                             alphaY, alphaX, ysim, xsim, alphaSimY, alphaSimX)
 !! LOADED VALUES/ARRAYS !!                        
  INTEGER           N
  INTEGER           loctype(N), stay(N)
  DOUBLE PRECISION  lly, llx, tau2y(N), tau2x(N), Qmat(N,5), Tmat(N, 4)
  DOUBLE PRECISION  x(N), y(N), lonadj(N), ysim(N), xsim(N)
  DOUBLE PRECISION  a1x(3), a1y(3), P1y(3,3), P1x(3,3)
  DOUBLE PRECISION  alphaY(3,1,N+1), alphaX(3,1,N+1)
  DOUBLE PRECISION  alphaSimY(N,3), alphaSimX(N,3)
                    
 !! DERIVED ARRAYS !!
  DOUBLE PRECISION vy(N), vx(N), Fy(N), Fx(N), vySim(N), vxSim(N)
  DOUBLE PRECISION Z(1,3)
  DOUBLE PRECISION Ky(3,1), Kx(3,1) 
  DOUBLE PRECISION Qy(3,3), Qx(3,3), T(3,3) 
  DOUBLE PRECISION Py(3,3,N+1), Px(3,3,N+1), Ly(3,3,N), Lx(3,3,N)
  DOUBLE PRECISION ay(3,1,N+1), ax(3,1,N+1), aySim(3,1,N+1), axSim(3,1,N+1)
  DOUBLE PRECISION ayHat(3), axHat(3), ayHatSim(3), axHatSim(3)
  DOUBLE PRECISION alphaMeanY(3,1), alphaMeanX(3,1)
  DOUBLE PRECISION ry(3,1), rySim(3,1), rx(3,1), rxSim(3,1)
  DOUBLE PRECISION c12y, c13y, c12x, c13x
  
 !! INTIAL CONDITIONS !! 
  ay(:,1,1) = a1y
  ax(:,1,1) = a1x
  aySim(:,1,1) = a1y
  axSim(:,1,1) = a1x
  Py(:,:,1) = P1y
  Px(:,:,1) = P1x
  alphaY(3,1,1) = a1y(3) + sqrt(P1y(3,3))*alphaY(3,1,1)
  alphaY(2,1,1) = a1y(2) + sqrt(P1y(2,2))*alphaY(2,1,1)
  IF(P1y(2,2)==0.0) THEN
   c12y = 0.0
  ELSE 
   c12y = P1y(1,2)/P1y(2,2)
  END IF
  IF(P1y(3,3) == 0.0) THEN
   c13y=0.0
  ELSE 
   c13y = P1y(1,3)/P1y(3,3)
  END IF  
  alphaY(1,1,1) = a1y(1) + c12y*(alphaY(2,1,1)-a1y(2)) &
                         + c13y*(alphaY(3,1,1)-a1y(3)) &
                         + sqrt(P1y(1,1) - P1y(1,2)*c12y - P1y(1,3)*c13y)*alphaY(1,1,1)
                         
  alphaX(3,1,1) = a1x(3) + sqrt(P1x(3,3))*alphaX(3,1,1)
  alphaX(2,1,1) = a1x(2) + sqrt(P1x(2,2))*alphaX(2,1,1)
  IF(P1x(2,2)==0.0) THEN
   c12x = 0.0
  ELSE 
   c12x = P1x(1,2)/P1x(2,2)
  END IF
  IF(P1x(3,3) == 0.0) THEN
   c13x=0.0
  ELSE 
   c13x = P1x(1,3)/P1x(3,3)
  END IF 
  alphaX(1,1,1) = a1x(1) + (P1x(1,2)/P1x(2,2))*(alphaX(2,1,1)-a1x(2)) &
                         + (P1x(1,3)/P1x(3,3))*(alphaX(3,1,1)-a1x(3)) &
                         + sqrt(P1x(1,1) - P1x(1,2)*c12x - P1x(1,3)*c13x)*alphaX(1,1,1)                   
  Z            = RESHAPE((/1.0, 0.0, 0.0/),(/1, 3/))
  Qy           = 0.0
  Qx           = 0.0
  T            = 0.0
  T(1,1)       = 1.0
  vy           = 0.0
  vx           = 0.0
  vySim        = 0.0
  vxSim        = 0.0
  Fy           = tau2y
  Fx           = tau2x
  ry           = 0.0
  rx           = 0.0
  rySim        = 0.0
  rxSim        = 0.0
    
 !! BEGIN FILTER LOOP !! 
  DO i=1,N
  
   !! GENERATE Q AND T MATRICES !!
    IF(stay(i)==1) THEN
      Qy = 0.0
      Qx = 0.0
      T(1,2) = 0.0
      T(2,2) = 0.0
     !! SIMULATE NEXT STATE VALUE !!
      alphaY(2,1,i+1) = 0.0
      alphaY(3,1,i+1) = 0.0
      alphaY(1,1,i+1) = alphaY(1,1,i)
      alphaX(2,1,i+1) = 0.0
      alphaX(3,1,i+1) = 0.0
      alphaX(1,1,i+1) = alphaX(1,1,i)

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
     !! SIMULATE NEXT STATE VALUE !!
      alphaMeanY = MATMUL(T,alphaY(:,:,i))
      alphaY(3,1,i+1) = alphaMeanY(3,1) + sqrt(Qy(3,3))*alphaY(3,1,i+1)
      alphaY(2,1,i+1) = alphaMeanY(2,1) + sqrt(Qy(2,2))*alphaY(2,1,i+1)
      alphaY(1,1,i+1) = alphaMeanY(1,1) + (Qy(1,2)/Qy(2,2))*(alphaY(2,1,i+1)-alphaMeanY(2,1)) &
                              + (Qy(1,3)/Qy(3,3))*(alphaY(3,1,i+1)-alphaMeanY(3,1)) &
                              + sqrt(Qy(1,1) - ((Qy(1,2)*Qy(1,2))/Qy(2,2)) &
                                      - ((Qy(1,3)*Qy(1,3))/Qy(3,3)))*alphaY(1,1,i+1)
      alphaMeanX = MATMUL(T,alphaX(:,:,i))
      alphaX(3,1,i+1) = alphaMeanX(3,1) + sqrt(Qx(3,3))*alphaX(3,1,i+1)
      alphaX(2,1,i+1) = alphaMeanX(2,1) + sqrt(Qx(2,2))*alphaX(2,1,i+1)
      alphaX(1,1,i+1) = alphaMeanX(1,1) + (Qx(1,2)/Qx(2,2))*(alphaX(2,1,i+1)-alphaMeanX(2,1)) &
                             + (Qx(1,3)/Qx(3,3))*(alphaX(3,1,i+1)-alphaMeanX(3,1)) &
                             + sqrt(Qx(1,1) - ((Qx(1,2)*Qx(1,2))/Qx(2,2)) &
                                     - ((Qx(1,3)*Qx(1,3))/Qx(3,3)))*alphaX(1,1,i+1) 
    END IF
        
   !! GENERAL KF FILTER !!
    Fy(i) = Py(1,1,i) + tau2y(i)
    Fx(i) = Px(1,1,i) + tau2x(i)/(lonadj(i)*lonadj(i))
    IF(loctype(i)==1 .OR. Fy(i)==0.0) THEN
      ay(:,:,i+1) = MATMUL(T,ay(:,:,i))
      aySim(:,:,i+1) = MATMUL(T, aySim(:,:,i))      
      Py(:,:,i+1) = MATMUL(MATMUL(T,Py(:,:,i)),TRANSPOSE(T)) + Qy      
      Ly(:,:,i) = T
      
    ELSE
      ysim(i) = alphaY(1,1,i) + sqrt(tau2y(i))*ysim(i)      
      vy(i) = y(i) - ay(1,1,i)      
      vySim(i) = ysim(i) - aySim(1,1,i)      
      Ky = MATMUL(MATMUL(T,Py(:,:,i)),TRANSPOSE(Z))/Fy(i)      
      Ly(:,:,i) = T - MATMUL(Ky,Z)      
      ay(:,:,i+1) = MATMUL(T,ay(:,:,i)) + Ky*vy(i)      
      aySim(:,:,i+1) = MATMUL(T,aySim(:,:,i)) + Ky*vySim(i)            
      Py(:,:,i+1) = MATMUL(MATMUL(T,Py(:,:,i)),TRANSPOSE(Ly(:,:,i))) + Qy      
      lly = lly - (log(Fy(i)) + vy(i)*vy(i)/Fy(i))/2      
    END IF
    IF(loctype(i)==1 .OR. Fx(i)==0.0) THEN
     ax(:,:,i+1) = MATMUL(T,ax(:,:,i))
     axSim(:,:,i+1) = MATMUL(T, axSim(:,:,i))
     Px(:,:,i+1) = MATMUL(MATMUL(T,Px(:,:,i)),TRANSPOSE(T)) + Qx
     Lx(:,:,i) = T
    ELSE
     xsim(i) = alphaX(1,1,i) + (sqrt(tau2x(i))/lonadj(i))*xsim(i)
     vx(i) = x(i) - ax(1,1,i)
     vxSim(i) = xsim(i) - axSim(1,1,i)
     Kx = MATMUL(MATMUL(T,Px(:,:,i)),TRANSPOSE(Z))/Fx(i)
     Lx(:,:,i) = T - MATMUL(Kx,Z)
     ax(:,:,i+1) = MATMUL(T,ax(:,:,i)) + Kx*vx(i)
     axSim(:,:,i+1) = MATMUL(T,axSim(:,:,i)) + Kx*vxSim(i)
     Px(:,:,i+1) = MATMUL(MATMUL(T,Px(:,:,i)),TRANSPOSE(Lx(:,:,i))) + Qx
     llx = llx - (log(Fx(i)) + vx(i)*vx(i)/Fx(i))/2
    END IF
  END DO

 !! BEGIN SMOOTHING LOOP!!  
  DO j=N,1,-1
    IF(loctype(j)==1 .OR. Fy(j)==0.0) THEN 
     ry = MATMUL(TRANSPOSE(Ly(:,:,j)),ry)      
     rySim = MATMUL(TRANSPOSE(Ly(:,:,j)),rySim)        
    ELSE
     ry = TRANSPOSE(Z)*vy(j)/Fy(j) + MATMUL(TRANSPOSE(Ly(:,:,j)),ry)          
     rySim = TRANSPOSE(Z)*vySim(j)/Fy(j) + MATMUL(TRANSPOSE(Ly(:,:,j)),rySim)      
    END IF
    IF(loctype(j)==1 .OR. Fx(j)==0.0) THEN
     rx = MATMUL(TRANSPOSE(Lx(:,:,j)),rx)
     rxSim = MATMUL(TRANSPOSE(Lx(:,:,j)),rxSim)
    ELSE
     rx = TRANSPOSE(Z)*vx(j)/Fx(j) + MATMUL(TRANSPOSE(Lx(:,:,j)),rx)
     rxSim = TRANSPOSE(Z)*vxSim(j)/Fx(j) + MATMUL(TRANSPOSE(Lx(:,:,j)),rxSim)
    END IF
    
    ayHat = RESHAPE(ay(:,:,j) + MATMUL(Py(:,:,j), ry), (/3/))
    axHat = RESHAPE(ax(:,:,j) + MATMUL(Px(:,:,j), rx), (/3/))
    ayHatSim = RESHAPE(aySim(:,:,j) + MATMUL(Py(:,:,j), rySim), (/3/))
    axHatSim = RESHAPE(axSim(:,:,j) + MATMUL(Px(:,:,j), rxSim), (/3/))
    alphaSimY(j,:) = ayHat + (RESHAPE(alphaY(:,:,j),(/3/)) - ayHatSim)
    alphaSimX(j,:) = axHat + (RESHAPE(alphaX(:,:,j),(/3/)) - axHatSim)    
  END DO  
END SUBROUTINE crwdrift_simulate

