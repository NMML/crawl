SUBROUTINE crw_simulate(tau2y, tau2x, Qmat, Tmat, x, y, loctype, &
                             stay, a1y, a1x, P1y, P1x, lonadj, N,  lly, llx, &
                             alphaY, alphaX, ysim, xsim, alphaSimY, alphaSimX)
 !! LOADED VALUES/ARRAYS !!                        
  INTEGER           N
  INTEGER           loctype(N), stay(N)
  DOUBLE PRECISION  lly, llx, tau2y(N),tau2x(N)
  DOUBLE PRECISION  x(N), y(N), lonadj(N), ysim(N), xsim(N)
  DOUBLE PRECISION  a1x(2), a1y(2), P1y(2,2), P1x(2,2)
  DOUBLE PRECISION  alphaY(2,1,N+1), alphaX(2,1,N+1)
  DOUBLE PRECISION  alphaSimY(N,2), alphaSimX(N,2), Qmat(N,3), Tmat(N,2)
                    
 !! DERIVED ARRAYS !!
  DOUBLE PRECISION vy(N), vx(N), Fy(N), Fx(N), vySim(N), vxSim(N)
  DOUBLE PRECISION Z(1,2)
  DOUBLE PRECISION Ky(2,1), Kx(2,1) 
  DOUBLE PRECISION Qy(2,2), Qx(2,2), T(2,2) 
  DOUBLE PRECISION Py(2,2,N+1), Px(2,2,N+1), Ly(2,2,N), Lx(2,2,N)
  DOUBLE PRECISION ay(2,1,N+1), ax(2,1,N+1), aySim(2,1,N+1), axSim(2,1,N+1)
  DOUBLE PRECISION ayHat(2), axHat(2), ayHatSim(2), axHatSim(2)
  DOUBLE PRECISION alphaMeanY(2,1), alphaMeanX(2,1)
  DOUBLE PRECISION ry(2,1), rySim(2,1), rx(2,1), rxSim(2,1)

  
 !! INTIAL CONDITIONS !! 
  ay(:,1,1) = a1y
  ax(:,1,1) = a1x
  aySim(:,1,1) = a1y
  axSim(:,1,1) = a1x
  Py(:,:,1) = P1y
  Px(:,:,1) = P1x
  
  alphaY(2,1,1) = a1y(2) + sqrt(P1y(2,2))*alphaY(2,1,1)
  alphaX(2,1,1) = a1x(2) + sqrt(P1x(2,2))*alphaX(2,1,1)
  IF(P1y(2,2)==0.0) THEN
  	alphaY(1,1,1) = a1y(1) + sqrt(P1y(1,1))*alphaY(1,1,1)
  ELSE
  	alphaY(1,1,1) = a1y(1) + (P1y(1,2)/P1y(2,2))*(alphaY(2,1,1)-a1y(2)) &
                           + sqrt(P1y(1,1) - ((P1y(1,2)*P1y(1,2))/P1y(2,2)))*alphaY(1,1,1)
  END IF
  IF(P1x(2,2)==0.0) THEN
  	alphaX(1,1,1) = a1x(1) + sqrt(P1x(1,1))*alphaX(1,1,1)
  
  ELSE
  	alphaX(1,1,1) = a1x(1) + (P1x(1,2)/P1x(2,2))*(alphaX(2,1,1)-a1x(2)) &
                           + sqrt(P1x(1,1) - P1x(1,2)*P1x(1,2)/P1x(2,2))*alphaX(1,1,1)
  END IF
  Z            = RESHAPE((/1.0, 0.0/),(/1, 2/))
!  Qy           = 0.0
!  Qx           = 0.0
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
      alphaY(1,1,i+1) = alphaY(1,1,i)
      alphaX(2,1,i+1) = 0.0
      alphaX(1,1,i+1) = alphaX(1,1,i)
    ELSE
      Qy(1,1) = Qmat(i,1)
      Qy(2,1) = Qmat(i,2)
      Qy(1,2) = Qy(2,1)
      Qy(2,2) = Qmat(i,3)
      Qx = Qy/(lonadj(i)*lonadj(i))      
      T(1,2) = Tmat(i,1)
      T(2,2) = Tmat(i,2)
      
     !! SIMULATE NEXT STATE VALUE !!
      alphaMeanY = MATMUL(T,alphaY(:,:,i))
      IF (Qy(2,2) == 0) THEN	
      	alphaY(2,1,i+1) = alphaMeanY(2,1)
      	alphaY(1,1,i+1) = alphaMeanY(1,1) + sqrt(Qy(1,1))*alphaY(1,1,i+1)
      ELSE 
      	alphaY(2,1,i+1) = alphaMeanY(2,1) + sqrt(Qy(2,2))*alphaY(2,1,i+1)
      	alphaY(1,1,i+1) = alphaMeanY(1,1) + (Qy(1,2)/Qy(2,2))*(alphaY(2,1,i+1)-alphaMeanY(2,1)) &
                        	+ sqrt(Qy(1,1)-((Qy(1,2)*Qy(1,2))/Qy(2,2)))*alphaY(1,1,i+1)
      END IF
      alphaMeanX = MATMUL(T,alphaX(:,:,i))
      IF (Qx(2,2) == 0) THEN	
      	alphaX(2,1,i+1) = alphaMeanX(2,1)
      	alphaX(1,1,i+1) = alphaMeanX(1,1) + sqrt(Qx(1,1))*alphaX(1,1,i+1)
      ELSE 
      	alphaX(2,1,i+1) = alphaMeanX(2,1) + sqrt(Qx(2,2))*alphaX(2,1,i+1)
      	alphaX(1,1,i+1) = alphaMeanX(1,1) + (Qx(1,2)/Qx(2,2))*(alphaX(2,1,i+1)-alphaMeanX(2,1)) &
                        	+ sqrt(Qx(1,1)-((Qx(1,2)*Qx(1,2))/Qx(2,2)))*alphaX(1,1,i+1)
      END IF
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
    
   !! STORE VALUES !!
    ayHat = RESHAPE(ay(:,:,j) + MATMUL(Py(:,:,j), ry), (/2/))
    axHat = RESHAPE(ax(:,:,j) + MATMUL(Px(:,:,j), rx), (/2/))
    ayHatSim = RESHAPE(aySim(:,:,j) + MATMUL(Py(:,:,j), rySim), (/2/))
    axHatSim = RESHAPE(axSim(:,:,j) + MATMUL(Px(:,:,j), rxSim), (/2/))
    alphaSimY(j,:) = ayHat + (RESHAPE(alphaY(:,:,j),(/2/)) - ayHatSim)
    alphaSimX(j,:) = axHat + (RESHAPE(alphaX(:,:,j),(/2/)) - axHatSim)    
  END DO  
END SUBROUTINE crw_simulate

