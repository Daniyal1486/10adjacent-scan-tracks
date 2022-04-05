      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1 TEMP,PRESS,SNAME)
C
      INCLUDE 'ABA_PARAM.INC'
C
      DIMENSION COORDS(3),FLUX(2),TIME(2)
      CHARACTER*80 SNAME
      REAL thermcon,alpha,vel,R0,R,pi,dist
     & x,y,z,power,lamda,abs,AI,M,Eff,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,u
      power=200
      vel=0.04
C ********************
      Eff=2. ! 11 this added t othis code
      pi=3.141592654
      abs=0.52
      R0=0.001
      u=0.006
      t1=u/vel
      t2=t1+u/vel
      t3=t2+u/vel
      t4=t3+u/vel
      t5=t4+u/vel
      t6=t5+u/vel
      t7=t6+u/vel
      t8=t7+u/vel
      t9=t8+u/vel
      t10=t9+u/vel
      AI=(Eff*power)/(pi*(R0*R0)) ! 11
      IF (Time(1) .GT. 0 .AND. TIME(1) .LE. t1) THEN          ! first pass
      dist=vel*Time(1)
      x=COORDS(1)-dist
      y=COORDS(3)
      z=COORDS(2)-0.021
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t1 .AND. TIME(1) .LE. t2) THEN    ! second pass
      dist=vel*t1+vel*(Time(1)-t1)
      x=COORDS(1)+0.006-dist
      y=COORDS(3)
      z=COORDS(2)-0.020
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t2 .AND. TIME(1) .LE. t3) THEN   ! third pass
      dist=vel*t1+vel*(t2-t1)+vel*(Time(1)-t2)
      x=COORDS(1)+0.012-dist
      y=COORDS(3)
      z=COORDS(2)-0.019
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t3 .AND. TIME(1) .LE. t4) THEN   ! fourth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(Time(1)-t3)
      x=COORDS(1)+0.018-dist
      y=COORDS(3)
      z=COORDS(2)-0.018
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t4 .AND. TIME(1) .LE. t5) THEN   ! fifth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(Time(1)-t4)
      x=COORDS(1)+0.024-dist
      y=COORDS(3)
      z=COORDS(2)-0.017
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t5 .AND. TIME(1) .LE. t6) THEN   ! sixth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(Time(1)-t5)
      x=COORDS(1)+0.030-dist
      y=COORDS(3)
      z=COORDS(2)-0.016
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t6 .AND. TIME(1) .LE. t7) THEN   ! seventh pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(Time(1)-t6)
      x=COORDS(1)+0.036-dist
      y=COORDS(3)
      z=COORDS(2)-0.015
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t7 .AND. TIME(1) .LE. t8) THEN   ! eighth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(Time(1)-t7)
      x=COORDS(1)+0.042-dist
      y=COORDS(3)
      z=COORDS(2)-0.014
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t8 .AND. TIME(1) .LE. t9) THEN   ! ninth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(t8-t7)+vel*(Time(1)-t8)
      x=COORDS(1)+0.048-dist
      y=COORDS(3)
      z=COORDS(2)-0.013
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      ELSE IF (Time(1) .GT. t9 .AND. TIME(1) .LE. t10) THEN   ! tenth pass
      dist=vel*t1+vel*(t2-t1)+vel*(t3-t2)+vel*(t4-t3)+vel*(t5-t4)+vel*(t6-t5)+vel*(t7-t6)+vel*(t8-t7)+vel*(t9-t8)+vel*(Time(1)-t9)
      x=COORDS(1)+0.054-dist
      y=COORDS(3)
      z=COORDS(2)-0.012
      R=sqrt(x*x+z*z)
      FLUX(1)=abs*AI*exp(-3*(R*R)/(R0*R0))
      FLUX(2)=0.
      JLTYP=0. 
      END IF
      IF (Time(1) .GT. t10) THEN
      FLUX(1)=0
      FLUX(2)=0.
      JLTYP=0.
      END IF
      RETURN
      END