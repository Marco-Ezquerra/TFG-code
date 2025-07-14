program lj  
! use modulo    !!!Importamos el modulo  
 implicit none

 integer i,j,nmax,bins,contar,k,binsradial, pasost,u
 integer npart,count,npartmax,nmezcla,pasos, offont,onb
 integer regular,binsmax,npromedio, pantalla, contador,check
 
 real*8 dargon,diametro,T1,T2,kb,margon,d,rnd
 real*8 sigma,v,suma,v1,v2,v3,dt
 real*8 ltot,r,lcaja,rcutoff,ecut,epsilon
 real*8 tsimulacion,a,b,s,vmax
 real*8 cineticatot,cineticac,cineticaf
 real*8 sumv2, sumvx,sumvy,sumvz,tt,pp
 real*8 energia,dr ,cinetica,potencial
 real*8 sx,sy,sz, t_target,q,zeta,t_targetfinal
 real*8 f,temp,qp,p_target,eta,xmax, ymax,zmax,virial
 real*8 beta, t_final, t_inicial, p0, pf, beta_p,rskin 
 real*8 ticpu, tfcpu, tiempocpu, temperatura
 
 parameter(npartmax=6**3) !numero total de part   !!elegir un numero de particulas que sea riaz cubica de 3
 parameter(dargon=0.84d0) !kg/m^3
 parameter(diametro=1.d0) !! diametro de las partculas de argon !!hay que fijarlo como 1
 parameter(t1=1.d0) !!K
 parameter(kb = 1.d0)
 parameter(margon = 1.d0)
 parameter(rskin=0.5d0)
 parameter(rcutoff=2.5d0) !!ya esta normalizado con respecto al diametro rcutoff=3*diametro
 parameter(epsilon=1.d0)
 parameter(ecut=4*epsilon*((1.d0/rcutoff**12.d0)+(1.d0/rcutoff**6.d0)))
 parameter(npromedio=100)
 parameter(bins=int(sqrt(real(npartmax))))
 parameter(binsmax=1000)


 !!--- VECTORES POSICION DE LAS PARTICULAS ---------------------------------------------------------------------------------
 
 real*8 :: x(1:npartmax) 
 real*8 :: y(1:npartmax)
 real*8 :: z(1:npartmax)
 real*8 :: xm(1:npartmax) 
 real*8 :: ym(1:npartmax)
 real*8 :: zm(1:npartmax)

 real*8 :: fx(1:npartmax)
 real*8 :: fy(1:npartmax)
 real*8 :: fz(1:npartmax)


!!-----VELOCIDIDADES DE LAS PARTICULAS--------------------------------------------------------------------------------------------
 real*8 :: vx(1:npartmax)  
 real*8 :: vy(1:npartmax)   
 real*8 :: vz(1:npartmax)     !!podria dimensionar de 1:nmax, hacer una matriz para ver la evolucion de v !!pero vamos a ahorrarnos esta memoria
                                  
 real*8 :: vmodu(1:npartmax)
 real*8 :: vini(1:npartmax) 
 real*8 :: px(1:bins)  !!vector de pesos de la distribucion
 real*8 :: gr(1:binsmax)
 real*8 :: rcontar(1:npartmax*npartmax)
 integer :: histo(1:1000)
!!---------------------------------------------------------------------------------------------------------------------------
 integer :: vecinos(1:npartmax,1:npartmax)
 integer :: nvecinos(1:npartmax)

 real*8 :: xacumu(1:npartmax)
 real*8 :: yacumu(1:npartmax)
 real*8 :: zacumu(1:npartmax)

 real*8 :: gr_temp(1:binsmax)
 real*8 :: grprom(1:binsmax)
 real*8 :: pxprom(1:bins)
 integer :: histo_temp(1:binsmax)

 
  call CPU_TIME(ticpu) 
  
  npart=6**3
  ltot=(npart*margon/(dargon))**(1.d0/3.d0) 
 print*, ltot ,"L"
  
 lcaja=ltot/(npart)**(1.d0/3.d0) 
 print*, lcaja
 
 call random_seed()

 
  
 !!UNIDADES REDUCIDAS-----------------------------------------------
 sigma=sqrt(kb*T1/margon)  !!velocidad cuadratica media para una componente
 sigma=sigma/sqrt(epsilon*kb/margon)    !!unidades reducidas
 dt=0.0001d0 !!para asegurar la estabilidad de verlet 
 ltot=ltot/diametro
 lcaja=lcaja/diametro
 tt=t1/epsilon
 contador=0

 
 print*, SIGMA, dt,ltot,lcaja
 print*, ltot, "longitud inical"
 !------------------------------------------------------------------------
 temp=t1
 !call initvelocidades(sigma,npart,vmodu,temp)
!  call histograma(vmodu, bins, npart, px, a, b)
!       open(100,file="distmb.txt")
!       do i=1,bins
!        write(100,*) a+(b-a)*i/bins, px(i)/(b-a)
!       enddo

  


 
  call posicioninicial(x,y,z,r,ltot,npart,lcaja,diametro,count) 
  !call inicializar_posiciones(x, y, z, npart, ltot, diametro)
  print*, "fin posicion inicial"
  npart=count

  
  !call boxmuller(npart,vx,vy,vz,sigma)
  ! vmax=abs(maxval(vx))
  ! open(21, file="mb.txt")
  !  do i=1,npart
  !   read(21,*) vini(i)
  !  enddo
  !  close(21)  
   !call maxwellboltzamann() !!esta es la subrutina dl modulo que ejecuta en python
   
   open(12,file="mbpy.txt")
    do i=1,npart
     read(12,*) vini(i)
     
     enddo
   
  
   call histograma(vini, bins, npart, px, a, b)
      open(100,file="distmb.txt")
      do i=1,bins
       write(100,*) a+(b-a)*i/bins, px(i)/(b-a)
      enddo

     do i=1,npart
    
    call direcciones(sx,sy,sz)
     vx(i)=vini(i)*sx
     vy(i)=vini(i)*sy
     vz(i)=vini(i)*sz
     enddo
  open(22,file="posini.txt")
  
  call reajuste(npart, vx,vy,vz)
   sumv2=0.d0
  do i=1,npart

  write(22,*) x(i),y(i), z(i)
 
  enddo
 





  xmax=maxval(x)
  print*, xmax
  
 open(12, file="virial.txt")
 !!!!Ya tenemos velocidades iniciales y posiciones iniciales ademas de tener vcm=0
 nmezcla=0
 open(33,file="energia.txt")
 !!!MUY IMPORTANTE PRIMERO HAY QUE TERMALIZAR EL SISTEMA PARA COMENZAR A TOMAR LAS MEDIDAS UTILIZAREMOS 1000 PASOS DE TERMALIZACION
 !!FASE DE PRECALENTAMIENTO, INCICIALMENTE EL SISTEMA EMPIEZA EN UNA CONFIGURACION CON ENERGIA POTENCIAL ALTA, VAMOS A TENER UN PICO DE ENERGIA CINETICA
 !!!UTILIZAREMOS EN LA FASE DE PRECALENTAMIENTO UN TERMOSTATO PARA LLEGAR A UNAS CONDICIONES INICALES  BUENAS PARA COMENZAR A MEDIR

 !!! Inicialmente el sistema empieza con un pico de energia, y luego ya se estabiliza asi pues la temeratura inicial ha de modicarse
 sumv2=0.d0
 do i=1,npart
 sumv2=sumv2+ vx(i)**2 +vy(i)**2 +vz(i)**2
  enddo
  
 t_target=(sumv2)/(3.d0*npart)   

 print*, t_target, "temperatura inicial"
  q=100
  pasos=400000
  offont=1
  onb=0
  t_target=0.5d0
  call listaverlet(npart, x,y,z, rcutoff,rskin, ltot, vecinos, nvecinos)
  call fuerzas(npart, x,y,z,fx,fy,fz, vecinos, nvecinos, ltot,potencial)
  grprom=0.d0
  do while(nmezcla .lt. pasos)

       if(mod(nmezcla,10) .eq. 0) then 
     call listaverlet(npart, x,y,z, rcutoff,rskin, ltot, vecinos, nvecinos) 
     endif

 call integrarmove(fx,fy,fz,x,y,z,vx,vy,vz,npart,dt,ltot,vecinos,nvecinos,potencial,t_target,q,zeta, offont,qp,p_target,eta,onb,u)
     if(mod(nmezcla, 10000).eq. 0) then
      sumv2=0.d0
      do i=1,npart
       sumv2=sumv2+vx(i)**2 +vz(i)**2 +vy(i)**2
       enddo
      print*, sumv2/(3.d0*npart), 100*real(nmezcla)/pasos 
      if(sumv2/(3.d0*npart) .lt. t_target) then 
        exit
        endif
      endif
     nmezcla=nmezcla+1

     if(mod(nmezcla, 1000) .eq. 0) then
        call reajuste(npart, vx,vy,vz)
       endif


  enddo

   sumv2=0.d0
   do i=1,npart
    sumv2=sumv2+vx(i)**2+vy(i)**2 +vz(i)**2
   enddo
   print*, sumv2/(3.d0*npart), "temeratura termalizacion"
   
    pasos=20000000
    dt=0.0001d0
    offont=0!!TERMOSTATO encendido
    q=30
    qp=500*npart
    onb=1
    t_target=0.55d0
    p_target=2.d0
 call listaverlet(npart, x,y,z, rcutoff,rskin, ltot, vecinos, nvecinos)
 call fuerzas(npart, x,y,z,fx,fy,fz, vecinos, nvecinos, ltot,potencial)

   onb=1!!!barostato encendido
   contador=0
   u=1
   print*, "empezamos enfriamiento y compresión"
   nmezcla=0
   cinetica=100000  !!PARA QUE NO SALGA DEL BUCLE
   
  do while(nmezcla .lt. pasos)
     
    if(mod(nmezcla,10) .eq. 0) then 
     call listaverlet(npart, x,y,z, rcutoff,rskin, ltot, vecinos, nvecinos) 
     endif

 call integrarmove(fx,fy,fz,x,y,z,vx,vy,vz,npart,dt,ltot,vecinos,nvecinos,potencial,t_target,q,zeta, offont,qp,p_target,eta,onb,u)
    nmezcla=nmezcla+1
     if(mod(nmezcla, 1000) .eq. 0) then
      print*, npart/ltot**3
      if(npart/(ltot**3) .gt. 1.2) then 
      print*, "se acabo la compresion"
      print*, "Ahora empezamos a medir con termostato activo t=0.33"
      exit
      endif
     endif
   if(mod(nmezcla, 10000) .eq. 0) then

     cinetica=0.d0
   do i=1,npart
  cinetica=cinetica+(vx(i)**2 +vy(i)**2 +vz(i)**2)
   enddo
   
     
    !print*, cinetica*0.5d0
  !stop "message"
    cinetica=0.5*cinetica/npart
    potencial=potencial/npart
  
    energia=cinetica+potencial
  !print*, cinetica, potencial, energia
    pantalla=pantalla+1
     write(33,*) nmezcla*dt, energia, potencial,cinetica
     if(mod(pantalla, 10) .eq. 0) then
     print*, real(nmezcla/real(pasos))*100, "%" , npart/ltot**3 , cinetica*(2.d0/3.d0)
     endif
     if(cinetica*(2.d0/3.d0) .gt. 0.7)then 
      offont=1
     elseif(cinetica*(2.d0/3.d0).lt.0.5d0 .and. offont .eq. 1)then
      offont=0
     endif
    endif
        
    
      !  if(t_target .gt. (cinetica*2.d0/(3.d0))) then
       
      !  print*, t_target, cinetica*(2.d0/3.d0), cinetica
      !  print*, "salmos del calentamiento"
      !  exit
       

      ! !  dt=0.0002d0
       
      ! !  !cinetica=0.d0
      ! !  do j=1, npart
      ! !  xacumu(j)=x(j)+xacumu(j)
      ! !  yacumu(j)=y(j)+yacumu(j)
      ! !  zacumu(j)=z(j)+zacumu(j)
      ! !  !cinetica=vx(j)**2+vy(j)**2+vz(j)**2 +cinetica
      ! !  enddo
      
      ! !  enddo
      
      ! !  endif
    
      ! !  endif
    
      ! !  if(mod(i, 2000) .eq. 0) then 
      ! !     call reajuste(npart, vx,vy,vz)
          
      !  endif
       if(mod(nmezcla, 1000) .eq. 0) then
        call reajuste(npart, vx,vy,vz)
       endif

      !  if(mod(nmezcla,pasos/npromedio) .eq. 0) then
        
      !   binsradial=150
      ! !  enddo
      
      ! !  endif
    
      !   call gradial(npart, x, y, z, ltot, binsradial,dr, gr,rcontar,HISTO) !!entrada nbins npart y las posiciones
      !   do i=1, binsradial
      !   grprom(i)=grprom(i)+gr(i)
      !   enddo
      !  endif


  enddo
  print*, "Enfriamiento terminado"   
  close(200)
  
   !!Desues del enfriamiento vamos a medir 
    offont=1
    pasos=400000
    nmezcla=0
    onb=0
    q=150
    t_target=0.66d0
    check=0
    do while(nmezcla .lt. pasos)
     
     if(mod(nmezcla,10) .eq. 0) then 
     call listaverlet(npart, x,y,z, rcutoff,rskin, ltot, vecinos, nvecinos) 
     endif

     if(mod(nmezcla, 20000) .eq. 0 .and. check .eq. 0) then
      temperatura=0.d0
       
        do i=1,npart
         temperatura=(vx(i)**2 + vy(i)**2 + vz(i)**2)+temperatura
         enddo
         temperatura=temperatura/(3*npart)
         print*, temperatura, "temperatura en la medida"
         endif
         

   call integrarmove(fx,fy,fz,x,y,z,vx,vy,vz,npart,dt,ltot,vecinos,nvecinos,potencial,t_target,q,zeta, offont,qp,p_target,eta,onb,u)
    nmezcla=nmezcla+1
    if(mod(nmezcla,pasos/npromedio) .eq. 0 .and. nmezcla .gt. pasos/2.d0) then
        
        !!vamos a medir tambien la distribucion de velocidades
        vmodu=0.d0
        do i=1,npart
         vmodu(i)=sqrt(vx(i)**2 + vy(i)**2 + vz(i)**2)
         enddo
         
         px=0.d0
        call histograma(vmodu, bins, npart, px, a, b)
        do i=1, bins
          pxprom(i)=pxprom(i)+ px(i)
          enddo

        binsradial=100
        call gradial(npart, x, y, z, ltot, binsradial,dr, gr,rcontar,HISTO) !!entrada nbins npart y las posiciones
        do i=1, binsradial
        grprom(i)=grprom(i)+gr(i)
        enddo
       endif
      

    enddo
      

      open(88,file="fotofinish.txt")
      do i=1,npart
 !print*,  x(i),y(i),z(i)
      write(88,*) x(i),y(i),z(i)
 
      enddo
        close(88)


        binsradial=100
    call gradial(npart, x, y, z, ltot, binsradial,dr, gr,rcontar,HISTO) !!entrada nbins npart y las posiciones
     open(89,file="grfinish_sinpromediar.txt")
     
     do i=1,binsradial
      write(89,*) (i-0.5d0)*dr, gr(i)
      enddo
      close(89)
       gr=0.d0
       rcontar=0.d0
       histo=0.d0
      print*, "distribucion gr calentamiento terminada"
     !!!ya hemos cogido los ultimos npromedio posiciones
     x=2*xacumu/npromedio
     y=2*yacumu/npromedio
     z=2*zacumu/npromedio
    


     grprom=2*grprom/npromedio

     open(899, file="gr_final.txt")
     do i=1,binsradial
      write(899,*) (i-0.5d0)*dr, grprom(i)
      enddo
      close(899)
       gr=0.d0
       rcontar=0.d0
       histo=0.d0
      print*, "distribucion gr promediada terminada"
     !!!ya hemos cogido los ultimos npromedio posiciones
     pxprom=(2*pxprom)/(npromedio)
     open(933,file="velocidades_barostato.txt")

     do i=1,bins
      write(933,*) a+(b-a)*i/bins, pxprom(i)/(b-a)
      enddo






end program

 subroutine posicioninicial(x,y,z,r,lcil,npart,lcaja,diametro,count)
  implicit none
  integer i, j, k, npart, nx, ny, nz, count
  real*8 r, lcil, celda,delta,xrnd,yrnd,zrnd,lcaja
  real*8 xf,yf,zf,diametro
  real*8 :: x(1:npart)
  real*8 :: y(1:npart)
  real*8 :: z(1:npart)
  
  ! Espaciado inicial
  print*, npart
  celda=lcaja
  delta= 0.d0*celda  !!con esto podemos generar algun solapamiento, veremos si nos puede dar problemas
  print*, lcaja, r, npart
  nx = nint(npart**(1.d0/3))
  ny = nint(npart**(1.d0/3))
  nz = nint(npart**(1.d0/3))
   
  count = 0
  do i = 0, nx-1
    do j = 0, ny-1
      do k = 0, nz-1
        if (count .ge. npart) then
        exit
        else
        x(count+1) = i*celda + r
        y(count+1) = j*celda + r
        z(count+1) = k*celda + r
        count = count + 1
        endif
      end do
    end do
  end do
  print*, count

  !Perturbamos las posiciones iniciales ya que las hemos colocado de manera ordenada
  do i = 1, count
    call random_number(xrnd)
    xf = (x(i) + delta*(xrnd - 0.5d0))
    if(xf .gt. 0.d0) then
    x(i)=xf-lcil*floor(x(i)/lcil)
    endif
    !x(i)=x(i)-lcil*floor(x(i)/lcil)
    call random_number(yrnd)
    yf = (y(i) + delta*(yrnd - 0.5d0))
    if(yf .gt. 0.d0) then
    y(i)=yf-lcil*floor(y(i)/lcil)
    endif
    !y(i)=y(i)-lcil*floor(y(i)/lcil)
    call random_number(zrnd)
    zf = (z(i) + delta*(zrnd - 0.5d0))
    if(zf .gt. 0.d0) then
    z(i)=zf-lcil*floor(z(i)/lcil)
    endif
    !z(i)=z(i)-lcil*floor(z(i)/lcil)
    
    
  end do

 end subroutine




subroutine listaverlet(n, x,y,z, rcutoff,rskin, l, vecinos, nvecinos)  !!esta lista se actualiza cada cierto tiempo
     implicit none
     integer i,j
     integer n
     real*8  rcutoff, l ,dx,dy,dz,r2,rc2, rskin
     
     real*8 :: x(1:n)
     real*8 :: y(1:n)
     real*8 :: z(1:n)
     integer :: vecinos(1:n,1:n) ! Lista de vecinos y número de vecinos por partícula
     integer :: nvecinos(1:n) 
     
    rc2 = (rcutoff+rskin)**2.d0  !!radio de corte/busqueda  
    nvecinos= 0  ! Inicializar el número de vecinos a 0
    vecinos=0
     
    do i = 1, n -1
      
        do j = i + 1, n
            dx = x(j)-x(i)
            dy = y(j)-y(i)
            dz = z(j)-z(i) 

            ! Condición de contorno periódica
            dx = dx -l*anint(dx/l)
            dy = dy -l*anint(dy/l)
            dz = dz -l*anint(dz/l)
            
            
            r2 = dx*dx + dy*dy + dz*dz  

            if (r2 .lt. rc2) then
                nvecinos(i) =nvecinos(i)+1
                nvecinos(j)=nvecinos(j)+1
                vecinos(i, nvecinos(i)) = j
                vecinos(j, nvecinos(j)) = i
            end if
        end do
    end do
end subroutine


subroutine fuerzas(n, x,y,z,fx,fy,fz, vecinos, nvecinos, L,potencial)  !!la longitud ya tiene que venir reducida 
    implicit none
    integer :: i, j, k,n
    integer :: vecinos(1:n,1:n)
    integer :: nvecinos(1:n)
    real*8 :: epsilon, sigma, L, rcutoff
    real*8 :: x(1:n)
    real*8 :: y(1:n)
    real*8 :: z(1:n)
    real*8 :: fx(1:n)
    real*8 :: fy(1:n)
    real*8 :: fz(1:n)
    real*8 :: rguardar(1:n,1:n)
    
    real*8 ecut,potencial,rmin,A,alpha
    real*8 :: dx, dy, dz, r2, inv_r2, inv_r6, inv_r12, flj
    parameter(rmin=0.88d0)  !!valor minimo recomendado para rmin
    parameter(rcutoff=2.5d0)
    parameter(A=1.84d0)
    parameter(alpha=14.d0)
    ecut=4.d0*(1.d0/(rcutoff)**6)*((1.d0/(rcutoff)**6)-1.d0)
    
    fx=0.d0
    fy=0.d0
    fz=0.d0
    potencial=0.d0
    do i = 1, n
        do k = 1, nvecinos(i)  !!!
            j = vecinos(i, k)  !!!
            
            dx = x(j)-x(i)
            dy = y(j)-y(i)
            dz = z(j)-z(i)
            
            dx= dx -L*anint(dx/L)
            dy= dy -L*anint(dy/L)
            dz= dz -L*anint(dz/L)
             !!!Aqui tendriamos que meter algo para evitar solapamientos
            r2 = dx*dx + dy*dy + dz*dz
            !!vamos a fijar un radio minimo r=0.88 por ejemplo para que no hay divergencias
            if(r2 .lt. 1.d0) then 
             if(r2 .lt. rmin) then 
             r2=rmin**2
             endif
             inv_r2 = 1.d0/r2
             inv_r6 = inv_r2**3.d0
             flj=alpha*A*exp(-alpha*sqrt(r2))/sqrt(r2) -24.d0*inv_r6*inv_r2
             

             fx(i)=fx(i)+flj*dx
             fy(i)=fy(i)+flj*dy
             fz(i)=fz(i)+flj*dz

             potencial=potencial-A*exp(-alpha*sqrt(r2))- 4.d0/(r2**3) -ecut
             else

             
             inv_r2 = 1.d0/r2
             inv_r6 = inv_r2**3.d0
            !inv_r12 = inv_r6*inv_r6
   
             flj = 48.d0*inv_r6*inv_r2*(inv_r6 -0.5d0) !! unidades reducidad F*=F*diametro/epsilon
                                                      !!r**-8 y r**-14 que multiplicando por elm vectro posicion r ya tenemso la 
                              
             fx(i)=fx(i)+flj*dx
             fy(i)=fy(i)+flj*dy
             fz(i)=fz(i)+flj*dz

            
             potencial=potencial+4.d0*inv_r6*(inv_r6-1.d0)-ecut    !!!actualizacion de la energia 
              endif                                        !! ecut es le energia en r=rcutoff, desplazamos el potencial para en que r=rc sea 0(para minimizar errores numericos)
        end do
    end do
    potencial=potencial/2.d0  !!se cuentas todos los pares i,j por separado, por lo que hay que dividir entre dos
end subroutine


subroutine integrarmove(fx,fy,fz,x,y,z,vx,vy,vz,npart,dt,l,vecinos,nvecinos,potencial,t_target,q,zeta,offont,qp,p_target,eta,onb,u)
implicit none 

integer i,j,npart,offont, onb,u, promediador
real*8 dt,sumav,sumav2,xx,vi,yy,zz,l
real*8 potencial, t_target,zeta,q,qp, p_target,eta
real*8 t_prom, virial_prom


real*8 :: x(1:npart)
real*8 :: vx(1:npart)
real*8 :: y(1:npart)
real*8 :: vy(1:npart)
real*8 :: z(1:npart)
real*8 :: vz(1:npart)
real*8 :: fx(1:npart)
real*8 :: fy(1:npart)
real*8 :: fz(1:npart)
real*8 :: fxn(1:npart)
real*8 :: fyn(1:npart)
real*8 :: fzn(1:npart)
integer :: vecinos(1:npart,1:npart)
integer :: nvecinos(1:npart)


sumav=0.d0
sumav2=0.d0
 
 !!ACTUALIZAMOS LAS NUEVAS POSICIONES
 do i=1,npart
  x(i)=x(i)+vx(i)*dt +0.5*fx(i)*dt**2
  y(i)=y(i)+vy(i)*dt +0.5*fy(i)*dt**2
  z(i)=z(i)+vz(i)*dt +0.5*fz(i)*dt**2

  x(i)=x(i)-l*floor(x(i)/l)
  y(i)=y(i)-l*floor(y(i)/l)
  z(i)=z(i)-l*floor(z(i)/l)

  
 enddo
 !!LLAMAMOS A LAS NUEVAS FUERZAS
  
  !!con las nuevas fuerzas y las anteriores actualizamos las velocidades
  
   do i=1,npart  !!Actualizamos paso intemredio con las fuerzas antiguas
   vx(i)=vx(i)+ 0.5d0*fx(i)*dt
   vy(i)=vy(i)+ 0.5d0*fy(i)*dt
   vz(i)=vz(i)+ 0.5d0*fz(i)*dt
    enddo
    !!EN ESTE PUNTO LLAMAMOS AL TERMOSTATO
    if(offont .eq. 1) then
    CALL nosehoover(npart,vx,vy,vz,dt,zeta,t_target,q)  !!Termo
    endif

   
   if(onb .eq. 1) then
    call barostato(npart,x,y,z,vx,vy,vz,fx,fy,fz,dt,qp,p_target,eta,l,u, promediador,vecinos, nvecinos,t_prom, virial_prom) !!Baro
    
    endif
   call fuerzas(npart, x,y,z,fxn,fyn,fzn, vecinos, nvecinos, L,potencial) !!Nuvea sfuerzas
  !!MOdificando temperatura.......
  do i=1,npart  !!Actualizamos velocidades
   vx(i)=vx(i)+ 0.5d0*fxn(i)*dt
   vy(i)=vy(i)+ 0.5d0*fyn(i)*dt
   vz(i)=vz(i)+ 0.5d0*fzn(i)*dt
    enddo

   do i=1,npart  
    fx(i)=fxn(i) 
    fy(i)=fyn(i) 
    fz(i)=fzn(i) 
    enddo
 return   !!nos devuelve las posiciones, velocidades y energia potencial
end subroutine 


subroutine gradial(npart, x, y, z, L, nbins,dr, gr,rcontar,histo) !!entrada nbins npart y las posiciones
  implicit none
  integer ::  npart, nbins
  real*8 :: x(1:npart), y(1:npart), z(1:npart), L
  real*8 :: gr(1:nbins),dr
  
  integer :: i, j, bin,contar,k
  real*8 :: dx, dy, dz, r, r2, vol_bin, norma, rho,pi
  integer :: histo(1:nbins)
  real*8 :: rcontar(1:npart**2)
  parameter(pi=4*atan(1.d0))
   contar=0
  histo = 0
  dr=(L/2.d0)/nbins  !!nbins lo metemos nosotros dr= tiene que ser entre 0.01
  k=1
  do i = 1, npart-1
    do j = i+1, npart
      dx = x(j)-x(i)
      dy = y(j)-y(i)
      dz = z(j)-z(i)

      dx = dx - L*nint(dx/L)
      dy = dy - L*nint(dy/L)
      dz = dz - L*nint(dz/L)

      r2 = dx**2 + dy**2 + dz**2

      if (r2 .lt. (L/2.d0)**2) then
       
        r = sqrt(r2)
        
        if(r .lt. 0.5d0) then !!!!!!!!
        print*, r, contar
        contar=contar+1
        elseif(r .ge. 0.88d0) then
        bin = int(r/dr) + 1
        rcontar(k)=r
        k=k+1
        if (bin .ge. 1.d0 .and. bin .le. nbins) then
          histo(bin) = histo(bin) + 2  ! Cada par cuenta para ambas partículas
        end if
      end if
      endif
      
    end do
  end do
  !print*, k , "pares de distancias para l/2"
  rho=npart/ l**3
  !Normalización de g(r)
  do bin = 1, nbins
    r = (bin - 0.5d0)*dr+0.12d0    !!nuestro 1 en la simulacion es 0.88 que es donde hacemos la f max, por o tanto desplazamos el r 0.12 
    vol_bin = 4.d0*pi* r**2 * dr  ! Volumen de la capa esférica
    norma = rho * npart * vol_bin
    if (norma .gt. 0.d0) then
      gr(bin) = dble(histo(bin)) / norma
    else
      gr(bin) = 0.0
    end if
  end do

   
end subroutine 



!!Inicializacion de direcciones aleorias------------------------------------------------------------------------------

subroutine direcciones(v1x,v1y,v1z)
 implicit none 

 integer i 
 real*8 r,v1x,v1y,v1z
 real*8 x,y,z,modus


 parameter(r=1.d0)
 modus=2.d0 
 do while(modus .gt. r)

 call random_number(x)
 v1x=(2.d0*x-1.d0)
 call random_number(y)
 v1y=(2.d0*y-1.d0)
 call random_number(z)
 v1z=(2.d0*z-1.d0)

 modus=dsqrt(v1x**2 +v1y**2 +v1z**2)

 enddo


 v1x=v1x/modus
 v1y=V1y/modus  !!Normalizacion
 v1z=v1z/modus

 return
end subroutine

subroutine reajuste(npart, vx,vy,vz)
integer i,npart
real*8 sumvx,sumvy,sumvz
real*8 :: vx(1:npart)
real*8 :: vy(1:npart)
real*8 :: vz(1:npart)

    sumvx=0.d0
    sumvy=0.d0
    sumvz=0.d0
  
  do i=1,npart
   sumvx=sumvx+vx(i)
   sumvy=sumvy+vy(i)
   sumvz=sumvz+vz(i)
   
   enddo
  sumvx=sumvx/npart
  sumvy=sumvy/npart
  sumvz=sumvz/npart
  do i=1,npart
   vx(i)=vx(i)-sumvx
   vy(i)=vy(i)-sumvy
   vz(i)=vz(i)-sumvz
   enddo
   return
end subroutine


 subroutine histograma(y, ncajas, n, px, a, b)
      implicit none

      integer ncajas, n,bin,i
      real*8  a, b,h
      real*8, intent(in) :: y(1:n)
      real*8, intent(out) :: px(1:ncajas)

         ! Inicializar px a cero
       px = 0.0
       b=maxval(y)
       a=minval(y)
  
       h = (b-a)/ncajas
       
       do i = 1, n
       if (y(i) .ge. a .and. y(i) .lt. b) then !
         bin = int((y(i)-a)/h) + 1
        if (bin .le. ncajas) then
        px(bin) = px(bin) + 1
         endif
       endif
       end do

         ! Normalizar px
        px = px/n

     end subroutine


     function f(v,t)
     real*8 f,v,pi,t

     parameter(pi=4.d0*atan(1.d0))

     f=8.d0*pi**(-0.5d0) *t*v**2 *exp(-(2.d0*t*v**2))

     return
     end  function



     subroutine nosehoover(npart,vx,vy,vz,dt,zeta,t_target,q)
     implicit none

     integer i,j,k,npart
     real*8 l,potencial,zeta,q,temperatura,t_target,dt
       !!!el valor de q es crucial, nos sale muy pequueño para lo que teoricamente deberia, pero si queremos 
    !!ver un cambio en la tempreratura "rapido" no nos queda otra.

     real*8 :: vx(1:npart)
     real*8 :: vy(1:npart)
     real*8 :: vz(1:npart)
 
     !!calcular la temperatura
     temperatura=0.d0
     
     do i=1,npart
     temperatura=(vx(i)**2 +vy(i)**2 +vz(i)**2)+temperatura
     enddo
     temperatura=temperatura/(3.d0*npart)
     !print*, temperatura, t_target
     
     !!actualizar zeta
     zeta=zeta+(dt/q)*(temperatura/t_target -1.d0)

     !!En q absorbemos el factor t_target ya que la expresion teorica dice que zeta=dt/Q *(t - t_target)
     !print*, zeta
     !!!aplicamos el efecto del termostato
      !print*, vx(1)
      do i=1,npart
      
      vx(i)=vx(i)*exp(-zeta*dt)
      vy(i)=vy(i)*exp(-zeta*dt)
      vz(i)=vz(i)*exp(-zeta*dt)

      enddo
     ! print*, vx(1)
     return
     
     end subroutine

    subroutine barostato(npart,x,y,z,vx,vy,vz,fx,fy,fz,dt,qp,p_target,eta,l,u, promediador,vecinos, nvecinos,t_prom,virial_prom)
     implicit none

     integer i,npart,u, promediador,k,j
     real*8 l,eta,qp,P_actual,p_target,dt, temperatura
     real*8 vol,xmax,ymax,zmax,virial
     real*8 rcutoff, rskin
     real*8 t_prom, virial_prom,pi, p_tail

      parameter(rskin=0.5d0)
      parameter(rcutoff=2.5d0)
      parameter(pi=4*atan(1.d0))

     real*8 :: x(1:npart)
     real*8 :: y(1:npart)
     real*8 :: z(1:npart)
     real*8 :: vx(1:npart)
     real*8 :: vy(1:npart)
     real*8 :: vz(1:npart)
     real*8 :: fx(1:npart)
     real*8 :: fy(1:npart)
     real*8 :: fz(1:npart)

     integer :: vecinos(1:npart,1:npart)
     integer :: nvecinos(1:npart)

       
        
        vol=l**3
       temperatura=0.d0
       do i=1, npart
      temperatura=(vx(i)**2 +vy(i)**2 +vz(i)**2)+temperatura
       enddo
       temperatura=temperatura/(3.d0*npart)
      call vvirial(npart, x,y,z,fx,fy,fz, vecinos, nvecinos, L, virial)
       !!!el virial sale normalizado virial/(3*V)
       
      
       p_tail=((16.d0*pi*npart/vol)/3.d0)*((1.d0/(3.d0*rcutoff**9)) -(1.d0/(rcutoff**3)))
       
       !print*, virial, temperatura*npart/vol, p_tail

       P_actual=npart*temperatura/vol + virial + p_tail
          !!!!!!!!!!!!!!!!!!!!
       

       !print*, P_actual, p_target, "presiones"
        
        
       
       P_target=2.d0
     !print*, P_actual, p_target, "presiones"
       eta = eta + (dt/QP)*(P_actual/P_target-1.d0)
     
     !print*, exp(eta*dt), eta, P_target, dt, qp, P_actual
     !print*, P_actual, "ppp", temperatura/vol , virial
     !print*, temperatura, P_actual, npart/vol, "densidad" 
     
      ! if(u .eq. 1) then
      !  eta=0.d0
      !  endif
      
      do i=1,npart
      x(i)=x(i)*exp(eta*dt)
      y(i)=y(i)*exp(eta*dt)
      z(i)=z(i)*exp(eta*dt)

      vx(i)=vx(i)*exp(-eta*dt)
      vy(i)=vy(i)*exp(-eta*dt)          !!!explicar en teoria por que 
      vz(i)=vz(i)*exp(-eta*dt)

      enddo

      l=l*exp(eta*dt)
      
     return
     end subroutine


     subroutine vvirial(n, x,y,z,fx,fy,fz, vecinos, nvecinos, L, virial)  !!la longitud ya tiene que venir reducida 
    implicit none
    integer :: i, j, k,n, contador
    integer :: vecinos(1:n,1:n)
    integer :: nvecinos(1:n)
    real*8 :: epsilon, sigma, L, rcutoff, virial
    real*8 :: x(1:n)
    real*8 :: y(1:n)
    real*8 :: z(1:n)
    real*8 :: fx(1:n)
    real*8 :: fy(1:n)
    real*8 :: fz(1:n)
    real*8 :: rguardar(1:n,1:n)
    
    
    real*8 ecut,potencial,rmin,A,alpha
    real*8 :: dx, dy, dz, r2, inv_r2, inv_r6, inv_r12, flj
    parameter(rmin=0.88d0)  !!valor minimo recomendado para rmin
    parameter(rcutoff=2.5d0)
    parameter(A=1.84d0)
    parameter(alpha=14.d0)
   
    
    fx=0.d0
    fy=0.d0
    fz=0.d0
    virial=0.d0
 

    
    do i = 1, n
        do k = 1, nvecinos(i)  !!!
            j = vecinos(i, k)  !!!
            !print*, j
            if(j .gt. i) then

            dx = x(j)-x(i)
            dy = y(j)-y(i)
            dz = z(j)-z(i)
            
            dx= dx -L*anint(dx/L)
            dy= dy -L*anint(dy/L)
            dz= dz -L*anint(dz/L)
             !!!Aqui tendriamos que meter algo para evitar solapamientos
            r2 = dx*dx + dy*dy + dz*dz
            !!vamos a fijar un radio minimo r=0.88 por ejemplo para que no hay divergencias
            if(r2 .lt. 1.d0 .and. sqrt(r2) .gt. 0.88d0) then    !!!No vamos a considerar para el calculo del virial
            !  if(r2 .lt. rmin) then 
            !  r2=rmin**2
            !  endif
            !  inv_r2 = 1.d0/r2
            !  inv_r6 = inv_r2**3.d0
            !  flj=alpha*A*exp(-alpha*sqrt(r2))/sqrt(r2) -24.d0*inv_r6*inv_r2
            !   !!!METER EL PARA IR GUARDANDO LOS TERMINOS DEL VIRIAL
            !   virial = virial - sqrt(r2)*flj
            !   !print*, flj, "flj"

             else

             inv_r2 = 1.d0/r2
             inv_r6 = inv_r2**3.d0
            !inv_r12 = inv_r6*inv_r6
   
             flj = 48.d0*inv_r6*inv_r2*(inv_r6 -0.5d0) 
             !!!METER EL PARA IR GUARDANDO LOS TERMINOS DEL VIRIAL
             virial = virial - sqrt(r2)*flj
             !print*, flj, "flj"
                              
              
              endif      
              endif                                  
        end do
    end do
    virial = virial/(3.d0* L**3)
    !print*, virial, "virial"
    
    return
    
end subroutine





      



    


!!------------RANDOM NUMBER-----------------------------------------------!!

subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    ! Inicializar el tamaño de la semilla
    call random_seed(size = n)
    allocate(seed(n))

    ! Usar el reloj del sistema si no hay acceso a /dev/urandom
    call system_clock(count)
    if (count /= 0) then
        t = transfer(count, t)
    else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
        t = transfer(tms, t)
    end if

    ! Calcular una semilla usando el tiempo y el PID
    s = ieor(t(1), t(2))
    pid = 1099279  ! Puedes usar una constante si no puedes obtener el PID
    s = ieor(s, pid)

    ! Asignar valores a la semilla
    if (n >= 3) then
        seed(1) = t(1) + 36269
        seed(2) = t(2) + 72551
        seed(3) = pid
        if (n > 3) then
            seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
    else
        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
    end if

    ! Inicializar el generador de números aleatorios con la semilla generada
    call random_seed(put=seed)

end subroutine init_random_seed
