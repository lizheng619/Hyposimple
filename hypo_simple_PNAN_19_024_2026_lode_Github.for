c-----------------------------------------------------------------------------
c Hyposimple model: A simple hypoplastic model for sand under drained condition
c and simple loading paths 
c-----------------------------------------------------------------------------
c UMAT subroutine developped by Zheng LI at University Gustave Eiffel
c email: zheng.li@univ-eiffel.fr
c------------------------------------------------------------------------------
c User subroutine for Abaqus 6.24
c-----------------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     &rpl,ddsddt,drplde,drpldt,
     &stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     &ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     &celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c------------------------------------------------------------------------------
      implicit none
c
      character*80 cmname
c
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt,
     &layer, kspt, kstep, kinc
c
      double precision stress(ntens), statev(nstatv),
     &ddsdde(ntens,ntens), ddsddt(ntens), drplde(ntens),
     &stran(ntens), dstran(ntens), time(2), predef(1), dpred(1),
     &props(nprops), coords(3), drot(3,3), dfgrd0(3,3), dfgrd1(3,3)
      double precision sse, spd, scd, rpl, drpldt, dtime, temp,
     &dtemp, pnewdt, celent
c------------------------------------------------------------------------------
c
      integer i,j
      integer nparms, counter
      integer flag, selection
c
      double precision stress_o(ntens),dstran_dt(ntens),dT(ntens)
      double precision cohesion,maxT,expo,parms(nprops)
      double precision DDe(6,6)
      double precision dev_sig_pre(6),sig_trace_pre,sig_p_pre
      double precision sig_e_pre,norm_sig_pre(ntens)
      double precision G_elas,nu_elas
      double precision theta,perturb
      double precision zero,one,two,three,six,pi
      double precision dot_vect
      double precision A,B,C,e_sand,sigma_cf,G_lit
      double precision sig_trial(ntens),dsig_trial(ntens)
      double precision norm_sig_trial(ntens),norm_sig(ntens)
      double precision dev_sig_ck(6),sig_trace_ck,sig_p_ck
      double precision sig_e_ck,norm_sig_ck(ntens)
      double precision M,M_d
      double precision fric_ang,fric_ang_cs,dila_ang
      double precision volume,S,volume_new
      double precision e_ini,coeff_CSL
      double precision p_atm,e0,lambda_c,xi_c
      double precision K_T(6,6)
      double precision f_y,d_lambda_trail,tol_fy
      double precision temp_vector(ntens),temp_temp
      double precision aCa,CaCa(ntens,ntens)
      double precision iter_max
      double precision load_factor
      double precision Lod,theta_lode,coef,cos3t,Lode_factor,c_factor
c
      parameter (zero = 0.0d0,one=1.0d0,two=2.0d0)
      parameter (three = 3.0d0,six=6.0d0,pi=3.1415926d0)
      parameter (perturb = 1.0d-10)
      parameter (tol_fy = 1.0d-9)
      parameter (iter_max = 500)
      parameter (selection = 1)
c-----------------------------------------------------------------------------------------------
      nparms=nprops
      call pzero(parms,nparms)
      call push(props,parms,nparms)
      G_elas      = parms(1)
      nu_elas     = parms(2)
      expo        = parms(3)
      cohesion    = parms(4)
      fric_ang_cs = parms(5)
      dila_ang    = parms(6)
      e_ini       = parms(7)
      p_atm       = parms(8)
      e0          = parms(9)
      lambda_c    = parms(10)
      xi_c        = parms(11)
c-----------------------------------------------------------------------------------------------
      e_sand=statev(7)
c-----------------------------------------------------------------------------------------------
      do i=1,ntens
         stress(i) = -stress(i)
      end do
      do i=1,ntens
         dstran(i) = -dstran(i)
      end do
c----------------------------------------------------------------------------------
      call lode_DM(stress,theta_lode,cos3t)
      theta_lode = anint(theta_lode*10.0**2)/10.0**2
c
      c_factor = 0.712d0
c
      Lode_factor = 2.0d0*c_factor/((1+c_factor)-(1-c_factor)*cos3t)
      Lode_factor = anint(Lode_factor*10.0**2)/10.0**2
      if (lode_factor .eq. zero) then
          lode_factor =  one
      endif
c-----------------------------------------------------------------------------------------------
      call pzero(stress_o,ntens)
      call push(stress,stress_o,ntens)
c----------------------------------------------------------------------------------------------
      sigma_cf=(stress(1)+stress(2)+stress(3))/3
      if ((stress(1) .le. zero) .or. (stress(2) .le. zero) .or. 
     &(stress(3) .le. zero)) then
      if (selection .eq. zero) then
         write(*,*) 'Tension stress !!! '
         write(*,*) 'Element number   = ',noel
         write(*,*) 'Step number      = ',kstep
         write(*,*) 'Increment number = ',kinc
         write(*,*) 'Current time     = ',time
         write(*,*) 'Time increment   = ',dtime
      endif
          fric_ang_cs = zero
          dila_ang    = zero
          M           = zero
          m_d         = zero
          maxT        = cohesion
          call el_stiff(parms,nparms,stress,ntens,statev,nstatv,DDe)
          call matrix_vect(DDe,ntens,ntens,dstran,dsig_trial)
          do i=1,ntens
             sig_trial(i) = stress(i)+dsig_trial(i)
          enddo
         call deviator(sig_trial,dev_sig_pre,sig_trace_pre,sig_p_pre)
         call f_yd(dev_sig_pre,ntens,maxT,M_d,sig_e_pre,norm_sig_trial,
     $f_y)
c
      if (f_y.le.zero) then
c
          flag = 0
c
          do i=1,ntens
              stress(i) = sig_trial(i)
          enddo
      else
c
         flag = 1
c
        counter = 1
        do while ((dabs(f_y) .gt. tol_fy) .and. (counter .le. iter_max))
           call matrix_vect(DDe,ntens,ntens,norm_sig_trial,temp_vector)
           temp_temp = dot_vect(3,norm_sig_trial,temp_vector,ntens)
           d_lambda_trail = f_y/temp_temp
           do i=1,ntens
             sig_trial(i) = sig_trial(i)-d_lambda_trail*temp_vector(i)
           enddo
          call deviator(sig_trial,dev_sig_pre,sig_trace_pre,sig_p_pre)
          call f_yd(dev_sig_pre,ntens,maxT,M_d,sig_e_pre,norm_sig_trial,
     $f_y)
         counter = counter + 1
       enddo
          do i=1,ntens
              stress(i) = sig_trial(i)
          enddo
      endif
c
      if (flag .eq. 0 ) then
         call el_stiff(parms,nparms,stress,ntens,statev,nstatv,DDe)
         do i=1,6
           do j=1,6
               DDsDDe(i,j)=DDe(i,j)
           enddo
         enddo
      else
cDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDE
         call el_stiff(parms,nparms,stress,ntens,statev,nstatv,DDe)
          do i=1,6
           do j=1,6
               DDsDDe(i,j)=zero
           enddo
          enddo
c------------------------------------------------------------------------------------------------
      call el_stiff(parms,nparms,stress,ntens,statev,nstatv,DDe)
      call deviator(stress,dev_sig_pre,sig_trace_pre,sig_p_pre)
      call f_yd(dev_sig_pre,ntens,maxT,M_d,sig_e_pre,norm_sig,
     $f_y)
      call matrix_vect(DDe,ntens,ntens,norm_sig,temp_vector)
      call cross_vect(temp_vector,temp_vector,ntens,CaCa)
      aCa = dot_vect(3,norm_sig,temp_vector,ntens)
      do i=1,6
        do j=1,6
         DDsDDe(i,j) = DDe(i,j)-CaCa(i,j)/aCa
        enddo
      enddo
cDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDEDDSDDE
      endif
c
      else !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
c
      call get_state(parms,nparms,stress,ntens,e_sand,e_ini,coeff_CSL)
c
      dila_ang = dila_ang*coeff_CSL
c
      if (coeff_CSL .gt. zero) then
         fric_ang = fric_ang_cs+0.8d0*dila_ang
      else
         fric_ang = fric_ang_cs+0.8d0*dila_ang
      endif
c--------------------------------------------------------------------
      M=six*sin(fric_ang*pi/180)/(three-sin(fric_ang*pi/180))
      maxT=sigma_cf*M
      maxT = maxT+cohesion
c
      M_d = six*sin(dila_ang*pi/180)/(three-sin(dila_ang*pi/180))
	M_d = M_d*lode_factor
c
c-----------------------------------------------------------------------------------------------
      call deviator(stress,dev_sig_pre,sig_trace_pre,sig_p_pre)
c-----------------------------------------------------------------------------------------------
      call f_yd(dev_sig_pre,ntens,maxT,M_d,sig_e_pre,norm_sig_pre,f_y)
c-----------------------------------------------------------------------------------------------
      call el_stiff(parms,nparms,stress,ntens,statev,nstatv,DDe)
      call pzero(dstran_dt,ntens)
c----------------------------------------------------------------------------------------------
      do i=1,ntens
         dstran_dt(i)=dstran(i)/dtime
      enddo
c
      call stress_rate(DDe,maxT,stress,dstran_dt,expo,
     &M,M_d,dtime,ntens,dT,K_T)
c-----------------------------------------------------------------------------------------------
      if (kstep .le. 1) then
         call el_stiff(parms,nparms,stress,ntens,statev,nstatv,DDSDDE)
      else
c-----------------------
         call pert_tang(DDe,stress_o,stress,dstran,dstran_dt,
     &maxT,M,M_d,dtime,expo,perturb,ntens,DDSDDE)
c-----------------------
c
      endif
c
      endif 
c-----------------------------------------------------------------------------------------------
      e_sand = e_sand-(1+e_sand)*((dstran(1)+dstran(2)+dstran(3)))
c
      statev(7) = e_sand
c-----------------------------------------------------------------------------------------------
      do i=1,ntens
         stress(i) = -stress(i)
      end do
      do i=1,ntens
         dstran(i) = -dstran(i)
      end do
c-----------------------------------------------------------------------------------------------
      return
      end
c------------------------------------------------------------------------------------------------
      subroutine stress_rate(DDe,maxT,stress,dstran_dt,expo,
     &M,M_d,dtime,ntens,dT,DDsDDe_temp)
c------------------------------------------------------------------------------------------------
c
      implicit none
c
      integer i,j,ntens
c
      double precision stress(ntens),dev_sig(ntens),sig_trace,sig_p
      double precision sig_e,norm_sig(ntens),maxT,dtime
      double precision norm_T(ntens),DDe(ntens,ntens)
      double precision dstran_dt(ntens),part_1(ntens),part_2(ntens)
      double precision part_22(ntens),dT(ntens)
      double precision expo,nonlinear,norm2,norm_d
      double precision nonlinear2
      double precision dot_vect
      double precision zero,one,two,three
      double precision M,M_d
      double precision tiny,factor_nonlinear
      double precision eta(6),NN(6),BB(6,6),DDsDDe_temp(6,6)
      double precision f_y
c
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter (tiny=1.0d-6)
c------------------------------------------------------------------------------------------------
      call deviator(stress,dev_sig,sig_trace,sig_p)
c------------------------------------------------------------------------------------------------
      call pzero(part_1,ntens)
      call pzero(part_2,ntens)
      call pzero(part_22,ntens)
c-----------------------------------------------------------------------------------------------
      call f_yd(dev_sig,ntens,maxT,M_d,sig_e,norm_sig,f_y)
c-----------------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------------
      nonlinear2 = dabs(sig_e/maxT)**expo
      nonlinear  = nonlinear2
c-----------------------------------------------------------------------------------------------
      norm2=dot_vect(6,dstran_dt,dstran_dt,6) ! 2013-10-18
      norm_d=dsqrt(norm2)
      do i=1,ntens
            part_22(i)=norm_d*nonlinear*norm_sig(i)
      enddo
c-----------------------------------------------------------------------------------------------
      call matmul_ZL(DDe,dstran_dt,part_1,6,6,1)
      call matmul_ZL(DDe,part_22,part_2,6,6,1)  
c-----------------------------------------------------------------------------------------------
      call pzero(dT,ntens)
c
      do i=1,6
         dT(i)=part_1(i)-part_2(i)
      enddo
c------------------------------------------------------------------------------------------------
      do i = 1,ntens
             stress(i) = stress(i)+dT(i)*dtime
      enddo
c------------------------------------------------------------------------------------------------
      if (norm_d.gt.zero) then
        do i = 1,ntens
             eta(i) = dstran_dt(i)/norm_d
        enddo
      else
        do i = 1,ntens
             eta(i) = zero
        enddo
      endif
c
      call matmul_ZL(DDe,norm_sig,NN,6,6,1)
      do i = 1,ntens
          NN(i) = -nonlinear*NN(i)
      enddo
c
      call cross_vect(NN,eta,6,BB)
c
      do i = 1,ntens
         do j = 1,ntens
          DDsDDe_temp(i,j) = DDe(i,j)+BB(i,j)
        enddo
      enddo
c------------------------------------------------------------------------------------------------
      return
      end
c------------------------------------------------------------------------------------------------
      subroutine RKF_stress_rate(DDe,maxT,stress,dstran_dt,expo,
     &M,M_d,dtime,ntens,dT)
c------------------------------------------------------------------------------------------------
c
      implicit none
c
      integer i,j,ntens
c
      double precision DDe(ntens,ntens),stress(ntens),DDsDDe(6,6)
      double precision maxT,dtime,dstran_dt(ntens)
      double precision dT(ntens),expo
      double precision stress_oo(ntens),stress_2(ntens),stress_3(ntens)
      double precision stress_4(ntens),stress_hat(ntens)
      double precision kRK_1(ntens),kRK_2(ntens),kRK_3(ntens)
      double precision kRK_4(ntens)
      double precision zero,one,two,three,four,half,six
      double precision one6,one3,two3
      double precision dot_vect
      double precision M,M_d
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(four=4.0d0,six=6.0d0,half=0.5d0)
c------------------------------------------------------------------------------------------------
      one6=one/six
      one3=one/three
      two3=two/three   
c------------------------------------------------------------------------------------------------
c    to be continue ...
c
      call push(stress,stress_oo,ntens)
      call stress_rate(DDe,maxT,stress_oo,dstran_dt,expo,
     &M,M_d,dtime,ntens,kRK_1,DDsDDe)
      call pzero(stress_2,ntens)
      do i=1,ntens
         stress_2(i)=stress(i)+half*kRK_1(i)*dtime
      enddo
c
      call stress_rate(DDe,maxT,stress_2,dstran_dt,expo,
     &M,M_d,dtime,ntens,kRK_2,DDsDDe)  
      call pzero(stress_3,ntens)
      do i=1,ntens
         stress_3(i)=stress(i)+half*kRK_2(i)*dtime
      enddo   
c 
      call stress_rate(DDe,maxT,stress_3,dstran_dt,expo,
     &M,M_d,dtime,ntens,kRK_3,DDsDDe)        
      call pzero(stress_4,ntens)
      do i=1,ntens
         stress_4(i)=stress(i)+kRK_3(i)*dtime
      enddo   
c   
      call stress_rate(DDe,maxT,stress_4,dstran_dt,expo,
     &M,M_d,dtime,ntens,kRK_4,DDsDDe)   
c              
      do i=1,ntens      
      stress_hat(i)=stress(i)+one6*kRK_1(i)*dtime+
     &one3*kRK_2(i)*dtime+one3*kRK_3(i)*dtime+
     &one6*kRK_4(i)*dtime
      enddo
c------------------------------------------------------------------------------------------------
      call push(stress_hat,stress,ntens)      
c------------------------------------------------------------------------------------------------
      return
      end
c------------------------------------------------------------------------------------------------
c------------------------------------------------------------------------------------------------
      subroutine el_stiff(parms,nparms,stress,ntens,statev,nstatv,
     &DDe)
c-----------------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------------
      implicit none
c
      integer i,j,nparms,nstatv,ntens
      double precision parms(nparms),stress(ntens),statev(nstatv)
      double precision Gt,nu,Et,LAMBDA,MU,p_atm,e_sand
      double precision zero,one,two,three
      double precision DDe(6,6),G0,p_mean
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
c------------------------------------------------------------------------------------------------
      G0   = parms(1)
      nu   = parms(2)
      p_atm    = parms(8)
      e_sand    = statev(7)
c
      p_mean    = dabs((one/three)*(stress(1)+stress(2)+stress(3)))
      if (p_mean .lt. 0.005*p_atm) then
         p_mean = 0.005*p_atm
      endif
c------------------------------------------------------------------------------------------------
c
      Gt = G0*p_atm*(((2.97-e_sand)**2.0d0)/(1+e_sand))*
     &((p_mean/p_atm)**(one/two))
c
      statev(6) = Gt
c
      do i=1,6
           do j=1,6
               DDe(i,j)=zero
           enddo
      enddo
c------------------------------------------------------------------------------------------------
      Et=Gt*(two*(one+nu))
c
      LAMBDA=Et*nu/((one+nu)*(one-two*nu))
      MU=Et/(2*(1+nu));
c
      do i=1,3
            do j=1,3
                 DDe(j,i)=LAMBDA
            enddo
            DDe(i,i)=LAMBDA+two*MU
      enddo
      do i=4,6
        DDe(i,i)=MU;
      enddo
c
      return
      end
c------------------------------------------------------------------------------------------------
      subroutine f_yd(dev_sig,ntens,maxT,M_d,sig_e,norm_sig,f_y)
c------------------------------------------------------------------------------------------------
      implicit none
      integer i, j, ntens
c
      double precision dev_sig_tensor(3,3),dev_sig(ntens),sig_e
      double precision norm_sig(ntens)
      double precision zero,one,two,three
      double precision dot_vect
      double precision norm_temp
      double precision tiny
      double precision M_d
      double precision maxT,f_y
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
      parameter(tiny=1.0d-12)
c
      dev_sig_tensor(1,1)=dev_sig(1)
      dev_sig_tensor(2,2)=dev_sig(2)
      dev_sig_tensor(3,3)=dev_sig(3)
      dev_sig_tensor(1,2)=dev_sig(4)
      dev_sig_tensor(1,3)=dev_sig(5)
      dev_sig_tensor(2,1)=dev_sig(4)
      dev_sig_tensor(2,3)=dev_sig(6)
      dev_sig_tensor(3,1)=dev_sig(5)
      dev_sig_tensor(3,2)=dev_sig(6)
c
      sig_e=0.0d0
      do i=1,3
           do j=1,3
                sig_e=sig_e+dev_sig_tensor(i,j)**2.0
           enddo
      enddo
      sig_e=((three/two)*sig_e)**(one/two)
c
      f_y = sig_e - maxT
c
      call pzero(norm_sig,ntens)
c
      if (sig_e.ne.zero) then 
                norm_sig(1)=(three/two)*dev_sig(1)/sig_e-M_d/3
                norm_sig(2)=(three/two)*dev_sig(2)/sig_e-M_d/3
                norm_sig(3)=(three/two)*dev_sig(3)/sig_e-M_d/3
c
                norm_sig(4)=(three)*dev_sig(4)/sig_e
                norm_sig(5)=(three)*dev_sig(5)/sig_e
                norm_sig(6)=(three)*dev_sig(6)/sig_e
          else
                norm_sig(1)=0
                norm_sig(2)=0
                norm_sig(3)=0
                norm_sig(4)=0
                norm_sig(5)=0
                norm_sig(6)=0
      endif
c
        norm_temp=dot_vect(6,norm_sig,norm_sig,ntens)
        norm_temp=dsqrt(norm_temp)
c
      if (norm_temp.le.zero) then
                norm_sig(1)=0
                norm_sig(2)=0
                norm_sig(3)=0
                norm_sig(4)=0
                norm_sig(5)=0
                norm_sig(6)=0
          else
                norm_sig(1)=norm_sig(1)/norm_temp
                norm_sig(2)=norm_sig(2)/norm_temp
                norm_sig(3)=norm_sig(3)/norm_temp
                norm_sig(4)=norm_sig(4)/norm_temp
                norm_sig(5)=norm_sig(5)/norm_temp
                norm_sig(6)=norm_sig(6)/norm_temp
       endif
c-------------------------------------------------------------------------------    
      return
      end
c-------------------------------------------------------------------------------
      subroutine pzero(v,nn)
c-------------------------------------------------------------------------------
      implicit  none

      integer n,nn
      double precision v(nn)
      save
c
      do n = 1,nn
        v(n) = 0.0d0
      enddo ! n
c
      return
      end
c------------------------------------------------------------------------------
      subroutine matmul_ZL(a,b,c,l,m,n)
c------------------------------------------------------------------------------
      implicit none
c
      integer i,j,k,l,m,n
c
      double precision a(l,m),b(m,n),c(l,n)
c
      do i=1,l
        do j=1,n
          c(i,j) = 0.0d0
          do k=1,m
            c(i,j) = c(i,j) + a(i,k)*b(k,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c------------------------------------------------------------------------------
      double precision function dot_vect(flag,a,b,n)
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      implicit none
      integer i,n,flag
      double precision a(n),b(n)
      double precision zero,half,one,two,coeff
c
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
c
      if(flag.eq.1) then
c
        coeff=two
c
      elseif(flag.eq.2) then
c
        coeff=half
c
      else
c
        coeff=one
c
      end if
c
      dot_vect=zero
c
      do i=1,n
        if(i.le.3) then
          dot_vect = dot_vect+a(i)*b(i)
        else
          dot_vect = dot_vect+coeff*a(i)*b(i)
        end if
      end do
c
      return
      end
c
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
      subroutine deviator(t,s,trace,mean)
c------------------------------------------------------------------------------
c------------------------------------------------------------------------------
c
      implicit none
c
      double precision t(6),s(6),trace,mean
      double precision one,three,onethird
c
      parameter (one=1.0d0,three=3.0d0)
c
      onethird=one/three
c
      trace=t(1)+t(2)+t(3)
      mean=onethird*trace
c
      s(1)=t(1)-mean
      s(2)=t(2)-mean
      s(3)=t(3)-mean
      s(4)=t(4)
      s(5)=t(5)
      s(6)=t(6)
c
      return
      end
c-----------------------------------------------------------------------------
      subroutine push(a,b,n)
c-----------------------------------------------------------------------------
      implicit none
      integer i,n
      double precision a(n),b(n)
c
      do i=1,n
        b(i)=a(i)
      enddo
c
      return
      end
c ----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
      subroutine pert_tang(DDe,stress_o,stress,dstran,dstran_dt,
     &maxT,M,M_d,dtime,expo,theta,ntens,DDSDDE)
c-----------------------------------------------------------------------------
c-----------------------------------------------------------------------------
      implicit none
c
      integer ntens,ii,jj,kk
      integer i,j
c
      double precision stress_o(ntens),stress_star(ntens)
      double precision stress(ntens)
      double precision maxT,expo
      double precision theta,dtime
      double precision dstran(6),dstran_dt(6),dstran_star_dt(6)
      double precision dev_sig_star,sig_trace_star,sig_p_star
      double precision dsig(6),dsig_star(6)
      double precision DD(6,6),DDe(6,6),DDSDDE(6,6)
      double precision zero,three
      double precision M,M_d
c
      parameter(zero=0.0d0,three=3.0d0)
c
      do i=1,ntens
           do j=1,ntens
                DD(i,j)=zero
           enddo
      enddo
      do i=1,ntens
           do j=1,ntens
                DDSDDE(i,j)=zero
           enddo
      enddo
      call pzero(stress_star,ntens)
c
	do jj=1,ntens
        call push(stress_o,stress_star,ntens)
        call push(dstran_dt,dstran_star_dt,ntens)
               dstran_star_dt(jj)=dstran_star_dt(jj)+theta/dtime
c
      call stress_rate(DDe,maxT,stress_star,dstran_star_dt,expo,
     &M,M_d,dtime,ntens,dsig_star,DDsDDe)     
c
          do kk=1,ntens
            dsig(kk)=stress_star(kk)-stress(kk)
            DD(kk,jj)=dsig(kk)/theta
          enddo !kk
c
        enddo !jj  
c----------------------------------------------------------------------------------
      do j=1,ntens
        do i=1,ntens
		if((i.le.3).and.(j.le.3)) then
          ddsdde(i,j) = DD(i,j)
		else
          ddsdde(i,j) = DD(i,j)
		end if
        end do
      enddo
c----------------------------------------------------------------------------------
      return
      end
c----------------------------------------------------------------------------------
c----------------------------------------------------------------------------------
      subroutine xit()
c      stop
      return
      end
c----------------------------------------------------------------------------------
      subroutine matrix_vect(matrix,row,col,vect,vect_out)
c------------------------------------------------------------------------------
      implicit none
      integer i,j,row,col
      double precision matrix(row,col),vect(col),vect_out(row)
      double precision zero,half,one,two,coeff
c
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
c
      call pzero(vect_out,row)
c
      do i=1,row
           do j=1,col
                vect_out(i)=vect_out(i)+matrix(i,j)*vect(j) 
           enddo
      enddo
c
      return
      end
c----------------------------------------------------------------------------------
      subroutine get_state(parms,nparms,stress,ntens,e_sand,e_ini,
     &coeff_CSL)
c------------------------------------------------------------------------------
      implicit none
      integer i,j,ntens,nparms
      double precision parms(nparms),stress(ntens)
      double precision G_elas,nu_elas,expo,cohesion,fric_ang
      double precision dila_ang,p_atm,e0,lambda_c,xi_c
      double precision p_mean,e_c,e_ini
      double precision big_psi_ini,big_psi,e_sand,coeff_CSL
      double precision zero,one,two,three
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0)
c
      G_elas   = parms(1)
      nu_elas  = parms(2)
      expo     = parms(3)
      cohesion = parms(4)
      fric_ang = parms(5)
      dila_ang = parms(6)
      e_ini    = parms(7)
      p_atm    = parms(8)
      e0       = parms(9)
      lambda_c = parms(10)
      xi_c     = parms(11)
c
      p_mean = (stress(1)+stress(2)+stress(3))/three
c
      e_c      = e0-lambda_c*(p_mean/p_atm)**xi_c
c
      big_psi_ini = dabs(e_c-e_ini)
      big_psi     = e_c-e_sand
c
      coeff_CSL   = 1.0d0*big_psi/big_psi_ini
      coeff_CSL   = dsign((dabs(coeff_CSL))**1.0d0,coeff_CSL)
c
      return
      end
c-------------------------------------------------------------------------------
      subroutine cross_vect(a,b,n,c)
c
      implicit none
      integer i,j,n
      double precision a(n),b(n),c(n,n)
      double precision zero,half,one,two,coeff
      parameter(zero = 0.0d0)
c
      do i=1,n
        do j=1,n
            c(i,j) = zero
        enddo
      enddo
c
      do i=1,n
        do j=1,n
            c(i,j) = a(i)*b(j)
        enddo
      enddo
c
      return
      end
c--------------------------------------------------------------------------------
      subroutine lode_DM(sig,theta,cos3t)
c------------------------------------------------------------------------------
c
      implicit none
c
      double precision sig(6),sig_trace,sig_p,sig_e
      double precision r(6),r2(6)
      double precision trr2,trr3,J2bar,J3bar,J2bar_sq
      double precision cM,n_VE,n_VEm1,numer,denom,cos3t
c
      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
      double precision alpha,beta,gth,dgdth
      double precision one,two,three
      double precision onethird,half,sqrt3,tiny
      double precision theta
c
      data one,two,three/1.0d0,2.0d0,3.0d0/
      data tiny,n_VE/1.0d-15,-0.25d0/
c
      onethird=one/three
      half=one/two
      sqrt3=dsqrt(three)
      call deviator(sig,r,sig_trace,sig_p)
c
      call dot_vect_ZL(1,r,r,6,trr2)
      J2bar=half*trr2
c
      r2(1)=r(1)*r(1)+r(4)*r(4)+r(5)*r(5)
      r2(2)=r(4)*r(4)+r(2)*r(2)+r(6)*r(6)
      r2(3)=r(6)*r(6)+r(5)*r(5)+r(3)*r(3)
      r2(4)=r(1)*r(4)+r(4)*r(2)+r(6)*r(5)
      r2(5)=r(5)*r(1)+r(6)*r(4)+r(3)*r(5)
      r2(6)=r(4)*r(5)+r(2)*r(6)+r(6)*r(3)
c
      if(trr2.lt.tiny) then
c
        cos3t=one
c
      else
c
        call dot_vect_ZL(1,r,r2,6,trr3)
c
        J3bar=onethird*trr3
        J2bar_sq=dsqrt(J2bar)
        numer=three*sqrt3*J3bar
        denom=two*(J2bar_sq**3)
        cos3t=numer/denom
        if(dabs(cos3t).gt.one) then
          cos3t=cos3t/dabs(cos3t)
        endif
c
      endif
c
      theta = (one/three)*acos(cos3t)
c
      return
      end
c------------------------------------------------------------------------------
      subroutine dot_vect_ZL(flag,a,b,n,dot_vect)
c------------------------------------------------------------------------------
      implicit none
      integer i,n,flag
c
      double precision dot_vect
      double precision a(n),b(n)
      double precision zero,half,one,two,coeff
c
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
c
      if(flag.eq.1) then
c
        coeff=two
c
      elseif(flag.eq.2) then
c
        coeff=half
c
      else
c
        coeff=one
c
      end if
c
      dot_vect=zero
c
      do i=1,n
        if(i.le.3) then
          dot_vect = dot_vect+a(i)*b(i)
        else
          dot_vect = dot_vect+coeff*a(i)*b(i)
        end if
      end do
c
      return
      end
c-------------------end of UMAT subroutines--------------------------