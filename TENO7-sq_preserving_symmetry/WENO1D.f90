module flow_array
    implicit none
    real(kind=8),allocatable,dimension(:) :: d,p,u,a
    real(kind=8),allocatable,dimension(:,:) :: u0
end
module constant
    real(kind=8),parameter :: pi=3.14159265358979323846264338327950288419716939937510_8
	real(kind=8),parameter :: gama=1.40_8
	integer,parameter :: nx=500
	real(kind=8),parameter :: hx=1.0_8/(nx-1)
	real(kind=8),parameter :: cfl=0.3_8
    integer,parameter :: mt=4 !2mt-1 eq order
    real(kind=8),parameter :: end_t=0.15_8
end
module coefficient
	real(kind=8),parameter :: c1=1.0_8/4.0_8,c2=3.0_8/4.0_8,c3=2.0_8/3.0_8,c4=1.0_8/3.0_8
	real(kind=8),parameter :: s1=-1.0/6.0_8,s2=5.0/6.0_8,s3=1.0/3.0_8,s4=-7.0/6.0_8,s5=11.0/6.0_8,&
		s6=1.0/4.0_8,s7=13.0/12.0_8,s8=-5.0/12.0_8,s9=1.0/12.0_8,s10=-1.0/4.0_8,&
		s11=-23.0/12.0_8,s12=25.0/12.0_8,s13=-11.0/6.0_8,s14=-3.0/2.0_8,s15=13.0/3.0_8,&
		s16=-5.0/2.0_8,s17=-1.0/2.0_8,s18=781.0/20.0_8,s19=1.0/2.0_8,s20=1.0/6.0_8,&
		s21=-1.0/3.0_8,s22=3.0/2.0_8,s23=-1.0/60.0_8,s24=3.0/20.0_8,s25=-3.0/4.0_8,&
		s26=3.0/4.0_8,s27=-3.0/20.0_8,s28=1.0/60.0_8,s29=1.0/90.0_8,s30=-49.0/18.0_8,&
		s31=1.0/8.0_8,s32=13.0/8.0_8,s33=-13.0/8.0_8,s34=-1.0/8.0_8,s35=-13.0/2.0_8,&
		s36=28.0/3.0_8,s37=5.0/2.0_8,s38=1680790021.0/1549497600.0_8,s39=7168762325689.0/6608785075200.0_8,s40=61.0/720.0_8,&
		s41=-43.0/2562.0_8,s42=949.0/11200.0_8,s43=-397.0/23652.0_8,s44=1.0/2520.0_8,s45=18.0/35.0_8,&
		s46=9.0/35.0_8,s47=3.0/35.0_8,s48=4.0/35.0_8,s49=1.0/35.0_8
end

program main
	use constant
    use flow_array
    real(kind=8) :: t1(2)
    call allocate_vars
    call cpu_time(t1(1))   
    call initialization
	
	open(unit=100,file='sym_err.dat')
	write(100,*) 'VARIABLES = t,L1,L2,Linf'
	
    call time_advance
    call cpu_time(t1(2))
    open(unit=11,file='time.dat')
    write(11,*) t1(2)-t1(1)
    call postprocess
	call deallocate_vars
end

subroutine allocate_vars
	use constant
    use flow_array
	allocate(d(1-mt:nx+mt),p(1-mt:nx+mt),u(1-mt:nx+mt),a(1-mt:nx+mt))
	allocate(u0(3,1-mt:nx+mt))
end

subroutine deallocate_vars
    use flow_array
	deallocate(d,p,u,a,u0)
end

subroutine initialization
	use constant
    use flow_array
    implicit none
	integer :: i
	do i=1,nx
		if ((i-1.d0)*hx<0.5d0) then
			d(i)=1.d0
			u(i)=-2.d0
			p(i)=0.4d0
        else
			d(i)=1.d0
			u(i)=2.d0
			p(i)=0.4d0
		end if
		a(i)=sqrt(gama*p(i)/d(i))
		u0(1,i)=d(i)
		u0(2,i)=d(i)*u(i)
		u0(3,i)=(p(i)/(gama-1.d0)+u(i)**2*0.5_8*d(i))
	end do
end
subroutine time_advance
    use constant
    implicit none
    integer :: i=0
    real(kind=8) :: dt,t=0.0_8,sym_err(3)

	call sym_out(t,sym_err)
    do 
        i=i+1
		call cal_dt(dt)
		if (t+dt>end_t) then
			dt=end_t-t
			t=end_t
            call rk3(dt)
			call sym_out(t,sym_err)
            exit       
		end if
        ! write(*,'(A,I5.5,A,F6.4)') 'Iter=',i,' Time=',t
        write(*,*) 'Iter=',i,' Time=',t,' Err=',sym_err(2)
        call rk3(dt)
        t=t+dt
		call sym_out(t,sym_err)
		
	end do
end

subroutine sym_out(t,sym_err)
	use constant
	use flow_array
    implicit none
    integer :: i
	real(kind=8) :: t,sym_err(3),temp_err(2)
	
	sym_err=0.d0
	do i=1,nx
		temp_err(1)=abs(d(i)-d(nx+1-i))
		temp_err(2)=(d(i)-d(nx+1-i))**2
		sym_err(1)=sym_err(1)+temp_err(1)
		sym_err(2)=sym_err(2)+temp_err(2)
		sym_err(3)=max(sym_err(3),temp_err(1))
	enddo
	sym_err(1)=sym_err(1)/(2.d0*nx+0.d0)
	sym_err(2)=sqrt(sym_err(2)/(2.d0*nx+0.d0))
	write(100,*) t,sym_err(1),sym_err(2),sym_err(3)
	
end

subroutine cal_dt(dt)
	use constant
    use flow_array
    implicit none
    real(kind=8) :: dt
    
    dt=cfl/(maxval(abs(u(1:nx))+abs(a(1:nx))))*hx
end

subroutine rk3(dt)
    use constant
    use flow_array
	use coefficient
    implicit none
    integer :: i
    real(kind=8),allocatable :: u_last(:,:)
	real(kind=8),allocatable :: f_hp(:,:),flx_hp(:,:)
    real(kind=8) :: dt
    allocate(u_last(3,1:nx))
	allocate( f_hp(3,0:nx),flx_hp(3,0:nx) )

    u_last(:,1:nx)=u0(:,1:nx)
    do i=1,3
        call boundarycondition
        call else_to_u0
		call restruct(f_hp,flx_hp,dt)
		
        if (i==1) then
			call pp(f_hp,flx_hp,u_last,dt)
        elseif (i==2) then
			call pp(c1*f_hp,c1*flx_hp,c2*u_last(:,1:nx)+c1*u0(:,1:nx),dt)
        elseif (i==3) then
			call pp(c3*f_hp,c3*flx_hp,c4*u_last(:,1:nx)+c3*u0(:,1:nx),dt)
        endif
        call u0_to_else
    enddo
	
	deallocate(u_last,f_hp,flx_hp)
	
end

subroutine pp(f_hp,flx_hp,u_last,dt)
    use constant
    use flow_array
    implicit none
    integer :: i,j
    real(kind=8) :: dt,ga_jm(nx),lamta(2),de(2),temp_p,temp_u0(3)
    real(kind=8) :: u_last(3,1:nx),flx_hp(3,0:nx),f_hp(3,0:nx)
    real(kind=8) :: eps_d,eps_p,eps_de
    eps_d=1.0E-13_8
    eps_p=1.0E-13_8
    do i=1,nx
		lamta(:)=1.0_8
        do j=1,1000
            temp_u0(:)=u_last(:,i)-dt/hx*( (lamta(2)*(f_hp(:,i)-flx_hp(:,i))+flx_hp(:,i)) - &
            & ( lamta(1)*(f_hp(:,i-1)-flx_hp(:,i-1))+flx_hp(:,i-1) ) )
            temp_p=(gama-1.0_8)*(temp_u0(3)-temp_u0(2)**2/temp_u0(1)*0.5_8)
            if (temp_u0(1)>=eps_d .and. temp_p>=eps_p) then
				u0(:,i)=temp_u0(:)
				exit
			else
				lamta(:)=0.90_8*lamta(:)
            endif
        enddo
    enddo
end

subroutine boundarycondition
    use constant
    use flow_array
    implicit none
    integer :: i
    do i=1,mt
        d(1-i)=1.d0
        p(1-i)=0.4d0
        u(1-i)=-2.d0
        d(nx+i)=1.d0
        p(nx+i)=0.4d0
        u(nx+i)=2.d0
    end do
end

subroutine u0_to_else
    use constant
    use flow_array
	integer i
	do i=1-mt,nx+mt
		d(i)=u0(1,i)
		u(i)=u0(2,i)/d(i)
		p(i)=(u0(3,i)/d(i)-u(i)*u(i)*0.5_8)*d(i)*(gama-1)
		a(i)=sqrt(gama*p(i)/d(i))
	end do
end

subroutine else_to_u0
    use constant
    use flow_array
	integer i
	do i=1-mt,nx+mt
        a(i)=sqrt(gama*p(i)/d(i))
        u0(1,i)=d(i)
        u0(2,i)=d(i)*u(i)
        u0(3,i)=d(i)*(p(i)/d(i)/(gama-1)+u(i)**2*0.5_8)
	end do
end

subroutine cal_f(f)
    use constant
    use flow_array
    implicit none
    integer :: i
    real(kind=8) :: f(3,1-mt:nx+mt)
    do i=1-mt,nx+mt
        f(1,i)=d(i)*u(i)
        f(2,i)=d(i)*u(i)*u(i)+p(i)
        f(3,i)=(p(i)/(gama-1.0_8)+d(i)*u(i)*u(i)*0.5_8+p(i))*u(i)
    end do 
end  

subroutine restruct(f_hp,flx_hp,dt)
    use constant
    use flow_array
    implicit none
    integer :: i,j,k,m
    real(kind=8) :: f_p(3,2*mt-1),f_n(3,2*mt-1),g_p(3,2*mt-1),g_n(3,2*mt-1)
    real(kind=8) :: dt,g_p_r1(3),g_n_r1(3),g_p_r(3),g_n_r(3),alpha,q(3,3),iq(3,3),g(3),g1(3)
    real(kind=8) :: f_hp(3,0:nx),flx_hp(3,0:nx)
    real(kind=8),allocatable :: f(:,:)
	real(kind=8),allocatable :: mid_u(:),mid_d(:),mid_a(:)
    
    allocate(mid_u(1-mt:nx+mt),mid_d(1-mt:nx+mt),mid_a(1-mt:nx+mt))
    allocate(f(3,1-mt:nx+mt))

    call roe(mid_u,mid_d,mid_a)

    call cal_f(f)
	
    do i=0,nx
        call cal_q(mid_u(i),mid_d(i),mid_a(i),q,iq)

        alpha=maxval(abs(u(1+i-mt:i+mt))+abs(a(i-mt+1:i+mt)))

        do j=i-mt+1,i+mt-1
            m=j-i+mt
            f_p(:,m)=(f(:,j)+alpha*u0(:,j))*0.5_8
            g_p(:,m)=matmul(iq,f_p(:,m)) 
        enddo

        do j=i+1-mt+1,i+1+mt-1
            m=j-i+mt-1
            f_n(:,m)=(f(:,j)-alpha*u0(:,j))*0.5_8
            g_n(:,m)=matmul(iq,f_n(:,m)) 
        enddo

        do j=1,3
            call weno_p(g_p(j,1:2*mt-1),g_p_r(j),g_p_r1(j))
            call weno_n(g_n(j,1:2*mt-1),g_n_r(j),g_n_r1(j))
        enddo

        ! f_hp(:,i)=matmul(q,(g_p_r+g_n_r))
        ! flx_hp(:,i)=matmul(q,(g_p_r1+g_n_r1))
		g=g_p_r+g_n_r
		f_hp(1,i)=q(1,1)*g(1)+(q(1,2)*g(2)+q(1,3)*g(3))
		f_hp(2,i)=q(2,1)*g(1)+(q(2,2)*g(2)+q(2,3)*g(3))
		f_hp(3,i)=q(3,1)*g(1)+(q(3,2)*g(2)+q(3,3)*g(3))
		g1=g_p_r+g_n_r
		flx_hp(1,i)=q(1,1)*g1(1)+(q(1,2)*g1(2)+q(1,3)*g1(3))
		flx_hp(2,i)=q(2,1)*g1(1)+(q(2,2)*g1(2)+q(2,3)*g1(3))
		flx_hp(3,i)=q(3,1)*g1(1)+(q(3,2)*g1(2)+q(3,3)*g1(3))
		
    enddo
	
	deallocate(f,mid_u,mid_d,mid_a)

end

subroutine roe(mid_u,mid_d,mid_a)
	use constant
    use flow_array
    implicit none
	integer :: i
    real(kind=8) :: mid_u(1-mt:nx+mt),mid_d(1-mt:nx+mt),mid_a(1-mt:nx+mt)
	real(kind=8) :: sq_d_i,sq_d_ip,sq_d_sum
    real(kind=8),allocatable :: h(:),mid_h(:)

    allocate(h(1-mt:nx+mt),mid_h(1-mt:nx+mt))
    
	do i=1-mt,nx+mt
		h(i)=gama*p(i)/d(i)/(gama-1.0_8)+u(i)**2*0.5_8
	end do
	do i=1-mt,nx+mt-1
		sq_d_i=sqrt(d(i))
		sq_d_ip=sqrt(d(i+1))
		sq_d_sum=sq_d_i+sq_d_ip
		mid_d(i)=(d(i)+d(i+1)+2*sq_d_i*sq_d_ip)*0.25_8
		mid_u(i)=(u(i)*sq_d_i+u(i+1)*sq_d_ip)/sq_d_sum
		mid_h(i)=(h(i)*sq_d_i+h(i+1)*sq_d_ip)/sq_d_sum
		mid_a(i)=sqrt((gama-1)*(mid_h(i)-mid_u(i)**2*0.5_8))
	end do
	mid_d(nx+mt)=mid_d(nx+mt-1)
	mid_u(nx+mt)=mid_u(nx+mt-1)
	mid_a(nx+mt)=mid_a(nx+mt-1)
	
	deallocate(mid_h,h)
	
end	

subroutine cal_q(mid_u,mid_d,mid_a,q,iq)
    use constant
    implicit none
    real(kind=8) :: q(3,3),iq(3,3)
    real(kind=8) :: mid_u,mid_d,mid_a
	iq(1,1)=mid_u**2/2.0_8-mid_a**2/(gama-1.0_8) ; iq(1,2)=-mid_u ; iq(1,3)=1.0_8 ;
	iq(2,1)=-mid_u-(gama-1.0_8)/mid_a*mid_u**2/2.0_8 ; iq(2,2)=1.0_8+(gama-1.0_8)/mid_a*mid_u ; iq(2,3)=-(gama-1.0_8)/mid_a;
	iq(3,1)=-mid_u+(gama-1.0_8)/mid_a*mid_u**2/2.0_8 ; iq(3,2)=1.0_8-(gama-1.0_8)/mid_a*mid_u ; iq(3,3)=(gama-1.0_8)/mid_a;

	q(1,1)=-(gama-1.0_8)/mid_a/mid_a ; q(1,2)=-1.0_8/2.0_8/mid_a ; q(1,3)=1.0_8/2.0_8/mid_a ;
	q(2,1)=-(gama-1.0_8)/mid_a/mid_a*mid_u ; q(2,2)=-(mid_u-mid_a)/2.0_8/mid_a ; q(2,3)=(mid_u+mid_a)/2.0_8/mid_a ; 
	q(3,1)=-(gama-1.0_8)/mid_a/mid_a*mid_u**2/2.0_8 ; q(3,2)=-(mid_u**2/2.0_8+mid_a**2/(gama-1.0_8)-mid_u*mid_a)/2.0_8/mid_a ;
													   q(3,3)=(mid_u**2/2.0_8+mid_a**2/(gama-1.0_8)+mid_u*mid_a)/2.0_8/mid_a ;
end

subroutine weno(v,x) 
    use constant
	use coefficient
    implicit none
	integer :: i
    real(kind=8) :: v(2*mt-1),b(2*mt-3),vr(2*mt-3),a(2*mt-3),a1(2*mt-3),x,tao,gammak(2*mt-3),gammak1(2*mt-3),delta(2*mt-3)
    real(kind=8) :: ct,vs(2*mt-2),gamma_s,a_s

	! real(kind=8),parameter :: s1=-1.0/6.0_8,s2=5.0/6.0_8,s3=1.0/3.0_8,s4=-7.0/6.0_8,s5=11.0/6.0_8,&
		! s6=1.0/4.0_8,s7=13.0/12.0_8,s8=-5.0/12.0_8,s9=1.0/12.0_8,s10=-1.0/4.0_8,&
		! s11=-23.0/12.0_8,s12=25.0/12.0_8,s13=-11.0/6.0_8,s14=-3.0/2.0_8,s15=13.0/3.0_8,&
		! s16=-5.0/2.0_8,s17=-1.0/2.0_8,s18=781.0/20.0_8,s19=1.0/2.0_8,s20=1.0/6.0_8,&
		! s21=-1.0/3.0_8,s22=3.0/2.0_8,s23=-1.0/60.0_8,s24=3.0/20.0_8,s25=-3.0/4.0_8,&
		! s26=3.0/4.0_8,s27=-3.0/20.0_8,s28=1.0/60.0_8,s29=1.0/90.0_8,s30=-49.0/18.0_8,&
		! s31=1.0/8.0_8,s32=13.0/8.0_8,s33=-13.0/8.0_8,s34=-1.0/8.0_8,s35=-13.0/2.0_8,&
		! s36=28.0/3.0_8,s37=5.0/2.0_8,s38=1680790021.0/1549497600.0_8,s39=7168762325689.0/6608785075200.0_8,s40=61.0/720.0_8,&
		! s41=-43.0/2562.0_8,s42=949.0/11200.0_8,s43=-397.0/23652.0_8,s44=1.0/2520.0_8,s45=18.0/35.0_8,&
		! s46=9.0/35.0_8,s47=3.0/35.0_8,s48=4.0/35.0_8,s49=1.0/35.0_8

	ct=1.0E-7_8

	vr(1)= s1*v(3) + s2*v(4) + s3*v(5)
	vr(2)= s3*v(4) + s2*v(5) + s1*v(6)
	vr(3)= s3*v(2) + s4*v(3) + s5*v(4)
	vr(4)= s6*v(4) + s7*v(5) + s8*v(6) + s9*v(7)
	vr(5)= s10*v(1)+ s7*v(2) + s11*v(3)+ s12*v(4)

	b(1)= s7*( v(3)-2*v(4)+v(5))**2 + s6*(v(3)-v(5))**2
	b(2)= s7*( v(4)-2*v(5)+v(6))**2 + s6*(3*v(4)-4*v(5)+v(6))**2
	b(3)= s7*( v(2)-2*v(3)+v(4))**2 + s6*(v(2)-4*v(3)+3*v(4))**2
	b(4)=         ( s13*v(4) +   3*v(5) + s14*v(6) + s3 *v(7))**2+&
		&    s15 *(     v(4) + s16*v(5) +   2*v(6) + s17*v(7))**2+&
		&    s18 *( s1 *v(4) + s19*v(5) + s17*v(6) + s20*v(7))**2
	b(5)=         ( s21*v(1) + s22*v(2)    -3*v(3) + s5 *v(4))**2+&
		&    s15 *( s17*v(1) +   2*v(2) + s16*v(3)      +v(4))**2+&
		& 	 s18 *( s1 *v(1) + s19*v(2) + s17*v(3) + s20*v(4))**2
	
	vs(1)= s23*v(1) + s24*v(2) + s25*v(3)            + s26*v(5) + s27*v(6) + s28*v(7)
	vs(2)= s29*v(1) + s27*v(2) + s22*v(3) + s30*v(4) + s22*v(5) + s27*v(6) + s29*v(7)
	vs(3)= s31*v(1) -     v(2) + s32*v(3)            + s33*v(5) +     v(6) + s34*v(7)
	vs(4)= s1 *v(1) +   2*v(2) + s35*v(3) + s36*v(4) + s35*v(5) +   2*v(6) + s1 *v(7)
	vs(5)= s17*v(1) +   2*v(2) + s16*v(3)            + s37*v(5) -   2*v(6) + s19*v(7)
	vs(6)=     v(1) -   6*v(2) +  15*v(3) -  20*v(4) +  15*v(5) -   6*v(6) +     v(7)

	tao= vs(1)**2 + vs(2)**2 + vs(3)**2 + vs(4)**2 + s38*vs(5)**2 + s39*vs(6)**2 + &
		& s40*( vs(3) + s41*vs(5) )**2 + s42*( vs(4) + s43*vs(6) )**2 + &
		&  s9*( vs(2) + s23*vs(4) + s44*vs(6) )**2
	
	tao=abs(tao + s1*( 4*b(1)+b(2)+b(3) ))

	gammak(1)=(1.0_8+tao/(1.0E-40_8+b(1)))**6
	gammak(2)=(1.0_8+tao/(1.0E-40_8+b(2)))**6
	gammak(3)=(1.0_8+tao/(1.0E-40_8+b(3)))**6
	gammak(4)=(1.0_8+tao/(1.0E-40_8+b(4)))**6
	gammak(5)=(1.0_8+tao/(1.0E-40_8+b(5)))**6
	gamma_s=sum(gammak)
	gammak1=gammak/gamma_s

	do i=1,2*mt-3
		if (gammak1(i)<ct) then
			delta(i)=0.0_8
		else
			delta(i)=1.0_8
		endif
	enddo    

	a(1)= s45*delta(1)
	a(2)= s46*delta(2)
	a(3)= s47*delta(3)
	a(4)= s48*delta(4)
	a(5)= s49*delta(5)
	a_s=sum(a)
	a1=a/a_s

	x=dot_product(vr,a1)

end

subroutine weno_p(v,x,x1) 
    use constant
    implicit none
    integer :: i
    real(kind=8) :: v(2*mt-1),x,x1
	call weno(v,x)
	x1=v(mt)
end

subroutine weno_n(v,x,x1) 
    use constant
    implicit none
    integer :: i
	real(kind=8) v_new(2*mt-1),v(2*mt-1),x,x1
	do i=1,7
		v_new(i)=v(8-i)
	enddo
    call weno_p(v_new,x,x1)
end

subroutine postprocess
    use constant
    use flow_array
    implicit none
    integer :: i
    real(kind=8),allocatable :: x(:)
    allocate(x(nx))
    do i=1,nx
        x(i)=(i-1)*hx
    end do
    open(unit=1,file='WENO_1D.dat')
    write(1,*) 'VARIABLES = x,d,p,u,a'
    do i=1,nx
        write(1,*) x(i),d(i),p(i),u(i),a(i)
    end do
	
	deallocate(x)
end
