	Program SwKeq
	implicit none

	integer i,l,j,k,cycleK,cycleC,ALLOC_ERR,mod_num
	character name*100,line*100,chara_mod*1,chara_si*50
	character chara_MW*4
	integer rn,ir,kk,rrr
	real*8 Ki,conct,conctM,c1,Swt,Rcorr_1,deltaRcorr
	real*8 Kim1,deltaK,ctm1,deltaC,Kmaxc,Kminc,deltaS,s1
	real*8 Kmin,Kmax,clow,chigh,MW,Rcorr,Kib,stand_err
	real*8, allocatable :: ct(:), Swe(:), Sw(:), Sw_1(:)



c	write(*,*)'Reaction numbers:'
c	read (*,*) rn
c	write(*,*)'Kmin [M^-1]:'
c	read (*,*) Kmin
c	write(*,*)'Kmax [M^-1]:'
c	read (*,*) Kmax
c	write(*,*) 'K_cycle_num:'
c	read (*,*) cycleK
c	write(*,*) 'Ct_cycle_num:'
c	read (*,*) cycleC

	!write(*,*) 'Ct [uM] vs. <S> [S] file'
	call getarg(1,chara_mod)
	call getarg(2,chara_si)
	call getarg(3,chara_MW)
	call getarg(4,name)
	if(name(:2).eq.'-H'.or.name(:2).eq.'  ')then
	 write(*,*)
	 write(*,*)'SwKeq calculates the equilibrium constant and the degr
     +ee of polymerization of a protein'
	 write(*,*)'using the concentration dependence of the sedimentatio
     +n coefficient to 20oC in water'
	 write(*,*)'and considering an isodesmic process.'
	 write(*,*)
	 write(*,*)
	 write(*,*)'please type: SwKeq.x models So Mw file_name_Ct[uM]_vs_
     +Sw20'
	 write(*,*)
	 write(*,*)
	 write(*,*)'models: In this option will choose the Hydrodinamic mo
     +del for protein polymers'
	 write(*,*)
	 write(*,*)'if type 1 the hydrodinamic model will be N^2/3'
	 write(*,*)'if type 2 the hydrodinamic model will be linear (spher
     +es)'
	 write(*,*)'if type 3 the hydrodinamic model will be linear (21 be
     +ads)'
	 write(*,*)'if type 4 the hydrodinamic model will be Helical (sphe
     +res)'
	 write(*,*)'if type 5 the hydrodinamic model will be Helical (21 b
     +eads)'
	 write(*,*)'if type 6 the hydrodinamic model will be Helical (42 b
     +eads)'
	 write(*,*)
	 write(*,*)'These models were tested by C.A. Sontag et al., Biophy
     +sical Chemistry 108 (2004) 215–230'
	 write(*,*)
	 write(*,*)'So is the Heterodimer Sedimentation Coefficient. For tu
     +bulin, authors recommend the following values:'
	 write(*,*)
	 write(*,*)'5.82 S obtained by J.G. de la Torre, J.M. Andreu, Hydr
     +odynamic analysis'
	 write(*,*)'of a,b-tubulin heterodimer and double rings, J.Mol. Bi 
     +ol. 238(1994) 223–225.'
	 write(*,*)
	 write(*,*)'Other value could be 5.36 S obtained by Alday and Corr
     +eia, Macromolecular' 
	 write(*,*)'Interaction of Halichondrin B Analogues Eribulin (E738
     +9) and ER-076349 with'
	 write(*,*)'Tubulin by Analytical Ultracentrifugation, Biochemistr
     +y 2009, 48, 7927–7938'
	 write(*,*)
	 write(*,*)'Mw is the molecular weight in kDa'
	 write(*,*)
	 write(*,*)'Version:2.101'
	 write(*,*)'Author: Javier Rodriguez-Salarichs'
	 write(*,*)
	stop
	endif

	 read(chara_si(:len_trim(chara_si)),*) s1

	open (1,file=name(:len_trim(name) ),status='old')

	read (chara_mod,'(i1)') mod_num
	if(mod_num.lt.1.or.mod_num.gt.6)stop('Only were considered 6 model
     +s')

	l=0
	ir=0

	do while (ir.eq.0)
	 read(1,*,iostat=ir)line
	if(ir.eq.0)l=l+1
	enddo

	rewind(1)

	allocate (ct(l),Swe(l),Sw(l),Sw_1(l),stat=ALLOC_ERR)

	do i=1,l
	 read(1,*)ct(i),Swe(i)
	enddo

	 read(chara_MW(:len_trim(chara_MW)),*) MW
	 MW=MW*1000.D0   !from kDa to Da
	 !MW=110000.D0   !g/mol Molecular weight of a,b-tubulin heterodimer

	do i=1,l
	 ct(i)=ct(i)*0.000001D0        ![M] [mol/L]
	enddo

	Kmin=1.D0
	Kmax=1000000000000.D0
	Kmaxc=Kmax
	Kminc=Kmin

	!do i=1,l

	rn=0
	Rcorr_1=-999999999999D0
	deltaRcorr=999999999999D0
	do while(deltaRcorr.gt.0.D0)
	rn=rn+1

	 Kmax=Kmaxc
	 Kmin=Kminc

c	 do k=1,cycleK
	 k=0
	 deltaK=9999999999.D0
	 do while(deltaK.gt.0.00000001)
	  k=k+1

	  Ki=(Kmax+Kmin)/2.D0

	  do i=1,l

	   clow=0.D0
	   chigh=ct(i)

	   j=0
	   deltaC=9999999999.D0
	   do while(deltaC.gt.0.00000001)
	    j=j+1

c	  do j=1,cycleC
	    c1=(clow+chigh)/2.D0

	    if (conctM(Ki,c1,rn).gt.ct(i))then    
	     chigh=(clow+chigh)/2.D0                          
	    elseif(conctM(Ki,c1,rn).lt.ct(i))then 
	     clow=(clow+chigh)/2.D0                           
	    endif

	    deltaC=abs((conctM(Ki,c1,rn)-ctm1)*1000000)
	   !write(*,*)deltaC,conctM(Ki,c1,rn),ctm1
	    ctm1=conctM(Ki,c1,rn)

	   enddo
	 ! write(*,*) 'ct[uM]=',conctM(Ki,c1,rn)*1000000
	 ! write(*,*) 'ci[uM]=',(((Ki**(kk-1))*(c1**kk)*1000000),kk=1,rn)
	 ! write(*,*) ''
	   Sw(i)=Swt(Ki,c1,rn,MW,mod_num,s1)
	  
	  enddo

	 ! write(*,*) (Swe(kk),Sw(kk),kk=1,l)
	 ! write(*,*)
	  !write(*,*)j,deltaC

	  ! write(*,*) deltaS(Sw,Swe,l)
	   if (deltaS(Sw,Swe,l).gt.1)then
	    Kmin=(Kmin+Kmax)/2.D0
	    !write(*,*)'r=',Rcorr(Sw,Swe,l)
	   elseif(deltaS(Sw,Swe,l).lt.1)then
	    Kmax=(Kmin+Kmax)/2.D0
	    !write(*,*)'r=',Rcorr(Sw,Swe,l)
	   endif

	  deltaK=abs(Ki-Kim1)
	  !write(*,*)deltaK,Ki,Kim1
	  Kim1=Ki

	 enddo

	 !write(*,*)rn,Ki,Rcorr(Sw,Swe,l),Rcorr_1,Rcorr(Sw,Swe,l)-Rcorr_1
	 deltaRcorr=Rcorr(Sw,Swe,l)-Rcorr_1
	 if(deltaRcorr.gt.0.D0)then
	  Rcorr_1=Rcorr(Sw,Swe,l)
	  Kib=Ki
	  do rrr=1,l
	   Sw_1(rrr)=Sw(rrr)
	  enddo
	 endif
	enddo

	do rrr=1,l
	 Sw(rrr)=Sw_1(rrr)
	enddo

	!write(*,*)k,deltaK
	write(*,*)
	write(*,'(a21,i3,a16,i3)')'Polimerization value=',rn,'; ki conside
     +red=',rn-1
	write(*,*)'The Ki [M^-1] is equal to: ',Kib
	write(*,*) 'correlation coefficient (r)=',Rcorr_1
	write(*,*) 'Standard error of the estimate=',stand_err(Sw,Swe,l)
	!write(*,*) 'ct[uM]=',conctM(Ki,c1,rn)*1000000
	!write(*,*) 'ci[uM]=',(((Ki**(kk-1))*(c1**kk)*1000000),kk=1,rn)
	write(*,*)'	ct[uM]			S20,w-teo		S20,w-exp'
	do rrr=1,l
	 write(*,*) ct(rrr)*1000000,Sw(rrr),Swe(rrr)
	enddo
	write(*,*) ''


	deallocate (ct,Swe,Sw,Sw_1)

	allocate (ct(100000),Sw(100000),stat=ALLOC_ERR)

	 ct(1)=0.001
	do i=2,1000
	 ct(i)=ct(i-1)+0.1D0
	enddo

	do i=1,1000
         ct(i)=ct(i)*0.000001D0        ![M] [mol/L]
        enddo

	open(5,file='SwKeq.log',status='new')

	do i=1,1000

           clow=0.D0
           chigh=ct(i)

           j=0
           deltaC=9999999999.D0
           do while(deltaC.gt.0.00000001)
            j=j+1

            c1=(clow+chigh)/2.D0

            if (conctM(Kib,c1,rn).gt.ct(i))then
             chigh=(clow+chigh)/2.D0
            elseif(conctM(Kib,c1,rn).lt.ct(i))then
             clow=(clow+chigh)/2.D0
            endif

            deltaC=abs((conctM(Kib,c1,rn)-ctm1)*1000000)
            ctm1=conctM(Kib,c1,rn)

	   enddo
           Sw(i)=Swt(Kib,c1,rn,MW,mod_num,s1)

	   write(5,*)ct(i)*1000000,Sw(i)
	enddo
	close(5)
	






	deallocate (ct,Sw)

c	write(*,*)'The Ki [M^-1] is equal to: ',Ki
c	write(*,*) 'ct[uM]=',conctM(Ki,c1,rn)*1000000
c	write(*,*) 'ci[uM]=',(((Ki**(i-1))*(c1**i)*1000000),i=1,rn)

	end




C FUNCTION BEGINING






	function conct(Ki,c1,rn,MW)
	integer i
	integer rn,ir
	real*8 Ki,conct,c1,MW


	conct=0

	do i=1,rn
	 conct=conct+(i*MW*(Ki**(i-1))*((c1)**(i)))   !conct [g/L]
	enddo

	return
	end



	function conctM(Ki,c1,rn)
	integer i
	integer rn,ir
	real*8 Ki,conctM,c1

	conctM=0

	do i=1,rn
	 conctM=conctM+(((Ki)**(i-1))*((c1)**(i)))          !conctM [mol/L]
	enddo

	return
	end



	function Swt(Ki,c1,rn,MW,mod_num,s1)
	integer i,mod_num
	integer rn,ir
	real*8 Ki,conct,c1,Swt,MW,si,s1

	Swt=0

	do i=1,rn
	 !si=5.09D0*(i**(2.D0/3.D0))       !5.82
	 Swt=Swt+((si(i,mod_num,s1)*(1.D0-0.018D0*c
     +onct(Ki,c1,rn,MW)))*((i*((Ki/MW)**(i-1))*((c1
     +*MW)**(i)))/(conct(Ki,c1,rn,MW))))

	enddo

	return
	end



	function si(kk,mod_num,s1)
	real*8 si,s1,i
	integer mod_num,kk

	  !s1=5.82D0  !5.82 T=20oC;  5.09 T<20oC

	i=dlog(kk*1.D0)
	s1=dlog(s1)

	if(mod_num.eq.1)then
	 si=s1+i*(2.D0/3.D0)
	elseif(mod_num.eq.2)then
	 si=s1+((0.4388D0*i)-(0.03172D0*(i**2.D0)))
	elseif(mod_num.eq.3)then
	 si=s1+((0.5932D0*i)-(0.0536D0*(i**2D0)))
	elseif(mod_num.eq.4)then
	 si=s1+((0.3664D0*i)-(0.01102D0*(i**2D0)))
	elseif(mod_num.eq.5)then
	 si=s1+((0.4911D0*i)-(0.0006354D0*(i**2D0)))
	elseif(mod_num.eq.6)then
	 si=s1+((0.4315D0*i)-(0.01306D0*(i**2D0)))
	endif

	si=dexp(si)
	s1=dexp(s1)

	end



	function deltaS(Sw,Swe,l)

c	Estoy calcluando la pendiente para n .eq. 0, antes estaba para n .ne. 0
	integer i,l
	real*8 Sw(l),Swe(l),deltaS,Swm,Swem,deltaS2

	deltaS=0
	!Swm=0
	!Swem=0

	!do i=1,l
	! Swm=Swm+Sw(i)
	! Swem=Swem+Swe(i)
	!enddo

	!Swm=Swm/l
	!Swem=Swem/l

	deltaS2=0
	do i=1,l
	! deltaS2=deltaS2+((Sw(i)-Swm)**2)
	 deltaS2=deltaS2+(Sw(i)**2)
	enddo

	do i=1,l
	 !deltaS=deltaS+(((Sw(i)-Swm)*(Swe(i)-Sewm))/(deltaS2))
	 deltaS=deltaS+((Sw(i)*Swe(i))/(deltaS2))
	enddo

	return
	end



	function Rcorr(Sw,Swe,l)

	integer i,l
	real*8 Sw(l),Swe(l),Rcorr,Swm,Swem,deltaS2,deltaSe2

	Rcorr=0
	Swm=0
	Swem=0

	do i=1,l
	 Swm=Swm+Sw(i)
	 Swem=Swem+Swe(i)
	enddo

	Swm=Swm/l
	Swem=Swem/l

	deltaS2=0
	deltaSe2=0
	do i=1,l
	 deltaS2=deltaS2+((Sw(i)-Swm)**2)
	 deltaSe2=deltaSe2+((Swe(i)-Swem)**2)
	enddo

	do i=1,l
	 Rcorr=Rcorr+(((Sw(i)-Swm)*(Swe(i)-Swem))/(dsqrt(deltaS2*
     +deltaSe2)))

	enddo

	return
	end



	function stand_err(Sw,Swe,l)
	integer i,l
	real*8 Sw(l),Swe(l),stand_err

	stand_err=0

	do i=1,l
	 stand_err=stand_err+(((Swe(i)-Sw(i))**2.d0)/l)
	enddo

	stand_err=dsqrt(stand_err)

	return
	end
