PROGRAM non_markov_cca

   IMPLICIT none
   integer :: i, k, j
   integer, parameter :: N_cav = 6000
   integer, parameter :: N_tot = N_cav + 1
   integer, parameter :: Nsob2 = N_tot/2
   integer, parameter :: Nsob2_plus1 = Nsob2 + 1
   real, dimension(N_tot, N_tot) :: H
   real, dimension(N_tot) :: AutValores
   integer :: omega_c, loop_amostra, aux_hpp
   real :: g
!  Lapack routine
   integer, parameter :: NP = N_tot, N = N_tot
   integer, parameter :: IMAX = 10 + 6*NP + 2*(NP**2), IMAX2 = 10 + 5*NP
   real, dimension(IMAX) :: WORK
   integer, dimension(IMAX2) :: IWORK
   integer :: INFO
   real :: ABSTOL 
   integer, parameter :: LWORK = 10 + 6*NP + 2*(NP**2), LIWORK = 10 + 5*NP
   CHARACTER (LEN=1) :: COMPZ="V", UPLO="U"
! Time evolution operator
   real, parameter :: dt = 0.05, tmax = 150.
   integer, parameter :: auxTMP = int(tmax/dt)
   real :: tt, norma
   complex, dimension(N_tot) :: f
   complex :: c1
   complex, parameter :: ic = cmplx(0.,1.)
   integer :: t
!  File variables
   character (len=200) :: arq, arq1, arq0
!  Variaveis para dinamica sistema 
   real, parameter :: ni = 0.1
! Energy variables
   double precision :: med1, med2, auxE, ir, kr, aux1, aux2
   double precision, parameter :: pi = 3.141592653589793238462643383279502884197169399375105820974944592307d0
   double precision, parameter :: two_pi = 2.0d0*pi
   double precision, parameter :: two_pi_N = two_pi/N
   double precision, dimension(N) :: phi, hpp
   integer :: semente
   real :: ran1
   double precision, parameter :: alpha = 0.0d0 ! grau de correlação
   double precision, parameter :: auxalpha = -alpha*0.5d0 

	DO loop_amostra = 1, 200

		H = 0.0
		AutValores = 0.
	
! Genarationg seed 
   		semente = 22 + loop_amostra
   		WRITE(arq, '("RTN_ni",f0.2,"alpha",f0.2,"amostra",i0,".dat")') ni, alpha, loop_amostra
    	WRITE(arq0, '("ener_ni",f0.2,"alpha",f0.2,"amostra",i0,".dat")') ni, alpha, loop_amostra
   		semente = -semente
  		phi = 0
  		
! Random number between [0, 1)*2*pi 
  		DO i = 1, (N/2)
      		phi(i) = two_pi*ran1(semente)
   		END DO    
! Correlated disorder
   		med1 = 0.0d0; med2 = 0.0d0; aux1 = 0.0d0; aux2 = 0.0d0
   		DO i = 1, N
      		ir = dble(i)
      		hpp(i) = 0.0E0
      		DO k = 1, Nsob2
         		kr = dble(k)
        		aux1 = kr**(auxalpha)
         		aux2 = cos( two_pi_N*ir*kr + phi(k) )
         		hpp(i) = hpp(i) + aux1*aux2
      		END DO
      		med1 = med1 + hpp(i)
      		med2 = med2 + hpp(i)*hpp(i)
   		END DO
! normalizing the Energy
   		med1 = med1/N
   		med2 = med2/N
! random energy
   		DO i = 1, N
      		hpp(i) = ( (hpp(i) - med1) / sqrt(med2 - med1*med1) ) 
      		OPEN(unit=1, file=arq0)
      		WRITE(1,*) i, hpp(i)
   		END DO

! filling the matrix
   		H = 0.; aux_hpp = 1.
   		g = ni*aux_hpp

   		DO i = 1, N_cav
      		H(i, i + 1) = hpp(i)
      		H(i + 1, i) = hpp(i)
   		END DO

   		DO i = 1, N_tot
      		H(i, i) = 0
   		END DO

  		H(Nsob2, Nsob2_plus1) = hpp(Nsob2)*g; H(Nsob2_plus1, Nsob2) = hpp(Nsob2)*g
   		H(Nsob2, Nsob2 + 2) = hpp(Nsob2 + 1); H(Nsob2 + 2, Nsob2) = hpp(Nsob2 + 1)
   		H(Nsob2_plus1, Nsob2_plus1) = 0
   		H(Nsob2_plus1, Nsob2 + 2) = 0; H(Nsob2 + 2, Nsob2_plus1) = 0

! printing matrix 
!   	DO i = 1, N_tot
!      		PRINT *, H(i, :)
!   	END DO

		PRINT *, "Escrevendo amostra", loop_amostra

  ! INFO = 0; ABSTOL = 1e-20
   		call SSYEVD(COMPZ,UPLO,N,H,NP,AutValores,WORK,LWORK,IWORK,LIWORK,INFO)

   		write(arq1, '("ENG_ni",f0.2,".dat")') ni
   		open(unit=16, file=arq1)          
   		do i = 1, N_tot
      		write(16, *) AutValores(i)
   		end do

   		do t = 1, auxTMP
   	   		norma = 0.
   	   		do i = 1, N
   	     		f(i) = cmplx(0.,0.)
   	     		do j = 1, N
   	        		c1 = ic*AutValores(j)*tt    !iEt
   	        		c1 = cexp(-c1)              !e^{-icEt}
   	        		f(i) = f(i) + H(i,j)*H(Nsob2_plus1,j)*c1   ! time evolution op
   	     		end do   
   	    		norma = norma + (cabs(f(i)))**2.
   	   		end do
   	   		open(unit=10, file=arq)
   	   		write(10, *) cabs(f(Nsob2_plus1))**2.
   	   		tt = tt + dt
   		end do
   
   END DO 


END PROGRAM non_markov_cca


FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836)
      PARAMETER (NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
      idum=max(-idum,1)
      do 11 j=NTAB+8,1,-1
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k

      if (idum.lt.0) idum=idum+IM
      if (j.le.NTAB) iv(j)=idum
11    continue
      iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
