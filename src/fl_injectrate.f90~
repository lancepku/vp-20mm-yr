!***********************
!inject rate depend on Axial valley depth
!*******************
subroutine fl_injectrate
use arrays
include 'precision.inc' ! vs BI(Behn and Ito 2008): the same
include 'params.inc'   ! vs BI: not the same
include 'arrays.inc'   ! vs BI: not the same
! (J-1,I-1)------(J-1,I)----(J-1,I+1)
!     |            |            |
!     |   a11      |   a12      |
!     | (j-1,i-1)  | (j,i-1)    |
!     |            |            |
!  (J,I-1)-------(J,I)--------(J,I+1)
!     |            |            |
!     |   a21      |    a22     |
!     | (j-1,i)    |   (j,i)    |
!     |            |            |
! (J+1,I-1)------(J+1,I)----(J+1,I+1)

!dimension njTinj(nz)!
D_axial = 0.25*(cord(1,nx/2,2)+cord(2,nx/2,2)+cord(1,nx/2+1,2)+cord(2,nx/2+1,2))
!print *,'D',D_axial
!define the depth of AMC
Z1100 = 1
do j=1,nz
if (temp(j,nx/2).lt.1100.)then
Z1100=j+1
else if (temp(j,nx/2).ge.1100)then
Z1100 = Z1100
end if
end do
Zamc = 0.25*abs(cord(int(Z1100),nx/2,2)+cord(int(Z1100)+1,nx/2+1,2)+cord(int(Z1100)+1,nx/2,2)+cord(int(Z1100),nx/2+1,2))+ D_axial
Z600 = 1
do j=1,nz
tempar = temp(j,nx/2)!+temp(j+1,nx/2)+temp(j,nx/2+1)+temp(j+1,nx/2+1))
!do j=1,nz
if (tempar.lt.600.)then
Z600=j+1
else if (tempar.ge.600)then
Z600 = Z600
end if
end do
!print *,Z1100
HL=0.25*abs(cord(int(Z600),nx/2,2)+cord(int(Z600)+1,nx/2+1,2)+cord(int(Z600)+1,nx/2,2)+cord(int(Z600),nx/2+1,2))+D_axial
orgh1 = 3000 ! density between magma and rock
orgh2 = 20000 ! density between water and rock
pd = 1.e7 ! driving pressure
D = -D_axial ! make sure axial depth for valley is positive
!print *,D
HG1 = -3.*(orgh1*HL+orgh2*D)
HG2 = sqrt(9*((orgh1*HL+orgh2*D)**2)+24*orgh1*pd*HL)
!print *,HG2
HG= (HG1+HG2)/(2.*orgh1) ! see Liu&Buck 2017 for detial
!HG=HL+2*D_axial  ! should redescribe the HG^C
Print *,'HL',HL,'HG',HG,D

if (Zamc.lt.6.e3) then
rate_inject_HG = Vp_
else if (Zamc.ge.6.e3)then
rate_inject_HG = Vp_*6.e3/(HL+0.5*HG)
endif
rate_inject= min(rate_inject_HG,Vp_)
!print *,'1', rate_inject
open (10,file='result.dat')
write(10,*),Zamc,HL,HG,D,rate_inject
close(10)
return
end
!*******************************
