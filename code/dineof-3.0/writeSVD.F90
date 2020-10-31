#include "ppdef.h"

subroutine writeSVD(DirOutput,s,P,u,n,m,v,sumVarIni,sumVarEnd)
 use ufileformat
 implicit none

 character(len=200),intent(in)     :: DirOutput
 real,intent(in)                   :: s(:,:),u(:,:),v(:,:),sumVarIni,sumVarEnd
 integer,intent(in)                :: P,n,m

 integer                           :: i,j,k
 real, parameter                   :: valex = 9999.
 real                              :: scvlsng  

 write(stdout,*) 'total variance of the initial matrix ',sumVarIni/(M*N)
 write(stdout,*) 'total variance of the reconstructed matrix ',sumVarEnd/(M*N)

 !   Variance detailed
 !   _______________

 open(1,file=trim(DirOutput)//'/outputEof.varEx')
 do i=1,P
   write(1,504) 'Mode ', i, ' = ',&
        &               100.0*s(i,1)*s(i,1)/sumVarEnd, ' %'
 enddo
 close(1)

 do i=1,P
   scvlsng=scvlsng+(s(i,1)*s(i,1))
 enddo
 write(stdout,*) 'Sum of the squares of the singular values of the ',P, 'eof retained', scvlsng

 !
 !   Singular values
 !   _______________
 open(1,file=trim(DirOutput)//'/outputEof.vlsng')
 do i=1,P
   write(1,501) s(i,1)
 enddo
 close(1)
 !
 !
 !   Spatial modes 
 !   _______________
 open(1,file=trim(DirOutput)//'/outputEof.lftvec')

 do j=1,m

   write(1,503) (u(j,i),i=1,P)
 enddo
 close(1)
 !
 !   Temporal modes 
 !   _______________
 open(1,file=trim(DirOutput)//'/outputEof.rghvec')
 do j=1,n

   write(1,503) (v(j,i),i=1,P)
 enddo
 close(1)


 !stop      

500 format(a72)
501 format(e10.4)
502 format(3e10.4)
503 format(300(e10.4,1X))
504 format(a5,i3,a3,f8.4,a2)

end subroutine writeSVD
