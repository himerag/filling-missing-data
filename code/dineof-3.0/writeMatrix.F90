#include "ppdef.h"

!---------------------------------------------------------
!----- subroutine writeMatrix ----------------------------
!---------------------------------------------------------

! reshapes matrix into 3D file, then writes results  
! into separate files
! 
!---------------------------------------------------------
                      
subroutine writeMatrix(X,xmean,valex,resultfnames,maskfile,norma,M,N,imax,jmax,first,DirOutput,valc,fileMean,fileStd)
use ufileformat
implicit none

real,intent(in)                  :: X(:,:),fileMean(:),fileStd(:)
integer,intent(in)               :: N,M,imax(:),jmax(:),first(:)
real,intent(in)                  :: xmean,valc(:)
real,intent(inout)               :: valex
integer,intent(in)               :: norma

character(len=200),intent(in)    :: resultfnames(:),maskfile(:)
character(len=200),intent(in)    :: DirOutput

real, allocatable                :: sst(:,:,:)
real,pointer                     :: mask(:,:)
real, parameter                  :: valexc = 9999.
real                             :: mean,stddev,var
integer                          :: i,j,s,t,q,nbvars,NN


nbvars=size(resultfnames)


!Reshape the matrix into a 3D matrix (two spatial dimensions and one temporal dimension)
!---------------------------------------------------------------------------------------
do q=1,nbvars
     
   allocate(sst(imax(q),jmax(q),N))
     
     call uload(maskfile(q),mask,valex)
     where (mask.eq.0) mask = valexc

     do t=1,N

        s = first(q)

        do i=1,imax(q)
           do j=1,jmax(q)

              if(mask(i,j).ne.valexc) then

                 sst(i,j,t)=X(s,t);

                  if(norma.eq.1) then
                    sst(i,j,t)=(sst(i,j,t)*fileStd(q))+fileMean(q)
                  end if

                 s = s+1
              else
                 sst(i,j,t)=valexc;
              endif
           enddo
        enddo
     enddo


    
  NN = count(sst.ne.valexc)
  mean = sum(sst,sst.ne.valexc)/NN
  var = sum((sst-mean)**2,sst.ne.valexc)/NN
  stddev = sqrt(var)

!   -----------------------------------------
!   some statistics about the filled matrices
!   -----------------------------------------

    write(stdout,*)
    write(stdout,*)'mean ',mean
    write(stdout,*)'Standard deviation',stddev,fileStd,fileMean
    write(stdout,*)
  
    call usave(trim(resultfnames(q)),sst(:,:,:),valexc)

deallocate(sst,mask)
enddo


!---------------------------------------------------------------------------------------

!   Valc: write expected error for each mode
!   ----------------
    call usave(trim(DirOutput)//'/valc.dat',valc,valexc)    
!---------------------------------------------------------------------------------------


end subroutine writeMatrix

   



   
