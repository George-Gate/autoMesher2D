program circle
    use mesh
    implicit none
    integer, parameter :: N=120
    real(8), parameter :: minR=1d-1, maxR=2d-1
    real(8) :: x(N),y(N), theta
    integer :: i, nKP, KP(3), smallEnd(2,1), bigEnd(2,2)
    
    nKP=1
    do i=1,N
        theta=2*pi/N*i
        x(i)=dcos(theta)
        y(i)=dsin(theta)
        if (theta<=2*pi*nKP/3 .and. 2*pi*nKP/3<2*pi/N*(i+1)) then
            KP(nKP)=i
            nKP=nKP+1
        end if
    end do
    
    smallEnd(:,1)=(/KP(1),KP(2)/)
    bigEnd(:,1)=(/KP(3),KP(1)/)
    bigEnd(:,2)=(/KP(2),KP(3)/)
    
    call initMesh(x,y,minR,maxR,keyPoints=KP,outputDetails=.true.)
    call generateVertexID(smallEndArea=smallEnd, bigEndArea=bigEnd)   
    open(81,file='finalMesh.txt')
    call outputMesh(81)
    close(81)
    
end program circle