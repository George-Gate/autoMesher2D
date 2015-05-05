program rectangle
    use mesh
    implicit none
    real(8), parameter :: minR1=8d-2, maxR1=1d-1
    real(8), parameter :: minR2=1d-2, maxR2=2d-2
    real(8) :: x(4),y(4)
    
    x=(/0d0,1d0,1d0 ,0d0/)
    y=(/0d0,0d0,1d-1,1d-1/)
    
    call initMesh(x,y,minR1,maxR1,step=minR2/2d0,outputDetails=.true.)
    pause
    call refineMesh(minR2,maxR2,outputDetails=.true.)
    pause
    call cleanMesh()
    call initMesh(x,y,minR2,maxR2,outputDetails=.true.)
    call generateVertexID()   
    open(81,file='finalMesh.txt')
    call outputMesh(81)
    close(81)
    
end program rectangle