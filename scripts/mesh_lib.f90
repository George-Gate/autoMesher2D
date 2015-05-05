!include "mesh.f90"

       
! initMesh(): 自动生成网格
!----------------------------------------------------------------------------------
!  1、调用makeGeometry()初始化求解域geometry
!  2、初始化meshInfo变量
!  3、建立最小网格
!  4、调用refineMesh()细化网格
!  5、调用validateMeshInfo()检查网格有效性
![Input]
!  x, y - 求解域边界的x、y坐标。求解域必须为单连通多边形区域，
!         以边界正(右手螺旋)方向按顺序给出多边形各顶点的xy坐标。
!         所谓的正方向是指：以正方向沿着边界走时，求解域在左手边。
! minR, maxR - 网格的最大/最小精度。
!              边界附近网格较密，中心区域网格较疏。 
! keyPoints - [可选]关键点在x、y中的下标，关键点一定会成为边界格点。
!             关键点的下标必须从小到大给出。
!             强烈建议将不同类型边界的分界点指定为关键点。
! outputDetails - [可选]指示是否需要将每一步迭代过后的网格状态输出至文件mesh.txt
!                 默认为否
!  step  - [可选]求解域边界采样时的几何精度，两点之间的距离不小于step/4, 不大于step，
!          默认取step=minR/2。
!----------------------------------------------------------------------------------
    subroutine initMesh(x, y, minR, maxR, keyPoints, outputDetails, step)
        use mesh
        implicit none
        real(8), intent(in) :: x(:), y(:)
        real(8), intent(in) :: minR, maxR
        integer, intent(in), optional :: keyPoints(:)
        logical, intent(in), optional :: outputDetails
        real(8), intent(in), optional :: step
        type(point_t)       :: p(0:3)
        integer             :: i, j, n, n0, iv(0:3), counter
        real(8) :: currentR
        logical :: succeed
        
        ! Make geometry
        if (present(step)) then
            call makeGeometry(x,y,step,keyPoints)
        else
            call makeGeometry(x,y,minR/2d0,keyPoints)
        end if
        
        ! Init meshInfo
        ! 分配内存
        n=max(100,int(4*maxval(x)*maxval(y)/(pi*minR**2)))   ! 参数待实测优化
        meshInfo.maxN=n
        meshInfo.nv=0
        meshInfo.ne=0
        meshInfo.ns=0
        meshInfo.nb=0
        allocate(meshInfo.vertex(n))
        allocate(meshInfo.edge(6*n))   ! 假设格点的平均度小于12
        allocate(meshInfo.surface(4*n))  ! 每条边与两个面元相邻，每个面元有3条边
        allocate(meshInfo.boundary(n-1))   ! 考虑极端情况，只有一个内点，其余都是边界
        
        ! 初始化最简网格
        currentR=min( maxval(x)-minval(x) , maxval(y)-minval(y) )  ! 以尽可能大的步长初始化网格，以提高网格质量
        succeed=.false.
        do while (currentR>minR)
            succeed=.false.
            ! 在边界上找凸角
            n=geometry.n
            i=1
            p(1)=geometry.b(i)
            iv(1)=i
            i=getNextBpoint(i, minR=currentR)
            if (i/=-1) then  ! 如果找不到点，则缩小currentR
                iv(2)=i
                p(2)=geometry.b(i)
                counter=geometry.n
                do while (counter>0)
                    counter=counter-1
                    i=getNextBpoint(i, minR=currentR)
                    if (i==-1) then  ! 找不到点则退出循环
                        exit
                    end if
                    iv(3)=i
                    p(3)=geometry.b(i)
                    p(0)=(p(1)+p(2)+p(3))/3d0
                    ! 判断找到的点是否符合条件：
                    ! 1. 为凸角
                    ! 2. p(0)在求解域内
                    ! 3. 三条内边完全在求解域内
                    ! 4. 三点的两两距离大于currentR
                    if ((p(2)-p(1))*(p(3)-p(2))>0 .and. &   ! 判断1
                        inArea(p(0))==1 .and. &             ! 判断2
                        inArea(p(1),p(0))==-1 .and. inArea(p(2),p(0))==-1 .and. inArea(p(3),p(0))==-1 .and. & ! 判断3  
                        norm(p(3)-p(1))>=currentR) then ! 判断4, 1->2和2->3的距离必定大于currentR
                        succeed=.true.
                        exit
                    end if
                    p(1:2)=p(2:3)
                    iv(1:2)=iv(2:3)
                    
                    ! 结束条件：iv(1)<1<=iv(2)
                    if (iv(1)>iv(2)) then
                        exit
                    end if
                end do
            end if
            currentR=currentR/1.1d0
        
            ! 如果前面找到满足条件的点，则初始化网格
            if (succeed) then
                ! 加格点
                iv(0)=addVertex(p(0))
                do j=1,3
                    iv(j)=addVertex(p(j),isB=.true.,bid=iv(j))
                end do
                ! 加边
                j=addEdge(iv(1),iv(2),isB=.true.)
                j=addEdge(iv(2),iv(3),isB=.true.)
                j=addEdge(iv(3),iv(1),isB=.true.)
                j=addEdge(iv(0),iv(1))
                j=addEdge(iv(0),iv(2))
                j=addEdge(iv(0),iv(3))
                ! 加面元
                j=addSurface(iv(0),iv(1),iv(2))
                j=addSurface(iv(0),iv(2),iv(3))
                j=addSurface(iv(0),iv(3),iv(1))
                exit
            end if
        end do
        
        if (.not. succeed) then
            write (*,*) "initMesh(): 无法初始化网格，请尝试输入更小的minR。"
            return
        end if
        
        ! 细化网格
        call refineMesh(minR, maxR, outputDetails)
        
        ! 检查meshInfo中的信息是否满足定义
        if (validateMeshInfo() /= 0) then
            write(*,"(A)") "initMesh(): meshInfo未通过有效性检查，请尝试改变网格划分条件。"
        end if
        return
    end subroutine initMesh

! cleanMesh(): 清空meshInfo和geometry的内容，释放内存
    subroutine cleanMesh()
        use mesh 
        implicit none
        ! clean geometry
        geometry.n=0
        deallocate(geometry.b)
        deallocate(geometry.l)
        deallocate(geometry.keyPointID)
        deallocate(geometry.bidOfKeyPoint)
        
        ! clean meshInfo
        meshInfo.maxN=0
        deallocate(meshInfo.vertex)
        deallocate(meshInfo.edge)
        deallocate(meshInfo.surface)
        deallocate(meshInfo.boundary)
    end subroutine cleanMesh
    
    
! validateMeshInfo(): 检查meshInfo中的信息是否满足定义的要求，并尝试修复错误。
!----------------------------------------------------------------------------------
!  具体的检查内容为：
!   1. 面元的面积、外接圆半径是否正确
!   2. 格点的bid与坐标信息是否相符
!   3. 是否有未标记为边界点的格点在边界上
!   4. 面元的1号顶点是否为非边界点
!   5. 面元中的顶点编号是否按照逆时针编号
!   6. 面元中的顶点与边端点是否匹配，顶点的对边的编号是否与该顶点相同
!   7. 边的长度是否正确
!   8. 边两侧的面元是否正确
!   9. 边的isB信息项是否正确
!   10.mehsInfo.boundary中的边是否都是边界边
![Return]
!  如果检查没有发现问题，则返回0；否则返回发现的问题在上述分类中的编号。
!  若存在多个问题，则会返回最后一个问题的编号。
!----------------------------------------------------------------------------------
    integer function validateMeshInfo() result(re)
        use mesh
        implicit none
        integer :: i, j, k, id(3), id2(3), counter
        real(8) :: R, S
        type(point_t) :: p(3)
        type(surface_t), pointer :: ps
        logical :: flag
        
        re=0
        ! 检查1, 如发现问题则自动修复
        do i=1,meshInfo.ns
            S=calcArea(i)
            R=calcRadius(i,S)
            if (abs(meshInfo.surface(i).S-S)>eps .or. abs(meshInfo.surface(i).R-R)>eps) then
                re=1
                meshInfo.surface(i).S=S
                meshInfo.surface(i).R=R
                write(*,"(A,I0,A)") "validateMeshInfo(): 面元 ",i," 的面积/半径信息错误，已修复。"
            end if
        end do
        
        ! 检查2
        do i=1,meshInfo.nv
            if (meshInfo.vertex(i).isB) then
                j=meshInfo.vertex(i).bid
                if (j<1 .or. j>geometry.n) then
                    re=2
                    meshInfo.vertex(i).bid=.false.
                    write(*,"(A,I0,A)") "validateMeshInfo(): 格点 ",i," 的bid信息越界，已修复。"
                    cycle
                end if
                if ( norm(geometry.b(j)-meshInfo.vertex(i).p)>eps ) then
                    re=2
                    write(*,"(A,I0,A)") "validateMeshInfo(): 格点 ",i," 的bid信息坐标不符，无法修复。"
                end if
            end if
        end do
        
        ! 检查3
        do i=1,meshInfo.nv
            if ( (.not. meshInfo.vertex(i).isB) .and. &
                 inArea(meshInfo.vertex(i).p)==2) then
                re=3
                write(*,"(A,I0,A)") "validateMeshInfo(): 格点 ",i," 在边界上，但未被标记为边界点，无法修复。"
            end if
        end do
        
        ! 检查4~6
        do i=1,meshInfo.ns
            ! 检查4
            k=meshInfo.surface(i).v(1)
            if (meshInfo.vertex(k).isB) then
                re=4
                id(2:3)=meshInfo.surface(i).v(2:3)
                if (.not. (meshInfo.vertex(id(2)).isB .and. meshInfo.vertex(id(3)).isB)) then
                    if (.not. meshInfo.vertex(id(2)).isB) then
                        meshInfo.surface(i).v(1)=meshInfo.surface(i).v(2)
                        meshInfo.surface(i).v(2)=k
                    else
                        meshInfo.surface(i).v(1)=meshInfo.surface(i).v(3)
                        meshInfo.surface(i).v(3)=k
                    end if
                    write(*,"(A,I0,A)") "validateMeshInfo(): 面元 ",i," 的v(1)为边界点，已修复。"
                else
                    write(*,"(A,I0,A)") "validateMeshInfo(): 面元 ",i," 的v(1)为边界点，无法修复。"
                end if
            end if
            
            ! 检查5
            id(:)=meshInfo.surface(i).v(:)
            p(:)=meshInfo.vertex(id(:)).p
            if ((p(1)-p(2))*(p(3)-p(2))>0) then
                re=5
                k=meshInfo.surface(i).v(2)
                meshInfo.surface(i).v(2)=meshInfo.surface(i).v(3)
                meshInfo.surface(i).v(3)=k
                write(*,"(A,I0,A)") "validateMeshInfo(): 面元 ",i," 的顶点没有按照逆时针编号，已修复。"
            end if
            
            ! 检查6
            id(:)=meshInfo.surface(i).v(:)
            id2(:)=meshInfo.surface(i).e(:)
            if (minval(id2)<1 .or. maxval(id2)>meshInfo.ne) then
                id2=1
            end if
            do j=1,3
                if (.not. ( (meshInfo.edge(id2(j)).v(1)==id(mod(j,3)+1) .and. meshInfo.edge(id2(j)).v(2)==id(mod(j+1,3)+1)) .or. &
                            (meshInfo.edge(id2(j)).v(2)==id(mod(j,3)+1) .and. meshInfo.edge(id2(j)).v(1)==id(mod(j+1,3)+1)) )) then
                    counter=0
                    do k=1,meshInfo.ne
                        if (meshInfo.edge(k).v(1)==id(1) .or. &
                            meshInfo.edge(k).v(1)==id(2) .or. &
                            meshInfo.edge(k).v(1)==id(3)) then
                            if      ( (meshInfo.edge(k).v(1)==id(2) .and. meshInfo.edge(k).v(2)==id(3)) .or. &
                                        (meshInfo.edge(k).v(2)==id(2) .and. meshInfo.edge(k).v(1)==id(3)) ) then
                                meshInfo.surface(i).e(1)=k
                                counter=counter+1
                            else if ( (meshInfo.edge(k).v(1)==id(3) .and. meshInfo.edge(k).v(2)==id(1)) .or. &
                                        (meshInfo.edge(k).v(2)==id(3) .and. meshInfo.edge(k).v(1)==id(1)) ) then
                                meshInfo.surface(i).e(2)=k
                                counter=counter+1
                            else if ( (meshInfo.edge(k).v(1)==id(1) .and. meshInfo.edge(k).v(2)==id(2)) .or. &
                                        (meshInfo.edge(k).v(2)==id(1) .and. meshInfo.edge(k).v(1)==id(2)) ) then
                                meshInfo.surface(i).e(3)=k
                                counter=counter+1
                            end if
                        end if
                    end do
                    re=6
                    if (counter==3) then
                        write(*,"(A,I0,A)") "validateMeshInfo(): 面元 ",i," 的面元信息有误，已修复。"
                    else
                        write(*,"(A,I0,A)") "validateMeshInfo(): 面元 ",i," 的面元信息有误，无法修复。"
                    end if
                    exit
                end if
            end do
        end do
        
        ! 检查7
        do i=1,meshInfo.ne
            R=norm(meshInfo.vertex(meshInfo.edge(i).v(1)).p-meshInfo.vertex(meshInfo.edge(i).v(2)).p)
            if (abs(R-meshInfo.edge(i).len)>eps) then
                re=7
                meshInfo.edge(i).len=R
                write(*,"(A,I0,A)") "validateMeshInfo(): 边 ",i," 的长度信息有误，已修复。"
            end if
        end do
        
        ! 检查8
        do i=1,meshInfo.ns
            id(:)=meshInfo.surface(i).v(:)
            id2(:)=meshInfo.surface(i).e(:)
            do j=1,3
                flag=.false.
                p(1:2)=meshInfo.vertex(meshInfo.edge(id2(j)).v(:)).p
                p(3)=meshInfo.vertex(id(j)).p
                if ((p(3)-p(1))*(p(2)-p(1))>0) then  ! 面元在边的右边
                    if ( meshInfo.edge(id2(j)).s(1)/=i ) then
                        flag=.true.
                        meshInfo.edge(id2(j)).s(1)=i
                    end if
                else ! 左边
                    if ( meshInfo.edge(id2(j)).s(2)/=i ) then
                        flag=.true.
                        meshInfo.edge(id2(j)).s(2)=i
                    end if
                end if
                if (flag) then
                    re=8
                    write(*,"(A,I0,A)") "validateMeshInfo(): 边 ",id2(j)," 与面元的相邻信息有误，已修复。"
                end if
            end do    
        end do
        
        ! 检查9
        do i=1,meshInfo.ne
            if (meshInfo.edge(i).isB) then
                if (meshInfo.edge(i).s(1)/=-1 .and. meshInfo.edge(i).s(2)/=-1) then
                    meshInfo.edge(i).isB=.false.
                    re=9
                    write(*,"(A,I0,A)") "validateMeshInfo(): 边 ",i," 被错误地设置为边界边，已修复。"
                end if
            elseif (meshInfo.edge(i).s(1)==-1 .or. meshInfo.edge(i).s(2)==-1) then
                if (.not. meshInfo.edge(i).isB) then
                    meshInfo.edge(i).isB=.true.
                    re=9
                    write(*,"(A,I0,A)") "validateMeshInfo(): 边 ",i," 被错误地设置为非边界边，已修复。"
                end if
            end if
        end do
        
        ! 检查10
        do i=1,meshInfo.nb
            if (.not. meshInfo.edge(meshInfo.boundary(i)).isB) then
                re=10
                write(*,"(A,I0,A)") "validateMeshInfo(): meshInfo.boundary(",i,")不是边界边......无法修复。"
            end if
        end do
        
        
    end function validateMeshInfo
    
! generateVertexID(): 建立格点编号,结果记录在meshInfo.vertex(:).id中
!----------------------------------------------------------------------------------
![Input]
! smallEndArea - [可选]希望最先编号的边界区域，位于这些边界区域的格点的编号将最小。
!                数组的尺寸为(2,n)，每一个区域给定起始关键点和结束关键点的编号，区
!                域的方向为边界正方向。(1,n)为起点，(2,n)为结束点。关键点编号指调用
!                initMesh()时该点在x、y数组中的下标。编号时将按照smallEndArea中区域
!                的给定顺序从起始点往结束点编号，先在smallEndArea中出现的区域编号较小。
!  bigEndArea  - [可选]希望最后编号的边界区域，位于这些边界区域的格点的编号将最大。
!                给定方式同smallEndArea。编号顺序则相反，将从最后一个区域的结束点开
!                始从大到小编号。故bigEndArea中最后出现的区域的结束点编号最大。
!
!  注意，以上的区域包括了区域边界的关键点。另外如果某个点同时出现在了多个不同的区域
!  中，那么该点的编号为该点第一次获得的编号。例如，某点同时包含在smallEndArea和bigEndArea
!  中，那么该点的编号会优先满足smallEndArea的要求；如果某点同时出现在了smallEndArea
!  的第一段边界和第三段边界中，那么该点的编号会在第一段边界中指定；再例如如果某点同
!  时是bigEndArea的最后一段边界的起点和结束点(并且该点不在smallEndArea中)，那么该点
!  的编号将会是所有点中最大的，因为当起点与结束点相同时，这段边界被认为是包含了整个
!  求解域边界。
!----------------------------------------------------------------------------------
    subroutine generateVertexID(smallEndArea, bigEndArea)
        use mesh
        implicit none
        integer, intent(in), optional :: smallEndArea(:,:), bigEndArea(:,:)
        integer, allocatable :: keyPointBID(:), bid2VertexSubscript(:)
        integer :: maxKeyPointID, i, j, k, ng, error, n
        integer :: bigEnd, smallEnd, b, e
        logical :: firstLoop
        
        ! 检查输入参数的形状
        if (present(smallEndArea)) then
            if (size(smallEndArea,1)/=2) then    
                write(*,"(A)") "generateVertexID(): 输入参数的尺寸不对，smallEndArea必须是2*n的整形数组。"
                return
            end if
        end if
        if (present(bigEndArea)) then
            if (size(bigEndArea,1)/=2) then    
                write(*,"(A)") "generateVertexID(): 输入参数的尺寸不对，bigEndArea必须是2*n的整形数组。"
                return
            end if
        end if
        
        ! 重置id
        meshInfo.vertex(:).id=-1
        
        ! 建立边界点、关键点与格点的对应关系
        if (present(smallEndArea) .or. present(bigEndArea)) then
            ng=geometry.n
            ! 分配内存
            maxKeyPointID=0
            do i=1,ng
                if (geometry.keyPointID(i)>maxKeyPointID) then
                    maxKeyPointID=geometry.keyPointID(i)
                end if
            end do
            allocate(keyPointBID(maxKeyPointID),stat=error)
            if (error/=0) then
                write(*,"(A)") 'generateVertexID(): 内存分配错误。'
            end if
            allocate(bid2VertexSubscript(ng),stat=error)
            if (error/=0) then
                write(*,"(A)") 'generateVertexID(): 内存分配错误。'
            end if
            
            ! 构建keyPointBID  记录keyPoint编号对应的边界点在geometry.b中的下标
            ! keyPointBID(i)=0 表示不存在编号为i的keyPoint
            keyPointBID=0
            do i=1,ng
                if (geometry.keyPointID(i)>0) then
                    keyPointBID(geometry.keyPointID(i))=i
                end if
            end do
            ! 构建bid2VertexSubscript 记录bid与vertex下标的对应关系
            ! bid2VertexSubscript(i)=0 表示i这个边界点不是格点
            bid2VertexSubscript=0
            do i=1,meshInfo.nv
                if (meshInfo.vertex(i).isB) then
                    bid2VertexSubscript(meshInfo.vertex(i).bid)=i
                end if
            end do
        end if
        
        ! 编号smallEndArea包含的边界点
        smallEnd=0
        if (present(smallEndArea)) then
            n=size(smallEndArea,2)
            do i=1,n
                b=smallEndArea(1,i)
                e=smallEndArea(2,i)
                ! 检查输入有效性
                if (max(b,e)>maxKeyPointID .or. min(b,e)<1) then
                    write(*,'("generateVertexID(): keyPoint ",I0," 或 ",I0," 不存在。")') b, e
                    cycle
                end if
                if (keyPointBID(b)==0 .or. keyPointBID(e)==0) then
                    write(*,'("generateVertexID(): keyPoint ",I0," 或 ",I0," 不存在。")') b, e
                    cycle
                end if
                ! 获取区域起/终点
                b=keyPointBID(b)
                e=keyPointBID(e)
                ! 编号
                j=b
                if (b==e .or. b==mod(e,ng)+1) then
                    firstLoop=.true.
                else
                    firstLoop=.false.
                end if
                e=mod(e,ng)+1
                do while (j/=e .or. firstLoop)
                    k=bid2VertexSubscript(j)
                    if (k>0) then
                        if (meshInfo.vertex(k).id==-1) then
                            smallEnd=smallEnd+1
                            meshInfo.vertex(k).id=smallEnd
                        end if
                    end if
                    if (j==e) then
                        firstLoop=.false.
                    end if
                    j=mod(j,ng)+1
                end do
            end do
        end if
        
        ! 编号bigEndArea包含的边界点
        bigEnd=meshInfo.nv+1
        if (present(bigEndArea)) then
            n=size(bigEndArea,2)
            do i=n,1,-1
                b=bigEndArea(1,i)
                e=bigEndArea(2,i)
                ! 检查输入有效性
                if (max(b,e)>maxKeyPointID .or. min(b,e)<1) then
                    write(*,'("generateVertexID(): keyPoint ",I0," 或 ",I0," 不存在。")') b, e
                    cycle
                end if
                if (keyPointBID(b)==0 .or. keyPointBID(e)==0) then
                    write(*,'("generateVertexID(): keyPoint ",I0," 或 ",I0," 不存在。")') b, e
                    cycle
                end if
                ! 获取区域起/终点
                b=keyPointBID(b)
                e=keyPointBID(e)
                ! 编号
                j=e
                if (b==e .or. b==mod(e,ng)+1) then
                    firstLoop=.true.
                else
                    firstLoop=.false.
                end if
                b=mod(b-2+ng,ng)+1
                do while (j/=b .or. firstLoop)
                    k=bid2VertexSubscript(j)
                    if (k>0) then
                        if (meshInfo.vertex(k).id==-1) then
                            bigEnd=bigEnd-1
                            meshInfo.vertex(k).id=bigEnd
                        end if
                    end if
                    if (j==b) then
                        firstLoop=.false.
                    end if
                    j=mod(j-2+ng,ng)+1
                end do
            end do
        end if
        
        ! 给剩余的点编号
        do i=1,meshInfo.nv
            if (meshInfo.vertex(i).id==-1) then
                smallEnd=smallEnd+1
                meshInfo.vertex(i).id=smallEnd
            end if
        end do
        
        ! 检查编号有效性
        if (smallEnd+1/=bigEnd) then
            write (*,"(A)") "generateVertexID(): 编号错误，总格点数不符。"
        end if
        
        return        
    end subroutine generateVertexID

! outputGeometry(): 输出geometry内的信息至文件
! 将geometry.b的内容输出到fileID所对应的文件中
    subroutine outputGeometry(fileID)
        use mesh 
        implicit none
        integer, intent(in) :: fileID
        integer :: i
        
        write(fileID,"('n=',I)") geometry.n
        write(fileID,"(A23,A23,A)") ' x ',' y ',' keyPointID '
        do i=1,geometry.n+1
             write(fileID,"(E23.16,E23.16,I)") geometry.b(i),geometry.keyPointID(i)
        end do
    end subroutine outputGeometry

! outputMesh(): 将meshInfo中的网格信息输出至fileID对应的文件
    subroutine outputMesh(fileID, tag)
        use mesh 
        implicit none
        integer, intent(in) :: fileID
        character(*), intent(in),optional :: tag
        integer :: i
        
        if (present(tag)) then
            write(fileID, "(A)") trim(tag)
        end if
        ! 输出vertex
        write(fileID,"('N_Vertex=',I11)") meshInfo.nv
        write(fileID,"(A23,A23,A)") '  x   ','  y   ', '  external_ID  '
        do i=1,meshInfo.nv
            write(fileID,"(E23.16,E23.16,I11)") meshInfo.vertex(i).p,meshInfo.vertex(i).id
        end do
        write(fileID,*)
        ! 输出edge
        write(fileID,"('N_Edge=',I11)") meshInfo.ne
        write(fileID,"(A11,A11)") '  v(1)  ', ' v(2)  '
        do i=1,meshInfo.ne
            write(fileID,"(I11,I11)") meshInfo.edge(i).v
        end do        
        write(fileID,*)
        ! 输出surface
        write(fileID,"('N_Surface=',I11)") meshInfo.ns
        write(fileID,"(6A11)") ' v(1) ',' v(2) ',' v(3) ',' e(1) ',' e(2) ',' e(3) '
        do i=1,meshInfo.ns
            write(fileID,"(6I11)") meshInfo.surface(i).v,meshInfo.surface(i).e
        end do
        write(fileID,*)
        ! 输出boundary
        write(fileID,"('N_Boundary=',I11)") meshInfo.nb
        write(fileID,"(A11)") ' edge_No '
        write(fileID,"(I11)") meshInfo.boundary(1:meshInfo.nb)
        write(fileID,*)
    end subroutine outputMesh
        
! arclen(): 返回边界边对应的实际边界长度
!--------------------------------------------------------------
![Input]
!  ib  - 边界边在meshInfo.boundary中的下标
![Return]
!    返回指定边界边对应的实际边界的曲线弧长
!--------------------------------------------------------------
    real(8) elemental function arclen(ib) result(l)
        use mesh
        implicit none
        integer, intent(in) :: ib
        integer             :: ie, s, t
        
        ie=meshInfo.boundary(ib)
        if ( meshInfo.edge(ie).s(1)==-1 ) then  ! 求解域在v(2)-v(1)的左手边
            s=meshInfo.edge(ie).v(1)
            t=meshInfo.edge(ie).v(2)
        else
            s=meshInfo.edge(ie).v(2)
            t=meshInfo.edge(ie).v(1)
        end if
        s=meshInfo.vertex(s).bid
        t=meshInfo.vertex(t).bid
        
        if (t>s) then
            l=geometry.l(t-1)-geometry.l(mod(s-2+geometry.n,geometry.n)+1)
        else
            l=geometry.l(geometry.n)-( geometry.l(s-1)-geometry.l(mod(t-2+geometry.n,geometry.n)+1) )
        end if
    end function arclen

! inArea_point(): 判断一个点是否在求解域内
!--------------------------------------------------------------
!  不在：      return 0
!  在边界上：  return 2
!  在求解域内：return 1
!--------------------------------------------------------------
    integer elemental function inArea_point(p0) result(re)
        use mesh
        implicit none
        type(point_t), intent(in) :: p0
        integer                   :: i, n
        real(8)                   :: angleCount
        type(point_t)             :: v0, v1
            
        angleCount=0d0
        n=geometry.n
        v0=geometry.b(1)-p0
        do i=2,n+1
            v1=geometry.b(i)-p0
            if (abs(v0*v1/norm(v0)/norm(v1))<eps) then
                re=2     ! On the boundary
                return
            end if
            if (v0*v1>0) then
                angleCount=angleCount+dacos((v0**v1)/norm(v0)/norm(v1))
            else
                angleCount=angleCount-dacos((v0**v1)/norm(v0)/norm(v1))
            end if
            v0=v1
        end do
            
        if (abs(angleCount-2*pi)>1d-6) then
            re=0     ! Outside the area
        else
            re=1     ! Inside the area
        end if
    end function inArea_point
        
    
! inArea_edge(): 判断一个线段(除两端点外)是否完全在求解域内
!--------------------------------------------------------------
![Input]
! p1, p2 - 线段的两个端点，其中p1必须保证在求解域内或边界上。
![Return]
!  完全在求解域内：      return -1
!  不完全在求解域内：    return 与线段相交的边界的端点中离p1最
!                        近的点在geometry.b中的下标。
!--------------------------------------------------------------
    integer elemental function inArea_edge(p1, p2) result(iminD)
        use mesh
        implicit none
        type(point_t), intent(in) :: p1, p2
        integer                   :: i, n
        real(8)                   :: minD, D, normv0_2
        type(point_t)             :: v0
        
        v0=p2-p1
        normv0_2=norm(v0)**2
        minD=2d0
        iminD=-1
        n=geometry.n
        do i=1,n
            if ( ((geometry.b( i )-p1)*v0) * ((geometry.b(i+1)-p1)*v0) &
                < -eps*normv0_2*norm(geometry.b( i )-p1)*norm(geometry.b(i+1)-p1) &
           .and. ((p1-geometry.b(i))*(geometry.b(i+1)-geometry.b(i))) * ((p2-geometry.b(i))*(geometry.b(i+1)-geometry.b(i))) &
                < -eps*norm(geometry.b(i+1)-geometry.b(i))**2*norm(p2-geometry.b(i))*norm(p1-geometry.b(i)))then
                D=min(norm(geometry.b(i)-p1),norm(geometry.b(i+1)-p1))/norm(v0)
                if (minD>D) then
                    minD=D
                    iminD=i
                end if
            end if
        end do
    end function inArea_edge

! getNextBpoint(): 在边界上取下一个点(沿边界正方向)
!----------------------------------------------------------------------------------
! 会优先返回两倍(最近点距离)范围内的keyPoint
!  [Input]
!    i0 - 当前点在geometry.b中的下标 (i0<=geometry.n)
!    minR, maxR 
!       - [可选]下一点与当前点的距离区间(直线距离)
!  iend - [可选]下一点不能超过这个点(geometry.b中的下标) (iend<=geometry.n)
!  [Return]
!      返回下一点在geometry.b中的下标(返回值<=geometry.n)，
!      若找不到满足条件的点，则返回-1。
!----------------------------------------------------------------------------------
    integer elemental function getNextBpoint(i0, minR, maxR, iend)  result(nexti)
        use mesh
        implicit none
        integer, intent(in)       :: i0
        real(8), intent(in), optional  :: minR, maxR
        integer, intent(in), optional  :: iend
        integer                   :: i,n,m
        real(8)                   :: Rub,Rlb,dis
        
        ! 根据输入参数设定筛选条件
        n=geometry.n
        if (present(iend)) then
            m=iend
        else
            m=i0
        end if
        if (present(minR)) then
            Rlb=minR
        else
            Rlb=0d0
        end if
        if (present(maxR)) then
            Rub=maxR
        else
            Rub=geometry.l(n)
        end if
        ! 寻找距离合适的点
        nexti=-1
        i=mod(i0,n)+1
        do while (i/=m)
            dis=norm(geometry.b(i)-geometry.b(i0))
            if (dis>Rub) then
                exit
            end if
            if (dis>Rlb) then
                if (geometry.keyPointID(i)>0) then
                    nexti=i  ! 遇到keyPoint直接返回
                    return
                else if (nexti==-1) then
                    nexti=i  ! 遇到非keyPoint则记下最先找到的
                    Rub=min(Rub,2d0*dis)
                end if
            end if
            i=i+1
            if (i==n+1) then
                i=1
            end if
        end do
        return
    end function getNextBpoint
        
!addVertex(): 添加格点
!---------------------------------------------------------------
![Input]
!  p0     - 需要添加的格点坐标
!  isB    - [可选]是否为边界格点
!  bid    - [可选]如果是边界格点，则为该点在geometry.b中的下标
![Return]
!    -1   - 已达到最大格点数，无法添加
!    -2   - 未提供bid
!格点下标 - 成功添加
!---------------------------------------------------------------
    integer function addVertex(p0, isB, bid)  result(re)
        use mesh
        implicit none
        type(point_t), intent(in)     :: p0
        logical, intent(in), optional :: isB
        integer, intent(in), optional :: bid
        logical                       :: isB0
        integer   :: nv
            
        if (.not. present(isB)) then
            isB0=.false.
        else
            isB0=isB
        end if
        if (isB0 .and. (.not. present(bid))) then
            re=-2
            return
        end if
        nv=meshInfo.nv+1
        if (nv > meshInfo.maxN) then
            call expandMeshInfo(mv=1.1d0)
            if (nv > meshInfo.maxN) then
                re=-1
                return
            end if
        end if
            
        meshInfo.vertex(nv).p=p0
        meshInfo.vertex(nv).isB=isB0
        if (isB0) then
            meshInfo.vertex(nv).bid=bid
        end if
        re=nv
        meshInfo.nv=nv
        return
    end function addVertex
        
!addEdge(): 添加边（自动识别并设置边界）
!---------------------------------------------------------------
![Input]
!  a,b    - 边的两个端点对应的vertex下标，a -> v(1), b -> v(2)
!  isB    - [可选]是否为边界边，默认为否
![Return]
!    -1   - 已达到最大边数，无法添加
!  边下标 - 成功添加
!---------------------------------------------------------------
    integer function addEdge(a, b, isB)  result(re)
        use mesh
        implicit none
        integer, intent(in)           :: a,b
        logical, intent(in), optional :: isB
        integer                       :: ne, nb
            
        ne=meshInfo.ne+1
        if (ne>size(meshInfo.edge)) then
            call expandMeshInfo(me=1.1d0)
            if (ne>size(meshInfo.edge)) then
                re=-1
                return
            end if
        end if
            
        if (present(isB) .and. isB) then
            nb=meshInfo.nb+1
            if (nb>size(meshInfo.boundary)) then
                call expandMeshInfo(mb=1.1d0)
                if (nb>size(meshInfo.boundary)) then
                    re=-1
                    return
                end if
            end if
            meshInfo.boundary(nb)=ne
            meshInfo.nb=nb
            meshInfo.edge(ne).isB=.true.
        end if
            
        meshInfo.edge(ne).v(1)=a
        meshInfo.edge(ne).v(2)=b
        meshInfo.edge(ne).len=norm(meshInfo.vertex(a).p-meshInfo.vertex(b).p)
        meshInfo.ne=ne
        re=ne
        return
    end function addEdge
        
!addSurface(): 添加面元
!---------------------------------------------------------------
![Input]
!  a,b,c  - 面元的三个顶点在vertex中的下标
!           a -> v(1), b -> v(2), c -> v(3)
!           e(1): b <--> c  e(2): c <--> a  e(3): a <--> b
!  x,y,z  - [可选]面元的三条边在edge中的下标，用于加速寻找边的过程。 
!           如果指定的边的端点不是a,b,c，程序将自动忽略该边。
![Return]
!    -1   - 已达到最面元数，无法添加
!面元下标 - 成功添加
!---------------------------------------------------------------
    integer function addSurface(a, b, c, x, y, z)  result(re)
        use mesh
        implicit none
        integer, intent(in)            :: a, b, c
        integer, intent(in), optional  :: x, y, z
        integer   :: ns, i, j, counter
        type(point_t) :: p(3)
            
        ns=meshInfo.ns+1
        if (ns>size(meshInfo.surface)) then
            call expandMeshInfo(ms=1.1d0)
            if (ns>size(meshInfo.surface)) then
                re=-1
                return
            end if
        end if
        ! 设置顶点
        meshInfo.surface(ns).v(1)=a
        meshInfo.surface(ns).v(2)=b
        meshInfo.surface(ns).v(3)=c
        ! 找/设置边
        counter=0
        if (present(x)) call testEdge(x,counter)
        if (present(y)) call testEdge(y,counter)
        if (present(z)) call testEdge(z,counter)
        if (counter<3) then
            do i=meshInfo.ne,1,-1
                if (meshInfo.edge(i).v(1)==a .or. &
                    meshInfo.edge(i).v(1)==b .or. &
                    meshInfo.edge(i).v(1)==c) then
                    call testEdge(i,counter)
                    if (counter==3) then
                        exit
                    end if
                end if
            end do
        end if
        ! 更新边两侧的面元 edge.s
        do i=1,3
            j=meshInfo.surface(ns).e(i)
            p(1:2)=meshInfo.vertex(meshInfo.edge(j).v(1:2)).p
            p(3)=meshInfo.vertex(meshInfo.surface(ns).v(i)).p
            if ( (p(3)-p(1))*(p(2)-p(1))>0 ) then  ! 面元在v(2)-v(1)的右边
                meshInfo.edge(j).s(1)=ns
            else
                meshInfo.edge(j).s(2)=ns
            end if
        end do
        ! 计算面积
        p=meshInfo.vertex(meshInfo.surface(ns).v(:)).p
        meshInfo.surface(ns).S=calcArea(p(1),p(2),p(3))
        ! 计算外接圆半径
        meshInfo.surface(ns).R=calcRadius(p(1),p(2),p(3),meshInfo.surface(ns).S)
        ! 返回
        meshInfo.ns=ns
        re=ns
        return
    contains
        subroutine testEdge(i, counter)
            implicit none
            integer, intent(in)    :: i
            integer, intent(inout) :: counter
            
            if      ( (meshInfo.edge(i).v(1)==b .and. meshInfo.edge(i).v(2)==c) .or. &
                        (meshInfo.edge(i).v(2)==b .and. meshInfo.edge(i).v(1)==c) ) then
                meshInfo.surface(ns).e(1)=i
                counter=counter+1
            else if ( (meshInfo.edge(i).v(1)==c .and. meshInfo.edge(i).v(2)==a) .or. &
                        (meshInfo.edge(i).v(2)==c .and. meshInfo.edge(i).v(1)==a) ) then
                meshInfo.surface(ns).e(2)=i
                counter=counter+1
            else if ( (meshInfo.edge(i).v(1)==a .and. meshInfo.edge(i).v(2)==b) .or. &
                        (meshInfo.edge(i).v(2)==a .and. meshInfo.edge(i).v(1)==b) ) then
                meshInfo.surface(ns).e(3)=i
                counter=counter+1
            end if
        end subroutine testEdge
    end function addSurface
        
! expandMeshInfo(): 增加meshInfo的大小，若增加失败，则会尝试回滚。
!---------------------------------------------------------------
![Input]
! mv, me, ms, mb  -  [可选]vertex, edge, surface, boundary数组
!                    尺寸的放大倍数，必须大于1。对于小于1的参数
!                    将自动忽略，不会报错。
!---------------------------------------------------------------
    subroutine expandMeshInfo(mv, me, ms, mb)
        use mesh
        implicit none
        real(8), intent(in), optional  :: mv, me, ms, mb
        type(meshInfo_t)              :: tmp
        integer                     :: n0, n, error
            
        if (present(mv) .and. mv>1d0) then
            n0=size(meshInfo.vertex)
            allocate(tmp.vertex(n0),stat=error)
            if (error==0) then
                tmp.vertex=meshInfo.vertex
                deallocate(meshInfo.vertex)
                n=int(n0*mv)
                allocate(meshInfo.vertex(n),stat=error)
                if (error/=0) then
                    allocate(meshInfo.vertex(n0),stat=error)
                end if
                meshInfo.vertex(1:n0)=tmp.vertex
                meshInfo.maxN=size(meshInfo.vertex)
                deallocate(tmp.vertex)
            end if
        end if
            
        if (present(me) .and. me>1d0) then
            n0=size(meshInfo.edge)
            allocate(tmp.edge(n0),stat=error)
            if (error==0) then
                tmp.edge=meshInfo.edge
                deallocate(meshInfo.edge)
                n=int(n0*me)
                allocate(meshInfo.edge(n),stat=error)
                if (error/=0) then
                    allocate(meshInfo.edge(n0),stat=error)
                end if
                meshInfo.edge(1:n0)=tmp.edge
                deallocate(tmp.edge)
            end if
        end if
                        
        if (present(ms) .and. ms>1d0) then
            n0=size(meshInfo.surface)
            allocate(tmp.surface(n0),stat=error)
            if (error==0) then
                tmp.surface=meshInfo.surface
                deallocate(meshInfo.surface)
                n=int(n0*ms)
                allocate(meshInfo.surface(n),stat=error)
                if (error/=0) then
                    allocate(meshInfo.surface(n0),stat=error)
                end if
                meshInfo.surface(1:n0)=tmp.surface
                deallocate(tmp.surface)
            end if
        end if
            
        if (present(mb) .and. mb>1d0) then
            n0=size(meshInfo.boundary)
            allocate(tmp.boundary(n0),stat=error)
            if (error==0) then
                tmp.boundary=meshInfo.boundary
                deallocate(meshInfo.boundary)
                n=int(n0*mb)
                allocate(meshInfo.boundary(n),stat=error)
                if (error/=0) then
                    allocate(meshInfo.boundary(n0),stat=error)
                end if
                meshInfo.boundary(1:n0)=tmp.boundary
                deallocate(tmp.boundary)
            end if
        end if
            
    end subroutine expandMeshInfo

! calcArea_p(): 返回由p1,p2,p3构成的三角形的面积
    real(8) elemental function calcArea_p(p1,p2,p3) result(S)
        use vector2D
        implicit none
        type(point_t), intent(in) :: p1,p2,p3
        
        S=abs((p2-p1)*(p3-p1)/2d0)
    end function calcArea_p
! calcArea_s(): 返回sID对应面元的面积
    real(8) elemental function calcArea_s(sID) result(S)
        use mesh 
        implicit none
        integer,intent(in) :: sID
        type(point_t) :: p(3)
        
        p=meshInfo.vertex(meshInfo.surface(sID).v(:)).p
        S=abs((p(2)-p(1))*(p(3)-p(1))/2d0)
    end function calcArea_s
    
! calcRadius_p(): 返回由p1,p2,p3构成的三角形的外接圆半径,可输入三角型面积加速计算
    real(8) elemental function calcRadius_p(p1,p2,p3,S) result(R)
        use vector2D
        implicit none
        type(point_t), intent(in)   :: p1,p2,p3
        real(8),intent(in),optional :: S
        real(8)                     :: S0
        
        if (present(S)) then
            S0=S
        else
            S0=abs((p2-p1)*(p3-p1)/2d0)
        end if
        R=norm(p1-p2)*norm(p2-p3)*norm(p3-p1)/(4d0*S0)
    end function calcRadius_p
! calcRadius_s(): 返回sID对应面元的外接圆半径,可输入三角型面积加速计算
    real(8) elemental function calcRadius_s(sID,S) result(R)
        use mesh 
        implicit none
        integer, intent(in) :: sID
        real(8),intent(in),optional :: S
        type(point_t)       :: p(3)
        
        p=meshInfo.vertex(meshInfo.surface(sID).v(:)).p
        R=calcRadius(p(1),p(2),p(3),S)
    end function calcRadius_s

! getShapeCoeff(): 计算形状因子 S/(pi*R^2)
    real(8) elemental function getShapeCoeff(p1,p2,p3) result(coeff)
        use mesh
        implicit none
        type(point_t),intent(in) :: p1, p2, p3
        real(8)       :: R, S
        S=calcArea(p1,p2,p3)
        R=calcRadius(p1,p2,p3,S)
        coeff=S/(pi*R**2)
    end function getShapeCoeff
    
! makeGeometry(): 离散化边界、初始化geometry变量
!--------------------------------------------------------------------
! geometry.keyPointID的取值：
!   不是keyPoint：0
!   是keyPoint：  该点在x、y数组中的下标
![Input]
!  step  - 几何精度，两点之间的距离不小于step/4, 不大于step。
! keyPoints - [可选]关键点在x、y中的下标，关键点一定会成为边界格点。
!             强烈建议将不同类型边界的分界点指定为关键点。
!--------------------------------------------------------------------
    subroutine makeGeometry(x,y,step,keyPoints)
        use mesh
        implicit none
        real(8), intent(in)  :: x(:), y(:), step
        integer, intent(in),optional  :: keyPoints(:)
        real(8), allocatable :: dtmp(:)
        integer, allocatable :: kPs(:), inputKP(:)
        integer              :: nKP
        type(point_t), allocatable :: tmpP(:), p(:)
        integer, allocatable :: isKey(:)
        type(point_t)        :: v0
        integer              :: i, j, k, n0, n, m
        real(8)              :: tmpLen
        
        geometry.steplb=step/4d0
        geometry.stepub=step
        
        ! 计算边界总长度
        n=size(x)
        allocate(dtmp(n))
        dtmp(1:n-1)=(x(2:n)-x(1:n-1))**2+(y(2:n)-y(1:n-1))**2
        dtmp(n)=(x(1)-x(n))**2+(y(1)-y(n))**2
        tmpLen=sum(sqrt(dtmp))
        deallocate(dtmp)
        
        ! 估算边界点数目
        n0=int(4.1*tmpLen/step)+100
        
        ! 分配内存
        allocate(p(n+1))
        allocate(kPs(n+1))
        allocate(tmpP(n0))
        allocate(isKey(n0))
        
        ! 复制点
        p(1:n).x=x(:)
        p(1:n).y=y(:)
        p(n+1).x=x(1)
        p(n+1).y=y(1)
        
        ! 构建keyPoints
        if (present(keyPoints)) then
            m=size(keyPoints)
            allocate(inputKP(m))
            inputKP=keyPoints
        else
            m=1
            allocate(inputKP(m))
            inputKP(1)=0
        end if
        nKP=0
        j=1
        if (inputKP(j)==1) then
            nKP=nKP+1
            kPs(nKP)=1
            if (m>j) j=j+1
        else if ( (p(n)-p(1))**(p(2)-p(1))/norm(p(n)-p(1))/norm(p(2)-p(1)) &
                  > dcos(pi*150/180)) then   ! 偏转角大于30度时设为关键点
            nKP=nKP+1
            kPs(nKP)=1
        end if
        do i=2,n
            if (inputKP(j)==i) then
                nKP=nKP+1
                kPs(nKP)=i
                if (m>j) j=j+1
            else if ( (p(i-1)-p(i))**(p(i+1)-p(i))/norm(p(i-1)-p(i))/norm(p(i+1)-p(i)) &
                      > dcos(pi*150/180)) then   ! 偏转角大于30度时设为关键点
                nKP=nKP+1
                kPs(nKP)=i
            end if
        end do
        if (kPs(1)==1) then
            nKP=nKP+1
            kPs(nKP)=n+1
        end if
        
        ! 采样
        isKey=0
        n=1
        k=1
        tmpP(n)=p(1)
        if (k<=nKP .and. kPs(k)==1) then
            isKey(n)=1
            if (k<nKP) k=k+1
        end if
        do i=2,size(p)
            if (norm(p(i)-tmpP(n)) > step) then
                v0=p(i)-p(i-1)
                v0=v0/norm(v0)
                n=n+1
                tmpP(n)=p(i-1)+step/4d0*v0
                
                tmpLen=norm(p(i)-tmpP(n))
                m=ceiling(tmpLen/step)
                v0=v0*(tmpLen/m)
                
                do j=1,m-1
                    n=n+1
                    tmpP(n)=tmpP(n-1)+v0
                end do
                n=n+1
                tmpP(n)=p(i)
                if (kPs(k)==i) then
                    isKey(n)=i
                    if (k<nKP) k=k+1
                end if
            else if (norm(p(i)-tmpP(n)) > step/4d0 ) then
                n=n+1
                tmpP(n)=p(i)
                if (kPs(k)==i) then
                    isKey(n)=i
                    if (k<nKP) k=k+1
                end if
            end if
        end do
        if (isKey(n)==size(p)) then
            isKey(n)=1
        end if
        
        ! 分配内存、复制点
        geometry.n=n-1
        allocate(geometry.b(n))
        allocate(geometry.l(n-1))
        allocate(geometry.keyPointID(n))
        geometry.b=tmpP(1:n)
        geometry.keyPointID=isKey(1:n)
        
        ! 构建geometry.bidOfKeyPoint
        allocate(geometry.bidOfKeyPoint(maxval(isKey)))
        geometry.bidOfKeyPoint=-1
        do i=1,geometry.n
            if (geometry.keyPointID(i)>0) then
                geometry.bidOfKeyPoint(geometry.keyPointID(i))=i
            end if
        end do
        
        ! 计算弧长
        geometry.l(1)=norm(geometry.b(2)-geometry.b(1))
        do i=2,n-1
            geometry.l(i)=geometry.l(i-1)+norm(geometry.b(i+1)-geometry.b(i))
        end do
        
        deallocate(tmpP)
        deallocate(p)
        deallocate(kPs)
        deallocate(isKey)
        
        return
    end subroutine makeGeometry
    
! refineMesh(): 细化网格
!----------------------------------------------------------------------------------
!  只能细化网格，不能使网格变粗
!  通过裂边的方式细化，不会移动现有格点
!  通过形状判定，限制钝角的出现
!  边界处面元的最大边长在 minR/2 ~ minR 之间
!  区域中间面元边长在 maxR/2 ~ maxR 之间
!  形状不好的面元尺寸会较小
![Input]
! minR, maxR - 面元的尺寸
! outputDetails - [可选]指示是否需要将每一步迭代过后的网格状态输出至文件mesh.txt
!                 默认为否
!----------------------------------------------------------------------------------
    subroutine refineMesh(minR, maxR, outputDetails)
        use mesh
        implicit none
        real(8), intent(in) :: minR, maxR
        logical, intent(in), optional :: outputDetails
        real(8)             :: actual_minR, actual_maxR, current_minR, current_maxR
        integer             :: i, j, s1, s2, ns, nb, ne
        logical             :: changed, relaxing
        character(100)      :: tag
        
        ! 输出初始网格信息
        if (present(outputDetails) .and. outputDetails) then
            open(80,file='geometry.txt')
            call outputGeometry(80)
            close(80)
        
        
            open(81,file='mesh.txt')
            j=0
            write(tag,"('[迭代=',I0,']')") j
            call outputMesh(81,tag)
        end if
        
        changed=.false.
        relaxing=.true.
        do while (changed .or. relaxing)
            if (changed) then
                relaxing=.true.
            else
                relaxing=.false.
            end if
            changed=.false.
            ne=meshInfo.ne
            actual_maxR=maxval(meshInfo.edge(1:ne).len)
            current_minR=minR
            current_maxR=(actual_maxR+3*maxR)/4d0
            
            ! 优先加密边界(目标: minR)
            ! relaxing时不加密边界
            if (.not. relaxing) then
                i=1
                nb=meshInfo.nb
                do while (i<=nb)
                    if (meshInfo.edge(meshInfo.boundary(i)).len > current_minR &
                  .or. arclen(i) > current_maxR) then
                        changed=changed .or. splitBorder(i, current_minR)
                    end if
                    i=i+1
                end do
            end if

            ! output mesh
            if (present(outputDetails) .and. outputDetails) then
                j=j+1
                write(tag,"('[iter=',I0,']')") j
                call outputMesh(81,tag)
            end if
            
            ! 分割形状不好的面元，边界边不分割(限制: minR)
            ! 边是否在边界的判定在splitSurface()内判定
            !ns=meshInfo.ns
            i=1
            do while (i<=meshInfo.ns)
                if (maxval(meshInfo.edge(meshInfo.surface(i).e).len) > current_minR &
              .and. meshInfo.surface(i).S/(pi*meshInfo.surface(i).R**2) < badShapeCoeff) then
                    changed=changed .or. splitSurface(i, current_minR/2d0, .true.)
                end if
                i=i+1
            end do
            
            ! output mesh
            if (present(outputDetails) .and. outputDetails) then
                write(tag,"('[iter=',I0,']')") j 
                call outputMesh(81,tag)
            end if
                
            ! 分割过长的边(目标: maxR)，边界边不分割
            !ne=meshInfo.ne
            i=1
            do while (i<=meshInfo.ne)
                if ((.not. meshInfo.edge(i).isB) &
                     .and. meshInfo.edge(i).len > current_maxR) then
                    s1=meshInfo.edge(i).s(1)
                    s2=meshInfo.edge(i).s(2)
                    ! 选择半径较大的面元进行分割
                    if (s1==-1) then
                        s1=s2
                    else if (s2/=-1 .and. &
                       meshInfo.surface(s2).R > meshInfo.surface(s1).R) then
                        s1=s2
                    end if
                    changed=changed .or. splitSurface(s1, current_maxR/2d0, .false.)
                end if
                i=i+1
            end do
            
            ! output mesh
            if (present(outputDetails) .and. outputDetails) then
                write(tag,"('[iter=',I0,']')") j
                call outputMesh(81,tag)
            end if
        end do
        
        ! 检查keyPoint是否全部设为格点
        call checkKeyPoint()
        
        ! output mesh
        if (present(outputDetails) .and. outputDetails) then
            write(tag,"(A)") "[Final]"
            call outputMesh(81,tag)
            close(81)
        end if
        
    end subroutine refineMesh    

! findBestBreakPoint(): 用三分法找出形状因子最好的分割点
!------------------------------------------------------------------ 
![Input]
! v - 长度为4的数组，记录四个顶点在meshInfo.vertex中的下标
!     v(1)-v(2)为待分割的边，v(3), v(4)为与该边相邻的面元的剩余顶点
! p - [与v选一个]长度为4的数组，记录四个顶点的坐标
! minlen - 分割后，新的边长不应小于这个值
!checkShape - 指示是否需要检查分割后的形状因子
![Return] 
! 返回分割点到v(1)点的距离
! 如果极大值点超出minR的约束，则取边界点
! 如果分割后的形状因子非常小，则返回一个小于零的值
!------------------------------------------------------------------
    real(8) pure function findBestBreakPoint_v(v,minlen,checkShape) result(pos)
        use mesh 
        implicit none
        integer, intent(in) :: v(4)
        real(8), intent(in) :: minlen
        logical, intent(in) :: checkShape
        type(point_t) :: p(4)
        
        p=meshInfo.vertex(v(:)).p
        pos=findBestBreakPoint(p,minLen,checkShape)
    end function findBestBreakPoint_v
    
    real(8) pure function findBestBreakPoint_p(p,minlen,checkShape) result(pos)
        use mesh 
        implicit none
        type(point_t), intent(in) :: p(4)
        real(8), intent(in) :: minlen
        logical, intent(in) :: checkShape
        real(8)  :: L,R,mL,mR,minL,minR, edgeLen
        type(point_t) :: bp
        
        ! iteration: find local maximum
        L=0d0
        R=1d0
        do while (R-L<1d-4)
            mL=(2d0*L+R)/3d0
            mR=(L+2d0*R)/3d0
            bp=p(1)+(p(2)-p(1))*mL
            minL=getShapeCoeff(bp,p(1),p(4))
            minL=min(minL,getShapeCoeff(bp,p(2),p(4)))
            minL=min(minL,getShapeCoeff(bp,p(1),p(3)))
            minL=min(minL,getShapeCoeff(bp,p(2),p(3)))
            bp=p(1)+(p(2)-p(1))*mR
            minR=getShapeCoeff(bp,p(1),p(4))
            minR=min(minR,getShapeCoeff(bp,p(2),p(4)))
            minR=min(minR,getShapeCoeff(bp,p(1),p(3)))
            minR=min(minR,getShapeCoeff(bp,p(2),p(3)))
            if (minL>minR) then
                R=mR
            else
                L=mL
            end if
        end do
        edgeLen=norm(p(1)-p(2))
        pos=(R+L)/2d0*edgeLen
        ! check bound
        if (edgeLen>2*minlen) then
            if (pos<minlen) then
                pos=minlen
            else if (pos>edgeLen-minlen) then
                pos=edgeLen-minlen
            end if
        end if
        ! check shape coefficient
        if (checkShape) then
            bp=p(1)+(p(2)-p(1))*pos/edgeLen
            minL=getShapeCoeff(bp,p(1),p(4))
            minL=min(minL,getShapeCoeff(bp,p(2),p(4)))
            minL=min(minL,getShapeCoeff(bp,p(1),p(3)))
            minL=min(minL,getShapeCoeff(bp,p(2),p(3)))
            if (minL<badShapeCoeff) then
                if (minL<min(getShapeCoeff(p(1),p(2),p(4)),getShapeCoeff(p(1),p(2),p(3)))) then
                    pos=-1d0
                    return
                end if
            end if
        end if
    end function findBestBreakPoint_p
    
    
! splitBorder(): 分割边界边
!--------------------------------------------------------------
! 在geometry中被定为keyPoint的点有很高概率会成为格点
! 可以直接调用此函数分割边界，这样做并不会破坏meshInfo的结构。
! 分割方式见splitBorder.jpg
![Input]
!   id - 要分割的边在meshInfo.boundary中的下标
! minR - 寻找分割点时采用的最小步长
![Return]
!  成功分割   - .true.
!  未成功分割 - .false.
!--------------------------------------------------------------
    logical function splitBorder(id, minR) result(re)
        use mesh
        implicit none
        integer, intent(in)  :: id
        real(8), intent(in)  :: minR
        type(point_t)        :: breakPoint, p(3)
        type(edge_t),pointer :: oldE
        type(surface_t),pointer :: oldS
        integer              :: i, j, k, ibestBreakPoint, v0, v20
        integer              :: e1, e2, ioldS, e30
        real(8)              :: bestShape, S1,S2,R1,R2
        
        ! 找出待分割的面元
        oldE=>meshInfo.edge(meshInfo.boundary(id))
        if (oldE.s(1)==-1) then
            ioldS=oldE.s(2)
        else
            ioldS=oldE.s(1)
        end if
        oldS=>meshInfo.surface(ioldS)
        
        ! 找分割点
        p(:)=meshInfo.vertex(oldS.v(:)).p
        i=meshInfo.vertex(oldS.v(2)).bid
        k=meshInfo.vertex(oldS.v(3)).bid
        ibestBreakPoint=-1
        bestShape=0d0
        i=getNextBpoint(i,minR=minR,iend=k)
        do while (i/=-1)
            ! 如果下一点距离边界的另一端距离过小，且不是keyPoint，则结束找点
            if (norm(geometry.b(i)-geometry.b(k))<minR .and. geometry.keyPointID(i)==0) then
                exit
            end if
            
            ! v(1)点与分割点的连线是否与边界相交
            j=inArea(p(1),geometry.b(i))
            ! 如果与边界相交，则选择交点为分割点
            if ( j==-1 ) then
                j=i
            end if
            breakPoint=geometry.b(j)
            
            ! 判断分割边是否在原三角形面元内部
            if ( ((p(2)-p(1))*(breakPoint-p(1)))*  &
                 ((p(3)-p(1))*(breakPoint-p(1))) < &
                  -eps*norm(p(2)-p(1))*norm(p(3)-p(1))*norm(breakPoint-p(1))**2 ) then
            
                if (geometry.keyPointID(j)>0) then   ! 遇到keyPoint直接分割
                    if ( norm(breakPoint-p(2))<3*minR .or. norm(breakPoint-p(3))<3*minR ) then
                        ibestBreakPoint=j
                        exit
                    end if
                end if
                ! 检查三角形的形状，并找出最好的分割方案
                S1=calcArea(p(1),p(2),breakPoint)
                R1=calcRadius(p(1),p(2),breakPoint,S1)
                S2=calcArea(p(1),p(3),breakPoint)
                R2=calcRadius(p(1),p(3),breakPoint,S2)
                S1=min(S1/(pi*R1**2),S2/(pi*R2**2))
                if (S1 < goodShapeCoeff) then
                    if (bestShape<S1) then
                        bestShape=S1
                        ibestBreakPoint=j
                    end if
                else
                    ! 在geometry.b中的下标
                    ibestBreakPoint=j
                    exit
                end if
            end if
            i=getNextBpoint(i,minR=minR,iend=k)
        end do
        
        if (ibestBreakPoint==-1) then
            ! 找不到分割点
            re=.false.
            return
        end if
        
        ! 划分面元
           ! 添加新的格点
        v0=addVertex(geometry.b(ibestBreakPoint),.true.,ibestBreakPoint)
        if (v0==-1) then
            write(*,*) "splitBorder(): 内存不足。"
            re=.false.
            return
        end if
           ! 修改旧的边
        if (oldE.s(1)==-1) then
            oldE.v(1)=v0
        else
            oldE.v(2)=v0
        end if
        oldE.len=norm(meshInfo.vertex(oldE.v(1)).p-meshInfo.vertex(oldE.v(2)).p)
            ! 添加两条边
        e1=addEdge(oldS.v(1),v0)     ! 分割边
        meshInfo.edge(e1).s(2)=ioldS
        e2=addEdge(oldS.v(2),v0,.true.)   ! 新边界边
        if (e1==-1 .or. e2==-1) then
            write(*,*) "splitBorder(): 内存不足。"
            re=.false.
            return
        end if
            ! 修改旧的面元
        v20=oldS.v(2)
        oldS.v(2)=v0
        e30=oldS.e(3)
        oldS.e(3)=e1
        oldS.S=calcArea(p(1),meshInfo.vertex(v0).p,p(3))
        oldS.R=calcRadius(p(1),meshInfo.vertex(v0).p,p(3),oldS.S)
            ! 添加新面元
        i=addSurface(oldS.v(1),v20,v0,e2,e1,e30)
        if (i==-1) then
            write(*,*) "splitBorder(): 内存不足。"
            re=.false.
            return
        end if
        re=.true.
        return
    end function splitBorder
    
! splitSurface(): 分割面元
!----------------------------------------------------------------------------------
! 1. 可以直接调用此函数分割面元，这样做并不会破坏meshInfo的结构。
! 2. 将面元分成两份，自动选择面元中最长的边进行分割。
! 3. 选择分割后最小形状因子{S/(pi*R^2)}最大的方案。
! 4. 在checkShape=.true.时，如果分割后形状因子从大于badShapeCoeff变成
!     小于badShapeCoeff，或者一开始就小于badShapeCoeff，但变得更小，则不分割。
! 5. 不分割边界边。
! 6. 分割方式见splitSurface.jpg
![Input]
!  id  - 要分割的面元在meshInfo.surface中的下标
! minR - 新的边长不应小于这个值
! checkShape - 是否要检查分割后的面元的形状因子
![Return]
!  成功分割   - .true.
!  未成功分割 - .false.
!----------------------------------------------------------------------------------
    logical function splitSurface(id, minR, checkShape) result(re)
        use mesh 
        implicit none
        integer, intent(in)  :: id
        real(8), intent(in)  :: minR
        logical, intent(in)  :: checkShape
        integer :: s(4), e(0:3), v(4), ibp, localIDe0(2)   ! 各元素的全局/局域下标
        integer :: e130, e140
        type(point_t)           :: breakpoint
        type(edge_t),pointer    :: e0
        type(edge_t)            :: oldE(3)
        integer                 :: i, j, imaxE
        real(8)                 :: bPos
        
        oldE=meshInfo.edge(meshInfo.surface(id).e(1:3))
        ! 找最长边
        if (oldE(1).len>=oldE(2).len .and. oldE(1).len>=oldE(3).len) then
            imaxE=1
        else if (oldE(2).len>=oldE(1).len .and. oldE(2).len>=oldE(3).len) then
            imaxE=2
        else
            imaxE=3  ! 最长边的局域下标
        end if
        ! 若最长边是边界边则不分割
        if (oldE(imaxE).isB) then
            re=.false.
            return
        end if
        
        ! 记录待分割的边和面元
        e(0)=meshInfo.surface(id).e(imaxE)
        e0=>meshInfo.edge(e(0))
        s(1:2)=e0.s
        !! 如果e0不是近似局域最长边，则不分割
        if (e0.s(1)==id) then
            i=e0.s(2)
        else
            i=e0.s(1)
        end if
        if (1.2d0*e0.len<maxval(meshInfo.edge( meshInfo.surface(i).e(:) ).len)) then
            re=.false.
            return
        end if
        
        ! 记录顶点
        v(1:2)=e0.v
        do i=1,2
            do j=1,3
                if (meshInfo.surface(s(i)).e(j) == e(0)) then
                    v(i+2)=meshInfo.surface(s(i)).v(j)
                    localIDe0(i)=j
                end if
            end do
        end do
        
        ! 找分割点, 坐标记录在breakpoint中
            ! 三分法找shape最好的分割点
        bPos=findBestBreakPoint(v,minR,checkShape)
        if (bPos<0d0) then
            re=.false.
            return
        end if
        breakpoint=bPos/e0.len*(meshInfo.vertex(e0.v(2)).p - meshInfo.vertex(e0.v(1)).p) + meshInfo.vertex(e0.v(1)).p
        
        ! 添加分割点
        ibp=addVertex(breakpoint)
        if (ibp==-1) then
            write(*,*) "splitSurface(): 内存不足。"
            re=.false.
            return
        end if
        
        ! 修改e0的端点及长度
        e0.v(1)=ibp
        e0.len=norm(meshInfo.vertex(e0.v(1)).p-meshInfo.vertex(e0.v(2)).p)
        
        ! 添加边
        e(1)=addEdge(v(1),ibp)
        e(2)=addEdge(v(3),ibp)
        e(3)=addEdge(ibp,v(4))
        if (e(1)==-1 .or. e(2)==-1 .or. e(3)==-1) then
            write(*,*) "splitSurface(): 内存不足。"
            re=.false.
            return
        end if
        
        ! 修改旧面元
            ! 修改顶点和边
        meshInfo.surface(s(1)).v( mod(localIDe0(1)+1,3)+1 )=ibp
        meshInfo.surface(s(2)).v( mod(localIDe0(2)  ,3)+1 )=ibp
        e130=meshInfo.surface(s(1)).e( mod(localIDe0(1)  ,3)+1 )
        e140=meshInfo.surface(s(2)).e( mod(localIDe0(2)+1,3)+1 )
        meshInfo.surface(s(1)).e( mod(localIDe0(1)  ,3)+1 )=e(2)
        meshInfo.surface(s(2)).e( mod(localIDe0(2)+1,3)+1 )=e(3)
            ! 更新面积和半径
        do i=1,2
            meshInfo.surface(s(i)).S=calcArea(s(i))
            meshInfo.surface(s(i)).R=calcRadius(s(i),meshInfo.surface(s(i)).S)
        end do
            ! 更新边与面元的相邻的信息
        meshInfo.edge(e(2)).s(1)=s(1)
        meshInfo.edge(e(3)).s(1)=s(2)
        
        ! 添加面元
        s(3)=addSurface(ibp,v(1),v(3),e(1),e(2),e130)
        s(4)=addSurface(ibp,v(4),v(1),e(1),e(3),e140)
        if (s(3)==-1 .or. s(4)==-1) then
            write(*,*) "splitSurface(): 内存不足。"
            re=.false.
            return
        end if
        
        re=.true.
        return
    end function splitSurface
    
! checkKeyPoint(): 检查是否所有关键点均为格点，如不是则分割边界
!---------------------------------------------------------------------
    subroutine checkKeyPoint()
        use mesh 
        implicit none
        integer :: i, j, k, ng, vk, e2, e3
        type(edge_t),pointer :: e
        type(surface_t),pointer :: s
        logical  :: added
        
        ng=geometry.n
        i=1
        do while (i<=meshInfo.nb)
            e=>meshInfo.edge(meshInfo.boundary(i))
            s=>meshInfo.surface(max(e.s(1),e.s(2)))
            j=mod(meshInfo.vertex(s.v(2)).bid,ng)+1
            k=meshInfo.vertex(s.v(3)).bid
            do while (j/=k)
                if (geometry.keyPointID(j)>0) then
                    ! 如果存在非格点的关键点，则分割边界
                    if (inArea(meshInfo.vertex(s.v(2)).p,meshInfo.vertex(s.v(3)).p)/=-1) then
                        write(*,'(A,2E10.3)') "checkKeyPoint(): 网格精度太低, 无法添加以下关键点为格点：", geometry.b(j)
                    else
                        added=addKeyPointToVertex(max(e.s(1),e.s(2)), j, i)
                        if (.not. added) then
                            write(*,'(A,2E10.3)') 'checkKeyPoint(): 内存不足, 无法添加以下关键点为格点：', geometry.b(j)
                        end if
                    end if
                end if
                j=mod(j,ng)+1
            end do
            i=i+1
        end do
        
    end subroutine checkKeyPoint
    
    
! addKeyPointToVertex(): 将keyPoint添加为格点
!---------------------------------------------------------------
! 分割方式见addKeyPointToVertex.jpg
![Input]
! s0   - 原始面元在meshInfo.surface中的下标
! ikp  - 需添加的keyPoint在geometry.b中的下标
! iedge - s0的边界边在meshInfo.boundary中的下标
![Return]
! re - 返回是否添加成功
!---------------------------------------------------------------
    logical function addKeyPointToVertex(s0, ikp, iedge) result(re)
        use mesh
        implicit none
        integer, intent(in)    :: s0, ikp, iedge
        integer :: v(3), vb, vk, itmp(4)
        integer :: e(5), e130, e230
        integer :: s(3)
        type(edge_t), pointer    :: oldE
        type(surface_t), pointer :: oldS
        type(point_t)  :: breakPoint, p1, p2
        real(8)        :: bPos
        
        ! 记录旧的信息
        oldS=>meshInfo.surface(s0)
        oldE=>meshInfo.edge(oldS.e(1))
        e130=oldS.e(2)
        e230=oldS.e(1)
        v(:)=oldS.v(:)
        
        ! 添加点vk
        vk=addVertex(geometry.b(ikp),.true.,ikp)
        if (vk==-1) then
            write(*,*) "addKeyPointToVertex(): 内存不足。"
            re=.false.
            return
        end if
        
        ! 计算breakPoint
        itmp(1)=v(2)
        itmp(2)=v(3)
        itmp(3)=v(1)
        itmp(4)=vk
        bPos=findBestBreakPoint(itmp(1:4),0d0,.false.)
        p1=meshInfo.vertex(v(2)).p
        p2=meshInfo.vertex(v(3)).p
        breakPoint=p1+(p2-p1)*bPos/norm(p2-p1)
        
        ! 添加点vb
        vb=addVertex(breakPoint)
        if (vb==-1) then
            write(*,*) "addKeyPointToVertex(): 内存不足。"
            re=.false.
            return
        end if
        
        ! 添加边
        e(1)=addEdge(v(2),vk,.false.)  ! 此处先设为false，后面手动修改
        e(2)=addEdge(vb,vk)
        e(3)=addEdge(vk,v(3),.true.)
        e(4)=addEdge(vb,v(3))
        e(5)=addEdge(vb,v(1))
        if (minval(e)==-1) then
            write(*,*) "addKeyPointToVertex(): 内存不足。"
            re=.false.
            return
        end if

        ! 修改旧边和旧面元
            ! 取消旧边的边界属性
        oldE.isB=.false.
        meshInfo.boundary(iedge)=e(1)  ! e(1)对应边界不可能含有关键点
        meshInfo.edge(e(1)).isB=.true.
            ! 修改旧边的端点
        if (oldE.v(1)==v(2)) then
            oldE.v(2)=vb
        else
            oldE.v(1)=vb
        end if
        oldE.len=norm(meshInfo.vertex(oldE.v(1)).p-meshInfo.vertex(oldE.v(2)).p)
            ! 修改旧面元
        oldS.v(3)=vb
        oldS.e(2)=e(5)
        meshInfo.edge(e(5)).s(2)=s0    
        oldS.S=calcArea(s0)
        oldS.R=calcRadius(s0,oldS.S)
        
        ! 添加面元
        s(1)=addSurface(vb,v(2),vk,e230,e(1),e(2))
        s(2)=addSurface(vb,vk,v(3),e(2),e(3),e(4))
        s(3)=addSurface(vb,v(3),v(1),e130,e(4),e(5))
        if (s(1)==-1 .or. s(2)==-1 .or. s(3)==-1) then
            write(*,*) "addKeyPointToVertex(): 内存不足。"
            re=.false.
            return
        end if
        re=.true.
        return
    end function addKeyPointToVertex
    

