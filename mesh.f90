!mesh.f90 
![Name]
!    2D Auto Mesher for Simply Connected Polygon Area
![Description]
!    按照给定精度，在单连通多边形区域内自动生成二维三角形网格。
!    生成好的网格记录在自定义类型变量meshInfo中，输入的求解域边界将被自动离散化并记
!录在geometry变量中。生成的网格满足如下性质：
!    1. 每个面元内部格点的编号逆时针递增
!    2. 面元的1号顶点不在边界上
!    3. 所有边界关键点都是格点
!    4. 边界格点的全局编号可由用户指定
!具体的数据结构说明，请参考module mesh_typedef中的定义及注释。
!
![Modules]
!  mesh           - 全局变量/常量的声明，use此module即可使用本mesher的全部功能。
!  mesh_typedef   - 数据类型的定义
!  mesh_interface - 函数的interface
!==================================================================================
![Important Functions]
!
!  initMesh(): 自动生成网格 
!----------------------------------------------------------------------------------
!  1、调用makeGeometry()初始化求解域geometry
!  2、初始化meshInfo变量
!  3、建立最小网格
!  4、调用refineMesh()细化网格
!  5、调用validateMeshInfo()检查网格有效性
![Syntax]
!   subroutine initMesh(x, y, minR, maxR, keyPoints, outputDetails, step)
!       real(8), intent(in) :: x(:), y(:)
!       real(8), intent(in) :: minR, maxR
!       integer, intent(in), optional :: keyPoints(:)
!       logical, intent(in), optional :: outputDetails
!       real(8), intent(in), optional :: step
!   end subroutine initMesh
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
!
! refineMesh(): 细化网格
!----------------------------------------------------------------------------------
!  只能细化网格，不能使网格变粗
!  通过裂边的方式细化，不会移动现有格点
!  通过形状判定，限制钝角的出现
!  边界处面元的最大边长在 minR/2 ~ minR 之间
!  区域中间面元边长在 maxR/2 ~ maxR 之间
!  形状不好的面元尺寸会较小
![Syntax]
!   subroutine refineMesh(minR, maxR, outputDetails)
!       real(8), intent(in) :: minR, maxR
!       logical, intent(in), optional :: outputDetails
!   end subroutine refineMesh
![Input]
! minR, maxR - 面元的尺寸
! outputDetails - [可选]指示是否需要将每一步迭代过后的网格状态输出至文件mesh.txt
!                 默认为否
!----------------------------------------------------------------------------------
!
!  validateMeshInfo(): 检查meshInfo中的信息是否满足定义的要求，并尝试修复错误。
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
![Syntax]
!   integer function validateMeshInfo()
!   end function validateMeshInfo
![Return]
!  如果检查没有发现问题，则返回0；否则返回发现的问题在上述分类中的编号。
!  若存在多个问题，则会返回最后一个问题的编号。
!----------------------------------------------------------------------------------
!
!  generateVertexID(): 建立格点编号,结果记录在meshInfo.vertex(:).id中。
!----------------------------------------------------------------------------------
![Syntax]
!   subroutine generateVertexID(smallEndArea, bigEndArea)
!       integer, intent(in), optional :: smallEndArea(:,:), bigEndArea(:,:)
!   end subroutine generateVertexID
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
!
!  outputGeometry(): 输出geometry内的信息至文件。
!----------------------------------------------------------------------------------
![Syntax]
!   subroutine outputGeometry(fileID)
!       integer, intent(in) :: fileID
!   end subroutine outputGeometry
!----------------------------------------------------------------------------------
!
!  outputMesh(): 将meshInfo中的网格信息输出至fileID对应的文件。
!----------------------------------------------------------------------------------
![Syntax]
!   subroutine outputMesh(fileID,tag)
!       integer, intent(in) :: fileID
!       character(*), intent(in),optional :: tag
!   end subroutine outputMesh
!----------------------------------------------------------------------------------
!
!==================================================================================
![Function List] 
!(警告：在直接调用下列函数之前，请详细了解该函数的功能，否则可能会损坏meshInfo的结构。
! 每个函数的详细说明包含在各函数前面的注释中，请直接到mesh_lib.f90中查看。)
!
! addEdge(): 添加边（自动识别并设置边界）
!            (警告：直接调用可能会破坏meshInfo结构)
! addKeyPointToVertex(): 将keyPoint添加为格点
!                        (警告：直接调用可能会破坏meshInfo结构)
! addSurface(): 添加面元(警告：直接调用可能会破坏meshInfo结构)
! addVertex(): 添加格点(警告：直接调用可能会破坏meshInfo结构)
! arclen(): 返回边界边对应的实际边界长度
! calcArea(): interface of
!   calcArea_p(): 返回由p1,p2,p3构成的三角形的面积
!   calcArea_s(): 返回sID对应面元的面积
! calcRadius(): interface of
!   calcRadius_p(): 返回由p1,p2,p3构成的三角形的外接圆半径,可输入三角型面积加速计算
!   calcRadius_s(): 返回sID对应面元的外接圆半径,可输入三角型面积加速计算
! checkKeyPoint(): 检查是否所有关键点均为格点，如不是则分割边界。
! cleanMesh(): 清空meshInfo和geometry的内容，释放内存
! expandMeshInfo(): 增加meshInfo的大小，若增加失败，则会尝试回滚。
!                   (警告：直接调用可能会破坏meshInfo结构)
! findBestBreakPoint(): interface of
!   findBestBreakPoint_v(): 用三分法找出形状因子最好的分割点。
!   findBestBreakPoint_p(): 用三分法找出形状因子最好的分割点。
! generateVertexID(): 建立格点编号,结果记录在meshInfo.vertex(:).id中。
! getNextBpoint(): 在边界上取下一个点(沿边界正方向)
! getShapeCoeff(): 计算形状因子 S/(pi*R^2)
! inArea(): interface of
!   inArea_point(): 判断一个点是否在求解域内
!   inArea_edge(): 判断一个线段(除两端点外)是否完全在求解域内
! initMesh(): 自动生成网格
! makeGeometry(): 离散化边界、初始化geometry变量。
!                 (警告：直接调用可能会破坏meshInfo结构)
! outputGeometry(): 输出geometry内的信息至文件
! outputMesh(): 将meshInfo中的网格信息输出至fileID对应的文件
! refineMesh(): 细化网格
! splitBorder(): 分割边界边
! splitSurface(): 分割面元
! validateMeshInfo(): 检查meshInfo中的信息是否满足定义的要求，并尝试修复错误。
!==================================================================================

!include "vector2D.f90"
    module mesh_typedef
        use vector2D
        implicit none
! <<< Type definition >>>
        type vertex_t    ! 记录格点信息
            type(point_t) :: p      ! 格点坐标
            integer       :: id=-1  ! 格点编号
            logical       :: isB    ! 记录这个点是否为边界点
            integer       :: bid    ! 如果格点在边界上，则记录该点在geometry.b中的下标
        end type vertex_t
    
        type edge_t    ! 记录边的信息
            integer :: v(2)       ! 边的端点
            integer :: s(2)=-1    ! 边两侧的面元, s(1)在v(2)-v(1)的右边
                                  ! 只有一边有面元的边为边界边
            real(8) :: len        ! 边的长度
            logical :: isB=.false.! 记录此边是否为边界边
        end type edge_t
    
        type surface_t   ! 记录面元的信息
            integer   :: v(3)   ! 格点, 按逆时针从小到大编号，如果面元有边界边，
                                ! 则边界边的两端必然是v(2)和v(3)，v(1)不允许在边界上。
            integer   :: e(3)   ! 边，v(i)的对边是e(i)
            real(8)   :: S=0d0  ! 面积
            real(8)   :: R=0d0  ! 外接圆半径
        end type surface_t
    
        type meshInfo_t   ! 记录整个网格的全部信息
            integer                      :: maxN               ! 当前允许的最大格点数
            integer                      :: nv, ne, ns, nb     ! 网格中的格点、边、面元、边界数
            type(vertex_t), allocatable  :: vertex(:)          ! 记录格点的信息
            type(edge_t), allocatable    :: edge(:)            ! 记录边的信息
            type(surface_t), allocatable :: surface(:)         ! 记录面元的信息
            integer, allocatable         :: boundary(:)        ! 记录边界包含哪些边
        end type meshInfo_t
    
        type geometry_t  ! 用于记录求解域的信息
            integer                     :: n         ! boundary中的点数
            type(point_t), allocatable  :: b(:)  
!-----------------------------------------------------------------------------------------------
!       记录一组点，将这些点按照1-->2-->...-->n+1的顺序连起来将得到一个闭合的单连通求解域的边界。
!       1-->2-->3-->...为该边界的正方向(右手螺旋方向)
!       b(1)=b(n+1)
!-----------------------------------------------------------------------------------------------
            real(8), allocatable :: l(:)             ! l(i)是从b(1)到b(i+1)的弧长，l(n)为整个边界的长度。
            real(8)              :: steplb, stepub   ! 初始化geometry时的采样步长的上下界,
                                                     ! 相当于边界上每条线段长度的上下界。    
            integer, allocatable :: keyPointID(:)    ! 指示一个点是否为求解域边界的关键点，关键点一定会成为格点。
                                                     ! (下标范围：1~n+1)曲率较大的点或不同类型边界的分界点将设为关键点。
                                                     ! keyPointID==0表示非关键点，大于零则是关键点。关于此数组的取值，
                                                     ! 请参考makeGeometry()的注释。
            integer, allocatable :: bidOfKeyPoint(:) ! 如果某个关键点的keyPointID为k，那么这个关键点在b数组中的下标为
                                                     ! bidOfKeyPoint(k), bidOfKeyPoint(i)=-1表示不存在keyPointID=i的关键点。
        end type geometry_t
    end module mesh_typedef

    module mesh_interface
        implicit none
! <<< Interfaces >>>
        interface validateMeshInfo
            integer function validateMeshInfo()
            end function validateMeshInfo
        end interface validateMeshInfo

        interface makeGeometry
            subroutine makeGeometry(x,y,step,keyPoints)
                real(8), intent(in) :: x(:), y(:), step
                integer, intent(in),optional  :: keyPoints(:)
            end subroutine makeGeometry
        end interface makeGeometry

        interface outputGeometry
            subroutine outputGeometry(fileID)
                integer, intent(in) :: fileID
            end subroutine outputGeometry
        end interface outputGeometry
        
        interface outputMesh
            subroutine outputMesh(fileID,tag)
                integer, intent(in) :: fileID
                character(*), intent(in),optional :: tag
            end subroutine outputMesh
        end interface outputMesh
        
        interface splitSurface
            logical function splitSurface(id, minR, checkShape)
                integer, intent(in)  :: id
                real(8), intent(in)  :: minR
                logical, intent(in)  :: checkShape
            end function splitSurface
        end interface splitSurface
        
        interface initMesh
            subroutine initMesh(x, y, minR, maxR, keyPoints, outputDetails, step)
                real(8), intent(in) :: x(:), y(:)
                real(8), intent(in) :: minR, maxR
                integer, intent(in), optional :: keyPoints(:)
                logical, intent(in), optional :: outputDetails
                real(8), intent(in), optional :: step
            end subroutine initMesh
        end interface initMesh
        
        interface splitBorder
            logical function splitBorder(id, minR)
                integer, intent(in) :: id
                real(8), intent(in) :: minR
            end function splitBorder
        end interface splitBorder
        
        interface inArea
            integer elemental function inArea_point(p0)
                use vector2D
                type(point_t), intent(in) :: p0
            end function inArea_point
            integer elemental function inArea_edge(p1, p2)
                use vector2D
                type(point_t), intent(in) :: p1, p2
            end function inArea_edge
        end interface inArea

        interface getNextBpoint
            integer elemental function getNextBpoint(i0, minR, maxR, iend)
                integer, intent(in)       :: i0
                real(8), intent(in), optional  :: minR, maxR
                integer, intent(in), optional  :: iend
            end function getNextBpoint
        end interface getNextBpoint
        
        interface addVertex
            integer function addVertex(p0, isB, bid)
                use vector2D
                type(point_t), intent(in)     :: p0
                logical, intent(in), optional :: isB
                integer, intent(in), optional :: bid 
            end function addVertex
        end interface addVertex
        
        interface addEdge
            integer function addEdge(a, b, isB)
                integer, intent(in)           :: a,b
                logical, intent(in), optional :: isB        
            end function addEdge
        end interface addEdge
        
        interface addSurface
            integer function addSurface(a, b, c, x, y, z)
                integer, intent(in)  :: a, b, c
                integer, intent(in), optional  :: x, y, z
            end function addSurface
        end interface addSurface

        interface expandMeshInfo
            subroutine expandMeshInfo(mv, me, ms, mb)
                real(8), intent(in), optional  :: mv, me, ms, mb
            end subroutine expandMeshInfo
        end interface expandMeshInfo
        
        interface arclen
            real(8) elemental function arclen(ib)
                integer, intent(in) :: ib
            end function arclen
        end interface arclen
        
        interface calcArea
            real(8) elemental function calcArea_p(p1,p2,p3)
                use vector2D
                type(point_t), intent(in) :: p1,p2,p3
            end function calcArea_p
            real(8) elemental function calcArea_s(sID)
                integer, intent(in) :: sID
            end function calcArea_s
        end interface calcArea
        
        interface calcRadius
            real(8) elemental function calcRadius_p(p1,p2,p3,S)
                use vector2D
                type(point_t), intent(in)   :: p1,p2,p3
                real(8),intent(in),optional :: S
            end function calcRadius_p
            real(8) elemental function calcRadius_s(sID,S)
                integer, intent(in)   :: sID
                real(8),intent(in),optional :: S
            end function calcRadius_s
        end interface calcRadius
        
        interface getShapeCoeff
            real(8) elemental function getShapeCoeff(p1,p2,p3)
                use vector2D
                type(point_t), intent(in)   :: p1,p2,p3
            end function getShapeCoeff
        end interface getShapeCoeff
        
        interface findBestBreakPoint
            real(8) pure function findBestBreakPoint_v(v,minlen,checkShape)
                    integer, intent(in) :: v(4)
                    real(8), intent(in) :: minlen
                    logical, intent(in) :: checkShape
            end function findBestBreakPoint_v
            real(8) pure function findBestBreakPoint_p(p,minlen,checkShape)
                    use vector2D
                    type(point_t), intent(in) :: p(4)
                    real(8), intent(in) :: minlen
                    logical, intent(in) :: checkShape
            end function findBestBreakPoint_p
        end interface findBestBreakPoint
        
        interface addKeyPointToVertex
            logical function addKeyPointToVertex(s0, ikp, iedge)
                integer, intent(in)    :: s0, ikp, iedge
            end function addKeyPointToVertex
        end interface addKeyPointToVertex
        
        interface generateVertexID
            subroutine generateVertexID(smallEndArea, bigEndArea)
                integer, intent(in), optional :: smallEndArea(:,:), bigEndArea(:,:)
            end subroutine generateVertexID
        end interface generateVertexID
        
        interface refineMesh
            subroutine refineMesh(minR, maxR, outputDetails)
                real(8), intent(in) :: minR, maxR
                logical, intent(in), optional :: outputDetails
            end subroutine refineMesh
        end interface refineMesh
        
        interface checkKeyPoint
            subroutine checkKeyPoint()
            end subroutine checkKeyPoint
        end interface checkKeyPoint
        
        interface cleanMesh
            subroutine cleanMesh()
            end subroutine cleanMesh
        end interface cleanMesh
        
    end module mesh_interface
    
    module mesh
        use mesh_typedef
        use mesh_interface
        implicit none
        
! <<< Var declaration >>>
        type(meshInfo_t), public, save, target  :: meshInfo    ! 网格信息
        type(geometry_t), public, save, target  :: geometry    ! 求解域信息
        real(8), parameter, public      :: pi=3.1415926535897932384626433832795d0, eps=1d-14
        real(8), parameter, public      :: badShapeCoeff=0.15  ! 判断三角形面元形状是否不好的阈值，
        real(8), parameter, public      :: goodShapeCoeff=0.35 ! 当S/(pi*R^2)<badShapeCoeff时，判定形状不好。
                                                               ! 当S/(pi*R^2)>goodShapeCoeff时，判定形状很好。
    end module mesh

