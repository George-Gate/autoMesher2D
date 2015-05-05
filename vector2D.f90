!vector2D.f90 
!
![Mudules]
!  vector2D - 2维矢量操作
!  
![Operators]
!    +   加法
!    -   减法
!    *   数乘或叉乘
!    **  数乘或内积
!    /   数乘(第二个操作数为实数)
!

    module vector2D
        implicit none
! <<< Type definition >>>
        type point_t    ! 矢量
            real(8)  :: x,y
        end type point_t
    
! <<< Interfaces >>>
        interface operator (+)
            module procedure padd
        end interface operator (+)

        interface operator (-)
            module procedure pminus
        end interface operator (-)
    
        interface operator (/)
            module procedure pdivide
        end interface operator (/)
    
        interface operator (*)
            module procedure pcross
            module procedure pntimesl
            module procedure pntimesr
        end interface operator (*)
        
        interface operator (**)
            module procedure ptimes
            module procedure pntimesl
            module procedure pntimesr
        end interface operator (**)
    
    contains
    
        type(point_t) elemental function padd(a,b)
            type(point_t), intent(in) :: a,b
            padd.x=a.x+b.x
            padd.y=a.y+b.y
        end function padd
    
        type(point_t) elemental function pminus(a,b)
            type(point_t), intent(in) :: a,b
            pminus%x=a%x-b%x
            pminus%y=a%y-b%y
        end function pminus
        
        type(point_t) elemental function pdivide(a,b)
            type(point_t), intent(in) :: a
            real(8), intent(in)       :: b
            pdivide%x=a%x / b
            pdivide%y=a%y / b
        end function pdivide
        
        type(point_t) elemental function pntimesr(a,b)
            type(point_t), intent(in) :: a
            real(8), intent(in)       :: b
            pntimesr%x=a%x * b
            pntimesr%y=a%y * b
        end function pntimesr
    
        type(point_t) elemental function pntimesl(a,b)
            type(point_t), intent(in) :: b
            real(8), intent(in)       :: a
            pntimesl%x=a * b%x
            pntimesl%y=a * b%y
        end function pntimesl
    
        real(8) elemental function pcross(a,b)
            type(point_t), intent(in) :: a,b
            pcross=a%x * b%y - b%x * a%y
        end function pcross

        real(8) elemental function ptimes(a,b)
            type(point_t), intent(in) :: a,b
            ptimes=a%x * b%x + a%y * b%y
        end function ptimes
    
        real(8) elemental function norm(a)
            type(point_t), intent(in) :: a
            norm=sqrt((a%x)**2+(a%y)**2)
        end function norm
    
        ! 返回两线段的交点，调用前请先确定两线段确实相交。
        type(point_t) elemental function insecPoint(a1,a2,b1,b2) result(re)
            type(point_t), intent(in) :: a1, a2, b1, b2
            type(point_t)             :: v0, v1, v2
            v0=a2-a1
            v1=b1-a1
            v2=b1-b2
            re=a1+v0*(v1.x*v2.y-v1.y*v2.x)/(v0.x*v2.y-v0.y*v2.x)
        end function insecPoint
        
    end module vector2D