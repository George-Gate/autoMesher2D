	  p1  �   k820309    V          13.1        <2=U                                                                                                           
       D:\文档\学校\大学\12物理\学科\2015-03 计算物理[贺彦章]\Homework_3\finiteElementAutoMeshing\mesh.f90 MESH_INTERFACE          @ "                                       u #VALIDATEMESHINFO             @  "                                      	   u #MAKEGEOMETRY             @  "                                      	   u #OUTPUTGEOMETRY             @  "                                      	   u #OUTPUTMESH             @ "                                       u #SPLITSURFACE             @  "                                      	   u #INITMESH             @ "                                       u #SPLITBORDER                                                           u #INAREA_POINT    #INAREA_EDGE 	            H "                                      u #GETNEXTBPOINT 
            @ "                                       u #ADDVERTEX             @ "                                       u #ADDEDGE             @ "                                       u #ADDSURFACE             @  "                                      	   u #EXPANDMESHINFO             H "                                      u #ARCLEN                                                           u #CALCAREA_P    #CALCAREA_S                                                           u #CALCRADIUS_P    #CALCRADIUS_S             H "                                      u #GETSHAPECOEFF                                                           u #FINDBESTBREAKPOINT_V    #FINDBESTBREAKPOINT_P             @ "                                       u #ADDKEYPOINTTOVERTEX             @  "                                      	   u #GENERATEVERTEXID             @  "                                      	   u #REFINEMESH             @  "                                      	   u #CHECKKEYPOINT             @  "                                      	   u #CLEANMESH    %         @ "                                                          #         @  "                                     	               #X    #Y    #STEP    #KEYPOINTS              
                                                   
              &                                                     
                                                   
              &                                                     
                                     
                
                                                                 &                                           #         @  "                                     	               #FILEID               
                                             #         @  "                                     	               #FILEID !   #TAG "             
                                 !                     
                               "                    1 %         @ "                                                         #ID #   #MINR $   #CHECKSHAPE %             
                                 #                     
                                $     
                
                                 %           #         @  "                                     	               #X &   #Y '   #MINR (   #MAXR )   #KEYPOINTS *   #OUTPUTDETAILS +   #STEP ,             
                                &                   
              &                                                     
                                '                   
              &                                                     
                                (     
                
                                )     
                
                                *                                 &                                                     
                                +                     
                               ,     
      %         @ "                                                         #ID -   #MINR .             
                                 -                     
                                .     
      %         H                                                         #INAREA_POINT%POINT_T /   #P0 2                     @                           /     '                    #X 0   #Y 1                �                              0                
                �                              1               
             
                                 2                   #INAREA_POINT%POINT_T /   %         H                               	                          #INAREA_EDGE%POINT_T 3   #P1 6   #P2 7                     @                           3     '                    #X 4   #Y 5                �                              4                
                �                              5               
             
                                 6                   #INAREA_EDGE%POINT_T 3             
                                 7                   #INAREA_EDGE%POINT_T 3   %         H "                             
                           #I0 8   #MINR 9   #MAXR :   #IEND ;             
                                 8                     
                               9     
                
                               :     
                
                                ;           %         @ "                                                        #ADDVERTEX%POINT_T <   #P0 ?   #ISB @   #BID A                     @                           <     '                    #X =   #Y >                �                              =                
                �                              >               
             
                                 ?                   #ADDVERTEX%POINT_T <             
                                @                     
                                A           %         @ "                                                         #A B   #B C   #ISB D             
                                 B                     
                                 C                     
                                D           %         @ "                                                         #A E   #B F   #C G   #X H   #Y I   #Z J             
                                 E                     
                                 F                     
                                 G                     
                                H                     
                                I                     
                                J           #         @  "                                     	               #MV K   #ME L   #MS M   #MB N             
                               K     
                
                               L     
                
                               M     
                
                               N     
      %         H "                                                 
       #IB O             
                                 O           %         H                                                  
       #CALCAREA_P%POINT_T P   #P1 S   #P2 T   #P3 U                     @                           P     '                    #X Q   #Y R                �                              Q                
                �                              R               
             
                                 S                   #CALCAREA_P%POINT_T P             
                                 T                   #CALCAREA_P%POINT_T P             
                                 U                   #CALCAREA_P%POINT_T P   %         H                                                   
       #SID V             
                                 V           %         H                                                  
       #CALCRADIUS_P%POINT_T W   #P1 Z   #P2 [   #P3 \   #S ]                     @                           W     '                    #X X   #Y Y                �                              X                
                �                              Y               
             
                                 Z                   #CALCRADIUS_P%POINT_T W             
                                 [                   #CALCRADIUS_P%POINT_T W             
                                 \                   #CALCRADIUS_P%POINT_T W             
                               ]     
      %         H                                                   
       #SID ^   #S _             
                                 ^                     
                               _     
      %         H "                                                
       #GETSHAPECOEFF%POINT_T `   #P1 c   #P2 d   #P3 e                     @                           `     '                    #X a   #Y b                �                              a                
                �                              b               
             
                                 c                   #GETSHAPECOEFF%POINT_T `             
                                 d                   #GETSHAPECOEFF%POINT_T `             
                                 e                   #GETSHAPECOEFF%POINT_T `   %         @                                                   
       #V f   #MINLEN g   #CHECKSHAPE h             
                                 f                       p          p            p                                    
                                g     
                
                                 h           %         @                                                  
       #FINDBESTBREAKPOINT_P%POINT_T i   #P l   #MINLEN m   #CHECKSHAPE n                     @                           i     '                    #X j   #Y k                �                              j                
                �                              k               
             
                                 l                        p          p            p                          #FINDBESTBREAKPOINT_P%POINT_T i             
                                m     
                
                                 n           %         @ "                                                         #S0 o   #IKP p   #IEDGE q             
                                 o                     
                                 p                     
                                 q           #         @  "                                     	               #SMALLENDAREA r   #BIGENDAREA s             
                                r                    	             &                   &                                                     
                                s                    
             &                   &                                           #         @  "                                     	               #MINR t   #MAXR u   #OUTPUTDETAILS v             
                                t     
                
                                u     
                
                                v           #         @  "                                     	                #         @  "                                     	                   �   �      fn#fn %   ,  V       gen@VALIDATEMESHINFO !   �  R       gen@MAKEGEOMETRY #   �  T       gen@OUTPUTGEOMETRY    (  P       gen@OUTPUTMESH !   x  R       gen@SPLITSURFACE    �  N       gen@INITMESH       Q       gen@SPLITBORDER    i  c       gen@INAREA "   �  S       gen@GETNEXTBPOINT      O       gen@ADDVERTEX    n  M       gen@ADDEDGE    �  P       gen@ADDSURFACE #     T       gen@EXPANDMESHINFO    _  L       gen@ARCLEN    �  `       gen@CALCAREA      d       gen@CALCRADIUS "   o  S       gen@GETSHAPECOEFF '   �  t       gen@FINDBESTBREAKPOINT (   6  Y       gen@ADDKEYPOINTTOVERTEX %   �  V       gen@GENERATEVERTEXID    �  P       gen@REFINEMESH "   5  S       gen@CHECKKEYPOINT    �  O       gen@CLEANMESH !   �  P       VALIDATEMESHINFO    '	  o       MAKEGEOMETRY    �	  �   a   MAKEGEOMETRY%X    "
  �   a   MAKEGEOMETRY%Y "   �
  @   a   MAKEGEOMETRY%STEP '   �
  �   a   MAKEGEOMETRY%KEYPOINTS    z  T       OUTPUTGEOMETRY &   �  @   a   OUTPUTGEOMETRY%FILEID      ]       OUTPUTMESH "   k  @   a   OUTPUTMESH%FILEID    �  L   a   OUTPUTMESH%TAG    �  r       SPLITSURFACE     i  @   a   SPLITSURFACE%ID "   �  @   a   SPLITSURFACE%MINR (   �  @   a   SPLITSURFACE%CHECKSHAPE    )  �       INITMESH    �  �   a   INITMESH%X    K  �   a   INITMESH%Y    �  @   a   INITMESH%MINR      @   a   INITMESH%MAXR #   W  �   a   INITMESH%KEYPOINTS '   �  @   a   INITMESH%OUTPUTDETAILS    #  @   a   INITMESH%STEP    c  b       SPLITBORDER    �  @   a   SPLITBORDER%ID !     @   a   SPLITBORDER%MINR    E  r       INAREA_POINT .   �  ^      INAREA_POINT%POINT_T+VECTOR2D 0     H   a   INAREA_POINT%POINT_T%X+VECTOR2D 0   ]  H   a   INAREA_POINT%POINT_T%Y+VECTOR2D     �  b   a   INAREA_POINT%P0      y       INAREA_EDGE -   �  ^      INAREA_EDGE%POINT_T+VECTOR2D /   �  H   a   INAREA_EDGE%POINT_T%X+VECTOR2D /   &  H   a   INAREA_EDGE%POINT_T%Y+VECTOR2D    n  a   a   INAREA_EDGE%P1    �  a   a   INAREA_EDGE%P2    0  v       GETNEXTBPOINT !   �  @   a   GETNEXTBPOINT%I0 #   �  @   a   GETNEXTBPOINT%MINR #   &  @   a   GETNEXTBPOINT%MAXR #   f  @   a   GETNEXTBPOINT%IEND    �  �       ADDVERTEX +   '  ^      ADDVERTEX%POINT_T+VECTOR2D -   �  H   a   ADDVERTEX%POINT_T%X+VECTOR2D -   �  H   a   ADDVERTEX%POINT_T%Y+VECTOR2D      _   a   ADDVERTEX%P0    t  @   a   ADDVERTEX%ISB    �  @   a   ADDVERTEX%BID    �  g       ADDEDGE    [  @   a   ADDEDGE%A    �  @   a   ADDEDGE%B    �  @   a   ADDEDGE%ISB      z       ADDSURFACE    �  @   a   ADDSURFACE%A    �  @   a   ADDSURFACE%B      @   a   ADDSURFACE%C    U  @   a   ADDSURFACE%X    �  @   a   ADDSURFACE%Y    �  @   a   ADDSURFACE%Z      h       EXPANDMESHINFO "   }  @   a   EXPANDMESHINFO%MV "   �  @   a   EXPANDMESHINFO%ME "   �  @   a   EXPANDMESHINFO%MS "   =  @   a   EXPANDMESHINFO%MB    }  X       ARCLEN    �  @   a   ARCLEN%IB      �       CALCAREA_P ,   �  ^      CALCAREA_P%POINT_T+VECTOR2D .   �  H   a   CALCAREA_P%POINT_T%X+VECTOR2D .   ;   H   a   CALCAREA_P%POINT_T%Y+VECTOR2D    �   `   a   CALCAREA_P%P1    �   `   a   CALCAREA_P%P2    C!  `   a   CALCAREA_P%P3    �!  Y       CALCAREA_S    �!  @   a   CALCAREA_S%SID    <"  �       CALCRADIUS_P .   �"  ^      CALCRADIUS_P%POINT_T+VECTOR2D 0   ##  H   a   CALCRADIUS_P%POINT_T%X+VECTOR2D 0   k#  H   a   CALCRADIUS_P%POINT_T%Y+VECTOR2D     �#  b   a   CALCRADIUS_P%P1     $  b   a   CALCRADIUS_P%P2     w$  b   a   CALCRADIUS_P%P3    �$  @   a   CALCRADIUS_P%S    %  `       CALCRADIUS_S !   y%  @   a   CALCRADIUS_S%SID    �%  @   a   CALCRADIUS_S%S    �%  �       GETSHAPECOEFF /   |&  ^      GETSHAPECOEFF%POINT_T+VECTOR2D 1   �&  H   a   GETSHAPECOEFF%POINT_T%X+VECTOR2D 1   "'  H   a   GETSHAPECOEFF%POINT_T%Y+VECTOR2D !   j'  c   a   GETSHAPECOEFF%P1 !   �'  c   a   GETSHAPECOEFF%P2 !   0(  c   a   GETSHAPECOEFF%P3 %   �(  s       FINDBESTBREAKPOINT_V '   )  �   a   FINDBESTBREAKPOINT_V%V ,   �)  @   a   FINDBESTBREAKPOINT_V%MINLEN 0   �)  @   a   FINDBESTBREAKPOINT_V%CHECKSHAPE %   *  �       FINDBESTBREAKPOINT_P 6   �*  ^      FINDBESTBREAKPOINT_P%POINT_T+VECTOR2D 8   +  H   a   FINDBESTBREAKPOINT_P%POINT_T%X+VECTOR2D 8   U+  H   a   FINDBESTBREAKPOINT_P%POINT_T%Y+VECTOR2D '   �+  �   a   FINDBESTBREAKPOINT_P%P ,   S,  @   a   FINDBESTBREAKPOINT_P%MINLEN 0   �,  @   a   FINDBESTBREAKPOINT_P%CHECKSHAPE $   �,  l       ADDKEYPOINTTOVERTEX '   ?-  @   a   ADDKEYPOINTTOVERTEX%S0 (   -  @   a   ADDKEYPOINTTOVERTEX%IKP *   �-  @   a   ADDKEYPOINTTOVERTEX%IEDGE !   �-  j       GENERATEVERTEXID .   i.  �   a   GENERATEVERTEXID%SMALLENDAREA ,   /  �   a   GENERATEVERTEXID%BIGENDAREA    �/  o       REFINEMESH      0  @   a   REFINEMESH%MINR     `0  @   a   REFINEMESH%MAXR )   �0  @   a   REFINEMESH%OUTPUTDETAILS    �0  H       CHECKKEYPOINT    (1  H       CLEANMESH 