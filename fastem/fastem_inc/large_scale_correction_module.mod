	  Q;  �   k820309    �          13.0        �wW                                                                                                           
       Large_Scale_Correction_Module.f90 LARGE_SCALE_CORRECTION_MODULE              IVAR_TYPE LARGE_SCALE_CORRECTION LARGE_SCALE_CORRECTION_TL LARGE_SCALE_CORRECTION_AD                      @                              
       FP                      @                              
       FITCOEFF_3D_TYPE                      @                              
  
     NPTS LPOLY_TYPE FIND_INDEX INTERP_4D INTERP_4D_TL INTERP_4D_AD CLEAR_LPOLY LPOLY LPOLY_TL LPOLY_AD                �                                       u #FIND_REGULAR_INDEX    #FIND_RANDOM_INDEX    #         @     @                                               #FIND_REGULAR_INDEX%SIZE    #FIND_REGULAR_INDEX%FLOOR    #FIND_REGULAR_INDEX%MIN    #FIND_REGULAR_INDEX%MAX    #X 	   #DX 
   #X_INT    #I1    #I2    #OUT_OF_BOUNDS                  @                                 SIZE               @                                 FLOOR               @                                 MIN               @                                 MAX           
  @                             	                   
 *             &                                                     
   @                             
     
                
  @                                  
                 
  @                                                    
  @                                                    
  @                                          #         @     @                                               #FIND_RANDOM_INDEX%SIZE    #FIND_RANDOM_INDEX%MIN    #FIND_RANDOM_INDEX%MAX    #X    #X_INT    #I1    #I2    #OUT_OF_BOUNDS                  @                                 SIZE               @                                 MIN               @                                 MAX           
  @                                                
 +             &                                                     
  @                                  
                 
  @                                                    
  @                                                    
  @                                                            !@               A                '�                    #IS_ALLOCATED    #RELEASE    #VERSION    #DIMENSIONS    #C                � $                                                                                                                                                 � $                                                                                                                            1                � $                                                                                                                            1               � $                                                           p          p            p                                                  �                                                          0           � $                                                          
            &                   &                   &                                                             !@                              '�              
      #ORDER    #NPTS     #LP_LEFT !   #LP_RIGHT "   #W_LEFT #   #W_RIGHT $   #DXI_LEFT %   #DXI_RIGHT &   #DX_LEFT '   #DX_RIGHT (               � $                                                                                                                              2                � $                                                                                                                                                � $                             !                            
  p          p            p                                                  �                          
                                 0.0            � $                             "                             
  p          p            p                                                  �                          
                                 0.0            � $                             #     8         
                                                 
                                 0.0                � $                             $     @         
                                                 
                                 0.0                � $                             %            H                
  p          p            p                                                  �                          
                                 0.0            � $                             &            `                
  p          p            p                                                  �                          
                                 0.0            � $                             '            x              	  
  p          p            p                                                  �                          
                                 0.0            � $                             (            �              
  
  p          p            p                                                  �                          
                                 0.0                 !                            )                                                                          !                            *                                                         #         @       !                            +                    #Z ,   #ULP -   #VLP .   #WLP /   #XLP 0   #Z_INT 1             
   @                             ,                   
              &                   &                   &                   &                                                     
   @                              -     �              #LPOLY_TYPE              
   @                              .     �              #LPOLY_TYPE              
   @                              /     �              #LPOLY_TYPE              
   @                              0     �              #LPOLY_TYPE              
  @                             1     
       #         @       !                            2                    #Z 3   #ULP 4   #VLP 5   #WLP 6   #XLP 7   #Z_TL 8   #ULP_TL 9   #VLP_TL :   #WLP_TL ;   #XLP_TL <   #Z_INT_TL =             
   @                             3                   
              &                   &                   &                   &                                                     
   @                              4     �              #LPOLY_TYPE              
   @                              5     �              #LPOLY_TYPE              
   @                              6     �              #LPOLY_TYPE              
   @                              7     �              #LPOLY_TYPE              
   @                             8                   
              &                   &                   &                   &                                                     
   @                              9     �              #LPOLY_TYPE              
   @                              :     �              #LPOLY_TYPE              
   @                              ;     �              #LPOLY_TYPE              
   @                              <     �              #LPOLY_TYPE              
  @                             =     
       #         @       !                            >                    #Z ?   #ULP @   #VLP A   #WLP B   #XLP C   #Z_INT_AD D   #Z_AD E   #ULP_AD F   #VLP_AD G   #WLP_AD H   #XLP_AD I             
   @                             ?                   
 &             &                   &                   &                   &                                                     
   @                              @     �              #LPOLY_TYPE              
   @                              A     �              #LPOLY_TYPE              
   @                              B     �              #LPOLY_TYPE              
   @                              C     �              #LPOLY_TYPE              
  @                             D     
                 
  @                             E                   
 '              &                   &                   &                   &                                                     
  @                              F     �               #LPOLY_TYPE              
  @                              G     �               #LPOLY_TYPE              
  @                              H     �               #LPOLY_TYPE              
  @                              I     �               #LPOLY_TYPE    #         @       !                            J                    #P K             
  @                              K     �               #LPOLY_TYPE    #         @       !                            L                    #X M   #X_INT N   #P O             
   @                             M                   
 ,             &                                                     
   @                             N     
                
  @                              O     �               #LPOLY_TYPE    #         @       !                            P                    #X Q   #X_INT R   #P S   #X_TL T   #X_INT_TL U   #P_TL V             
   @                             Q                   
 -             &                                                     
   @                             R     
                
   @                              S     �              #LPOLY_TYPE              
   @                             T                   
 .             &                                                     
   @                             U     
                
  @                              V     �               #LPOLY_TYPE    #         @       !                            W                    #X X   #X_INT Y   #P Z   #P_AD [   #X_AD \   #X_INT_AD ]             
   @                             X                   
 /             &                                                     
   @                             Y     
                
   @                              Z     �              #LPOLY_TYPE              
  @                              [     �               #LPOLY_TYPE              
  @                             \                   
 0              &                                                     
  @                             ]     
                      �  @                          ^     '�                   #WIND_SPEED _   #ZCOEFF_INVALID `   #SEC_Z a   #ZCOEFF_V b   #ZCOEFF_H c   #LSCI d               � D                             _               
                                                 
                                 0.0                � D                              `                                                                                  ��������                        � D                             a              
                                                 
                                 0.0                � D                             b                            
  p          p            p                                                  �                          
                                 0.0            � D                             c            H                
  p          p            p                                                  �                          
                                 0.0             � D                              d            x       �             #IINTERP_TYPE e   p          p            p                                         @  @                         e     '�                    #LP f   #I1 g   #I2 h   #OUTBOUND i   #XINT j   #X k                �                               f     �                      #LPOLY_TYPE                 �                               g     �                          �                               h     �                          �                               i     �                          �                              j     �          
                �                              k            �                 
  p          p            p                          #         @                                   l                    #LSCCOEFF m   #FREQUENCY n   #COS_Z o   #WIND_SPEED p   #RV_LARGE q   #RH_LARGE r   #IVAR s             
   @                              m     �              #FITCOEFF_3D_TYPE              
  @@                             n     
                
   @                             o     
                
   @                             p     
                D  @                             q     
                 D  @                             r     
                 
D @@                              s     �              #IVAR_TYPE ^   #         @                                   t                    #WIND_SPEED_TL u   #RV_LARGE_TL v   #RH_LARGE_TL w   #IVAR x             
   @                             u     
                D  @                             v     
                 D  @                             w     
                 
   @                              x     �             #IVAR_TYPE ^   #         @                                   y                    #RV_LARGE_AD z   #RH_LARGE_AD {   #WIND_SPEED_AD |   #IVAR }             
D  @                             z     
                 
D  @                             {     
                 
D  @                             |     
                 
   @                              }     �             #IVAR_TYPE ^      �   H      fn#fn 3   �   e   b   uapp(LARGE_SCALE_CORRECTION_MODULE    M  C   J  TYPE_KINDS     �  Q   J  FITCOEFF_DEFINE #   �  �   J  CRTM_INTERPOLATION 2   �  o       gen@FIND_INDEX+CRTM_INTERPOLATION 6   �  �      FIND_REGULAR_INDEX+CRTM_INTERPOLATION @   �  =      FIND_REGULAR_INDEX%SIZE+CRTM_INTERPOLATION=SIZE B   (  >      FIND_REGULAR_INDEX%FLOOR+CRTM_INTERPOLATION=FLOOR >   f  <      FIND_REGULAR_INDEX%MIN+CRTM_INTERPOLATION=MIN >   �  <      FIND_REGULAR_INDEX%MAX+CRTM_INTERPOLATION=MAX 8   �  �   e   FIND_REGULAR_INDEX%X+CRTM_INTERPOLATION 9   j  @   e   FIND_REGULAR_INDEX%DX+CRTM_INTERPOLATION <   �  @   e   FIND_REGULAR_INDEX%X_INT+CRTM_INTERPOLATION 9   �  @   e   FIND_REGULAR_INDEX%I1+CRTM_INTERPOLATION 9   *  @   e   FIND_REGULAR_INDEX%I2+CRTM_INTERPOLATION D   j  @   e   FIND_REGULAR_INDEX%OUT_OF_BOUNDS+CRTM_INTERPOLATION 5   �  �      FIND_RANDOM_INDEX+CRTM_INTERPOLATION ?   y  =      FIND_RANDOM_INDEX%SIZE+CRTM_INTERPOLATION=SIZE =   �  <      FIND_RANDOM_INDEX%MIN+CRTM_INTERPOLATION=MIN =   �  <      FIND_RANDOM_INDEX%MAX+CRTM_INTERPOLATION=MAX 7   .  �   e   FIND_RANDOM_INDEX%X+CRTM_INTERPOLATION ;   �  @   e   FIND_RANDOM_INDEX%X_INT+CRTM_INTERPOLATION 8   �  @   e   FIND_RANDOM_INDEX%I1+CRTM_INTERPOLATION 8   :	  @   e   FIND_RANDOM_INDEX%I2+CRTM_INTERPOLATION C   z	  @   e   FIND_RANDOM_INDEX%OUT_OF_BOUNDS+CRTM_INTERPOLATION 1   �	  �       FITCOEFF_3D_TYPE+FITCOEFF_DEFINE >   M
  �   a   FITCOEFF_3D_TYPE%IS_ALLOCATED+FITCOEFF_DEFINE 9   �
  �   a   FITCOEFF_3D_TYPE%RELEASE+FITCOEFF_DEFINE 9   �  �   a   FITCOEFF_3D_TYPE%VERSION+FITCOEFF_DEFINE <   ;  �   a   FITCOEFF_3D_TYPE%DIMENSIONS+FITCOEFF_DEFINE 3   4  �   a   FITCOEFF_3D_TYPE%C+FITCOEFF_DEFINE .   �  �       LPOLY_TYPE+CRTM_INTERPOLATION 4   �  �   a   LPOLY_TYPE%ORDER+CRTM_INTERPOLATION 3   n  �   a   LPOLY_TYPE%NPTS+CRTM_INTERPOLATION 6     �   a   LPOLY_TYPE%LP_LEFT+CRTM_INTERPOLATION 7     �   a   LPOLY_TYPE%LP_RIGHT+CRTM_INTERPOLATION 5     �   a   LPOLY_TYPE%W_LEFT+CRTM_INTERPOLATION 6   �  �   a   LPOLY_TYPE%W_RIGHT+CRTM_INTERPOLATION 7   V  �   a   LPOLY_TYPE%DXI_LEFT+CRTM_INTERPOLATION 8   Q  �   a   LPOLY_TYPE%DXI_RIGHT+CRTM_INTERPOLATION 6   L  �   a   LPOLY_TYPE%DX_LEFT+CRTM_INTERPOLATION 7   G  �   a   LPOLY_TYPE%DX_RIGHT+CRTM_INTERPOLATION    B  p       FP+TYPE_KINDS (   �  p       NPTS+CRTM_INTERPOLATION -   "  ~       INTERP_4D+CRTM_INTERPOLATION /   �  �   e   INTERP_4D%Z+CRTM_INTERPOLATION 1   t  X   e   INTERP_4D%ULP+CRTM_INTERPOLATION 1   �  X   e   INTERP_4D%VLP+CRTM_INTERPOLATION 1   $  X   e   INTERP_4D%WLP+CRTM_INTERPOLATION 1   |  X   e   INTERP_4D%XLP+CRTM_INTERPOLATION 3   �  @   e   INTERP_4D%Z_INT+CRTM_INTERPOLATION 0     �       INTERP_4D_TL+CRTM_INTERPOLATION 2   �  �   e   INTERP_4D_TL%Z+CRTM_INTERPOLATION 4   �  X   e   INTERP_4D_TL%ULP+CRTM_INTERPOLATION 4   �  X   e   INTERP_4D_TL%VLP+CRTM_INTERPOLATION 4   S  X   e   INTERP_4D_TL%WLP+CRTM_INTERPOLATION 4   �  X   e   INTERP_4D_TL%XLP+CRTM_INTERPOLATION 5     �   e   INTERP_4D_TL%Z_TL+CRTM_INTERPOLATION 7   �  X   e   INTERP_4D_TL%ULP_TL+CRTM_INTERPOLATION 7   /  X   e   INTERP_4D_TL%VLP_TL+CRTM_INTERPOLATION 7   �  X   e   INTERP_4D_TL%WLP_TL+CRTM_INTERPOLATION 7   �  X   e   INTERP_4D_TL%XLP_TL+CRTM_INTERPOLATION 9   7   @   e   INTERP_4D_TL%Z_INT_TL+CRTM_INTERPOLATION 0   w   �       INTERP_4D_AD+CRTM_INTERPOLATION 2   2!  �   e   INTERP_4D_AD%Z+CRTM_INTERPOLATION 4   "  X   e   INTERP_4D_AD%ULP+CRTM_INTERPOLATION 4   ^"  X   e   INTERP_4D_AD%VLP+CRTM_INTERPOLATION 4   �"  X   e   INTERP_4D_AD%WLP+CRTM_INTERPOLATION 4   #  X   e   INTERP_4D_AD%XLP+CRTM_INTERPOLATION 9   f#  @   e   INTERP_4D_AD%Z_INT_AD+CRTM_INTERPOLATION 5   �#  �   e   INTERP_4D_AD%Z_AD+CRTM_INTERPOLATION 7   z$  X   e   INTERP_4D_AD%ULP_AD+CRTM_INTERPOLATION 7   �$  X   e   INTERP_4D_AD%VLP_AD+CRTM_INTERPOLATION 7   *%  X   e   INTERP_4D_AD%WLP_AD+CRTM_INTERPOLATION 7   �%  X   e   INTERP_4D_AD%XLP_AD+CRTM_INTERPOLATION /   �%  O       CLEAR_LPOLY+CRTM_INTERPOLATION 1   )&  X   e   CLEAR_LPOLY%P+CRTM_INTERPOLATION )   �&  a       LPOLY+CRTM_INTERPOLATION +   �&  �   e   LPOLY%X+CRTM_INTERPOLATION /   n'  @   e   LPOLY%X_INT+CRTM_INTERPOLATION +   �'  X   e   LPOLY%P+CRTM_INTERPOLATION ,   (  �       LPOLY_TL+CRTM_INTERPOLATION .   �(  �   e   LPOLY_TL%X+CRTM_INTERPOLATION 2   )  @   e   LPOLY_TL%X_INT+CRTM_INTERPOLATION .   U)  X   e   LPOLY_TL%P+CRTM_INTERPOLATION 1   �)  �   e   LPOLY_TL%X_TL+CRTM_INTERPOLATION 5   9*  @   e   LPOLY_TL%X_INT_TL+CRTM_INTERPOLATION 1   y*  X   e   LPOLY_TL%P_TL+CRTM_INTERPOLATION ,   �*  �       LPOLY_AD+CRTM_INTERPOLATION .   T+  �   e   LPOLY_AD%X+CRTM_INTERPOLATION 2   �+  @   e   LPOLY_AD%X_INT+CRTM_INTERPOLATION .    ,  X   e   LPOLY_AD%P+CRTM_INTERPOLATION 1   x,  X   e   LPOLY_AD%P_AD+CRTM_INTERPOLATION 1   �,  �   e   LPOLY_AD%X_AD+CRTM_INTERPOLATION 5   \-  @   e   LPOLY_AD%X_INT_AD+CRTM_INTERPOLATION    �-  �       IVAR_TYPE %   A.  �   !   IVAR_TYPE%WIND_SPEED )   �.  �   !   IVAR_TYPE%ZCOEFF_INVALID     �/  �   !   IVAR_TYPE%SEC_Z #   30  �   !   IVAR_TYPE%ZCOEFF_V #   .1  �   !   IVAR_TYPE%ZCOEFF_H    )2  �   !   IVAR_TYPE%LSCI    �2  �      IINTERP_TYPE     ^3  `   a   IINTERP_TYPE%LP     �3  H   a   IINTERP_TYPE%I1     4  H   a   IINTERP_TYPE%I2 &   N4  H   a   IINTERP_TYPE%OUTBOUND "   �4  H   a   IINTERP_TYPE%XINT    �4  �   a   IINTERP_TYPE%X '   z5  �       LARGE_SCALE_CORRECTION 0    6  ^   a   LARGE_SCALE_CORRECTION%LSCCOEFF 1   ~6  @   a   LARGE_SCALE_CORRECTION%FREQUENCY -   �6  @   a   LARGE_SCALE_CORRECTION%COS_Z 2   �6  @   a   LARGE_SCALE_CORRECTION%WIND_SPEED 0   >7  @   a   LARGE_SCALE_CORRECTION%RV_LARGE 0   ~7  @   a   LARGE_SCALE_CORRECTION%RH_LARGE ,   �7  W   a   LARGE_SCALE_CORRECTION%IVAR *   8  �       LARGE_SCALE_CORRECTION_TL 8   �8  @   a   LARGE_SCALE_CORRECTION_TL%WIND_SPEED_TL 6   �8  @   a   LARGE_SCALE_CORRECTION_TL%RV_LARGE_TL 6   9  @   a   LARGE_SCALE_CORRECTION_TL%RH_LARGE_TL /   \9  W   a   LARGE_SCALE_CORRECTION_TL%IVAR *   �9  �       LARGE_SCALE_CORRECTION_AD 6   ::  @   a   LARGE_SCALE_CORRECTION_AD%RV_LARGE_AD 6   z:  @   a   LARGE_SCALE_CORRECTION_AD%RH_LARGE_AD 8   �:  @   a   LARGE_SCALE_CORRECTION_AD%WIND_SPEED_AD /   �:  W   a   LARGE_SCALE_CORRECTION_AD%IVAR 