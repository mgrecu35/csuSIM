	  	&  _   k820309    �          13.0        �wW                                                                                                           
       ACCoeff_Binary_IO.f90 ACCOEFF_BINARY_IO              ACCOEFF_BINARY_INQUIREFILE ACCOEFF_BINARY_READFILE ACCOEFF_BINARY_WRITEFILE ACCOEFF_BINARY_IOVERSION                      @                              
       LONG DOUBLE                      @                              
       FILE_OPEN FILE_EXISTS                      @                              
       SUCCESS FAILURE INFORMATION DISPLAY_MESSAGE                      @                              
       OPEN_BINARY_FILE                      @                              
       ACCOEFF_TYPE ACCOEFF_ASSOCIATED ACCOEFF_DESTROY ACCOEFF_CREATE ACCOEFF_VALIDRELEASE ACCOEFF_INFO                �                                       u #FILE_UNIT_EXISTS    #FILE_NAME_EXISTS    %         @    @                                                      #FILEID              
   @                                         %         @    @                                                     #FILENAME 	             
   @                             	                    1                �                                       u #FILE_OPEN_BY_UNIT 
   #FILE_OPEN_BY_NAME    %         @    @                          
                           #FILEID              
   @                                         %         @    @                                                     #FILENAME              
   @                                                 1                   !@               A                '�                   #IS_ALLOCATED    #RELEASE    #VERSION    #N_FOVS    #N_CHANNELS    #SENSOR_ID    #WMO_SATELLITE_ID    #WMO_SENSOR_ID    #SENSOR_CHANNEL    #A_EARTH    #A_SPACE    #A_PLATFORM                � $                                                                                                                                                 � $                                                                                                                            1                � $                                                                                                                            1                � $                                                                                                                            0                � $                                                                                                                            0                � $                                                                                                                     &              C                                                    � $                                  (                                                                           �              1023                � $                                  ,                                                                           �              2047               � $                                         0              	               &                                                      � $                                         x              
   
            &                   &                                                      � $                                         �                 
            &                   &                                                      � $                                         8                
            &                   &                                                            !                                                                                                      !                                                                                                      !                                                                                   0                 !                                                                                  3                 !                                                                                  1#         @ �     !                                               #DISPLAY_MESSAGE%TRIM !   #DISPLAY_MESSAGE%PRESENT "   #ROUTINE_NAME #   #MESSAGE $   #ERROR_STATE %   #MESSAGE_LOG &                 @                            !     TRIM               @                            "     PRESENT           
   @                             #                    1           
   @                             $                    1           
   @                              %                     
  @                             &                    1 %         @      !                          '                          #OPEN_BINARY_FILE%PRESENT (   #OPEN_BINARY_FILE%TRIM )   #FILENAME *   #FILEID +   #FOR_OUTPUT ,   #NO_CHECK -                 @                            (     PRESENT               @                            )     TRIM           
   @                             *                    1             @                              +                      
  @                              ,                     
  @                              -           %         H      !                         .                           #ACCOEFF /             
   @                              /     �             #ACCOEFF_TYPE    #         H       !                          0                    #ACCOEFF 1               @                              1     �              #ACCOEFF_TYPE    #         H       !                          2                    #ACCOEFF 3   #N_FOVS 4   #N_CHANNELS 5               @                              3     �              #ACCOEFF_TYPE              
   @                              4                     
   @                              5           %         @      !                          6                           #ACCOEFF 7             
   @                              7     �             #ACCOEFF_TYPE    #         @       !                           8                   #ACCOEFF_INFO%ACHAR 9   #ACCOEFF_INFO%MIN :   #ACCOEFF_INFO%LEN ;   #ACCOEFF_INFO%LEN_TRIM <   #ACCOEFF =   #INFO >                 @                            9     ACHAR               @                            :     MIN               @                            ;     LEN               @                            <     LEN_TRIM           
   @                              =     �             #ACCOEFF_TYPE                @                             >                     1 %         @                                 ?                          #ACCOEFF_BINARY_INQUIREFILE%PRESENT @   #ACCOEFF_BINARY_INQUIREFILE%TRIM A   #FILENAME B   #N_FOVS C   #N_CHANNELS D   #RELEASE E   #VERSION F   #SENSOR_ID G   #WMO_SATELLITE_ID H   #WMO_SENSOR_ID I                 @                            @     PRESENT               @                            A     TRIM           
@ @@                             B                    1           F @@                              C                      F @@                              D                      F @@                              E                      F @@                              F                      F @@                             G                     1           F @@                              H                      F @@                              I            %         @                                 J                          #ACCOEFF_BINARY_READFILE%PRESENT K   #ACCOEFF_BINARY_READFILE%TRIM L   #FILENAME M   #ACCOEFF N   #NO_CLOSE O   #QUIET P   #DEBUG Q                 @                            K     PRESENT               @                            L     TRIM           
@ @@                             M                    1           D @@                              N     �              #ACCOEFF_TYPE              
 @@                              O                     
 @@                              P                     
 @@                              Q           %         @                                 R                          #ACCOEFF_BINARY_WRITEFILE%PRESENT S   #ACCOEFF_BINARY_WRITEFILE%TRIM T   #FILENAME U   #ACCOEFF V   #NO_CLOSE W   #QUIET X   #DEBUG Y                 @                            S     PRESENT               @                            T     TRIM           
@ @@                             U                    1           
  @@                              V     �             #ACCOEFF_TYPE              
 @@                              W                     
 @@                              X                     
 @@                              Y           #         @                                   Z                    #ID [             D  @                             [                     1    �   0      fn#fn '   �   u   b   uapp(ACCOEFF_BINARY_IO    E  L   J  TYPE_KINDS    �  V   J  FILE_UTILITY     �  l   J  MESSAGE_HANDLER $   S  Q   J  BINARY_FILE_UTILITY    �  �   J  ACCOEFF_DEFINE -   E  l       gen@FILE_EXISTS+FILE_UTILITY .   �  \      FILE_UNIT_EXISTS+FILE_UTILITY 5     @   e   FILE_UNIT_EXISTS%FILEID+FILE_UTILITY .   M  ^      FILE_NAME_EXISTS+FILE_UTILITY 7   �  L   e   FILE_NAME_EXISTS%FILENAME+FILE_UTILITY +   �  n       gen@FILE_OPEN+FILE_UTILITY /   e  \      FILE_OPEN_BY_UNIT+FILE_UTILITY 6   �  @   e   FILE_OPEN_BY_UNIT%FILEID+FILE_UTILITY /     ^      FILE_OPEN_BY_NAME+FILE_UTILITY 8   _  L   e   FILE_OPEN_BY_NAME%FILENAME+FILE_UTILITY ,   �        ACCOEFF_TYPE+ACCOEFF_DEFINE 9   �  �   a   ACCOEFF_TYPE%IS_ALLOCATED+ACCOEFF_DEFINE 4   ]  �   a   ACCOEFF_TYPE%RELEASE+ACCOEFF_DEFINE 4   	  �   a   ACCOEFF_TYPE%VERSION+ACCOEFF_DEFINE 3   �	  �   a   ACCOEFF_TYPE%N_FOVS+ACCOEFF_DEFINE 7   L
  �   a   ACCOEFF_TYPE%N_CHANNELS+ACCOEFF_DEFINE 6   �
  �   a   ACCOEFF_TYPE%SENSOR_ID+ACCOEFF_DEFINE =   �  �   a   ACCOEFF_TYPE%WMO_SATELLITE_ID+ACCOEFF_DEFINE :   j  �   a   ACCOEFF_TYPE%WMO_SENSOR_ID+ACCOEFF_DEFINE ;     �   a   ACCOEFF_TYPE%SENSOR_CHANNEL+ACCOEFF_DEFINE 4   �  �   a   ACCOEFF_TYPE%A_EARTH+ACCOEFF_DEFINE 4   R  �   a   ACCOEFF_TYPE%A_SPACE+ACCOEFF_DEFINE 7   �  �   a   ACCOEFF_TYPE%A_PLATFORM+ACCOEFF_DEFINE "   �  p       DOUBLE+TYPE_KINDS       p       LONG+TYPE_KINDS (   �  q       SUCCESS+MESSAGE_HANDLER (   �  q       FAILURE+MESSAGE_HANDLER ,   l  q       INFORMATION+MESSAGE_HANDLER 0   �  �       DISPLAY_MESSAGE+MESSAGE_HANDLER :   �  =      DISPLAY_MESSAGE%TRIM+MESSAGE_HANDLER=TRIM @   �  @      DISPLAY_MESSAGE%PRESENT+MESSAGE_HANDLER=PRESENT =     L   e   DISPLAY_MESSAGE%ROUTINE_NAME+MESSAGE_HANDLER 8   f  L   e   DISPLAY_MESSAGE%MESSAGE+MESSAGE_HANDLER <   �  @   e   DISPLAY_MESSAGE%ERROR_STATE+MESSAGE_HANDLER <   �  L   e   DISPLAY_MESSAGE%MESSAGE_LOG+MESSAGE_HANDLER 5   >  �       OPEN_BINARY_FILE+BINARY_FILE_UTILITY E   �  @      OPEN_BINARY_FILE%PRESENT+BINARY_FILE_UTILITY=PRESENT ?   ?  =      OPEN_BINARY_FILE%TRIM+BINARY_FILE_UTILITY=TRIM >   |  L   e   OPEN_BINARY_FILE%FILENAME+BINARY_FILE_UTILITY <   �  @   e   OPEN_BINARY_FILE%FILEID+BINARY_FILE_UTILITY @     @   e   OPEN_BINARY_FILE%FOR_OUTPUT+BINARY_FILE_UTILITY >   H  @   e   OPEN_BINARY_FILE%NO_CHECK+BINARY_FILE_UTILITY 2   �  ]       ACCOEFF_ASSOCIATED+ACCOEFF_DEFINE :   �  Z   e   ACCOEFF_ASSOCIATED%ACCOEFF+ACCOEFF_DEFINE /   ?  U       ACCOEFF_DESTROY+ACCOEFF_DEFINE 7   �  Z   e   ACCOEFF_DESTROY%ACCOEFF+ACCOEFF_DEFINE .   �  q       ACCOEFF_CREATE+ACCOEFF_DEFINE 6   _  Z   e   ACCOEFF_CREATE%ACCOEFF+ACCOEFF_DEFINE 5   �  @   e   ACCOEFF_CREATE%N_FOVS+ACCOEFF_DEFINE 9   �  @   e   ACCOEFF_CREATE%N_CHANNELS+ACCOEFF_DEFINE 4   9  ]       ACCOEFF_VALIDRELEASE+ACCOEFF_DEFINE <   �  Z   e   ACCOEFF_VALIDRELEASE%ACCOEFF+ACCOEFF_DEFINE ,   �  �       ACCOEFF_INFO+ACCOEFF_DEFINE 8   �  >      ACCOEFF_INFO%ACHAR+ACCOEFF_DEFINE=ACHAR 4   �  <      ACCOEFF_INFO%MIN+ACCOEFF_DEFINE=MIN 4   (  <      ACCOEFF_INFO%LEN+ACCOEFF_DEFINE=LEN >   d  A      ACCOEFF_INFO%LEN_TRIM+ACCOEFF_DEFINE=LEN_TRIM 4   �  Z   e   ACCOEFF_INFO%ACCOEFF+ACCOEFF_DEFINE 1   �  L   e   ACCOEFF_INFO%INFO+ACCOEFF_DEFINE +   K        ACCOEFF_BINARY_INQUIREFILE 3   d  @      ACCOEFF_BINARY_INQUIREFILE%PRESENT 0   �  =      ACCOEFF_BINARY_INQUIREFILE%TRIM 4   �  L   a   ACCOEFF_BINARY_INQUIREFILE%FILENAME 2   -  @   a   ACCOEFF_BINARY_INQUIREFILE%N_FOVS 6   m  @   a   ACCOEFF_BINARY_INQUIREFILE%N_CHANNELS 3   �  @   a   ACCOEFF_BINARY_INQUIREFILE%RELEASE 3   �  @   a   ACCOEFF_BINARY_INQUIREFILE%VERSION 5   -  L   a   ACCOEFF_BINARY_INQUIREFILE%SENSOR_ID <   y  @   a   ACCOEFF_BINARY_INQUIREFILE%WMO_SATELLITE_ID 9   �  @   a   ACCOEFF_BINARY_INQUIREFILE%WMO_SENSOR_ID (   �  �       ACCOEFF_BINARY_READFILE 0   �   @      ACCOEFF_BINARY_READFILE%PRESENT -   !  =      ACCOEFF_BINARY_READFILE%TRIM 1   L!  L   a   ACCOEFF_BINARY_READFILE%FILENAME 0   �!  Z   a   ACCOEFF_BINARY_READFILE%ACCOEFF 1   �!  @   a   ACCOEFF_BINARY_READFILE%NO_CLOSE .   2"  @   a   ACCOEFF_BINARY_READFILE%QUIET .   r"  @   a   ACCOEFF_BINARY_READFILE%DEBUG )   �"  �       ACCOEFF_BINARY_WRITEFILE 1   �#  @      ACCOEFF_BINARY_WRITEFILE%PRESENT .   �#  =      ACCOEFF_BINARY_WRITEFILE%TRIM 2   $  L   a   ACCOEFF_BINARY_WRITEFILE%FILENAME 1   S$  Z   a   ACCOEFF_BINARY_WRITEFILE%ACCOEFF 2   �$  @   a   ACCOEFF_BINARY_WRITEFILE%NO_CLOSE /   �$  @   a   ACCOEFF_BINARY_WRITEFILE%QUIET /   -%  @   a   ACCOEFF_BINARY_WRITEFILE%DEBUG )   m%  P       ACCOEFF_BINARY_IOVERSION ,   �%  L   a   ACCOEFF_BINARY_IOVERSION%ID 