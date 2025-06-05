MODULE OVL_MODULE
CHARACTER(256) :: FILE
LOGICAL(4) :: TINIT=.FALSE.
LOGICAL(4) :: TWAKE=.FALSE.
END MODULE OVL_MODULE
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OVL$SETL4(ID,VAL)
!      *************************************************************************
       USE OVL_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4),INTENT(IN) :: VAL
!      *************************************************************************
       IF(ID.EQ.'INIT') THEN
         TINIT=VAL
       ELSE IF(ID.EQ.'WAKE') THEN
         TWAKE=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OVL$SETL4')
       END IF
       RETURN
       END SUBROUTINE OVL$SETL4
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OVL$GETL4(ID,VAL)
!      *************************************************************************
       USE OVL_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       LOGICAL(4),INTENT(OUT) :: VAL
!      *************************************************************************
       IF(ID.EQ.'INIT') THEN
         VAL=TINIT
       ELSE IF(ID.EQ.'WAKE') THEN
         VAL=TWAKE
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OVL$GETL4')
       END IF
       RETURN
       END SUBROUTINE OVL$GETL4
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OVL$SETCH(ID,VAL)
!      *************************************************************************
       USE OVL_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       CHARACTER(*),INTENT(IN) :: VAL
!      *************************************************************************
       IF(ID.EQ.'FILE') THEN
         FILE=VAL
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OVL$SETCH')
       END IF
       RETURN
       END SUBROUTINE OVL$SETCH
!
!      ..1.........2.........3.........4.........5.........6.........7.........8
       SUBROUTINE OVL$GETCH(ID,VAL)
!      *************************************************************************
       USE OVL_MODULE
       IMPLICIT NONE
       CHARACTER(*),INTENT(IN) :: ID
       CHARACTER(*),INTENT(OUT) :: VAL
!      *************************************************************************
       IF(ID.EQ.'FILE') THEN
         VAL=FILE
       ELSE
         CALL ERROR$MSG('ID NOT RECOGNIZED')
         CALL ERROR$CHVAL('ID',ID)
         CALL ERROR$STOP('OVL$GETCH')
       END IF
       RETURN
       END SUBROUTINE OVL$GETCH