! *************************************************************************************
! MOD_MATRIX_DOUBLING contains the algorithms that compute the matrix doubling among
! vegetation sublayers and between vegetation and soil. Both the procedure are called
! throught MATRIX_DOUBLING interface
!
! MATRIX_DOUBLING_VEGETATION
! MATRIX_DOUBLING_SOIL
! *************************************************************************************
MODULE MOD_MATRIX_DOUBLING

USE MOD_GLOBAL_PARAMETERS

IMPLICIT NONE

PRIVATE

INTERFACE MATRIX_DOUBLING
  MODULE PROCEDURE MATRIX_DOUBLING_MULTIPLE_VEGETATION ,&
									 MATRIX_DOUBLING_DOUBLE_VEGETATION ,&
                   MATRIX_DOUBLING_SOIL
END INTERFACE 

PUBLIC MATRIX_DOUBLING 

CONTAINS

	SUBROUTINE MATRIX_DOUBLING_MULTIPLE_VEGETATION(sv,tv,N_iter)
	! Dummy variables declaration
	REAL, DIMENSION(:,:), INTENT(INOUT) :: sv,tv
	INTEGER, INTENT(IN)									:: N_iter
	! Local variables declaration
	REAL, DIMENSION(NN,NN)	:: id,ami,am1,am2,sv1,tv1
	REAL, DIMENSION(NN)			:: work
	INTEGER	:: info,i

		id(:,:) = 0.
		FORALL(i = 1:NN) id(i,i) = 1.

		DO i = 1,N_iter
			am1(:,:) = id(:,:) - MATMUL(sv(:,:),sv(:,:))
			ami(:,:) = am1(:,:)
			CALL RINV(NN,ami,NN,work,info)
			IF(info /= 0) STOP 'Error during matrix inversion in MATRIX_DOUBLING_VEGETATION function'
			am1(:,:) = MATMUL(tv(:,:),sv(:,:))
			am2(:,:) = MATMUL(am1(:,:),ami(:,:))
			sv1(:,:) = MATMUL(am2(:,:),tv(:,:))
			am1(:,:) = MATMUL(tv(:,:),ami(:,:))
			tv1(:,:) = MATMUL(am1(:,:),tv(:,:))
			sv(:,:) = sv1(:,:) + sv(:,:)
			tv(:,:) = tv1(:,:)
		END DO

		RETURN
	END SUBROUTINE

	SUBROUTINE MATRIX_DOUBLING_DOUBLE_VEGETATION(s_up,t_up,s_down,t_down)
	! Dummy variables declaration
	REAL, DIMENSION(:,:), INTENT(INOUT) :: s_up,t_up
	REAL, DIMENSION(:,:), INTENT(IN)	  :: s_down,t_down
	! Local variables declaration
	REAL, DIMENSION(NN,NN)	:: id,ami,am1,am2,sv1,tv1
	REAL, DIMENSION(NN)			:: work
	INTEGER	:: info,i

		id(:,:) = 0.
		FORALL(i = 1:NN) id(i,i) = 1.
		am1(:,:) = id(:,:) - MATMUL(s_up(:,:),s_down(:,:))
		ami(:,:) = am1(:,:)
		CALL RINV(NN,ami,NN,work,info)
		IF(info /= 0) STOP 'Error during matrix inversion in MATRIX_DOUBLING_VEGETATION function'
		am1(:,:) = MATMUL(t_up(:,:),s_down(:,:))
		am2(:,:) = MATMUL(am1(:,:),ami(:,:))
		sv1(:,:) = MATMUL(am2(:,:),t_up(:,:))
		am1(:,:) = MATMUL(t_down(:,:),ami(:,:))
		tv1(:,:) = MATMUL(am1(:,:),t_up(:,:))
		s_up(:,:) = sv1(:,:) + s_up(:,:)
		t_up(:,:) = tv1(:,:)

		RETURN
	END SUBROUTINE

	SUBROUTINE MATRIX_DOUBLING_SOIL(sv_star,tv_star,tv,sg0,svg1)
	! Dummy variables declaration
	REAL, DIMENSION(:,:), INTENT(INOUT) :: svg1
	REAL, DIMENSION(:,:), INTENT(IN)	  :: sv_star,tv_star,tv,sg0

	! Local variables declaration
	REAL, DIMENSION(NN,NN)	:: id,ami,am1,am2
	REAL, DIMENSION(NN)			:: work
	INTEGER	:: info,i

		id(:,:) = 0.
		FORALL(i = 1:NN) id(i,i) = 1.
		am1(:,:) = id(:,:) - MATMUL(sv_star(:,:),sg0(:,:))
		ami(:,:) = am1(:,:)
		CALL RINV(NN,ami,NN,work,info)
		IF(info /= 0) STOP 'Error during matrix inversion in MATRIX_DOUBLING_VEGETATION function'
		am1(:,:) = MATMUL(sg0(:,:),ami(:,:))
		am2(:,:) = MATMUL(tv_star(:,:),am1(:,:))
		svg1(:,:) = MATMUL(am2(:,:),tv(:,:))            ! Interaction between vegetation and soil

		RETURN
	END SUBROUTINE MATRIX_DOUBLING_SOIL
END MODULE MOD_MATRIX_DOUBLING