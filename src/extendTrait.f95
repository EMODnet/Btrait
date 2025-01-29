  SUBROUTINE extendtraitf( nRowTaxon, nColTaxon, Taxonomy,                      &
                           nRowTrait, nColTrait, Trait, Tnames,                 &
                           nMax, nNew, newTrait, newName, numLevs, Levs)

  IMPLICIT NONE
  ! dimensions of taxonomy
  INTEGER, INTENT (IN)   :: nColTaxon, nRowTaxon
  INTEGER, INTENT(INOUT) :: Taxonomy(nRowTaxon, nColTaxon) 
  
  
  ! dimensions of trait
  INTEGER, INTENT (IN)   :: nColTrait, nRowTrait
  INTEGER, INTENT (IN)   :: Tnames(nRowTrait)
  INTEGER, INTENT (OUT)  :: nNew
  DOUBLE PRECISION, INTENT(IN)  :: Trait(nRowTrait, nColTrait) 
  
  ! outputted traits
  INTEGER, INTENT (IN)   :: nMax
  DOUBLE PRECISION, INTENT(OUT) :: newTrait(nMax, nColTrait) 
  INTEGER, INTENT(OUT) :: newName(nMax), numLevs(nMax), Levs(nMax)
  
  ! locals
  DOUBLE PRECISION :: sumT(nColTrait) 
  INTEGER          :: TAX, I, J, K, L, M, nT, TinTrait, II
  LOGICAL          :: Known

! -----------------------------------------------------------------------------
! initialisation

  DO I = 1, nMax
    DO J = 1, nColTrait
      newtrait(I,J) = 0.D0
    END DO
  END DO 
  
! looping  
  
  nNew = 0                 ! number of traits that will be recorded in TraitNew
  
  DO J = 2, nColTaxon      ! columns in taxonomy except the first one CHNGED FROM 2 TO 1

    DO I = 1, nRowTaxon    ! loop over all rows of taxonomy
      TAX = Taxonomy(I,J)  ! taxon to extend with traits
      
      Known = .TRUE.
      IF (TAX .NE. 0) THEN ! it is a regular taxon that has not yet been processed
        Known = .FALSE.
        
        ! check if traits were not there to start with
        DO K = 1, nRowTrait  
           IF (TAX == Tnames(K)) THEN
              Known = .TRUE.
              EXIT  
           ENDIF
        END DO 
        
      END IF 
    
      IF (.NOT. known) THEN
        DO K = 1, nColtrait
          SumT (K) = 0.D0
        END DO
        
        nT = 0   ! Counter on number of traits that will be used 
        
        DO K = I, nRowtaxon  ! loop over all rows to find lower level taxa    
          
          IF (TAX == Taxonomy(K,J)) THEN   ! same taxon
          !!! CHECK: (K,1) or (K,I-1)?????
            TinTrait = Taxonomy(K,1)  ! taxon name to be found in trait data
            
            DO L = 1, nRowTrait
              IF (Tnames(L) == TinTrait) THEN  ! lower level taxon present
                Taxonomy(K,J) = 0              ! next time it won t be reused
                nT = nT + 1
                DO M = 1, nColTrait
                   SumT(M) = SumT(M)+Trait(L,M)
                END DO
                EXIT
              END IF
            END DO  ! with L
          END IF  ! same taxon
        
        END DO  ! with K
        
        IF (nT > 0) THEN                     ! traits have been estimated
           
          II = 0        ! if it has been estimated on lower level - overwrite
          DO M = 1, nNew
              IF (newName(M) == TAX) THEN
                 II = M
                 EXIT
              END IF 
          ENDDO
          IF (II == 0) THEN
            nNew          = nNew + 1
            newName(nNew) = TAX                ! Name of the taxon  
            II = nNew
          END IF

          NumLevs(II) = nT                 ! number of low-level traits on which based 
          Levs(II)    = J                  ! taxon level (column)
          
          DO M = 1, nColTrait
            newTrait(II, M) = SumT(M)/nT
          END DO
        END IF
            
      END IF ! not known          
    END DO ! with I
  END DO !  with JJ
  
END SUBROUTINE extendtraitf  
