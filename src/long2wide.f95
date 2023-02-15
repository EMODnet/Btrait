! ===============================================================
! Converts descriptor, taxa, value data to wide format,
! summing value data for one taxon and descriptor, and
! taking the averages over nAvg replicates 
! ===============================================================

SUBROUTINE long2wide(ndat, ntaxa, nstat,                   &   ! dimensions
                     nAvg,                                 &   ! input
                     descriptor, taxon, value,             &   ! input
				          	 wide)                                     ! output
IMPLICIT NONE

INTEGER, INTENT(IN) :: ndat, nstat, ntaxa  ! number of data, stations, taxa
INTEGER, INTENT(IN) :: nAvg(nstat)         ! Number of replicates for each station
INTEGER, INTENT(IN) :: taxon(ndat)         ! taxonomic identity
INTEGER, INTENT(IN) :: descriptor(ndat)    ! descriptor identity
DOUBLE PRECISION, INTENT(IN)  :: value(ndat)  ! the actual value
DOUBLE PRECISION, INTENT(OUT) :: wide(nstat,ntaxa)

INTEGER :: I, itax, idesc

wide(:,:) = 0.D0

! Calculate the summed value for a descriptor-taxon pair  
DO I = 1, ndat
  idesc = descriptor(I)
  itax  = taxon(I)
  wide(idesc, itax) = wide(idesc, itax) + value(I)
END DO

! Divide by the number of replicates over which to average
DO I = 1, nstat
  wide(I,:) = wide(I, :)/nAvg(I)
END DO

END SUBROUTINE long2wide 

! ===============================================================
! Converts descriptor, taxa, value data by,
! *summing* value data for same taxon and descriptors, and
! taking the averages over nAvg replicates 
! ===============================================================

SUBROUTINE long2long(ndat, ntaxa, nstat, nAvg,             &   ! input
                     descriptor, taxon, value,             &   ! in/output
				          	 wide, nout)                               ! output
IMPLICIT NONE

INTEGER, INTENT(IN)    :: ndat, nstat, ntaxa
INTEGER, INTENT(IN)    :: nAvg(nstat)          ! Number of replicates 
INTEGER, INTENT(INOUT) :: taxon(ndat)          ! taxonomic identity
INTEGER, INTENT(INOUT) :: descriptor(ndat)     ! descriptor identity (sorted)
DOUBLE PRECISION, INTENT(INOUT) :: value(ndat) ! the actual value
INTEGER, INTENT(INOUT)          :: nout        ! the number of relevant data
DOUBLE PRECISION, INTENT(INOUT) :: wide(ntaxa)

INTEGER :: I, J, itax, idesc

idesc   = descriptor(1)
DO J = 1, ntaxa               
  wide(J) = 0.D0 
END DO
nout    = 0

DO I = 1, ndat+1
  
  ! do this at the end of each descriptor
  IF (descriptor(I) .NE. idesc .OR. I .EQ. ndat+1) THEN  

  	DO J = 1, ntaxa                ! save data (nonzeros) in vectors
	   IF (wide(J) .GT. 0.D0) THEN 
	     nout = nout+1
		   descriptor(nout) = idesc
		   taxon(nout)      = J
		   value(nout)      = wide(J)/nAvg(idesc)  ! standardize for number replicates
	   END IF  
	  END DO  
	  IF (I .GT. ndat) EXIT
    DO J = 1, ntaxa                
      wide(J) = 0.D0 
    END DO
    idesc = descriptor(I)
  ENDIF

  itax = taxon(I)
  wide(itax) = wide(itax) + value(I)
END DO

END SUBROUTINE long2long 

! ===============================================================
! Converts descriptor, taxa, value data by,
! *averaging* value/taxon data for descriptors, and
! taking the averages over nAvg replicates TO UPDATE!!
! ===============================================================

SUBROUTINE long2longMean(ndat, ntaxa, nstat, nAvg,         &   ! input
                     descriptor, taxon, toAverage, value,  &   ! in/output
				          	 replicate, wide, nout)                    ! output
IMPLICIT NONE

INTEGER, INTENT(IN)    :: ndat, nstat, ntaxa
INTEGER, INTENT(IN)    :: nAvg(nstat)         ! Number of replicates 
INTEGER, INTENT(IN)    :: replicate(ndat)     ! Replicate ID
INTEGER, INTENT(INOUT) :: taxon(ndat)         ! taxonomic identity
INTEGER, INTENT(INOUT) :: descriptor(ndat)    ! descriptor identity (sorted)
INTEGER, INTENT(IN)    :: toAverage(ndat)     ! replicates 

DOUBLE PRECISION, INTENT(INOUT) :: value(ndat)  ! the actual value
INTEGER,          INTENT(INOUT) :: nout         ! the number of relevant data
DOUBLE PRECISION, INTENT(INOUT) :: wide(ntaxa)

INTEGER :: I, J, itax, idesc
INTEGER :: ntax(ntaxa)

idesc   = descriptor(1)
DO J = 1, ntaxa               
  wide(J) = 0.D0 
  ntax(J) = 0 
END DO

nout    = 0
DO I = 1, ndat+1
  
  ! do this at the end of each descriptor
  IF (descriptor(I) .NE. idesc .OR. I .EQ. ndat+1) THEN  

  	DO J = 1, ntaxa                ! save data (nonzeros) in vectors
	   IF (wide(J) .GT. 0.D0) THEN 
	     nout = nout+1
		   descriptor(nout) = idesc
		   taxon(nout)      = J
		    
		   ! standardize for number replicates and take mean per taxon
		   value(nout)      = wide(J)/nAvg(idesc)/ntax(j)  
	   END IF  
	  END DO  
	  IF (I .GT. ndat) EXIT
    DO J = 1, ntaxa                
      wide(J) = 0.D0 
      ntax(J) = 0 
    END DO
    idesc = descriptor(I)
  ENDIF

  itax = taxon(I)
  wide(itax) = wide(itax) + value(I)
  ntax(itax) = ntax(itax) + 1 

END DO

END SUBROUTINE long2longMean 

