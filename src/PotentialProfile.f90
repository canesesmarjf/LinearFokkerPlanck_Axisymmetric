subroutine PotentialProfile(iPotential)
!
!     This subroutine is called to assign
!     a model (pre-determined) potential profile. In the future this
!     may be calulated more self-consistently
!
!*******************************************************************
!
USE local
USE plasma_params
USE spline_fits     
implicit none

INTEGER :: i
LOGICAL :: iPotential

if (iPotential) then
    do i=1,nz
        if(z_Ref(i) .ge. 0. .and. z_Ref(i) .lt. s1) then
            Phi_Ref(i) = phi1
        else if(z_Ref(i) .ge. s1 .and. z_Ref(i) .lt. s2) then
            Phi_Ref(i) = phi1 - 2.*(phi2 - phi1)*((z_Ref(i) - s1)**3)/((s2 - s1)**3) &
            + 3.*(phi2 - phi1)*((z_Ref(i) - s1)**2)/((s2 - s1)**2)
        else if(z_Ref(i) .ge. s2 .and. z_Ref(i) .lt. s3) then
            Phi_Ref(i) = phi3 - 2.*(phi2 - phi3)*((z_Ref(i) - s3)**3)/((s2 - s3)**3) &
            + 3.*(phi2 - phi1)*((z_Ref(i) - s3)**2)/((s2 - s3)**2)
        else if(z_Ref(i) .ge. s3) then
            Phi_Ref(i) = phi3
        endif
        !	write(*,*) z_Ref(i), Phi_Ref(i)
    end do
else
    do i = 1,nz
        Phi_Ref(i) = 0
    end do
end if

return
end