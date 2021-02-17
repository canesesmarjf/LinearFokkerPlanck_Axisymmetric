subroutine PotentialProfile(fieldspline,params)
USE local
USE spline_fits
USE dataTYP

IMPLICIT NONE
TYPE(fieldSplineTYP) :: fieldspline
TYPE(paramsTYP) :: params
REAL(r8) :: x1, x2, x3
REAL(r8) :: p1, p2, p3
INTEGER :: i

x1 = params%s1
x2 = params%s2
x3 = params%s3
p1 = params%phi1
p2 = params%phi2
p3 = params%phi3

do i=1,params%nz
    if(fieldspline%V%x(i) .ge. 0. .and. fieldspline%V%x(i) .lt. x1) then
        fieldspline%V%y(i) = p1
    else if(fieldspline%V%x(i) .ge. x1 .and. fieldspline%V%x(i) .lt. x2) then
        fieldspline%V%y(i) = p1 - 2.*(p2 - p1)*((fieldspline%V%x(i) - x1)**3)/((x2 - x1)**3) &
        + 3.*(p2 - p1)*((fieldspline%V%x(i) - x1)**2)/((x2 - x1)**2)
    else if(fieldspline%V%x(i) .ge. x2 .and. fieldspline%V%x(i) .lt. x3) then
        fieldspline%V%y(i) = p3 - 2.*(p2 - p3)*((fieldspline%V%x(i) - x3)**3)/((x2 - x3)**3) &
        + 3.*(p2 - p1)*((fieldspline%V%x(i) - x3)**2)/((x2 - x3)**2)
    else if(fieldspline%V%x(i) .ge. x3) then
        fieldspline%V%y(i) = p3
    endif
end do

return
end
