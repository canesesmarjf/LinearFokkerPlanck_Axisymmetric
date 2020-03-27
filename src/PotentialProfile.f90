subroutine PotentialProfile(spline0,in0)
USE local
USE spline_fits
USE dataTYP

IMPLICIT NONE
TYPE(splTYP) :: spline0
TYPE(inTYP) :: in0
REAL(r8) :: x1, x2, x3
REAL(r8) :: p1, p2, p3
INTEGER :: i

x1 = in0%s1
x2 = in0%s2
x3 = in0%s3
p1 = in0%phi1
p2 = in0%phi2
p3 = in0%phi3

do i=1,in0%nz
    if(spline0%x(i) .ge. 0. .and. spline0%x(i) .lt. x1) then
        spline0%y(i) = p1
    else if(spline0%x(i) .ge. x1 .and. spline0%x(i) .lt. x2) then
        spline0%y(i) = p1 - 2.*(p2 - p1)*((spline0%x(i) - x1)**3)/((x2 - x1)**3) &
        + 3.*(p2 - p1)*((spline0%x(i) - x1)**2)/((x2 - x1)**2)
    else if(spline0%x(i) .ge. x2 .and. spline0%x(i) .lt. x3) then
        spline0%y(i) = p3 - 2.*(p2 - p3)*((spline0%x(i) - x3)**3)/((x2 - x3)**3) &
        + 3.*(p2 - p1)*((spline0%x(i) - x3)**2)/((x2 - x3)**2)
    else if(spline0%x(i) .ge. x3) then
        spline0%y(i) = p3
    endif
end do

return
end
