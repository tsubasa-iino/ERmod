#define assert(x) if(.not. x) then; \
print *, ("Assertion "//#x//" failed!");\
stop;\
end if
#define asserteq(x,y) if(x /= y) then; \
print *, "Failure: "//#x//" /= "//#y;\
print *, "Expected: ", x;\
print *, "but was: ", y;\
stop;\
end if
#define assertnear_t(x,y,THRESHOLD) if(abs(x - y) > abs(x) * THRESHOLD .or. abs(x - y) > abs(y) * THRESHOLD) then; \
print *, "Failure: "//#x//" !~= "//#y//" at line: ", __LINE__;\
print *, "Expected: ", x;\
print *, "but was: ", y;\
stop;\
end if
#define assertnear(x,y) assertnear_t(x,y,1e-4)

program spline_test
  use spline

  call spline_init(2)
  assertnear(0.0, spline_value(0.0d0))
  assertnear(0.0, spline_value(2.0d0))
  assertnear(1.0, spline_value(1.0d0))
  assertnear(0.5, spline_value(0.5d0))
  assertnear(0.5, spline_value(1.5d0))
  call spline_cleanup()

  call spline_init(6)
  assertnear(0.0, spline_value(0.0d0))
  assertnear(0.0, spline_value(6.0d0))
  assertnear(2.604166666666667E-004, spline_value(0.5d0))
  assertnear(6.171875000000000E-002, spline_value(1.5d0))
  assertnear(0.438020833333333E+000, spline_value(2.5d0))
  assertnear(0.438020833333333E+000, spline_value(3.5d0))
  assertnear(6.171875000000000E-002, spline_value(4.5d0))
  assertnear(2.604166666666667E-004, spline_value(5.5d0))
  call spline_cleanup()
end program spline_test
