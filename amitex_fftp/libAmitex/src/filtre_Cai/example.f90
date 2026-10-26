program fortranfilter

  use caifilters_mod

  implicit none

  double precision :: a,k
  
  write(*,*) ' w(a,k) examples: '
  a=0
  k=5
  write(*,*) ' wk(0.0, 5.0) = ',wk(a,k)
  
  write(*,*) ' w(a,k) examples: '
  a=1
  k=5
  write(*,*) ' wk(1.0, 5.0) = ',wk(a,k)
  k=10
  write(*,*) ' wk(1.0,10.0) = ',wk(a,k)
  k=50
  write(*,*) ' wk(1.0,50.0) = ',wk(a,k)
  a=3
  k=5
  write(*,*) ' wk(3.0, 5.0) = ',wk(a,k)
  k=10
  write(*,*) ' wk(3.0,10.0) = ',wk(a,k)
  k=50
  write(*,*) ' wk(3.0,50.0) = ',wk(a,k)  
  write(*,*)
  
  write(*,*) ' wtilde (a,k) examples: '
  a=1
  k=5
  write(*,*) ' wtildek(1.0, 5.0) = ',wtildek(a,k)
  k=10
  write(*,*) ' wtildek(1.0,10.0) = ',wtildek(a,k)
  k=50
  write(*,*) ' wtildek(1.0,50.0) = ',wtildek(a,k)
  a=3
  k=5
  write(*,*) ' wtildek(3.0, 5.0) = ',wtildek(a,k)
  k=10
  write(*,*) ' wtildek(3.0,10.0) = ',wtildek(a,k)
  k=50
  write(*,*) ' wtildek(3.0,50.0) = ',wtildek(a,k)
  
end program fortranfilter
