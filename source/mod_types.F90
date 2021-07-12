!> module for variable declaration
module mod_types
  implicit none
  integer, parameter :: sp = kind(1.0)
  !! Single precision
  integer, parameter :: dp = selected_real_kind(2*precision(1.0_sp))
  !! Double precision
  integer, parameter :: qp = selected_real_kind(2*precision(1.0_dp))
  !! Quadruple precision
end module mod_types


