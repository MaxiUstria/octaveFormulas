%# Copyright (C) 2013, Roberto Porcu' <roberto.porcu@polimi.it>
%# OdePkg - A package for solving ordinary differential equations and more
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%#
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with this program; If not, see <http://www.gnu.org/licenses/>.

%# -*- texinfo -*-
%# @deftypefn {Command} {[@var{J}] =} approx_Constr_grad (@var{@@fun}, @var{x})
%# 
%# This function file can be used to approximate the gradient of the constraints in a constrained Hamiltonian system solved with @command{odeRATTLE}.
%# 
%# If no explicit expression for the gradient is given to @command{odeRATTLE}, it calls this function to approximate it.
%# 
%# The output argument is a matrix where the first dimension refers to each different constraint.
%# 
%# The first input argument must be a function_handle or an inline function and must define the constraints.
%# 
%# The second input argument is just the point at which the gradient will be evaluated.
%#
%# @seealso{odepkg}
%# @end deftypefn

function J = approx_Constr_grad (f, q)

  %# evaluating the gradient
  f_q = f(q);

  %# determining the number of constraints
  constraints_nb = length(f_q);

  %# determining the size of the variable
  dim_q = length(q);

  %# initializing the output matrix
  J = zeros(constraints_nb,dim_q);

  %# step for approximation of the derivatives
  delta = sqrt(eps);

  %# approximating the derivatives
  for j=1:1:dim_q
    J(:,j) = (1/delta)*(f(q + delta.*[zeros(j-1,1);1;zeros(dim_q-j,1)]) - f_q);
  end

end

%!function [g] = constraint (vy)
%! g = (vy(:)')*vy(:) - 25.0;
%!end
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  res = approx_Constr_grad (@constraint, [5.0*cos(-pi/4);5.0*sin(-pi/4)]);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
