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
%# @deftypefn {Command} {[@var{J}] =} approx_Constr_Hess (@var{@@fun}, @var{x})
%# 
%# This function file can be used to approximate the Hessian of the constraints in a constrained Hamiltonian system solved with @command{odeRATTLE}.
%# 
%# If no explicit expression for the Hessian is given to @command{odeRATTLE}, it calls this function to approximate it.
%# 
%# The output argument is a three-dimensional matrix where the last dimension refers to each different constraint.
%# 
%# The first input argument must be a function_handle or an inline function and must define the gradient of the constraints.
%# 
%# The second input argument is just the point at which the Hessian will be evaluated.
%#
%# @seealso{odepkg}
%# @end deftypefn

function H = approx_Constr_Hess(f,q)

  %# evaluating the gradient
  f_q = f(q);

  %# determining the number of constraints
  [constraints_nb,dim_f] = size(f_q);

  %# determining the size of the variable
  dim_q = length(q);

  %# if the gradient has wrong column dimension, it is an error
  if( dim_f ~= dim_q )
    error('OdePkg:InvalidArgument', ...
      'gradient of constrain function has wrong dimensions');
  end

  %# initializing the three-dimensional matrix
  H = zeros(dim_f,dim_q,constraints_nb);

  %# step for approximation of the derivatives
  delta = sqrt(eps);

  %# approximating the derivatives
  for k = 1:1:constraints_nb
    for j = 1:1:dim_q
      H(:,j,k) = (1/delta)*((f(q + delta.*[zeros(j-1,1);1;zeros(dim_q-j,1)])(k,:)) - f_q(k,:))';
    end
  end
end

%!function [G] = constraint_grad (vy)
%!  G = 2*vy(:)';
%!end
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  res = approx_Constr_Hess (@constraint_grad, [5.0*cos(-pi/4);5.0*sin(-pi/4)]);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
