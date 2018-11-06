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
%# @deftypefn {Command} {[@var{J}] =} approx_Hamilt_Hess (@var{@@fun}, @var{t}, @var{x})
%# 
%# This function file can be used to approximate the Hessian of the Hamiltonian in a generic Hamiltonian system solved with @command{odeSE}, @command{odeSV}, @command{odeSPVI}, @command{odeRATTLE}.
%# 
%# If no explicit expression for the Hessian is given to the solver, it calls this function to approximate it.
%# 
%# The output argument is a square matrix. 
%# 
%# The first input argument must be a function_handle or an inline function and must define the Hamilton's equations of motion:
%# @ifhtml
%# @example
%# @math{q' = dH/dp (t,[q;p])}
%# @math{p' = - dH/dq (t,[q;p])}
%# @end example
%# @end ifhtml
%# @ifnothtml
%# @math{q' = dH/dp (t,[q;p])}
%# @math{p' = - dH/dq (t,[q;p])}.
%# @end ifnothtml
%# 
%# The second and third input arguments are just the time and the point at which the Hessian will be evaluated.
%#
%# @seealso{odepkg}
%# @end deftypefn

function H = approx_Hamilt_Hess (f, t, x)

  %# evaluating the Hamilton's equations
  f_x = f(t,x);

  %# determining the size of the unknowns
  dim_f_half = length(f_x)/2;
  dim_x = length(x);

  %# if Hamilton's equations have wrong dimensions, it is an error
  if( 2*dim_f_half ~= dim_x )
    error('OdePkg:InvalidArgument', ...
      'first argument returns a wrong number of equations');
  end

  %# initializing the output matrix
  H = zeros(2*dim_f_half,dim_x);

  %# step for approximation of the derivatives
  delta = sqrt(eps);

  %# approximating the derivatives
  for j = 1:1:dim_x
    H(1:dim_f_half,j) = -(1/delta)*(f(t,x+delta*[zeros(j-1,1);1;zeros(dim_x-j,1)])(dim_f_half+1:end) - f_x(dim_f_half+1:end));
    H(dim_f_half+1:end,j) = (1/delta)*(f(t,x+delta*[zeros(j-1,1);1;zeros(dim_x-j,1)])(1:dim_f_half) - f_x(1:dim_f_half));
  end
end

%!function [ydot] = armonic_oscillator (vt,vy)
%!  ydot = [vy(length(vy)/2+1:end);-vy(1:length(vy)/2)];
%!end
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  res = approx_Hamilt_Hess (@armonic_oscillator, 0, [1;0]);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
