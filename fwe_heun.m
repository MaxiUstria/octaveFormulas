%# Copyright (C) 2015, Jacopo Corno <jacopo.corno@gmail.com>
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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} fwe_heun (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{err}] =} fwe_heun (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%#
%# This function can be used to integrate a set of non--stiff ordinary
%# differential equations (non--stiff ODEs) with a given initial condition
%# @var{x} from @var{t} to @var{t+dt}, with the Forward-Euler method using the
%# Heun-Euler method to estimate the error. For the definition of this method
%# see @url{http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods}.
%#
%# First output argument is the final integration time value.
%#
%# Second output parameter is the computed solution at time @var{t_next}.
%#
%# Third output parameter is a higher order solution for the estimation of the
%# error.
%#
%# First input argument is the function describing the system of ODEs to be
%# integrated.
%#
%# Second input parameter is the first extreme of integration interval.
%#
%# Third input argument is the initial condition of the system.
%#
%# Fourth input argument is the timestep, that is the length of the integration
%# interval.
%#
%# Fifth input parameter is optional and describes a set of options useful to
%# adapt the computation to what is needed.
%#
%# @seealso{odepkg, ode12}
%# @end deftypefn

function [t_next, x_next, x_est, k] = fwe_heun (f, t, x, dt, options = [],
                                                k_vals = [], t_next = t + dt)
                                                  
  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end
  
  k = zeros (size (x), 2);

  %# one-step scheme as defined in <http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods>
  k(:,1) = feval (f, t, x, args);
  k(:,2) = feval (f, t + dt, x + dt*k(:,1), args);

  %# computing new time and new values for the unkwnowns
  t_next = t + dt;
  x_next = x + dt*k(:,1);

  %# if the estimation of the error is required
  if(nargout >= 3)
    %# new solution to be compared with the previous one
    x_est = x + dt*((1/2)*k(:,1) + (1/2)*k(:,2));
  end

end

%! # We are using the "Van der Pol" implementation.
%!function [ydot] = fpol (vt, vy) %# The Van der Pol
%!  ydot = [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!end
%!
%!shared opts
%!  opts = odeset;
%!  opts.funarguments = {};
%!
%!test
%!  [t,y] = fwe_heun (@fpol, 0, [2;0], 0.05, opts);
%!test
%!  [t,y,x] = fwe_heun (@fpol, 0, [2;0], 0.05, opts);

%# Local Variables: ***
%# mode: octave ***
%# End: ***
