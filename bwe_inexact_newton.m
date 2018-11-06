%# Copyright (C) 2013, Simone Panza
%# Copyright (C) 2013, Fabio Cisaria
%# Copyright (C) 2013, Carlo de Falco
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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} bwe_inexact_newton (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} bwe_inexact_newton (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# 
%# This function can be used to integrate a set of non--stiff ordinary differential equations (non--stiff ODEs) with a given initial condition @var{x} from @var{t} to @var{t+dt}, with the Backward-Euler method using an inexact Newton method to solve the nonlinear system and estimating the error comparing the solution to that one given in two substeps by the Crank-Nicolson method.
%# 
%# First output argument is the final integration time value. 
%# 
%# Second output parameter is the computed solution at time @var{t_next}.
%# 
%# Third output parameter is a solution of higher order for the estimation of the error.
%# 
%# First input argument is the function describing the system of ODEs to be integrated.
%# 
%# Second input parameter is the first extreme of integration interval.
%# 
%# Third input argument is the initial condition of the system.
%# 
%# Fourth input argument is the timestep, that is the length of the integration interval.
%# 
%# Fifth input parameter is optional and describes a set of options useful to adapt the computation to what is needed. Possible fields are
%# @table @option
%# @item  @code{UseJacobian}
%# is a switch variable @var{['yes','no']} that tells the program how to set the @var{'Jacobian'} option within optimset function,
%# @item  @code{NewtonTol}
%# tolerance for termination test on residual,
%# @item  @code{MaxNewtonIterations}
%# maximum number of nonlinear iterations,
%# @item  @code{Eta}
%# initial forcing term (must be in the interval [0,1)). For details see [1],
%# @item  @code{Choice}
%# formula to use to select the forcing term may be 1 or 2, default value is 1. For details see [1],
%# @item  @code{Algorithm}
%# iterative method to solve the linearized system (default is @code{'GMRES'}),
%# @item  @code{Restart}
%# restart parameter for the GMRES solver (ignored for other solvers).
%# @end table
%# 
%# References:
%# [1]  S.C. Eisenstat and H.F. Walker, "Choosing the Forcing Terms in an Inexact Newton Method." SIAM Journal on Scientific Computing, 17(1), pp. 16-16, 1996.
%# @end deftypefn
%#
%# @seealso{odepkg}

function [t_next, x_next, x_est, k] = bwe_inexact_newton (f, t, x, dt,
                                                          options = [],
                                                          k_vals = [],
                                                          t_next = t + dt)

  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end

  if (strcmp(options.UseJacobian, 'yes'))
    H = options.Jacobian;
    opts = optimset ('Jacobian', 'on', 'TolFun', options.NewtonTol, 'MaxIter', options.MaxNewtonIterations);
  else
    H = [];
    opts = optimset ('Jacobian', 'off', 'TolFun', options.NewtonTol, 'MaxIter', options.MaxNewtonIterations);
  end

  dim = length (x);
  havejachandle = is_function_handle (options.Jacobian);
  BWE = @(y)bwe_system (y, dim, f, H, t, x, dt, havejachandle, args);

  if (strcmp (options.InexactSolver, 'inexact_newton'))

    opts.Eta       = options.Eta;
    opts.Choice    = options.Choice;
    opts.Algorithm = options.Algorithm;
    opts.Restart   = options.Restart;

    x_next = inexact_newton (BWE, x, opts);
  else

    x_next = fsolve (BWE, x, opts);

  end
  
  #t_next = t + dt;

  if (nargout >= 3)
    %# second order approximation with CRANK-NICOLSON
    x_est = x + (dt/2)*(feval (f, t, x, args{:}) + ...
        feval (f, t+dt, x_next, args{:}));
    k = []; % no runge-kutta evaluations
  end

end


function [F,J] = bwe_system (y, dim, f, H, t, x0, dt, jac_handle, vargs)

  F = y - x0 - dt*feval (f, t+dt, y, vargs{:});

  if (nargout == 2)
    if (jac_handle)
      J = eye (dim) - dt*feval (H, t+dt, y, vargs{:});
    else
      J = eye (dim) - dt*H;
    end
  end

end

%! # We are using the "Van der Pol" implementation.
%!function [ydot] = fpol (vt, vy) %# The Van der Pol
%!  ydot = [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!end
%!function [vjac] = fjac (vt, vy, varargin) %# its Jacobian
%!  vjac = [0, 1; -1 - 2 * vy(1) * vy(2), 1 - vy(1)^2];
%!end
%!
%!shared opts
%!  opts = odeset;
%!  opts.funarguments = {};
%!  opts.UseJacobian = 'no';
%!  opts.havejachandle = true;
%!  opts.Jacobian = @fjac;
%!  opts.InexactSolver = 'inexact_newton';
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  [t,y] = bwe_inexact_newton (@fpol, 0, [2;0], 0.05, opts);
%!test
%!  opts.UseJacobian = 'yes';
%!  opts.InexactSolver = 'fsolve';
%!  [t,y,x] = bwe_inexact_newton (@fpol, 0, [2;0], 0.05, opts);
%!test
%!  opts.UseJacobian = 'no';
%!  opts.InexactSolver = 'inexact_newton';
%!  opts.NewtonTol = 1.e-10; opts.MaxNewtonIterations = 1000;
%!  [t,y,x] = bwe_inexact_newton (@fpol, 0, [2;0], 0.05, opts);
%!test
%!  opts.NewtonTol = 1.e-10; opts.MaxNewtonIterations = 1000;
%!  opts.Choice = 2;
%!  [t,y,x] = bwe_inexact_newton (@fpol, 0, [2;0], 0.05, opts);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
