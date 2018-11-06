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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} symplectic_euler (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} symplectic_euler (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%#
%# This function can be used to integrate a Hamiltonian system with a given
%# initial condition @var{x} from @var{t} to @var{t+dt}, with the Symplectic
%# Euler method (also known as Semi--implicit Euler) using the Richardson
%# extrapolation method to estimate the error. For details about this method see
%# E. Hairer and C. Lubich and G. Wanner, "Geometric Numerical Integration:
%# Structure-Preserving Algorithms for Ordinary Differential Equations",
%# Springer.
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
%# adapt the computation to what is needed. The possible fields are
%# @table @option
%# @item  @code{Explicit}
%# a switch @var{['yes','no']} which tells the solver if the system to be solved
%# is explicit or not,
%# @item  @code{HamiltonianHessFcn}
%# a function_handle that represents the Hessian of the Hamiltonian of the
%# system.
%# @end table
%#
%# References: [1] E. Hairer, C. Lubich and G. Wanner, "Geometric Numerical
%# Integration: Structure-Preserving Algorithms for Ordinary Differential
%# Equations." Second Edition, Springer.
%#
%# @seealso{odepkg}
%# @end deftypefn

function [t_next, x_next, x_est, k] = symplectic_euler (f, t, x, dt,
                                                        options = [],
                                                        k_vals = [],
                                                        t_next = t + dt)

  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end
  
  %# setting initial conditions
  dim = length (x) / 2;
  q0 = x(1:dim);
  p0 = x(dim+1:end);

  %# if the system is explicit the computation will be faster
  if (strcmp (options.Explicit, "yes"))

    %# one-step scheme as defined in [1]
    p_next = p0 + dt*(feval (f, t, [q0; p0], args)(dim+1:end));
    q_next = q0 + dt*feval (f, t, [q0; p_next], args)(1:dim);

    %# defining new time and new values for the unkwnowns
    t_next = t+dt;
    x_next = [q_next; p_next];

    %# if the estimation of the error is required
    if (nargout >= 3)
      %# new solution computed with 2 substeps
      dt_new = dt/2;

      %# one-step scheme from t to t+(dt/2)
      p_next = p0 + dt_new*(feval (f, t, [q0; p0], args)(dim+1:end));
      q_next = q0 + dt_new*feval (f, t, [q0; p_next], args)(1:dim);

      %# updating partial solutions
      q0 = q_next;
      p0 = p_next;
      t  = t + dt_new;

      %# one-step scheme from t+(dt/2) to t+dt
      p_next = p0 + dt_new*(feval (f, t, [q0; p0], args)(dim+1:end));
      q_next = q0 + dt_new*feval (f, t, [q0; p_next], args)(1:dim);

      %# new solution to be compared with the previous one
      x_est = [q_next;p_next];
      k = []; % no runge-kutta values
      
    end

  else   %# if the system is implicit

    %# getting information about the Hessian of the Hamiltonian
    H = options.HamiltonianHessFcn;
    if (! isempty (H))
      opts = optimset ('Jacobian', 'on', 'TolFun', options.NewtonTol, 'MaxIter', options.MaxNewtonIterations);
    else
      opts = optimset ('Jacobian', 'off', 'TolFun', options.NewtonTol, 'MaxIter', options.MaxNewtonIterations);
    end

    if (strcmp(options.InexactSolver, 'inexact_newton'))
    opts.Eta       = options.Eta;
    opts.Choice    = options.Choice;
    opts.Algorithm = options.Algorithm;
    opts.Restart   = options.Restart;
    end

    %# defining the nonlinear system to be solved
    SE = @(y) se_system (y, dim, f, H, t, q0, p0, dt, args);

    %# setting initial conditions
    y0 = [q0;p0];

    %# solving the system
    y0 = options.solver (SE, y0, opts);

    %# computing new time and new values for the unkwnowns
    t_next = t + dt;
    x_next = y0;

    %# if the estimation of the error is required
    if (nargout >= 3)
      %# new solution computed with 2 substeps
      dt_new = dt/2;

      %# defining the nonlinear system from t to t+(dt/2)
      SE = @(y) se_system (y, dim, f, H, t, q0, p0, dt_new, args);

      %# setting initial conditions
      y0 = [q0;p0];

      %# solving the system
      y0 = options.solver (SE, y0, opts);

      %# updating partial solutions
      q0 = y0(1:dim);
      p0 = y0(dim+1:end);
      t  = t + dt_new;

      %# defining the nonlinear system from t to t+(dt/2)
      SE = @(y) se_system (y, dim, f, H, t, q0, p0, dt_new, args);

      %# setting initial conditions
      y0 = [q0;p0];

      %# solving the system
      y0 = options.solver (SE, y0, opts);

      %# new solution to be compared with the previous one
      x_est = y0;
      k = []; %no runge-kutta values

    end
  end

end


%# nonlinear system to be solved. For its definition I refer to [1]
function [F, J] = se_system (y, dim, f, H, t, q0, p0, dt, vargs)

  F = zeros (2*dim, 1);

  F(1:dim) = y(1:dim) - q0 - dt*feval (f, t, [q0;y(dim+1:end)], vargs{:})(1:dim);
  F(dim+1:end) = y(dim+1:end) - p0 - dt*feval (f, t, [q0;y(dim+1:end)], vargs{:})(dim+1:end);

  %# explicit expression for the Jacobian of the nonlinear system
  if (nargout == 2)
    J = zeros (2*dim, 2*dim);

    J(1:dim,1:dim) = eye (dim);
    J(1:dim,dim+1:end) = -dt*feval (H, t, [q0;y(dim+1:end)], vargs{:})(dim+1:end,dim+1:end);

    J(dim+1:end,dim+1:end) = eye (dim) + dt*feval (H, t, [q0;y(dim+1:end)], vargs{:})(1:dim,dim+1:end);

  end

end

%! Armonic-oscillator
%!function [ydot] = hamilt (vt, vy)
%!  ydot = [vy(length(vy)/2+1:end); -vy(1:length(vy)/2)];
%!end
%!
%!shared opts
%!  opts = odeset;
%!  opts.vfunarguments = {};
%!  opts.solver = @fsolve;
%!  opts.Explicit = 'no';
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  [t,y] = symplectic_euler (@hamilt, 0, [1;0], 0.05, opts);
%!test
%!  [t,y,x] = symplectic_euler (@hamilt, 0, [1;1], 0.1, opts);
%!test
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
