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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} rattle (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} rattle (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%#
%# This function can be used to integrate a constrained Hamiltonian system with
%# a given initial condition @var{x} from @var{t} to @var{t+dt}, with the RATTLE
%# method using the Richardson extrapolation method to estimate the error. For
%# details about this method see [1].
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
%# adapt the computation to what is needed. Possible fields are
%# @table @option
%# @item  @code{HamiltonianHessFcn}
%# a function_handle that represents the Hessian of the Hamiltonian of the system,
%# @item  @code{ConstraintFcn}
%# a function_handle that represents the constraints of the system,
%# @item  @code{ConstraintsNb}
%# a positive integer equal to the number of constraints of the system,
%# @item  @code{ConstraintGradFcn}
%# a function_handle that represents the gradient of the constraints of the system,
%# @item  @code{ConstraintHessFcn}
%# a function_handle that represents the Hessian of the constraints of the system.
%# @end table
%#
%# References:
%# [1]  E. Hairer, C. Lubich and G. Wanner, "Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations." Second Edition, Springer.
%#
%# @seealso{odepkg}
%# @end deftypefn

function [t_next, x_next, x_est, k] = rattle (f, t, x, dt, options = [],
                                                      k_vals = [],
                                                      t_next = t + dt)

  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end

  %# getting information about the Hessian of the Hamiltonian and the Hessian
  %# of the constraint function
  H  = options.HamiltonianHessFcn;
  GG = options.ConstraintHessFcn;

  %# if they are not empty it means that it is possible to use the explicit
  %# expression of te Jacobian of the nonlinear system
  if ( ~isempty(H) && ~isempty(GG) )
    opts = optimset ('Jacobian', 'on', 'TolFun', options.NewtonTol,
                     'MaxIter', options.MaxNewtonIterations);
  else
    opts = optimset ('Jacobian', 'off', 'TolFun', options.NewtonTol,
                     'MaxIter', options.MaxNewtonIterations);
  end

  if (strcmp(options.InexactSolver, 'inexact_newton'))
    opts.Eta       = options.Eta;
    opts.Choice    = options.Choice;
    opts.Algorithm = options.Algorithm;
    opts.Restart   = options.Restart;
  end

  %# getting information about the constraints function and its gradient
  g = options.ConstraintFcn;
  G = options.ConstraintGradFcn;

  %# getting the number of constraints
  constr_nb = options.ConstraintsNb;

  %# setting initial conditions
  dim = length (x) / 2;
  q0 = x(1:dim);
  p0 = x(dim+1:end);

  %# defining the nonlinear system to be solved
  RATTLE = @(y) constrained_system (y, dim, constr_nb, f, H, g, G, GG,
                                    t, q0, p0, dt, args);

  %# setting initial conditions for the nonlinear system
  y0 = [q0;p0;p0;zeros(2*constr_nb, 1)];

  %# solving the system
  y0 = options.solver (RATTLE, y0, opts);

  %# computing new time and new values for the unkwnowns
  t_next = t + dt;
  x_next = [y0(1:dim);y0(2*dim+1:3*dim)];

  %# if the estimation of the error is required
  if (nargout >= 3)
    %# new solution computed with 2 substeps
    dt_new = dt/2;

    %# defining the nonlinear system to be solved
    RATTLE = @(y) constrained_system (y, dim, constr_nb, f, H, g, G, GG,
                                      t, q0, p0, dt_new, args);

    %# setting initial conditions
    y0 = [q0;p0;p0;zeros(2*constr_nb, 1)];

    %# solving the system from t to t+(dt/2)
    y0 = options.solver (RATTLE, y0, opts);

    %# updating partial solutions
    q0 = y0(1:dim);
    p0 = y0(2*dim+1:3*dim);
    t = t + dt_new;

    %# defining the nonlinear system to be solved
    RATTLE = @(y) constrained_system (y, dim, constr_nb, f, H, g, G, GG,
                                      t, q0, p0, dt_new, args);

    %# setting initial conditions
    y0 = [q0;p0;p0;zeros(2*constr_nb, 1)];

    %# solving the system from t+(dt/2) to t+dt
    y0 = options.solver (RATTLE, y0, opts);

    %# new solution to be compared with the previous one
    x_est = [y0(1:dim);y0(2*dim+1:3*dim)];
    k = []; % no runge-kutta values
  end

end


%# nonlinear system to be solved. For its definition I refer to [1]
function [F, J] = constrained_system (y, dim, c_nb, f ,H, g, G, GG, t, q0, p0, dt, vargs)

  F = zeros (3*dim+2*c_nb, 1);

  F(1:dim) = y(1:dim) - q0 - (dt/2).*(feval (f, t, [q0;y(dim+1:2*dim)], vargs{:})(1:dim) + feval (f, t+dt, y(1:2*dim), vargs{:})(1:dim));
  F(dim+1:2*dim) = y(dim+1:2*dim) - p0 - (dt/2).*(feval (f, t, [q0;y(dim+1:2*dim)], vargs{:})(dim+1:end) - G(q0)'*y(3*dim+1:3*dim+c_nb));
  F(2*dim+1:3*dim) = y(2*dim+1:3*dim) - y(dim+1:2*dim) - (dt/2)*(feval (f, t+dt, y(1:2*dim), vargs{:})(dim+1:end) - G(y(1:dim))'*y(3*dim+c_nb+1:end));
  F(3*dim+1:3*dim+c_nb) = g(y(1:dim));
  F(3*dim+c_nb+1:end) = G(y(1:dim))*(feval (f, t+dt, [y(1:dim);y(2*dim+1:3*dim)], vargs{:})(1:dim));

  %# explicit expression for the Jacobian of the nonlinear system
  if (nargout == 2)

    J = zeros (3*dim+2*c_nb, 3*dim+2*c_nb);

    J(1:dim,1:dim) = eye(dim) - (dt/2)*(feval (H, t+dt, y(1:2*dim), vargs{:})(dim+1:end,1:dim));
    J(1:dim,dim+1:2*dim) = -(dt/2)*(feval (H, t, [q0;y(dim+1:2*dim)], vargs{:})(dim+1:end,dim+1:end) + feval (H, t+dt, y(1:2*dim), vargs{:})(dim+1:end,dim+1:end));

    J(dim+1:2*dim,dim+1:2*dim) = eye (dim) + (dt/2)*(feval (H, t, [q0;y(dim+1:2*dim)], vargs{:})(1:dim,dim+1:end));
    J(dim+1:2*dim,3*dim+1:3*dim+c_nb) = (dt/2)*G(q0)';

    J(2*dim+1:3*dim,1:dim) = (dt/2)*(feval (H, t+dt, y(1:2*dim), vargs{:})(1:dim,1:dim));
    for k = 1:1:c_nb
      J(2*dim+1:3*dim,1:dim) = J(2*dim+1:3*dim,1:dim) - (dt/2)*(y(3*dim+c_nb+k)*(GG(y(1:dim))(:,:,k)));
    end

    J(2*dim+1:3*dim,dim+1:2*dim) = -eye (dim) + (dt/2)*(feval (H, t+dt, y(1:2*dim), vargs{:})(1:dim,dim+1:end));
    J(2*dim+1:3*dim,2*dim+1:3*dim) = eye (dim) + (dt/2)*(feval (H, t+dt, y(1:2*dim), vargs{:})(1:dim,dim+1:end));
    J(2*dim+1:3*dim,3*dim+c_nb+1:end) = (dt/2)*G(y(1:dim))';

    J(3*dim+1:3*dim+c_nb,1:dim) = G(y(1:dim));

    J(3*dim+c_nb+1:end,1:dim) = G(y(1:dim))*(feval (H, t+dt, [y(1:dim);y(2*dim+1:3*dim)], vargs{:})(dim+1:end,1:dim));
    for k = 1:1:c_nb
      J(3*dim+c_nb+k,1:dim) = J(3*dim+c_nb+k,1:dim) + ((GG(y(1:dim))(:,:,k))*(feval (f, t+dt, [y(1:dim);y(2*dim+1:3*dim)], vargs{:})(1:dim)))';
    end

    J(3*dim+c_nb+1:end,2*dim+1:3*dim) = G(y(1:dim))*(feval (H, t+dt, [y(1:dim);y(2*dim+1:3*dim)], vargs{:})(dim+1:end,dim+1:end));

  end

end

%! Simple-pendulum
%!
%!function [ydot] = pendulum (vt, vy)
%!  ydot = [0.5*vy(length(vy)/2+1:end); -(2*5*9.806)*[zeros(length(vy)/4);ones(length(vy)/4)]];
%!end
%!function [g] = constraint (vy)
%!  g = (vy(:)')*vy(:) - 25.0;
%!end
%!function [G] = grad_constr (vy)
%!  G = 2.0*vy(:)';
%!end
%!
%!shared opts
%!  opts = odeset;
%!  opts.vfunarguments = {};
%!  opts.solver = @fsolve;
%!  opts.ConstraintFcn = @constraint;
%!  opts.ConstraintsNb = 1;
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  opts.ConstraintGradFcn = @grad_constr;
%!  [t,y] = rattle (@pendulum, 0, [5.0*cos(-pi/4);5.0*sin(-pi/4);0;0], 0.05, opts);
%!test
%!  opts.ConstraintGradFcn = @grad_constr;
%!  [t,y,x] = rattle (@pendulum, 0, [5.0*cos(-pi/4);5.0*sin(-pi/4);0;0], 0.1, opts);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
