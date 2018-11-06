%# Copyright (C) 2009-2012, Sebastian Schoeps <schoeps@math.uni-wuppertal.de>
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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} bwe_newton_raphson (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} bwe_newton_raphson (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# 
%# This function can be used to integrate a set of non--stiff ordinary differential equations (non--stiff ODEs) with a given initial condition @var{x} from @var{t} to @var{t+dt}, with the Backward-Euler method using the Newton-Raphson method to solve the nonlinear system and estimating the error comparing the solution to that one given in two substeps by the same method.
%# 
%# First output argument is the final integration time value. 
%# 
%# Second output parameter is the computed solution at time @var{t_next}.
%# 
%# Third output parameter is a higher order solution for the estimation of the error.
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
%# @item  @code{NewtonTol}
%# tolerance for termination test on residual,
%# @item  @code{MaxNewtonIterations}
%# maximum number of nonlinear iterations,
%# @item  @code{Mass}
%# a function_handle or a constant matrix defining the Mass matrix,
%# @item  @code{Jacobian}
%# a function_handle defining the Jacobian,
%# @item  @code{havemasshandle}
%# whether or not Mass is a function_handle,
%# @item  @code{massdependence}
%# whether or not Mass depends on the solution.
%# @end table
%# @end deftypefn
%#
%# @seealso{odepkg}

function [t_next, x_next, x_est, k] = bwe_newton_raphson (f, t, x, dt,
                                                          options = [],
                                                          k_vals = [],
                                                          t_next = t + dt)

  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end

  jacfun = false;
  jacmat = false;
  if (! isempty (options.Jacobian))
    if (ischar (options.Jacobian))
      jacfun = true;
      jac = str2fun (options.Jacobian);
    elseif (is_function_handle (options.Jacobian))
      jacfun = true;
      jac = options.Jacobian;
    elseif (ismatrix (options.Jacobian))
      jacmat = true;
      jac = options.Jacobian;
    else
      error ("ode23s: the jacobian should be passed as a matrix, a string or a function handle")
    end
  end

  jacpat = false;
  if (! isempty (options.JPattern))
    jacpat = true;
    [ii, jj] = options.Jacobian;
    pattern = sparse (ii, jj, true);
  end

  %# Jacobian matrix, dfxpdp
  if (jacmat)
    J = jac;
  elseif (jacfun)
    J = jac (t, x);
  elseif (! jacpat)
    J = __dfxpdp__ (x, @(a) feval (f, t, a, args{:}), options.RelTol);
  elseif (jacpat)
    J = __dfxpdp__ (x, @(a) feval (f, t, a, args{:}), options.RelTol, pattern);
  end

  y = zeros (size (x, 1), (nargout-1));

  for j=1:1:(nargout-1)
    %# Initial value (result of the previous timestep)
    y0 = x;

    %# Initial guess for Newton-Raphson
    y(:,j) = x;

    %# We do not use a higher order approximation for the
    %# comparsion, but two steps by the Backward Euler
    %# method
    for i=1:1:j
      %# initialize the time stepping parameters
      step = dt/j;
      time = t + i*step;
      newton_iteration = 1;
      residual_norm = inf (1, options.MaxNewtonIterations);

      %# Start the Newton iteration
      while ( (newton_iteration < options.MaxNewtonIterations) && (residual_norm(newton_iteration) > options.NewtonTol) )

        %# Compute the Jacobian of the non-linear equation,
        %# that is the matrix pencil of the mass matrix and
        %# the right-hand-side's Jacobian. Perform a (sparse)
        %# LU-Decomposition afterwards.
        if ( (newton_iteration==1) || (~options.issimplified) )
          %# Get the mass matrix from the left-hand-side
          if (isa (options.Mass, "function_handle"))      %# Handle only the dynamic mass matrix,
            if (! strcmp (options.MStateDependence, "none"))
              mass = options.Mass (t, x, args{:});
            else      %# if (vmassdependence == false)
              mass = options.Mass (t, options.funarguments{:});
            end
          else
            mass = options.Mass;
          end

          full_jacobian = mass/step - J;
          %# one could do a matrix decomposition of vfulljac here, 
          %# but the choice of decomposition depends on the problem
          %# and therefore we use the backslash-operator in row 105 
        end

        %# Compute the residual
        residual = mass/step * (y(:,j)-y0) - feval (f, time, y(:,j), args{:});
        residual_norm(newton_iteration+1) = norm (residual, inf);

        %# Solve the linear system
        y(:,j) = full_jacobian\(-residual + full_jacobian*y(:,j));

        %# Prepare next iteration
        newton_iteration = newton_iteration+1;
      end   %# while Newton

      %# Leave inner loop if Newton diverged
      if (residual_norm(newton_iteration) > options.NewtonTol)
        break;
      end

      %# Save intermediate solution as initial value
      %# for the next intermediate step
      y0 = y(:,j);
    end   %# for steps

    %# Leave outer loop if Newton diverged
    if (residual_norm(newton_iteration) > options.NewtonTol)
      break;
    end

  end %# for estimators

  % if all Newton iterations converged
  if (residual_norm(newton_iteration) < options.NewtonTol)
    %# First order approximation using step size h
    y1 = y(:,1);
    %# If adaptive: first order approximation using step
    %# size h/2, if fixed: y1=y2=y3
    y2 = y(:,(nargout-1));
    %# Second order approximation by ("Richardson")
    %# extrapolation using h and h/2
    y3 = y2 + (y2-y1);
  end

  %# defining new time and new values for the unkwnowns
  #t_next = t + dt;
  x_next = y2;

  %# computing the estimation of the error
  if (residual_norm(newton_iteration) <= options.NewtonTol && (nargout>=3))
    x_est = y1;
    k = []; % no runge-kutta evaluations
  end

end

%! # We are using the "Van der Pol" implementation.
%!function [ydot] = fpol (vt, vy) %# The Van der Pol
%!  ydot = [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!end
%!function [vjac] = fjac (vt, vy) %# its Jacobian
%!  vjac = [0, 1; -1 - 2 * vy(1) * vy(2), 1 - vy(1)^2];
%!end
%!
%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  opts.issimplified = false; 
%!  opts.MaxNewtonIterations = 100; 
%!  opts.NewtonTol = 1.e-6; 
%!  [t,y] = bwe_newton_raphson (@fpol, 0, [2;0], 0.05, opts);
%!test
%!  opts.issimplified = true; 
%!  opts.MaxNewtonIterations = 1000; 
%!  opts.NewtonTol = 1.e-10; 
%!  [t,y,x] = bwe_newton_raphson (@fpol, 0, [2;0], 0.05, opts);
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
