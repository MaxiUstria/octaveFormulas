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
%# @deftypefn {Command} {[@var{x}, @var{fval}, @var{exitflag}, @var{output}, @var{jacobian}] =} inexact_newton (@var{@@fun}, @var{x0}, [@var{options}])
%# 
%# Solve the nonlinear equation @code{f(x) = 0} using the inexact Newton method, starting from the initial guess @var{x0}. 
%# 
%# If an exact jacobian is not provided it is approximated by finite differences, if the method for solving linear systems is not specified it is set as defaults to @var{'gmres'}, other possible values are @var{'pcg'} or @var{'bicgstab'}. 
%# 
%# First output parameter @var{x} is the computed solution.
%# 
%# Second output parameter is @var{@@fval} that is the value of the function at computed solution.
%# 
%# Third output parameter is exitflag and it has not yet been implemented.
%# 
%# Fourth output parameter is output which is a struct containing the field @var{'iterations'} that is the number if iterations required to converge and @var{'funcCount'} which contains the number of function evaluations computed.
%# 
%# Last output argument is jacobian and is the value of the Jacobian at the computed solution.
%# 
%# The first input parameter @var{@@fun} must be a function_handle.
%# 
%# The second input argument @var{x0} must be the initial guess.
%# 
%# The third input argument is optional and describes a set of options useful to adapt the computation to what is needed. Possible fields are
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
%# @seealso{gmres, pcg, bicgstab, odebwe}

function [x, fval, exitflag, output, jacobian] = inexact_newton (fun, x0, options) 
  
  %# check if first input argument is a valid one
  if ~(isa(fun,'function_handle') || isa(fun,'inline') )
    error('OdePkg:InvalidArgument', ...
      'function %s must be a function_handle or inline function',func2str(fun));
  end

  %# checking if the user has requested to use the Jacobian
  set_jacob = 'off';
  jcb = optimget(options,'Jacobian');
  if( strcmp(jcb,'on') )
    set_jacob = 'on';
  end

  %# when the Jacobian is requested, check if the given function returns a second output argument
  have_jacobian = false;
  flag = 0;
  if( strcmp(set_jacob,'on') )
    try
      [~,jcb] = fun(x0);
    catch
      flag = 1;
      warning('OdePkg:InvalidArgument', ...
        'request of using Jacobian cannot be satisfied');
    end
    if(flag==0)
      have_jacobian = true;
    end
  end

  %# get solving options
  %# Newton Tolerance
  tol = options.TolFun;
  %# Maximum number of nonlinear iterations
  maxit = options.MaxIter;

  %# Eta parameter
  eta = options.Eta;

  %# choice on the error estimation
  choice = options.Choice; % default, see [1] for details

  %# getting the iterative method for linearized system
  iterative_method = options.Algorithm; % default

  switch iterative_method
    case 'gmres'
      %# restart parameter for gmres solver
      restart = options.Restart;

      solver = @(A,b,rtol) gmres (A,b,restart,rtol);
    case 'pcg'
      solver = @(A,b,rtol) pcg (A,b,rtol);
    case 'bicgstab'
      solver = @(A,b,rtol) bicgstab (A,b,rtol);
    otherwise
      error('OdePkg:InvalidArgument', ...
        'invalid iterative solver %s',iterative_method);
  end

  backtr_max = 10;
  %# initializing iterations and function evaluations
  iter = 0;
  fevals = 0;

  x = x0;
  f_x = fun(x);

  fevals++;
  res = norm(f_x,2);
  
  %# compute delta with the empirical formula given in [2]
  ord = 1; %# 1st order

  delta = eps^(1/(ord+1)); 

  while ((res > tol) && (iter++ < maxit))

    %# approximating the Jacobian if this is not given
    if (!have_jacobian)     
      A = @(v) (1/delta)*(fun(x + delta*v) - f_x);
    else
      [~,A] = fun(x); %# the Jacobian is given as second output argument
      fevals++;
    end   

    x_old = x;
    f_x_old = f_x;

    %# solve the linearized system
    [step,~,~,inner_iter] = solver(A,-f_x,eta); % use the chosen solver

    if( !have_jacobian )
      fevals = fevals + 2*inner_iter;
    end

    %# backtracking: reduce step size if nonlinear residual is bigger
    %# than the residual at the end of the past newton iteration 
    backtr_iter = 0;
    
    if 0
      factor = 1/2; %# this is an arbitrary choice
      while ((norm(fun(x+step),2) > res) && (backtr_iter++ < backtr_max))
        fevals++;
        step *= factor;
        eta = eta/2 + 1/2; %# see [1]
      end
      if (backtr_iter == backtr_max)
        error('OdePkg:InvalidArgument', ...
          'residual cannot be reduced')
      end
    end

    %# update x
    x = x + step;
    f_x = fun(x);
    fevals++;
    
    resold = res;
    res = norm(f_x,2);
    
    %# compute the forcing term
    if (choice == 1) % see [1]
      
      %# eta = norm(F(x) - F(x_old) - F'(x_old)*step_old)/norm(F(x_old)) 
      %# F'(x_old)*step_old   can be approximated with first order
      %# finite differences formula. see [2] for details 
      eta_old = eta^((1+sqrt(5))/2);
      
      if (have_jacobian)
        eta = norm(f_x-f_x_old-A*step,2)/resold;
      else
        eta = norm(f_x-f_x_old-A(step),2)/resold;
        fevals = fevals + 2;
      end
      
    else %# choice 2

      alpha = (1+sqrt(5))/2; %# see [1]
      gamma = 1;
      eta_old = gamma * eta ^ alpha;
      eta = gamma * (res / resold) ^ alpha;
      
    end
      
    if (eta_old > .1) %# safeguard
      eta = max(eta,eta_old);
    end
    eta = min(eta,.99); %# eta must belong to interval [0,1) 	
    
  end

  fval = f_x;

  exitflag = [];

  if (!have_jacobian)
    jacobian = A(x);
    fevals = fevals + 2;
  else
    [~,jacobian] = fun(x);
    fevals++;
  end

  output.iterations = iter;
  output.funcCount = fevals;
  
end

%! %# Turn off output of warning messages for all tests, turn them on
%! %# again if the last test is called
%!test
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  f = @(x) [((x(1)+1)^2 + 3*x(2)^7 + 7*(x(1)+1)); sum(x)+1];
%!  opts.NewtonTol = 1.e-10; opts.MaxNewtonIterations = 100;
%!  sol = inexact_newton (f, [.1;.1], opts);
%!  assert (sol, [-1; 0], 1e-9)
%!  assert (f (sol), [0; 0], 1e-10)
%!test
%!function [f,jac] = func(x)
%!  f = [((x(1)+1)^2 + 3*x(2)^7 + 7*(x(1)+1)); sum(x)+1];
%!  jac = [2*(x(1)+1)+7, 3*7*x(2)^6; 1, 1];
%!end
%!  opts.NewtonTol = 1.e-10; opts.MaxNewtonIterations = 100;
%!  opts.Jacobian = 'on';
%!  [sol, it] = inexact_newton (@func, [.1;.1], opts);
%!  assert (sol, [-1; 0], 1e-9)
%!  assert (f (sol), [0; 0], 1e-10)
%!test 
%!function [f,jac] = func(x)
%!  f = [((x(1)+1)^2 + 3*x(2)^7 + 7*(x(1)+1)); sum(x)+1];
%!  jac = [2*(x(1)+1)+7, 3*7*x(2)^6; 1, 1];
%!end
%!
%!  opts.NewtonTol = 1.e-10; opts.MaxNewtonIterations = 100;
%!  opts.Jacobian = 'on'; opts.Eta = 0.1; opts.Choice = 2;
%!  opts.Algorithm = 'bicgstab';
%!  [sol, it] = inexact_newton (@func, [.1;.1], opts);
%!  assert (sol, [-1; 0], 1e-9)
%!  assert (func (sol), [0; 0], 1e-10)
%!
%!  warning ('on', 'OdePkg:InvalidArgument');

%# Local Variables: ***
%# mode: octave ***
%# End: ***
