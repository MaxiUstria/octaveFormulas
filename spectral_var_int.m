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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} spectral_var_int (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} spectral_var_int (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# 
%# This function can be used to integrate a Hamiltonian system with a given initial condition @var{x} from @var{t} to @var{t+dt}, with the spectral variational integrators method. The error is estimated comparing the solution with that one obtained with the same method but one order higher polynomials and one degree higher quadrature rule. For details about the theory see [1].
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
%# Fifth input parameter is optional and describes a set of options useful to adapt the computation to what is needed. The possible fields are
%# @table @option
%# @item  @code{HamiltonianHessFcn}
%# a function_handle that represents the Hessian of the Hamiltonian of the system,
%# @item  @code{Q_DoFs}
%# a positive integer equal to Gauss quadrature rule degrees of freedom,
%# @item  @code{P_DoFs}
%# a positive integer equal to Legendre polynomial degrees of freedom,
%# @item  @code{Nodes}
%# a vector containing Gauss quadrature nodes,
%# @item  @code{Weights}
%# a vector containing Gauss quadrature weights,
%# @item  @code{Legendre}
%# a matrix containing Legendre polynomials values at quadrature nodes,
%# @item  @code{Derivatives}
%# a matrix containing derivatives of Legendre polynomials evaluated at quadrature nodes,
%# @item  @code{Extremes}
%# a matrix containing Legendre polynomials values at extremes,
%# @item  @code{Q_DoFs_err}
%# a positive integer equal to Gauss quadrature rule degrees of freedom for the solution used to estimate the error,
%# @item  @code{P_DoFs_err}
%# a positive integer equal to Legendre polynomial degrees of freedom for the solution used to estimate the error,
%# @item  @code{Nodes_err}
%# a vector containing Gauss quadrature nodes for the solution used to estimate the error,
%# @item  @code{Weights_err}
%# a vector containing Gauss quadrature weights for the solution used to estimate the error,
%# @item  @code{Legendre_err}
%# a matrix containing Legendre polynomials values at quadrature nodes for the solution used to estimate the error,
%# @item  @code{Derivatives_err}
%# a matrix containing derivatives of Legendre polynomials evaluated at quadrature nodes for the solution used to estimate the error,
%# @item  @code{Extremes_err}
%# a matrix containing Legendre polynomials values at extremes for the solution used to estimate the error.
%# @end table
%#
%# References:
%# [1]  J.E. Marsden and M. West, "Discrete Mechanics and Variational Integrators." Acta Numerica (2001), pp 1-158, Cambridge University Press.
%#
%# @seealso{odepkg}
%# @end deftypefn

function [t_next, x_next, x_est, k] = spectral_var_int (f, t, x, dt,
                                                        options = [],
                                                        k_vals = [],
                                                        t_next = t + dt)

  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end

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

  %# getting information about Gauss quadrature rule and Legendre polynomials
  q_dofs         = options.Q_DoFs;
  p_dofs         = options.P_DoFs;
  nodes          = options.Nodes;
  weights        = options.Weights;
  Legendre       = options.Legendre;
  derivatives    = options.Derivatives;
  extreme_values = options.Extremes;

  %# setting parameters and initial conditions
  dim = length (x) / 2;
  N   = q_dofs+p_dofs+2;
  q0  = x(1:dim)';
  p0  = x(dim+1:end)';

  %# defining the nonlinear system to be solved
  SPVI = @(y) svi_system (y, q_dofs, p_dofs, weights, Legendre, derivatives, ...
      extreme_values, dim, N, f, H, t, q0, p0, dt, args);

  %# setting initial conditions for the nonlinear system
  y0 = [q0;zeros(q_dofs-1,dim);ones(p_dofs,1)*p0;q0;p0];
  y0 = reshape (y0, dim*N, 1);

  %# solving the system
  z = options.solver (SPVI, y0, opts);
  z = reshape (z, N, dim);

  %# temporary variable
  y = [Legendre*z(1:q_dofs,1:dim),z((q_dofs+1):(end-2),1:dim)];

  %# computing new times and new values for the unkwnowns
  x_next = [y;z(end-1,:),z(end,:)]'; %x_next

  comp = odeget (options, 'comp', 0.0, 'fast_not_empty'); #FIXME: this was never introduced before...
  local_comp = comp*ones(1,p_dofs);
  t_next = t*ones(1,p_dofs);
  for i=1:1:p_dofs
    [t_next(i), local_comp(i)] = kahan (t_next(i), local_comp(i), (dt/2)*(nodes(i)+1));
  end
  t_next = [t_next,t+dt]';

  %# if the estimation of the error is required
  if (nargout >= 3)
    %# new solution computed with one higher order polynomials and one higher degree of the quadrature rule
    q_dofs_err         = options.Q_DoFs_err;
    p_dofs_err         = options.P_DoFs_err;
    nodes_err          = options.Nodes_err;
    weights_err        = options.Weights_err;
    Legendre_err       = options.Legendre_err;
    derivatives_err    = options.Derivatives_err;
    extreme_values_err = options.Extremes_err;

    %# setting parameters and initial conditions
    N_err = q_dofs_err+p_dofs_err+2;
    q1 = x_next(1:dim,end)';
    p1 = x_next(dim+1:end,end)';

    %# defining the nonlinear system to be solved
    SPVI = @(y) svi_system (y, q_dofs_err, p_dofs_err, weights_err, Legendre_err, ...
        derivatives_err, extreme_values_err, dim, N_err, f, H, t, q0 ,p0 , dt, args);

    %# initializing initial conditions
    p_interp = zeros (p_dofs_err, dim);

    %# setting initial nodal values for momenta: they are just a simple mean of those of lower order, because doing more sophisticated initializations doesn't change solver's iterations number
    p_lower_order = [p0;z(q_dofs+1:q_dofs+p_dofs,:);p1];
    for i=1:1:p_dofs_err
      p_interp(i,:) = .5*(p_lower_order(i,:) + p_lower_order(i+1,:));
    end

    %# setting initial conditions
    y0 = [z(1:q_dofs,:);zeros(1,dim);p_interp;q1;p1];
    y0 = reshape (y0, dim*N_err, 1);

    %# solving the system
    z = options.solver (SPVI, y0, opts);
    z = reshape (z, N_err, dim);

    %# new solution to be compared with the previous one
    x_est = [z(end-1,:),z(end,:)]';
    k = []; % no runge-kutta values

  end

t_next = t_next(end);
x_next = x_next(:,end);

end

%# nonlinear system to be solved. For its definition I refer to [1]
function [F, Jac] = svi_system (y, q_dofs, p_dofs, w, L, D, C, dim, N, f, H, t, q0, p0, dt, vargs)

  F = zeros (N*dim, 1);
  V = zeros (p_dofs, dim*2);
  X = zeros (dim*2, 1);
  W = reshape (y, N, dim);
  
  for i = 1:1:p_dofs
    X = [L(i,:)*W(1:q_dofs,:),W(i+q_dofs,:)]';
    V(i,:) = feval (f, t, X, vargs{:});
  end

  for i = 1:1:dim
    F((1+N*(i-1)):(q_dofs+N*(i-1)),1) = (ones (1, p_dofs)*(((w.*y((q_dofs+1+N*(i-1)):(q_dofs+p_dofs+N*(i-1))))*ones (1, q_dofs)).*D + ...
      (((dt/2).*w.*V(:,i+dim))*ones (1, q_dofs)).*L) + (p0(i)*ones (1, q_dofs)).*C(1,:) - ...
      (y(N*i)*ones (1, q_dofs)).*C(2,:))';
    F(1+N*(i-1)+q_dofs:N*(i-1)+q_dofs+p_dofs,1) = V(:,i) - (2.0/dt)*(D*y((1+N*(i-1)):(q_dofs+N*(i-1))));
    F(N*i-1) = C(2,:)*y((1+N*(i-1)):(q_dofs+N*(i-1))) - y(N*i-1);
    F(N*i) = C(1,:)*y((1+N*(i-1)):(q_dofs+N*(i-1))) - q0(i);
  end

  %# explicit expression for the Jacobian of the nonlinear system
  if (nargout == 2)
    warning ('off', 'Octave:broadcast');

    DV = zeros ((dim*2)^2, p_dofs);

    for i = 1:1:p_dofs
      X = [L(i,:)*W(1:q_dofs,:),W(i+q_dofs,:)]';
      DV(:,i) = transform (feval (H, t, X, vargs{:}));
    end

    DV  = DV';
    Jac = zeros (N*dim, N*dim);

    for u = 1:1:dim
      for i = 1:1:q_dofs
        Jac(i+N*(u-1),1+N*(u-1):q_dofs+N*(u-1)) = -(dt/2).*(ones (1, p_dofs)*((w.*L(:,i).*DV(:,u+(u-1)*(2*dim))).*L));
      end
      Jac(1+N*(u-1):q_dofs+N*(u-1),1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1)) = (w.*(D - (dt/2).*DV(:,u+(u+dim-1)*(2*dim)).*L))';
      Jac(1+N*(u-1):q_dofs+N*(u-1),N*u) = -C(2,:)';
      Jac(1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1),1+N*(u-1):q_dofs+N*(u-1)) = DV(:,u+dim+(u-1)*(2*dim)).*L - (2.0/dt).*D;
      Jac(1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1),1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1)) = diag(DV(:,u+dim+(u+dim-1)*(2*dim)));  
      Jac(N*u-1,1+N*(u-1):q_dofs+N*(u-1)) = C(2,:);
      Jac(N*u,1+N*(u-1):q_dofs+N*(u-1)) = C(1,:);
      Jac(N*u-1,N*u-1) = -1;
      for v=u+1:1:dim
        for i=1:1:q_dofs
          Jac(i+N*(u-1),1+N*(v-1):q_dofs+N*(v-1)) = -(dt/2).*(ones (1, p_dofs)*((w.*L(:,i).*DV(:,u+(v-1)*(2*dim))).*L));
        end
        Jac(1+N*(v-1):q_dofs+N*(v-1),1+N*(u-1):q_dofs+N*(u-1)) = Jac(1+N*(u-1):q_dofs+N*(u-1),1+N*(v-1):q_dofs+N*(v-1));
        Jac(1+N*(u-1):q_dofs+N*(u-1),1+q_dofs+N*(v-1):p_dofs+q_dofs+N*(v-1)) = -(dt/2).*((w.*DV(:,u+(v+dim-1)*(2*dim))).*L)';
        Jac(1+N*(v-1):q_dofs+N*(v-1),1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1)) = Jac(1+N*(u-1):q_dofs+N*(u-1),1+q_dofs+N*(v-1):p_dofs+q_dofs+N*(v-1));
        Jac(1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1),1+N*(v-1):q_dofs+N*(v-1)) = DV(:,u+dim+(v-1)*(2*dim)).*L;
        Jac(1+q_dofs+N*(v-1):p_dofs+q_dofs+N*(v-1),1+N*(u-1):q_dofs+N*(u-1)) = Jac(1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1),1+N*(v-1):q_dofs+N*(v-1));
        Jac(1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1),1+q_dofs+N*(v-1):p_dofs+q_dofs+N*(v-1)) = diag(DV(:,u+dim+(v+dim-1)*(2*dim)));
        Jac(1+q_dofs+N*(v-1):p_dofs+q_dofs+N*(v-1),1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1)) = Jac(1+q_dofs+N*(u-1):p_dofs+q_dofs+N*(u-1),1+q_dofs+N*(v-1):p_dofs+q_dofs+N*(v-1));   
      end
    end
  end

end

%# function useful to vectorize a matrix
function v = transform (M)
  [r,c] = size (M);
  v = zeros (r*c, 1);
  for i = 1:1:r
    for j = 1:1:c
      v(i+(j-1)*r) = M(i,j);
    end
  end
end

%! Armonic-oscillator
%!function [ydot] = hamilt (vt, vy)
%!  ydot = [vy(length(vy)/2+1:end); -vy(1:length(vy)/2)];
%!end
%!
%!test
%!  vopt = odeset;
%!  vopt.vfunarguments = {};
%!  vopt.solver = @fsolve;
%!  [t,y] = spectral_var_int (@hamilt, 0, [1;0], 0.05, vopt);
%!test
%!  vopt = odeset;
%!  vopt.vfunarguments = {};
%!  vopt.solver = @fsolve;
%!  [t,y,x] = spectral_var_int (@hamilt, 0, [1;1], 0.1, vopt);

%# Local Variables: ***
%# mode: octave ***
%# End: ***
