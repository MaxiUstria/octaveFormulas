%# Copyright (C) 2014-2016, Jacopo Corno <jacopo.corno@gmail.com>
%# Copyright (C) 2006-2012, Thomas Treichl <treichl@users.sourceforge.net>
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
%# @deftypefn  {Function File} {[@var{t}, @var{y}] =} odebwe (@var{fun}, @var{trange}, @var{init})
%# @deftypefnx {Function File} {[@var{t}, @var{y}] =} odebwe (@var{fun}, @var{trange}, @var{init}, @var{ode_opt})
%# @deftypefnx {Function File} {[@var{t}, @var{y}] =} odebwe (@dots{}, @var{par1}, @var{par2}, @dots{})
%# @deftypefnx {Function File} {[@var{t}, @var{y}, @var{te}, @var{ye}, @var{ie}] =} odebwe (@dots{})
%# @deftypefnx {Function File} {@var{solution} =} odebwe (@dots{})
%#
%# Solve a set of Ordinary Differential Equations with the well known implicit
%# Backward-Euler method of order 1. The method is indicated both for stiff and
%# non--stiff problems.
%#
%# @var{fun} is a function handle, inline function, or string containing the
%# name of the function that defines the ODE: @code{y' = f(t,y)}.  The function
%# must accept two inputs where the first is time @var{t} and the second is a
%# column vector of unknowns @var{y}.
%#
%# @var{trange} specifies the time interval over which the ODE will be
%# evaluated.  Typically, it is a two-element vector specifying the initial and
%# final times (@code{[tinit, tfinal]}).  If there are more than two elements
%# then the solution will also be evaluated at these intermediate time
%# instances using linear interpolation.
%#
%# @var{init} contains the initial value for the unknowns.  If it is a row
%# vector then the solution @var{y} will be a matrix in which each column is
%# the solution for the corresponding initial value in @var{init}.
%#
%# The optional fourth argument @var{ode_opt} specifies non-default options to
%# the ODE solver.  It is a structure generated by @code{odeset}.  @code{odebwe}
%# will ignore the following options: "BDF", "InitialSlope", "MassSingular",
%# "MaxOrder", "MvPattern", "NonNegative".  The user can set extra options by
%# adding fields to the structure.  @code{odebwe} allows for:
%#
%# @table @asis
%# @item  NewtonTol
%# the tolerance for the Newton algorithm [Default = 1e-7].
%#
%# @item  MaxNewtonIterations
%# maximum number of Newton iterations [Default = 100].
%#
%# @item  InexactSolver
%# the Newton solver. It can be set to either ["newton_raphson"],
%# "inexact_newton" or "fsolve".
%# @end table
%#
%# If the "inexact_newton" solver is set, a further set of options can be
%# provided:
%#
%# @table @asis
%# @item  UseJacobian
%# specify if the Jacobian information should be used. It can be set to either
%# "yes" or ["no].
%#
%# @item  Eta
%# initial forcing term (must be in the interval [0,1)).  The dafault value is
%# set to 0.5.  For details see [1].
%#
%# @item  Choice
%# select the forcing term. It can be 1 (Default) or 2.  For details see [1].
%#
%# @item  Algorithm
%# iterative method to solve the linearized system.  It can be ["gmres"], "pcg"
%# or "bicgstab".
%#
%# @item  Restart
%# restart parameter for the GMRES solver (ignored for other solvers).
%# [Default = 20].
%# @end table
%#
%# The function typically returns two outputs.  Variable @var{t} is a
%# column vector and contains the times where the solution was found.  The
%# output @var{y} is a matrix in which each column refers to a different
%# unknown of the problem and each row corresponds to a time in @var{t}.
%#
%# The output can also be returned as a structure @var{solution} which
%# has field @var{x} containing the time where the solution was evaluated and
%# field @var{y} containing the solution matrix for the times in @var{x}.
%# Use @code{fieldnames (@var{solution})} to see the other fields and additional
%# information returned.
%#
%# If using the @code{"Events"} option then three additional outputs may
%# be returned.  @var{te} holds the time when an Event function returned a
%# zero.  @var{ye} holds the value of the solution at time @var{te}.  @var{ie}
%# contains an index indicating which Event function was triggered in the case
%# of multiple Event functions.
%#
%# For example, solve an anonymous implementation of the Van der Pol equation
%#
%# @example
%# @group
%# fvdp = @@(@var{t},@var{y}) [@var{y}(2); (1 - @var{y}(1)^2) * @var{y}(2) - @var{y}(1)];
%# [@var{t},@var{y}] = odebwe (fvdp, [0, 2], [2, 0]);
%# @end group
%# @end example
%#
%# References:
%#
%# [1] S.C. Eisenstat and H.F. Walker, "Choosing the Forcing Terms in an Inexact
%# Newton Method." SIAM Journal on Scientific Computing, 17(1), pp. 16-16, 1996.
%# @seealso{odeset, odefwe, inexact_newton}
%# @end deftypefn

function varargout = odebwe (fun, trange, init, varargin)

  if (nargin < 3)
    print_usage ();
  end
  
  order = 1; # we don't use local extrapolation
  solver = "odebwe";

  if (nargin >= 4)
    if (! isstruct (varargin{1}))
      %# varargin{1:len} are parameters for fun
      odeopts = odeset ();
      funarguments = varargin;
    elseif (length (varargin) > 1)
      %# varargin{1} is an ODE options structure opt
      odeopts = varargin{1};
      funarguments = {varargin{2:length(varargin)}};
    else  # if (isstruct (varargin{1}))
      odeopts = varargin{1};
      funarguments = {};
    end
  else  # nargin == 3
    odeopts = odeset ();
    funarguments = {};
  end

  if (! isnumeric (trange) || ! isvector (trange))
    error ("Octave:invalid-input-arg",
           "odebwe: TRANGE must be a numeric vector");
  end

  if (length (trange) < 2)
    error ("Octave:invalid-input-arg",
           "odebwe: TRANGE must contain at least 2 elements");
  elseif (trange(1) == trange(2))
    error ("Octave:invalid-input-arg",
           "odebwe: invalid time span, TRANGE(1) == TRANGE(2)");
  else
    direction = sign (trange(2) - trange(1));
  end
  trange = trange(:);

  if (! isnumeric (init) || ! isvector (init))
    error ("Octave:invalid-input-arg",
           "odebwe: INIT must be a numeric vector");
  end
  init = init(:);

  if (ischar (fun))
    try
      fun = str2func (fun);
    catch
      warning (lasterr);
    end
  end
  if (! isa (fun, "function_handle"))
    error ("Octave:invalid-input-arg",
           "odebwe: FUN must be a valid function handle");
  end

  %# Start preprocessing, have a look which options are set in odeopts,
  %# check if an invalid or unused option is set
  [defaults, classes, attributes] = odedefaults (numel (init), trange(1),
                                                 trange(end));
  persistent odebwe_ignore_options = ...
    {"BDF", "InitialSlope", "MassSingular", ...
     "MaxOrder", "MvPattern", "NonNegative"};

  defaults   = rmfield (defaults, odebwe_ignore_options);
  classes    = rmfield (classes, odebwe_ignore_options);
  attributes = rmfield (attributes, odebwe_ignore_options);

  % Specific options for odebwe
  defaults.NewtonTol   = 1e-7;
  classes.NewtonTol    = {"float"};
  attributes.NewtonTol = {"scalar", "positive"};
  defaults.MaxNewtonIterations   = 100;
  classes.MaxNewtonIterations    = {"float"};
  attributes.MaxNewtonIterations = {"scalar", "integer", "positive"};
  defaults.InexactSolver   = "newton_raphson";
  classes.InexactSolver    = {"char"};
  attributes.InexactSolver = {"newton_raphson", "inexact_newton", "fsolve"};
  # Specific options for Inexact Newton solver (see inexact_newton.m)
  defaults.UseJacobian   = "no";
  classes.UseJacobian    = {"char"};
  attributes.UseJacobian = {"yes", "no"};
  defaults.Eta   = 0.5;
  classes.Eta    = {"float"};
  attributes.Eta = {"scalar", ">=", 0, "<", 1};
  defaults.Choice   = 1;
  classes.Choice    = {"float"};
  attributes.Choice = {">=", 1, "<=", 2, "integer"};
  defaults.Algorithm   = "gmres";
  classes.Algorithm    = {"char"};
  attributes.Algorithm = {"gmres", "pcg", "bicgstab"};
  defaults.Restart   = 20;
  classes.Restart    = {"float"};
  attributes.Restart = {"integer"};

  odeopts = odemergeopts ("odebwe", odeopts, defaults, classes, attributes);

  odeopts.funarguments = funarguments;
  odeopts.direction    = direction;

  %# the following option enables the simplified Newton method
  %# which evaluates the Jacobian only once instead of the
  %# standard method that updates the Jacobian in each iteration
  odeopts.issimplified = false;
  odeopts.havenonnegative = false;

  if (isempty (odeopts.OutputFcn) && nargout == 0)
    odeopts.OutputFcn = @odeplot;
    odeopts.haveoutputfunction = true;
  else
    odeopts.haveoutputfunction = ! isempty (odeopts.OutputFcn);
  end

  if (isempty (odeopts.InitialStep))
    odeopts.InitialStep = odeopts.direction * ...
                          starting_stepsize (order, fun, trange(1),
                                             init, odeopts.AbsTol,
                                             odeopts.RelTol,
                                             strcmp (odeopts.NormControl,
                                             "on"), odeopts.funarguments);
  end

  if (! isempty (odeopts.Mass) && isnumeric (odeopts.Mass))
    havemasshandle = false;
    mass = odeopts.Mass;  # constant mass
  elseif (isa (odeopts.Mass, "function_handle"))
    havemasshandle = true;    # mass defined by a function handle
  else  # no mass matrix - creating a diag-matrix of ones for mass
    havemasshandle = false;
    odeopts.Mass = eye (length (init));
  end

  if (havemasshandle)   # Handle only the dynamic mass matrix,
    if (! strcmp (odeopts.MStateDependence, "none")) # constant mass matrices have already
      mass = @(t,x) odeopts.Mass (t, x, odeopts.funarguments{:});
      fun = @(t,x) mass (t, x, odeopts.funarguments{:}) ...
             \ fun (t, x, odeopts.funarguments{:});
    else           # if ((! strcmp (odeopts.MStateDependence, "none")) == false)
      mass = @(t) odeopts.Mass (t, odeopts.funarguments{:});
      fun = @(t,x) mass (t, odeopts.funarguments{:}) ...
             \ fun (t, x, odeopts.funarguments{:});
    end
  end

  if (nargout == 1)
    %# Single output requires auto-selected intermediate times,
    %# which is obtained by NOT specifying specific solution times.
    trange = [trange(1); trange(end)];
    odeopts.Refine = [];  # disable Refine when single output requested
  elseif (numel (trange) > 2)
    odeopts.Refine = [];  # disable Refine when specific times requested
  end

  if (strcmp (odeopts.InexactSolver, "newton_raphson"))
    solution = integrate_adaptive (@bwe_newton_raphson,
                                  order, fun, trange, init, odeopts);
  else
    solution = integrate_adaptive (@bwe_inexact_newton,
                                   order, fun, trange, init, odeopts);
  end

  %# Postprocessing, do whatever when terminating integration algorithm
  if (odeopts.haveoutputfunction)  # Cleanup plotter
    feval (odeopts.OutputFcn, [], [], "done", odeopts.funarguments{:});
  end
  if (! isempty (odeopts.Events))   # Cleanup event function handling
    ode_event_handler (odeopts.Events, solution.t(end), ...
                       solution.x(end,:).', "done", odeopts.funarguments{:});
  end

  %# Print additional information if option Stats is set
  if (strcmp (odeopts.Stats, "on"))
    havestats = true;
    nsteps    = solution.cntloop;                    %# vcntloop from 2..end
    nfailed   = solution.cntcycles - nsteps; %# vcntcycl from 1..end
    ndecomps  = 0;                             %# number of LU decompositions
    npds      = 0;                             %# number of partial derivatives
    nlinsols  = 0;                             %# no. of solutions of linear systems
    %# Print cost statistics if no output argument is given
    if (nargout == 0)
      vmsg = fprintf (1, "Number of successful steps: %d\n", nsteps);
      vmsg = fprintf (1, "Number of failed attempts:  %d\n", nfailed);
    end
  else
    havestats = false;
  end

  if (nargout == 1) %# Sort output variables, depends on nargout
    varargout{1}.x = solution.t.';  %# Time stamps are saved in field x
    varargout{1}.y = solution.x.';  %# Results are saved in field y
    varargout{1}.solver = solver;   %# Solver name is saved in field solver
    if (! isempty (odeopts.Events))
      varargout{1}.ie = solution.event{2};  %# Index info which event occurred
      varargout{1}.xe = solution.event{3};  %# Time info when an event occurred
      varargout{1}.ye = solution.event{4};  %# Results when an event occurred
    end
    if (havestats)
      varargout{1}.stats = struct;
      varargout{1}.stats.nsteps   = nsteps;
      varargout{1}.stats.nfailed  = nfailed;
      varargout{1}.stats.npds     = npds;
      varargout{1}.stats.ndecomps = ndecomps;
      varargout{1}.stats.nlinsols = nlinsols;
    end
  elseif (nargout == 2)
    varargout{1} = solution.t;   %# Time stamps are first output argument
    varargout{2} = solution.x;   %# Results are second output argument
  elseif (nargout == 5)
    varargout{1} = solution.t;   %# Same as (nargout == 2)
    varargout{2} = solution.x;   %# Same as (nargout == 2)
    varargout{3} = [];              %# LabMat doesn't accept lines like
    varargout{4} = [];              %# varargout{3} = varargout{4} = [];
    varargout{5} = [];
    if (! isempty (odeopts.Events))
      varargout{3} = solution.event{3};     %# Time info when an event occurred
      varargout{4} = solution.event{4};     %# Results when an event occurred
      varargout{5} = solution.event{2};     %# Index info which event occurred
    end
  end

end


%!demo
%! %# Demonstrate convergence order for odebwe
%! tol = 1e-4 ./ 10.^[0:5];
%! for i = 1 : numel (tol)
%!   opt = odeset ("RelTol", tol(i), "AbsTol", realmin);
%!   [t, y] = odebwe (@(t, y) -y, [0, 1], 1, opt);
%!   h(i) = 1 / (numel (t) - 1);
%!   err(i) = norm (y .* exp (t) - 1, Inf);
%! end
%!
%! %# Estimate order visually
%! loglog (h, tol, "-ob",
%!         h, err, "-b",
%!         h, (h/h(end)) .^ 1 .* tol(end), "k-");
%! axis tight
%! xlabel ("h");
%! ylabel ("err(h)");
%! title ("Convergence plot for odebwe");
%! legend ("imposed tolerance", "odebwe (relative) error",
%!         "order 1", "location", "northwest");
%!
%! %# Estimate order numerically
%! p = diff (log (err)) ./ diff (log (h))

%# We are using the Van der Pol equation for all tests that are done
%# for this function.
%# For further tests we also define a reference solution (computed at high
%# accuracy)
%!function ydot = fpol (t, y)  # The Van der Pol ODE
%!  ydot = [y(2); (1 - y(1)^2) * y(2) - y(1)];
%!function ref = fref ()       # The computed reference solution
%!  ref = [0.32331666704577, -1.83297456798624];
%!function [vjac] = fjac (vt, vy, varargin) %# its Jacobian
%!  vjac = [0, 1; -1 - 2 * vy(1) * vy(2), 1 - vy(1)^2];
%!function [vjac] = fjcc (vt, vy, varargin) %# sparse type
%!  vjac = sparse ([0, 1; -1 - 2 * vy(1) * vy(2), 1 - vy(1)^2]);
%!function [val, trm, dir] = feve (t, y, varargin)
%!  val = fpol (t, y, varargin);    # We use the derivatives
%!  trm = zeros (2,1);              # that's why component 2
%!  dir = ones (2,1);               # does not seem to be exact
%!function [val, trm, dir] = fevn (t, y, varargin)
%!  val = fpol (t, y, varargin);    # We use the derivatives
%!  trm = ones (2,1);               # that's why component 2
%!  dir = ones (2,1);               # does not seem to be exact
%!function mas = fmas (t, y, varargin)
%!  mas = [1, 0; 0, 1];            # Dummy mass matrix for tests
%!function mas = fmsa (t, y, varargin)
%!  mas = sparse ([1, 0; 0, 1]);   # A sparse dummy matrix
%!function out = fout (t, y, flag, varargin)
%!  out = false;
%!  if (strcmp (flag, "init"))
%!    if (! isequal (size (t), [2, 1]))
%!      error ('fout: step "init"');
%!    end
%!  elseif (isempty (flag))
%!    if (! isequal (size (t), [1, 1]))
%!      error ('fout: step "calc"');
%!    end
%!  elseif (strcmp (flag, "done"))
%!    if (! isempty (t))
%!      warning ('fout: step "done"');
%!    end
%!  else
%!    error ("fout: invalid flag <%s>", flag);
%!  end
%!
%!test  # two output arguments
%! [t, y] = odebwe (@fpol, [0 2], [2 0]);
%! assert ([t(end), y(end,:)], [2, fref], 5e-2);
%!test  # not too many steps
%! [t, y] = odebwe (@fpol, [0 2], [2 0]);
%! assert (size (t) < 120);
%!test  # anonymous function instead of real function
%! fvdp = @(t,y) [y(2); (1 - y(1)^2) * y(2) - y(1)];
%! [t, y] = odebwe (fvdp, [0 2], [2 0]);
%! assert ([t(end), y(end,:)], [2, fref], 5e-2);
%!test  # string instead of function
%! [t, y] = odebwe ("fpol", [0 2], [2 0]);
%! assert ([t(end), y(end,:)], [2, fref], 5e-2);
%!test  # extra input arguments passed through
%! [t, y] = odebwe (@fpol, [0 2], [2 0], 12, 13, "KL");
%! assert ([t(end), y(end,:)], [2, fref], 5e-2);
%!test  # empty ODEOPT structure *but* extra input arguments
%! opt = odeset;
%! [t, y] = odebwe (@fpol, [0 2], [2 0], opt, 12, 13, "KL");
%! assert ([t(end), y(end,:)], [2, fref], 5e-2);
%!test  # Solve another anonymous function below zero
%! vref = [0, 14.77810590694212];
%! [t, y] = odebwe (@(t,y) y, [-2 0], 2);
%! assert ([t(end), y(end,:)], vref, 5e-1);
%!test  # InitialStep option
%! opt = odeset ("InitialStep", 1e-8);
%! [t, y] = odebwe (@fpol, [0 0.2], [2 0], opt);
%! assert ([t(2)-t(1)], [1e-8], 1e-9);
%!test  # MaxStep option
%! opt = odeset ("MaxStep", 1e-3);
%! sol = odebwe (@fpol, [0 0.2], [2 0], opt);
%! assert ([sol.x(5)-sol.x(4)], [1e-3], 1e-2);
%!test  # Solve with intermediate step
%! [t, y] = odebwe (@fpol, [0 1 2], [2 0]);
%! assert (any((t-1) == 0));
%! assert ([t(end), y(end,:)], [2, fref], 5e-2);
%!test  # Solve in backward direction starting at t=0
%! vref = [-1.205364552835178, 0.951542399860817];
%! sol = odebwe (@fpol, [0 -2], [2 0]);
%! assert ([sol.x(end); sol.y(:,end)], [-2; vref'], 1e-1);
%!test  # Solve in backward direction starting at t=2
%! vref = [-1.205364552835178, 0.951542399860817];
%! sol = odebwe (@fpol, [2 -2], fref);
%! assert ([sol.x(end); sol.y(:,end)], [-2; vref'], 1e-1);
%!test  # Solve in backward direction starting at t=2, with intermediate step
%! vref = [-1.205364552835178, 0.951542399860817];
%! [t, y] = odebwe (@fpol, [2 0 -2], fref);
%! idx = find(y < 0, 1, "first") - 1;
%! assert ([t(idx), y(idx,:)], [0,2,0], 5e-1);
%! assert ([t(end), y(end,:)], [-2, vref], 5e-2);
%!test  # Solve another anonymous function in backward direction
%! vref = [-1, 0.367879437558975];
%! sol = odebwe (@(t,y) y, [0 -1], 1);
%! assert ([sol.x(end); sol.y(:,end)], vref', 1e-2);
%!test  # Solve another anonymous function below zero
%! vref = [0, 14.77810590694212];
%! sol = odebwe (@(t,y) y, [-2 0], 2);
%! assert ([sol.x(end); sol.y(:,end)], vref', 5e-1);
%!test  # Solve in backward direction starting at t=0 with MaxStep option
%! vref = [-1.205364552835178, 0.951542399860817];
%! opt = odeset ("MaxStep", 1e-3);
%! sol = odebwe (@fpol, [0 -2], [2 0], opt);
%! assert ([abs(sol.x(8)-sol.x(7))], [1e-3], 1e-2);
%! assert ([sol.x(end); sol.y(:,end)], [-2; vref'], 1e-2);
%!test  # AbsTol option
%! opt = odeset ("AbsTol", 1e-5);
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test  # AbsTol and RelTol option
%! opt = odeset ("AbsTol", 1e-8, "RelTol", 1e-8);
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-2);
%!test  # RelTol and NormControl option -- higher accuracy
%! opt = odeset ("RelTol", 1e-8, "NormControl", "on");
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 2e-3);
%!test  # Details of OutputSel and Refine can't be tested
%! opt = odeset ("OutputFcn", @fout, "OutputSel", 1, "Refine", 5);
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%!test  # Stats must add further elements in sol
%! opt = odeset ("Stats", "on");
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert (isfield (sol, "stats"));
%! assert (isfield (sol.stats, "nsteps"));
%!test  # Events option add further elements in sol
%! opt = odeset ("Events", @feve);
%! sol = odebwe (@fpol, [0 10], [2 0], opt);
%! assert (isfield (sol, "ie"));
%! assert (sol.ie(1), 2);
%! assert (isfield (sol, "xe"));
%! assert (isfield (sol, "ye"));
%!test  # Events option, now stop integration
%! warning ("off", "integrate_adaptive:unexpected_termination", "local");
%! opt = odeset ("Events", @fevn, "NormControl", "on");
%! sol = odebwe (@fpol, [0 10], [2 0], opt);
%! assert ([sol.ie, sol.xe, sol.ye],
%!         [2.0, 2.496110, -0.830550, -2.677589], 6e-1);
%!test  # Events option, five output arguments
%! warning ("off", "integrate_adaptive:unexpected_termination", "local");
%! opt = odeset ("Events", @fevn, "NormControl", "on");
%! [t, y, vxe, ye, vie] = odebwe (@fpol, [0 10], [2 0], opt);
%! assert ([vie, vxe, ye],
%!         [2.0, 2.496110, -0.830550, -2.677589], 6e-1);
%!test  # Mass option as function
%! opt = odeset ("Mass", @fmas);
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test  # Mass option as matrix
%! opt = odeset ("Mass", eye (2,2));
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test  # Mass option as sparse matrix
%! opt = odeset ("Mass", sparse (eye (2,2)));
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test  # Mass option as function and sparse matrix
%! opt = odeset ("Mass", @fmsa);
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test  # Mass option as function and MStateDependence
%! opt = odeset ("Mass", @fmas, "MStateDependence", "strong");
%! sol = odebwe (@fpol, [0 2], [2 0], opt);
%! assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test %# Jacobian option
%!  opt = odeset ('Jacobian', @fjac);
%!  sol = odebwe (@fpol, [0 2], [2 0], opt);
%!  assert ([sol.x(end); sol.y(:,end)], [2; fref'], 5e-2);
%!test %# Jacobian option and sparse return value
%!  opt = odeset ('Jacobian', @fjcc);
%!  sol = odebwe (@fpol, [0 2], [2 0], opt);
%!  assert ([sol.x(end); sol.y(:,end)], [2; fref'], 1e-1);
%!test %# inexact_newton
%!  opt = odeset ('InexactSolver', 'inexact_newton', 'NewtonTol', 1e-8, 'Algorithm', 'gmres');
%!  sol = odebwe (@fpol, [0 2], [2 0], opt);
%!  assert ([sol.x(end); sol.y(:,end)], [2; fref'], 5e-2);
%!test %# inexact_newton and pcg
%!  opt = odeset ('InexactSolver', 'inexact_newton', 'NewtonTol', 1e-8, 'Algorithm', 'gmres');
%!  sol = odebwe (@fpol, [0 2], [2 0], opt);
%!  assert ([sol.x(end); sol.y(:,end)], [2; fref'], 5e-2);
%!test %# inexact_newton and UseJacobian
%!  opt = odeset ('InexactSolver', 'inexact_newton', 'Jacobian', @fjac, 'UseJacobian', 'yes');
%!  sol = odebwe (@fpol, [0 2], [2 0], opt);
%!  assert ([sol.x(end); sol.y(:,end)], [2; fref'], 5e-2);
%!
%! %# test for JPattern option is missing
%! %# test for Vectorized option is missing
%!

%# Note: The following options have no effect on this solver
%#       therefore it makes no sense to test them here:
%#
%# "BDF"
%# "InitialSlope"
%# "MassSingular"
%# "MaxOrder"
%# "MvPattern"
%# "NonNegative"

%!test # Check that imaginary part of solution does not get inverted
%! sol = odebwe (@(x,y) 1, [0 1], 1i);
%! assert (imag (sol.y), ones (size (sol.y)))
%! [x, y] = odebwe (@(x,y) 1, [0 1], 1i);
%! assert (imag (y), ones (size (y)))

%!error odebwe ()
%!error odebwe (1)
%!error odebwe (1,2)
%!error <TRANGE must be a numeric> odebwe (@fpol, {[0 25]}, [3 15 1])
%!error <TRANGE must be a .* vector> odebwe (@fpol, [0 25; 25 0], [3 15 1])
%!error <TRANGE must contain at least 2 elements> odebwe (@fpol, [1], [3 15 1])
%!error <invalid time span> odebwe (@fpol, [1 1], [3 15 1])
%!error <INIT must be a numeric> odebwe (@fpol, [0 25], {[3 15 1]})
%!error <INIT must be a .* vector> odebwe (@fpol, [0 25], [3 15 1; 3 15 1])
%!error <FUN must be a valid function handle> odebwe (1, [0 25], [3 15 1])