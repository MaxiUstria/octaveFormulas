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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} runge_kutta_45_fehlberg (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} runge_kutta_45_fehlberg (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%#
%# This function can be used to integrate a system of ODEs with a given initial condition @var{x} from @var{t} to @var{t+dt}, with the Fehlberg method. For the definition of this method see @url{http://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method}.
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
%# Fifth input parameter is optional and describes a set of options useful to adapt the computation to what is needed.
%#
%# @seealso{odepkg}
%# @end deftypefn

function varargout = runge_kutta_45_fehlberg (f, t, x, dt, varargin)

  options = varargin{1};
  k = zeros (size (x, 1), 6);
  
  %# one-step scheme as defined in <http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods>
  k(:,1) = feval (f, t, x, options.vfunarguments{:});
  k(:,2) = feval (f, t +   (1/4)*dt, x + dt*       (1/4)*k(:,1), options.vfunarguments{:});
  k(:,3) = feval (f, t +   (3/8)*dt, x + dt*(     (3/32)*k(:,1) +      (9/32)*k(:,2)), options.vfunarguments{:});
  k(:,4) = feval (f, t + (12/13)*dt, x + dt*((1932/2197)*k(:,1) - (7200/2197)*k(:,2) +  (7296/2197)*k(:,3)), options.vfunarguments{:});
  k(:,5) = feval (f, t +         dt, x + dt*(  (439/216)*k(:,1) -         (8)*k(:,2) +   (3680/513)*k(:,3) -    (845/4104)*k(:,4)), options.vfunarguments{:});
  k(:,6) = feval (f, t +   (1/2)*dt, x + dt*(    -(8/27)*k(:,1) +         (2)*k(:,2) -  (3544/2565)*k(:,3) +   (1859/4104)*k(:,4) - (11/40)*k(:,5)), options.vfunarguments{:});

  %# computing new time and new values for the unkwnowns
  varargout{1} = t + dt; %t_next
  varargout{2} = x + dt*((25/216)*k(:,1) + (1408/2565)*k(:,3) + (2197/4104)*k(:,4) - (1/5)*k(:,5)); %x_next

  %# if the estimation of the error is required
  if(nargout >= 3)
    %# new solution to be compared with the previous one
    varargout{3} = x + dt*((16/135)*k(:,1) + (6656/12825)*k(:,3) + (28561/56430)*k(:,4) - (9/50)*k(:,5) + (2/55)*k(:,6)); %x_est
    varargout{4} = k;
  end

end

%! # We are using the "Van der Pol" implementation.
%!function [ydot] = fpol (vt, vy) %# The Van der Pol
%!  ydot = [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!end
%!
%!shared opts
%!  opts = odeset;
%!  opts.vfunarguments = {};
%!
%!test
%!  [t,y] = runge_kutta_45_fehlberg (@fpol, 0, [2;0], 0.05, opts);
%!test
%!  [t,y,x] = runge_kutta_45_fehlberg (@fpol, 0, [2;0], 0.1, opts);

%# Local Variables: ***
%# mode: octave ***
%# End: ***
