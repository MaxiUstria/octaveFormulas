%# Copyright (C) 2014, Jacopo Corno <jacopo.corno@gmail.com>
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
%# @deftypefn {Command} {[@var{t_next}, @var{x_next}] =} runge_kutta_78 (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%# @deftypefnx {Command} {[@var{t_next}, @var{x_next}, @var{x_est}] =} runge_kutta_78 (@var{@@fun}, @var{t}, @var{x}, @var{dt}, [@var{options}])
%#
%# This function can be used to integrate a system of ODEs with a given initial condition @var{x} from @var{t} to @var{t+dt}. For the definition of this method see p.91 in Ascher & Petzold.
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

function [t_next, x_next, x_est, k] = runge_kutta_78 (f, t, x, dt,
                                                      options = [],
                                                      k_vals = [],
                                                      t_next = t + dt)

  persistent va = [                    0,    0,      0,                        0,    0, 0, 0, 0, 0, 0, 0, 0;
                                    1/18,    0,      0,                        0,    0, 0, 0, 0, 0, 0, 0, 0;
                                    1/48, 1/16,      0,                        0,    0, 0, 0, 0, 0, 0, 0, 0;
                                    1/32,    0,   3/32,                        0,    0, 0, 0, 0, 0, 0, 0, 0;
                                    5/16,    0, -75/64,                    75/64,    0, 0, 0, 0, 0, 0, 0, 0;
                                    3/80,    0,      0,                     3/16, 3/20, 0, 0, 0, 0, 0, 0, 0;
                      29443841/614563906,    0,      0,       77736538/692538347, -28693883/1125000000, 23124283/1800000000, 0, 0, 0, 0, 0, 0;
                      16016141/946692911,    0,      0,       61564180/158732637, 22789713/633445777, 545815736/2771057229, -180193667/1043307555, 0, 0, 0, 0, 0;
                      39632708/573591083,    0,      0,     -433636366/683701615, -421739975/2616292301, 100302831/723423059, 790204164/839813087, 800635310/3783071287, 0, 0, 0, 0;
                    246121993/1340847787,    0,      0, -37695042795/15268766246, -309121744/1061227803, -12992083/490766935, 6005943493/2108947869, 393006217/1396673457, 123872331/1001029789, 0, 0, 0;
                   -1028468189/846180014,    0,      0,     8478235783/508512852, 1311729495/1432422823, -10304129995/1701304382, -48777925059/3047939560, 15336726248/1032824649, -45442868181/3398467696, 3065993473/597172653, 0, 0;
                     185892177/718116043,    0,      0,    -3185094517/667107341, -477755414/1098053517, -703635378/230739211, 5731566787/1027545527, 5232866602/850066563, -4093664535/808688257, 3962137247/1805957418, 65686358/487910083, 0;
                     403863854/491063109,    0,      0,    -5068492393/434740067, -411421997/543043805, 652783627/914296604, 11173962825/925320556, -13158990841/6184727034, 3936647629/1978049680, -160528059/685178525, 248638103/1413531060, 0];
  persistent vb7 = [13451932/455176623; 0; 0; 0; 0; -808719846/976000145; ...
                    1757004468/5645159321; 656045339/265891186; -3867574721/1518517206; ...
                       465885868/322736535; 53011238/667516719; 2/45; 0];
  persistent vb8 = [14005451/335480064; 0; 0; 0; 0; -59238493/1068277825; 181606767/758867731; ...
                    561292985/797845732; -1041891430/1371343529; 760417239/1151165299; ...
                    118820643/751138087; -528747749/2220607170; 1/4];
  persistent vc = sum (va, 2);

  if (! isempty (options))  # extra arguments for function evaluator
    args = options.funarguments;
  else
    args = {};
  end

  # one-step scheme as defined in <http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods>
  k(:,1)  = feval (f, t            , x                                                                                                                                                                                                                                  , args{:});
  k(:,2)  = feval (f, t +  vc(2)*dt, x + dt*  va(2,1)*k(:,1)                                                                                                                                                                                                            , args{:});
  k(:,3)  = feval (f, t +  vc(3)*dt, x + dt*( va(3,1)*k(:,1) + va(3,2)*k(:,2))                                                                                                                                                                                          , args{:});
  k(:,4)  = feval (f, t +  vc(4)*dt, x + dt*( va(4,1)*k(:,1)                   + va(4,3)*k(:,3))                                                                                                                                                                        , args{:});
  k(:,5)  = feval (f, t +  vc(5)*dt, x + dt*( va(5,1)*k(:,1)                   + va(5,3)*k(:,3) +  va(5,4)*k(:,4))                                                                                                                                                      , args{:});
  k(:,6)  = feval (f, t +  vc(6)*dt, x + dt*( va(6,1)*k(:,1)                                    +  va(6,4)*k(:,4) +  va(6,5)*k(:,5))                                                                                                                                    , args{:});
  k(:,7)  = feval (f, t +  vc(7)*dt, x + dt*( va(7,1)*k(:,1)                                    +  va(7,4)*k(:,4) +  va(7,5)*k(:,5) +  va(7,6)*k(:,6))                                                                                                                  , args{:});
  k(:,8)  = feval (f, t +  vc(8)*dt, x + dt*( va(8,1)*k(:,1)                                    +  va(8,4)*k(:,4) +  va(8,5)*k(:,5) +  va(8,6)*k(:,6) +  va(8,7)*k(:,7))                                                                                                , args{:});
  k(:,9)  = feval (f, t +  vc(9)*dt, x + dt*( va(9,1)*k(:,1)                                    +  va(9,4)*k(:,4) +  va(9,5)*k(:,5) +  va(9,6)*k(:,6) +  va(9,7)*k(:,7) +  va(9,8)*k(:,8))                                                                              , args{:});
  k(:,10) = feval (f, t + vc(10)*dt, x + dt*(va(10,1)*k(:,1)                                    + va(10,4)*k(:,4) + va(10,5)*k(:,5) + va(10,6)*k(:,6) + va(10,7)*k(:,7) + va(10,8)*k(:,8) + va(10,9)*k(:,9))                                                            , args{:});
  k(:,11) = feval (f, t + vc(11)*dt, x + dt*(va(11,1)*k(:,1)                                    + va(11,4)*k(:,4) + va(11,5)*k(:,5) + va(11,6)*k(:,6) + va(11,7)*k(:,7) + va(11,8)*k(:,8) + va(11,9)*k(:,9) + va(11,10)*k(:,10))                                        , args{:});
  k(:,12) = feval (f, t + vc(12)*dt, x + dt*(va(12,1)*k(:,1)                                    + va(12,4)*k(:,4) + va(12,5)*k(:,5) + va(12,6)*k(:,6) + va(12,7)*k(:,7) + va(12,8)*k(:,8) + va(12,9)*k(:,9) + va(12,10)*k(:,10) + va(12,11)*k(:,11))                    , args{:});
  k(:,13) = feval (f, t + vc(13)*dt, x + dt*(va(13,1)*k(:,1)                                    + va(13,4)*k(:,4) + va(13,5)*k(:,5) + va(13,6)*k(:,6) + va(13,7)*k(:,7) + va(13,8)*k(:,8) + va(13,9)*k(:,9) + va(13,10)*k(:,10) + va(13,11)*k(:,11) + va(13,12)*k(:,12)), args{:});

  %# computing new time and new values for the unkwnowns
  # t_next = t + dt; %t_next
  x_next = x + dt*(vb7(1)*k(:,1) + vb7(6)*k(:,6) + vb7(7)*k(:,7) + vb7(8)*k(:,8) + vb7(9)*k(:,9) + vb7(10)*k(:,10) + vb7(11)*k(:,11) + vb7(12)*k(:,12));
  
  %# if the estimation of the error is required
  if (nargout >= 3)
    %# new solution to be compared with the previous one
    x_est= x + dt*(vb8(1)*k(:,1) + vb8(6)*k(:,6) + vb8(7)*k(:,7) + vb8(8)*k(:,8) + vb8(9)*k(:,9) + vb8(10)*k(:,10) + vb8(11)*k(:,11) + vb8(12)*k(:,12) + vb8(13)*k(:,13));
  end

end

%! # We are using the "Van der Pol" implementation.
%!function [ydot] = fpol (vt, vy) %# The Van der Pol
%!  ydot = [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%!end
%!
%!test
%!  [t,y] = runge_kutta_78 (@fpol, 0, [2;0], 0.05);
%!test
%!  [t,y,x] = runge_kutta_78 (@fpol, 0, [2;0], 0.1);

%# Local Variables: ***
%# mode: octave ***
%# End: ***
