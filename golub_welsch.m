%# Copyright (C) 2013, Mattia Penati
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
%# @deftypefn {Command} {[@var{nodes}] =} golub_welsch (@var{m},@var{n})
%# @deftypefnx {Command} {[@var{nodes}, @var{weights}] =} golub_welsch (@var{m},@var{n})
%# @deftypefnx {Command} {[@var{nodes}, @var{weights}, @var{values}] =} golub_welsch (@var{m},@var{n})
%# @deftypefnx {Command} {[@var{nodes}, @var{weights}, @var{values}, @var{der}] =} golub_welsch (@var{m},@var{n})
%# 
%# This function can be used to compute Gauss quadrature nodes and weights and values of Legendre polynomials (and their derivatives) at Gauss quadrature nodes just in one shot, whitout using an iterative method. For all the theory about this topic see [1].
%# 
%# First output parameter contains Gauss quadrature nodes in @math{(-1,1)}.
%# 
%# Second output argument contains computed Gauss quadrature weights.
%# 
%# Third output parameter contains the values of Legendre polynomials at computed Gauss quadrature nodes.
%# 
%# Fourth output argument contains the values of the derivatives of Legendre polynomials at Gauss quadrature nodes.
%# 
%# First input argument must be a positive integer scalar and represents the maximum degree of Legendre polynomials.
%# 
%# Second input argument must be a positive integer scalar and represents the Gauss quadrature order.
%# 
%# References:
%# [1] G.H. Golub and J.H. Welsch, "Calculation of Gauss Quadrature Rules." Mathematics of computation, Vol. 23, No. 106 (Apr. 1969), pp. 221-230+s1-s10, American Mathematical Society.
%#
%# @seealso{odepkg}
%# @end deftypefn

function [nodes, weights, values, der] = golub_welsch (m, n)

  %# m is equal to the polynomial degree
  %# n is equal to the quadrature order

  Nvec = [1:n-1];

  %# definition of subdiagonal and overdiagonal elements
  beta = sqrt((Nvec./(2.*Nvec-1)).*(Nvec./(2.*Nvec+1)));

  %# building up the matrix A
  A = diag(beta,-1)+diag(beta,+1);

  [v,D] = eig(A);

  %# Gauss quadrature nodes are the eigenvalues of A
  nodes = diag(D)';

  %# computing Gauss quadrature weights
  weights = (2.*(v(1,:).^2))';
  
  Nvec = [0,Nvec]';
  temp = sqrt(2.*Nvec+1) * v(1,:);

  %# normalizing eigenvectors
  v = v./temp;

  %# initializing values of Legendre polynomials evaluated at Gauss quadrature nodes
  values = zeros(m,n);

  %# computing values of Legendre polynomials evaluated at Gauss quadrature nodes, as described in [1]
  values(1:min(m,n),:) = v(1:min(m,n),:);
  if(m > n)
    if(n==1)
      values(2,:) = nodes.*values(1,:);
    end
    for n = max(2,n):(m-1)
      m = 1.0 / n;
      values(n+1,:) = (2-m).*nodes.*values(n,:) - (1-m).*values(n-1,:);
    end
  end

  %# computing values of derivative Legendre polynomials evaluated at Gauss quadrature nodes, as described in [1]
  der = zeros(size(values));
  if(m>1)
    der(2,:) = (nodes.*values(2,:) -1) ./ (nodes.^2 -1);
    for n=2:(m-1)
      der(n+1,:) = der(n-1,:) + (2*n-1).*values(n,:);
    end
  end

  values = values';
  der = der';
end

%!test
%!  [nodes,w,L,D] = golub_welsch(2,2);
%!test
%!  [nodes,w,L,D] = golub_welsch(2,4);
%!test
%!  [nodes,w,L,D] = golub_welsch(3,2);

%# Local Variables: ***
%# mode: octave ***
%# End: ***
