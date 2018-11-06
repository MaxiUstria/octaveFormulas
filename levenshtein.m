## Copyright (C) 2013, Roberto Porcu' <roberto.porcu@polimi.it>
## OdePkg - A package for solving ordinary differential equations and more
##
## This file is part of Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Command} {[@var{dist}] =} levenshtein (@var{"string1"}, @var{"string2"}, [@var{upper_bound}])
## @deftypefnx {Command} {[@var{dist}, @var{d}] =} levenshtein (@var{"string1"}, @var{"string2"}, [@var{upper_bound}])
##
## This function file can be used to compare two strings and it uses the
## Levenshtein distance as definition of metric between strings.
## More details at @url{http://en.wikipedia.org/wiki/Levenshtein_distance}.
##
## This function must be called with two output arguments: @var{dist} is the
## distance between the two strings and @var{d} is the matrix computed by
## Levenshtein algorithm.
##
## The first and the second input arguments are the two strings to be compared.
## This comparison is case-sensitive.
##
## The third argument is optional and fixes an upper bound for the distance.
## If the distance is greater than this limit then the function ends and
## returns a value equal to Inf.
##
## @seealso{odeset, odeget, fuzzy_compare}
## @end deftypefn

function [dist, d] = levenshtein(string1, string2, upper_bound)

  ## check on output arguments
  if( nargout>2 )
    error('OdePkg:InvalidArgument',...
      'too many output arguments');
  end

  ## check on input arguments
  if( nargin<2 || nargin >3 )
    error('OdePkg:InvalidArgument',...
      'input arguments number must be equal to 2 or 3');
  end
  if( ~ischar(string1) || ~ischar(string2) )
    error('OdePkg:InvalidArgument',...
      'input arguments must be strings');
  end

  ## lengths of the two strings
  m = length(string1);
  n = length(string2);

  ## taking lower cases of the two strings if no case sensitive comparison is required
  if( nargin!=3 )
    upper_bound = inf;
  elseif( isempty(upper_bound) )
    upper_bound = inf;
  end

  ## initializing the matrix
  d = zeros(m+1,n+1);
  d(2:m+1,1) = [1:1:m];
  d(1,2:n+1) = [1:1:n];

  ## the following algorithm is described in <http://en.wikipedia.org/wiki/Levenshtein_distance>
  for j=2:1:n+1
    for k=2:1:m+1
      if(string1(k-1)==string2(j-1))
        d(k,j)=d(k-1,j-1);
      else
        d(k,j)=min(d(k-1,j)+1,min(d(k,j-1)+1,d(k-1,j-1)+1));
      end
    end

    if( d(min(j,m+1),j) > upper_bound )
      dist = inf;
      return;
    end

  end

  ## the distance between the strings is given by the last element of the matrix
  dist = d(m+1,n+1);
end

%! # A simple test which is taken from <http://en.wikipedia.org/wiki/Levenshtein_distance>.
%!
%! ## Turn off output of warning messages for all tests, turn them on
%! ## again if the last test is called
%!test ## one output argument
%!  warning ('off', 'OdePkg:InvalidArgument');
%!  res=levenshtein('sitting','kitten');
%!  assert (res == 3);
%!test ## two output arguments
%!  [res,matrix] = levenshtein('sitting','kitten');
%!  ex_matrix = [0 1 2 3 4 5 6;
%!               1 1 2 3 4 5 6;
%!               2 2 1 2 3 4 5;
%!               3 3 2 1 2 3 4;
%!               4 4 3 2 1 2 3;
%!               5 5 4 3 2 2 3;
%!               6 6 5 4 3 3 2;
%!               7 7 6 5 4 4 3];
%!  assert (res == 3);
%!  assert (matrix == ex_matrix);
%!test ## two output arguments
%!  [res,matrix]=levenshtein('Sunday','Saturday');
%!  ex_matrix = [0 1 2 3 4 5 6 7 8;
%!               1 0 1 2 3 4 5 6 7;
%!               2 1 1 2 2 3 4 5 6;
%!               3 2 2 2 3 3 4 5 6;
%!               4 3 3 3 3 4 3 4 5;
%!               5 4 3 4 4 4 4 3 4;
%!               6 5 4 4 5 5 5 4 3];
%!  assert (res == 3);
%!  assert (matrix == ex_matrix);
%!test ## two output arguments case sensitiveness
%!  [res,matrix]=levenshtein('Sunday','saturday');
%!  assert (res == 4);
%!test ## two output arguments case no-sensitiveness
%!  [res,matrix]=levenshtein('Sunday','saturday',1);
%!  assert (res > 1e6);
%!  [res,matrix]=levenshtein('Sunday','saturday',10);
%!  assert (res == 4);
%!  warning ('on', 'OdePkg:InvalidArgument');

## Local Variables: ***
## mode: octave ***
## End: ***
