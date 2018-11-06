%% The following function is copied from the optim
%% package of Octave-Forge
%% Copyright (C) 1992-1994 Richard Shrager
%% Copyright (C) 1992-1994 Arthur Jutan
%% Copyright (C) 1992-1994 Ray Muzic
%% Copyright (C) 2010, 2011 Olaf Till <olaf.till@uni-jena.de>

function prt = __dfxpdp__ (p, func, rtol, pattern)

  f = func (p)(:);
  m = numel (f);
  n = numel (p);

  diffp = rtol .* ones (n, 1);
  sparse_jac = false;
  if (nargin > 3 && issparse (pattern))
    sparse_jac = true;
  end
  %% initialise Jacobian to Zero

  if (sparse_jac)
    prt = pattern;
  else
    prt = zeros (m, n);
  end

  del = ifelse (p == 0, diffp, diffp .* p);
  absdel = abs (del);

  p2 = p1 = zeros (n, 1);

  %% double sided interval
  p1 = p + absdel/2;
  p2 = p - absdel/2;

  ps = p;
  if (! sparse_jac)

    for j = 1:n
      ps(j) = p1(j);
      tp1 = func (ps);
      ps(j) = p2(j);
      tp2 = func (ps);
      prt(:, j) = (tp1(:) - tp2(:)) / absdel(j);
      ps(j) = p(j);
    end

  else

    for j = find (any (pattern, 1))
      ps(j) = p1(j);
      tp1 = func (ps);
      ps(j) = p2(j);
      tp2 = func (ps);
      nnz = find (pattern(:, j));
      prt(nnz, j) = (tp1(nnz) - tp2(nnz)) / absdel(j);
      ps(j) = p(j);
    end

  end

end
