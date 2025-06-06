function b=dst2(arg1,mrows,ncols)

[m, n] = size(arg1);
% Basic algorithm.
if (nargin == 1),
  if (m > 1) & (n > 1),
    b = dst(dst(arg1).').';
    return;
  else
    mrows = m;
    ncols = n;
  end
end

% Padding for vector input.
a = arg1;
if nargin==2, ncols = mrows(2); mrows = mrows(1); end
mpad = mrows; npad = ncols;
if m == 1 & mpad > m, a(2, 1) = 0; m = 2; end
if n == 1 & npad > n, a(1, 2) = 0; n = 2; end
if m == 1, mpad = npad; npad = 1; end   % For row vector.

% Transform.

b = dst(a, mpad);
if m > 1 & n > 1, b = dst(b.', npad).'; end
