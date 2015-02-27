function [d,beta] = stat_dprime(pHit,pFA)
%-- Convert to Z scores, no error checking
zHit = norminv(pHit,0,1) ;
zFA  = norminv(pFA,0,1) ;
%-- Calculate d-prime
d = zHit - zFA ;
%-- If requested, calculate BETA
if (nargout > 1)
  yHit = normpdf(zHit) ;
  yFA  = normpdf(zFA) ;
  beta = yHit ./ yFA ;
end
%end dprime()

function pdf = stdnormal_pdf (x)
%http://www.dynare.org/dynare-matlab-m2html/matlab/missing/stats/stdnormal_pdf.html
if (nargin ~= 1)
 error('stdnormal_pdf: you should provide one argument');
end
sz = size(x);
pdf = zeros (sz);
k = find (isnan (x));
if (any (k))
 pdf(k) = NaN;
end
k = find (~isinf (x));
if (any (k))
 pdf (k) = (2 * pi)^(- 1/2) * exp (- x(k) .^ 2 / 2);
end
%end stdnormal_pdf
 

function pdf = normpdf (x, m, s)
%http://www.dynare.org/dynare-matlab-m2html/matlab/missing/stats/normpdf.html
if (nargin ~= 1 && nargin ~= 3)
 error('normpdf: you must give one or three arguments');
end
if (nargin == 1)
 m = 0;
 s = 1;
end
if (~isscalar (m) || ~isscalar (s))
 [retval, x, m, s] = common_size (x, m, s);
 if (retval > 0)
     error ('normpdf: x, m and s must be of common size or scalars');
 end
end
sz = size (x);
pdf = zeros (sz);
if (isscalar (m) && isscalar (s))
 if (find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf)))
     pdf = NaN * ones (sz);
 else
     pdf = stdnormal_pdf ((x - m) ./ s) ./ s;
 end
else
 k = find (isinf (m) | isnan (m) | ~(s >= 0) | ~(s < Inf));
 if (any (k))
     pdf(k) = NaN;
 end
 k = find (~isinf (m) & ~isnan (m) & (s >= 0) & (s < Inf));
 if (any (k))
     pdf(k) = stdnormal_pdf ((x(k) - m(k)) ./ s(k)) ./ s(k);
 end
end
pdf((s == 0) & (x == m)) = Inf;
pdf((s == 0) & ((x < m) | (x > m))) = 0;
%end normpdf()
 
function x = norminv(p,m,s)
%http://fast-toolbox.googlecode.com/svn-history/r2/trunk/innards/norminv.m
% allocate output memory and check size of arguments
x = sqrt(2)*erfinv(2*p - 1).*s + m;  % if this line causes an error, input arguments do not fit.
x((p>1) | (p<0) | isnan(p) | isnan(m) | isnan(s) | (s<0)) = nan;
k = (s==0) & ~isnan(m);		% temporary variable, reduces number of tests.
x((p==0) & k) = -inf;
x((p==1) & k) = +inf;
k = (p>0) & (p<1) & k;
if prod(size(m))==1,
    x(k) = m;
else
    x(k) = m(k);
end
%end norminv()
