function logFacs = GetLogFactorials(n)
% GetLogFactorials
%   Calculates the logarithm of x! for all x <= n and places those in a vector. 
%   To get around that you cannot index a vector with 0, you should always 
%   index this vector with the value + 1.
%   The purpose of this function is that x! cannot be calculated for 
%   x > 170, while that is not a problem for log(x!), but you cannot use
%   the built-in factorial function in matlab and do log afterwards!
%   
% Input:
%   n               The max value needed
%
% Usage: logFacs = GetLogFactorials(n);
%        val = logFacs(k+1); % gets log(k!)
%
% Johan Gustafsson, 2019-11-12
%

a = 1:n;
loga = log(a);
logFacs = [0 cumsum(loga)].';

end


