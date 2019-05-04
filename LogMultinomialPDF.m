function ll = LogMultinomialPDF(obs, prob)
n = sum(obs,1);

%formula:
%x = obs;
%p = prob;
% PDF = n!/(x1! ... xn!) * p1^x1 * .... * pn^xn
% =>
% log(PDF) = log(n!) - log(x1!) ... - log(xn!) + log(p1^x1) + .... + log(pn*xn)


%generate list of n, n-1, ..1, cannot calculate n! for large n otherwise
facs = n:-1:1;

s1 = sum(log2(facs),2);
s2 = log2(factorial(obs));%don't sum yet
s3 = nansum(log2(prob).*obs);

% Some of the s2s may have become inf due to that matlab cannot calculate
% n! for large n (> 170). Replace those values by values calculated by the
% n method
sel = obs > 150;%take all large values; not sure about cancellations 
%below is for test
%sel = obs > 3;%add this line when running test below
ind = (1:size(obs,1)).';
selInd = ind(sel);
for i = 1:size(selInd)
    facs = obs(selInd(i)):-1:1;
    s2(selInd(i)) = sum(log2(facs),2);
end

s2 = -sum(s2,1);

ll = s1 + s2 + s3;



end

%{
test code: run the code below to test that this function gives the same
results as the built-in matlab function (which will not work for large
number of bins)

obs = [1;1;4;2;0;1];
prob = [0.2;0.2;0.1;0.3;0.1;0.1];
y = mnpdf(obs,prob);
log2(y)
ll = LogMultinomialPDF(obs, prob)

%log2(y) and ll should be the same, which seems to work!


%}

