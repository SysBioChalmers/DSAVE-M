%Reads housekeeping genes from file
%set forward to false if you want to delog
%numToAdd - If you don't want to do TPM + 1 since this ruins the lowly
%expressed genes, send in another number that corresponds to something
%a little bit lower than the lowest signal
function result = LogTrans(M, forward, numToAdd)
if nargin < 3
    numToAdd = 1;
end
if (forward)
    result = log(M+numToAdd);
else
    result = exp(1).^M - numToAdd;
end