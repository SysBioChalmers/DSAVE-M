%converts the dataset to TPM
function outObj = TPM(obj)
    if isobject(obj)
        outObj = obj;
        currSum = sum(obj.data);
        outObj.data = outObj.data * 1000000 ./ currSum;
        outObj.data(isnan(outObj.data)) = 0;%problems with division by zero
    else
        outObj = obj;
        currSum = sum(obj);
        outObj = outObj * 1000000 ./ currSum;
        outObj(isnan(outObj)) = 0;
    end
end