%reads a file typically exported from R
function s = ImportTabSepSamples(filename)
    f = readtable(filename, 'ReadVariableNames',true, 'ReadRowNames', true, 'Delimiter', '\t');
    s = Samples;
    s.data = table2array(f);
    s.genes = f.Properties.RowNames;
    s.sampleIds = f.Properties.VariableNames;
end


