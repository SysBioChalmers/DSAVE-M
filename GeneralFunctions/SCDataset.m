classdef SCDataset
   properties
      name %the name of the dataset
      data %genes as rows, samples as columns
      genes % one column with all gene names. Use name convention "GAPDH" etc and not "ENS..." or "NM..."
      cellIds % one row of all cell ids
      sampleIds % one row of all tumors/patients/some other group. Rename this?
      cellType % one row describing cell type, the source is what was published. Shall match the Celltype enum
      subCellType % one row which can be used for dividing the cells within a celltype into more categories
      extraCellInfo % one row (cell array of char arrays) that can be used to describe the cell in more detail
   end
   methods
       
       %usage:
       %l is a logical vector or a vector of indices
       %for example
       %ds2 = ds.cellSubsetLog([1,3,8:23]);
       %ds2 = ds.cellSubsetLog(strcmp(classification,'Macrophages/Monocytes'));
       %classification in the example above is expected to be a cell array
       %with the cell classifications
       function ds = cellSubset(this, l) %logical vector
           ds = SCDataset;
           ds.name = strcat(this.name, '_subset');
           ds.data = this.data(:, l);
           ds.genes = this.genes;
           ds.cellIds = this.cellIds(1, l);
           ds.sampleIds = this.sampleIds(1, l);
           if (~isempty(this.cellType))
            ds.cellType = this.cellType(1, l);
           end
           if (~isempty(this.subCellType))
               ds.subCellType = this.subCellType(1, l);
           end
           if (~isempty(this.extraCellInfo))
               ds.extraCellInfo = this.extraCellInfo(1, l);
           end
       end
       
       %randomly selects numCells cells from the dataset, without
       %replacement
       function ds = randSample(this, numCells)
           indToKeep = randsample(size(this.data,2),numCells);
           ds = this.cellSubset(indToKeep);
       end

       function [ds,ia] = geneSubset(this, genesToKeep) %genesToKeep should be a vertical cell array, logical array or array of indices
           ds = this;
           if isnumeric(genesToKeep) || islogical(genesToKeep)
               ds.genes = this.genes(genesToKeep);
               ds.data = ds.data(genesToKeep,:);
           else
               [ds.genes,ia,~] = intersect(this.genes, genesToKeep, 'stable');
               ds.data = ds.data(ia,:);
           end
       end
       
       function ds = setCellType(this, cellIds, values)
           ds = this;
           [~,ia,ib] = intersect(ds.cellIds, cellIds);
           ds.cellType(ia) = values(ib);
       end
       
       function ds = setSubCellType(this, cellIds, values)
           ds = this;
           [~,ia,ib] = intersect(ds.cellIds, cellIds);
           ds.subCellType(ia) = values(ib);
       end

       
       %removes any genes that do not exist in both datasets
       function dsRes = innerJoin(this, ds)
           dsRes = SCDataset;
           dsRes.name = strcat('inner join (', this.name,', ',ds.name, ')');
           [dsRes.genes, ia, ib] = intersect(this.genes, ds.genes);
           dsRes.data = [this.data(ia,:) ds.data(ib,:)];
           dsRes.cellIds = [this.cellIds ds.cellIds];
           dsRes.sampleIds = [this.sampleIds ds.sampleIds];
           dsRes.cellType = [this.cellType ds.cellType];
           dsRes.subCellType = [this.subCellType ds.subCellType];
           dsRes.extraCellInfo = [this.extraCellInfo ds.extraCellInfo];
       end
       
       %keeps all genes that exist in any dataset and sets them to zero for
       %cells where there is no data
       function dsRes = fullOuterJoin(this, ds)
           dsRes = SCDataset;
           dsRes.name = strcat('full outer join (', this.name,', ',ds.name, ')');
           dsRes.genes = union(this.genes, ds.genes);
           %join the data matrices and adapt to the new gene set
           [~,ires,ids] = intersect(dsRes.genes,ds.genes);
           [numGenes,~] = size(dsRes.genes);
           [~, numCells] = size(ds.cellIds);
           dsData = zeros(numGenes, numCells);
           dsData(ires,:) = ds.data(ids,:);
           [~,ires,ithis] = intersect(dsRes.genes,this.genes);
           [~, numCells] = size(this.cellIds);
           thisData = zeros(numGenes, numCells);
           thisData(ires,:) = this.data(ithis,:);
           dsRes.data = [thisData dsData];
           
           dsRes.cellIds = [this.cellIds ds.cellIds];
           dsRes.sampleIds = [this.sampleIds ds.sampleIds];
           dsRes.cellType = [this.cellType ds.cellType];
           dsRes.subCellType = [this.subCellType ds.subCellType];
           dsRes.extraCellInfo = [this.extraCellInfo ds.extraCellInfo];
       end
       
       %We could also implement left join, which would keep the genes of
       %this and fill in with zeros for ds if missing, and discard the 
       %extra genes from the ds set. If we need it...
       
       %will save the data, genes and cellIds as a table, the rest is not
       %saved
       function saveDataTable(this, filename)
           t = array2table(this.data);
           t.Properties.RowNames = this.genes;
           t.Properties.VariableNames = this.cellIds;
           writetable(t,filename,'Delimiter','\t','WriteRowNames',true);
       end
       
       function saveInfoTable(this, filename)
           t = table(this.cellType.', CelltypeId2CelltypeName(this.cellType).', this.sampleIds.','RowNames', this.cellIds.','VariableNames',{'celltype', 'celltype_text', 'sample_id'});
           %t.Properties.RowNames = this.genes;
           %t.Properties.VariableNames = this.cellIds;
           writetable(t,filename,'Delimiter','\t','WriteRowNames',true);
       end
       
       %if any of the non-critical member vectors are empty, they
       %will be filled with appropriate data.
       % cellIds - generated from name and a number series
       % sampleIds - all set to 'unknown'
       % cellType - all set to Celltype.Unknown
	   % subCellType - all set to 0
	   % extraCellInfo - all set to ''
       function ds = fillEmpties(this)
           %must copy and return a new val unless we make this a handle
           %class, which we do not want. In that case objects are 
           %normally not deep copied when copied, which will cause
           %a lot of confusion
           ds = this;
           [r,c] = size(ds.data);
           if isempty(ds.cellIds)
               formatString = strcat(ds.name, '_%d');
               ds.cellIds = sprintfc(formatString,(1:size(ds.data,2)));
           end
           if isempty(ds.sampleIds)
               [ds.sampleIds{1:c}] = deal('Unknown');
           end
           if isempty(ds.cellType)
               ds.cellType = ones(1,c);%one happens to be Celltype.Unknown.
           end
           if isempty(ds.subCellType)
               ds.subCellType = zeros(1,c);
           end
           if isempty(ds.extraCellInfo)
               [ds.extraCellInfo{1:c}] = deal('');
           end
           
       end
       
       function samples = splitIntoRandomGroupSamples(this, groupSize) %logical vector
           totNumCells = size(this.cellIds,2);
           indices = randperm(totNumCells);
           
           numSamp = floor(totNumCells/groupSize);
           numGenes = size(this.genes,1);
           
           samples = Samples();
           samples.genes = this.genes;
           samples.sampleIds = cell(1,numSamp);
           samples.data = zeros(numGenes,numSamp);
           
           for i = 1:numSamp
               st = 1 + (i-1)*numSamp;
               en = i*numSamp;
               ind = indices(st:en);
               subds = this.cellSubset(ind);
               samples.sampleIds{1,i} = strcat(this.name,'_',num2str(i));
               samples.data(:,i) = mean(subds.data,2);
           end 
       end
   end
end