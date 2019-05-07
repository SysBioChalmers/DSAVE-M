classdef Samples
   properties
      data %columns are samples, rows are genes
      genes % one column with all gene names. Use name convention "GAPDH" etc and not "ENS..." or "NM..."
      sampleIds % one row with sample ids
   end
   methods
       %usage:
       %l is a logical vector or a vector of indices
       function s = cellSubset(this, l) %logical vector
           s = Samples;
           s.data = this.data(:, l);
           s.genes = this.genes;
           s.sampleIds = this.sampleIds(1, l);
       end

       function [s,ia] = geneSubset(this, genesToKeep, sortGenes) %genesToKeep should be a vertical cell array
           if nargin < 3
               sortGenes = 0;
           end
           s = this;
           if sortGenes
             [s.genes,ia,~] = intersect(this.genes, genesToKeep, 'sorted');
           else
             [s.genes,ia,~] = intersect(this.genes, genesToKeep, 'stable');
           end
           s.data = s.data(ia,:);
       end
       
       %removes any genes that do not exist in both datasets
       function dsRes = innerJoin(this, ds)
           dsRes = Samples;
           %dsRes.name = strcat('inner join (', this.name,', ',ds.name, ')');
           [dsRes.genes, ia, ib] = intersect(this.genes, ds.genes);
           dsRes.data = [this.data(ia,:) ds.data(ib,:)];
           dsRes.sampleIds = [this.sampleIds ds.sampleIds];
           dsRes = TPM(dsRes);
       end
       
       %keeps all genes that exist in any dataset and sets them to zero for
       %cells where there is no data
       function sRes = fullOuterJoin(this, s)
           %not tested at all
           sRes = Samples;
           sRes.genes = union(this.genes, s.genes);
           %join the data matrices and adapt to the new gene set
           [~,ires,ids] = intersect(sRes.genes,s.genes);
           [numGenes,~] = size(sRes.genes);
           [~, numSamples] = size(s.sampleIds);
           sData = zeros(numGenes, numSamples);
           sData(ires,:) = s.data(ids,:);
           [~,ires,ithis] = intersect(sRes.genes,this.genes);
           [~, numSamples] = size(this.sampleIds);
           thisData = zeros(numGenes, numSamples);
           thisData(ires,:) = this.data(ithis,:);
           sRes.data = [thisData sData];
           sRes.sampleIds = [this.sampleIds s.sampleIds];
       end
       
       %We could also implement left join, which would keep the genes of
       %this and fill in with zeros for ds if missing, and discard the 
       %extra genes from the ds set. If we need it...
       
       %if any of the non-critical member vectors are empty, they
       %will be filled with appropriate data.
       % sampleIds - generated from name and a number series
       function s = fillEmpties(this)
           s = this;
           [r,c] = size(s.data);
           if isempty(s.sampleIds)
               formatString = strcat(s.name, '_%d');
               s.sampleIds = sprintfc(formatString,(1:size(s.data,2)));
           end
       end
       
       function writeToTextFile(this, filename)
           fileID = fopen(filename,'w');
           fprintf(fileID,'\t%s', this.sampleIds{:,:});
           fprintf(fileID,'\n');
           rows = size(this.data, 1);
           for i = 1:rows
               fprintf(fileID,'%s', this.genes{i});
               fprintf(fileID,'\t%f', this.data(i,:));
               fprintf(fileID,'\n');
           end
           fclose(fileID);
       end
   end
end