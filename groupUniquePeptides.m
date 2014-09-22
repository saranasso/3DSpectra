% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
% This software is an implementation of the method described in	
% "3DSpectra: A 3-Dimensional Quantification Algorithm For LC-MS Labeled Profile Data" 	
% Copyright (C) 2012 Sara Nasso	
%     	
% 3DSpectra is free software; you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation; either version 2 of the License, or (at	
% your option) any later version.	
% 	
% 3DSpectra is distributed in the hope that it will be useful, but	
% WITHOUT ANY WARRANTY; without even the implied warranty of	
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU	
% General Public License for more details.	
% 	
% You should have received a copy of the GNU General Public License	
% along with this program; if not, write to the Free Software	
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307	
% USA	
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
function [resultsTmp]=groupUniquePeptides(results)
% it estimates one ratio taking into account all peptides occurences that have
% been quantified

pepSequences={results.pep_seq}';
[pepNames pepZeros]=strtok(pepSequences,'0');
replicates={results(:).replicate}';

for i=1:length(pepNames)
    pepNames{i,1}=[pepNames{i},replicates{i}];
    [pepname{i,1} pepLabel{i,1}]=getlabel(pepSequences{i});
end

resultsTmpIter=1;
analyzedIdxs=[];
for i=1:length(pepNames)
    if ~ismember(i,analyzedIdxs)
        clear ratioWeight pepRatio goodPepId
        currentPeptide=pepNames{i};
        idxPeptide = strmatch(currentPeptide, pepNames, 'exact');
        numPepOcc=numel(idxPeptide);
        if numPepOcc>1
            goodPepId=[];
            for j=1:numPepOcc
                    goodPepId=[goodPepId idxPeptide(j)];
            end
            pepResult.pep_charge=[results(goodPepId).pep_charge];
            pepResult.pep_idx=[results(goodPepId).pep_idx];
            resultsTmp(resultsTmpIter)=pepResult;
            resultsTmpIter=resultsTmpIter+1;
        else
            resultsTmp(resultsTmpIter)=results(idxPeptide);
            resultsTmpIter=resultsTmpIter+1;
        end
        analyzedIdxs=[analyzedIdxs;idxPeptide];
    end
end


	
