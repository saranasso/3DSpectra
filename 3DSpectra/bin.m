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
function peptide=bin(peptide, tosum)
% it performs the binning operation on the m/z dimension summing up as many
% m/z values as tosum.

j=1;
rows=ceil(size(peptide,1)/tosum);
pep_tmp=zeros(rows,size(peptide,2));

for i=1:size(peptide,1)
    if mod(i,tosum)~=0
        pep_tmp(j,:)=pep_tmp(j,:)+peptide(i,:);
    elseif mod(i,tosum)==0
        pep_tmp(j,:)=pep_tmp(j,:)+peptide(i,:);
        j=j+1;
    end
end
peptide.int=pep_tmp;	
