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
function peptide=mzRTreeAccess(filePath,mzi,mzf,rti,rtf,varargin)
% it accesses the mzRTree and returns peptide data

warning off all
mzRTree = javaObject('mzRTree.MzRTree', filePath);

peptide=double(mzRTree.range_query(rti-1, rtf, mzi+0.001, mzf+0.001)');

if ~isempty(varargin)
    
   tosum=varargin{1};
   
   peptide=bin(peptide,tosum);
   
elseif isempty(varargin)
    
    peptide.int=peptide;
  
end

load([filePath 'information.mat'])
peptide.info=mzRTreeInfo;
	
