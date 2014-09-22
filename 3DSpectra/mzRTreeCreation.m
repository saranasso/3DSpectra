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
function [mzRTreeInfo]=mzRTreeCreation(file_mzXML,file_mzrtree, number_of_strips, mz_resolution)
% it allows to efficiently create mzRTree data structure throughout Matlab
% if the mz_resolution is given as [] (hence, empty) the function will read
% the whole dataset file scan by scan to look for its mz_resolution (be careful: it can take quite a long time)

flag=0;
if isempty(mz_resolution)
    queryBounds=defineQueryBounds([],[],[],[],[],[],0,1);
    mzRTreeInfo=jrapAccess(file_mzXML,queryBounds,[],1);
    mz_resolution=mzRTreeInfo.mz_res;
    mzRTreeInfo.mzi=mzRTreeInfo.mz_min;
    mzRTreeInfo.mzf=mzRTreeInfo.mz_max;
    mzRTreeInfo.rti=mzRTreeInfo.rtt_min;
    mzRTreeInfo.rtf=mzRTreeInfo.rtt_max;
    mzRTreeInfo.folderName=fliplr(strtok(fliplr(mzRTreePath),filesep));
end

% mz_resolution_for_mzRTree=ceil(inv(mz_resolution));
mz_resolution_for_mzRTree=round(inv(mz_resolution));
disp(mz_resolution_for_mzRTree)
% return

while ~flag
    try
        mzRTree = javaObject('mzRTree.MzRTree', file_mzXML,file_mzrtree, mz_resolution_for_mzRTree,number_of_strips);
    catch
        disp('You need to increase the number of strips!')
        while ~isnumeric(number_of_strips)
            number_of_strips=input('Type the increased (integer) number of strips:');
        end
    end
    flag=1;
end

mzRTreeInfo.mz_resolution=mz_resolution;
mzRTreeInfo.time_resolution='scan representation: it depends on the scan sampling rate'; % we are using a scan representation of the data along the temporal dimension
save([file_mzrtree 'information.mat'])	
