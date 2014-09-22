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
function queryBounds=defineQueryBounds(mzi,mzf,rti,rtf,startScan,endScan,MS_boolean,scan_representation)
% it constructs the object defining the bounds for the range query we would
% like to perform on our dataset stored as mzXML. We can read only the file info
% just passing the first 6 arguments as [].

if ~isempty(mzi) && ~isempty(mzf)
    queryBounds.mzi=mzi;
    queryBounds.mzf=mzf;
else
    queryBounds.mzi=[];
    queryBounds.mzf=[];
end

if ~isempty(rti) && ~isempty(rtf)
    queryBounds.rti=rti;
    queryBounds.rtf=rtf;
else
    queryBounds.rti=[];
    queryBounds.rtf=[];
end

if ~isempty(startScan) && ~isempty(endScan)
    queryBounds.startScan=startScan;
    queryBounds.endScan=endScan;
else
    queryBounds.startScan=[];
    queryBounds.endScan=[];
end

if ~isempty(MS_boolean)
    queryBounds.MS1=MS_boolean;
else
    queryBounds.MS1=[];
end

if ~isempty(scan_representation)
    queryBounds.scan_representation=scan_representation;
else
    queryBounds.scan_representation=0;
end	
