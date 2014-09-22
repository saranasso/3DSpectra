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
function [MSall_to_MS1 tot_scan]=ms_ruler(filePath)
% association between MS>1 and level 1 data

% filePath='C:\Drosophila data\exp_026\20050412_KcNR_1.mzXML';
parser = javaObject('org.systemsbiology.jrap.stax.MSXMLParser', filePath);
tot_scan=parser.getScanCount();
j=0;
MSall_to_MS1=[];
% for i=start_scan:end_scan
for i=1:tot_scan
    current_scan= parser.rap(i);
    msLevel(i)=current_scan.getHeader().getMsLevel();
    retentionTime = current_scan.getHeader().getRetentionTime();
    RT{i} = java.lang.Float.parseFloat(retentionTime.substring(2, retentionTime.length() - 1));
%     peakList{i,1} = current_scan.getMassIntensityList();
    if msLevel(i)==1
        j=j+1;
    end
    MSall_to_MS1(i)=j;
end
	
