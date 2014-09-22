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
function setYaxis(mz_peaks_pep, old_rti,old_rtf,az,el,y_loci)
% it enables the visualization of correct labels

axis tight
hold on;
ylabel('m/z')
xlabel('rtt')

x_label=(old_rti:10:old_rtf);
set(gca,'XtickLabel',x_label)
ylabel('m/z')
peak_dist_step=abs(mz_peaks_pep(1)-mz_peaks_pep(2));


mzi=mz_peaks_pep(1);
mzf=mz_peaks_pep(end);

set(gca,'Ytick',y_loci)
y_label=(mzi:peak_dist_step:mzf);
set(gca,'YtickLabel',y_label)

set(gcf,'Color',[1,1,1])
pause(1)
view(az, el);
pause(1)	
