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
function visualInspection(Chroms_pep,isot_distr,peak_distance,k,pep_scan,mz_peaks_pep,npeaks,titolo,peakpicking)
% it provides a rough visualization

if (pep_scan-200)>30 && (pep_scan+200)<(size(Chroms_pep.int,2)-30)
    int_max=max(max(Chroms_pep.int(:,(pep_scan-200):(pep_scan+200))));
elseif (pep_scan-200)<=30 && (pep_scan+200)<(size(Chroms_pep.int,2)-30)
    int_max=max(max(Chroms_pep.int(:,1:(pep_scan+200))));
elseif (pep_scan-200)<=30 && (pep_scan+200)>=(size(Chroms_pep.int,2)-30)
    int_max=max(max(Chroms_pep.int(:,1:(size(Chroms_pep.int,2)))));
elseif (pep_scan-200)>30 && (pep_scan+200)>=(size(Chroms_pep.int,2)-30)
    int_max=max(max(Chroms_pep.int(:,(pep_scan-200):(size(Chroms_pep.int,2)))));
end

surf(Chroms_pep.int);

shading interp;
axis tight
axis manual
hold on;
ylabel('m/z')
xlabel('rtt')
% x_label=(1:1:numberofscans);
% set(gca,'XtickLabel',x_label)
y_label=mz_peaks_pep;
set(gca,'YtickLabel',y_label)
title(titolo)

if peakpicking
    int_max=max(max(Chroms_pep.int));
    numberofscans=size(Chroms_pep.int,2);
    peak_dist=round(peak_distance*k);
    
    for j=1:npeaks
        stem3(pep_scan,k+(j-1)*peak_dist,int_max*isot_distr.i_peaks(j),'g','filled');
    end
end
pause(1)
az = 0;
el = 90;
view(az, el);
pause(1)
	
