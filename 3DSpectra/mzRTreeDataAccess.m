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
function [Chroms_partner Chroms_pep mz_peaks_partner mz_peaks_pep peak_rtt peak_scan]=mzRTreeDataAccess(fileHeader,partner,mzpartner,idx_pep,pep_library,npeaks,settings)
% this function accesses peptide and its partner data using mzRTree

%% some useful metadata
pep_seq=pep_library.peptidesSequence{idx_pep};
pep_charge=pep_library.charges(idx_pep);
peak_distance=1/pep_charge;
peak_scan=pep_library.MS_scan(idx_pep);
peak_rtt=peak_scan;
mzpeptide=((pep_library.massH1(idx_pep)+1)/pep_charge); 
mzpartner=mzpartner(1); % they're all the same value
numberofscans=numel(fileHeader.rtts);

% Partner isotopic distribution
isot_distr_partner=distribution(partner,npeaks);
indx=find(max(isot_distr_partner.i_peaks)); % since monoisotopic peak is not the most probable after 800Da 
mz_peaks_partner(1:indx)=mzpartner-(peak_distance*[0:(indx-1)]);
mz_peaks_partner((indx+1):npeaks)=mzpartner+(peak_distance*[1:(npeaks-indx)]);

% Peptide isotopic distribution
isot_distr_pep=distribution(pep_seq,npeaks);
peak_mz=(pep_library.massH1(idx_pep)+1)/pep_charge;
indx_pep=find(max(isot_distr_pep.i_peaks)); % since monoisotopic peak is not the most probable after 800Da 
mz_peaks_pep(1:indx_pep)=peak_mz-(peak_distance*[0:(indx_pep-1)]);
mz_peaks_pep((indx_pep+1):npeaks)=peak_mz+(peak_distance*[1:(npeaks-indx_pep)]);


% actual data access
Chroms_partner=mzRTreeAccess(fileHeader.localFullPath,mzpartner-1,mz_peaks_partner(end)+1,1,numberofscans); 
Chroms_pep=mzRTreeAccess(fileHeader.localFullPath,mz_peaks_pep(1)-1,mz_peaks_pep(end)+1,1,numberofscans); 
% 1 Da additional at end and beginning.

% visual inspection 
if settings.visualize
    hf=figure;
    subplot 211
    [pepname pepzero]=strtok(pep_seq,'0');
    title=['Peptide: ' pepname '. Charge: ' num2str(pep_charge)];
    visualInspection(Chroms_pep,[],[],[],peak_scan,mz_peaks_pep,npeaks,title,false)
    subplot 212
    [pepname pepzero]=strtok(partner,'0');
    title=['Partner: ' pepname '. Charge: ' num2str(pep_charge)];
    visualInspection(Chroms_partner,[],[],[],peak_scan,mz_peaks_partner,npeaks,title,false)
    attachPlotToFile(hf,'3DSpectra')
    close all
end


	
