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
function [ results]=quantitation(results,npeaks,pep_info,j,k,replicate)
% it performs the quantitation on filtered data 

pep_seq=pep_info.pep_seq;
partner_seq=pep_info.partner_seq;
pep_charge=pep_info.pep_charge;
label=pep_info.label;
pep_idx=pep_info.pep_idx;
isot_distr_pep=distribution(pep_seq,npeaks);
peak_distance=1/pep_charge;
isot_distr_partner=distribution(partner_seq,npeaks);
mz_pep=pep_info.mz_pep;
mz_partner=pep_info.mz_partner;
Chroms_pep=pep_info.Chroms_pep;
Chroms_partner=pep_info.Chroms_partner;

%% spectra processing
[Spectra_pep]=getSpectra(Chroms_pep,mz_pep,isot_distr_pep,peak_distance,k);
[Spectra_partner]=getSpectra(Chroms_partner,mz_partner,isot_distr_partner,peak_distance,k);

% %% max spectra/chromatograms (double point of view)
Max_Spectra_pep=quant_max(Spectra_pep);
Max_Spectra_partner=quant_max(Spectra_partner);

%% compute correlation  peptide-partner
if isequal(size(Chroms_pep.int),size(Chroms_partner.int))
    corr(1,1)=corr2(Chroms_pep.int,Chroms_partner.int);
else
    corr(1,1)=NaN;
end

%% quantitation
[spectrum_pep_m flag_pep]=fitModel(Max_Spectra_pep,isot_distr_pep,pep_seq);
[spectrum_partner_m flag_partner]=fitModel(Max_Spectra_partner,isot_distr_partner,partner_seq);
if flag_pep==true && flag_partner==true
    [results_max_spec]=getresults(spectrum_pep_m,spectrum_partner_m,pep_seq,pep_charge,pep_idx,label,replicate);
    results_max_spec.corr=corr;
    if ~isempty(results_max_spec.ratios)
        results(j)=results_max_spec;
    end
    clear results_max_spec mean_ratio_max_spec std_ratio_max_spec
end	
