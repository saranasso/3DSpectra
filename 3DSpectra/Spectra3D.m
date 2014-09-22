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
function [results]=Spectra3D(results,library,fileHeader,settings,pep_idx,j)
% it implements the 3DSpectra algorithm

NlogL_pep=[];
noisy_gaussian_pep=[];

replicate=settings.replicate;
savename=[settings.ratio replicate];
npeaks=settings.npeaks;
diff_mass=settings.diff_mass;
k=settings.k;

peptides=[];
upperbound=length(library.peptidesSequence);
while pep_idx<=upperbound % on all peptides
    try
        pep_info.pep_seq=library.peptidesSequence{pep_idx}
        pepname=strtok(pep_info.pep_seq,'0');
        pep_info.pep_charge=library.charges(pep_idx);
        pep_info.pep_idx=pep_idx;
        indxs=strmatch([pepname num2str(pep_info.pep_charge)],peptides,'exact');
        
        if isempty(indxs) % peptide never analyzed before for this dataset
            
            % peptide metadata retrieving
            [pep_info.pep_seq, flag, idx, pep_info.label]=findPeptide(pep_info.pep_seq,pep_idx,library,pep_info.pep_charge);
            pep_info.pep_idx=idx;
            if flag==1 % il peptide è valido
                results.IDs=results.IDs+1;
                results.seqs{results.IDs}=pepname;
                [idx]=checkShift(library,idx);
                warning off all
                %% check in library whether also the isotopic partner was identified
                [flagPartner pep_info.partner_seq mzpartner]=findPartner(pep_info.pep_seq, library, idx, pep_info.label, diff_mass);
                
                if flagPartner
                    % access to peptide data using mzRTree
                    [pep_info.Chroms_partner pep_info.Chroms_pep mz_peaks_partner mz_peaks_pep peak_rtt peak_scan]=mzRTreeDataAccess(fileHeader,pep_info.partner_seq,mzpartner,idx,library,npeaks,settings);
                    %% isotopic distribution features
                    isot_distr_pep=distribution(pep_info.pep_seq,npeaks);
                    peak_distance=1/pep_info.pep_charge;
                    peak_dist=round(peak_distance*k);
                    distr_width=round(peak_distance*k*npeaks);
                    isot_distr_partner=distribution(pep_info.partner_seq,npeaks);
                    
                    %% Peak detection and GMM fit
                    mz_shifted_pep=library.massToSearch(idx);
                    mz_shift=mz_shifted_pep-mz_peaks_pep(1);
                    disp(pepname)
                    [ pep_info.scan_pep mz_pep w_mzi w_mzf w_rti w_rtf]=peakDetection(pep_info.Chroms_pep,settings,peak_scan,mz_shift,distr_width,isot_distr_pep,peak_dist);
                    [pep_info.Chroms_pep gmm_pep tmp_pep s_pep pep_info.mz_pep NlogL_pep(j) noisy_gaussian_pep(j) ]=getGmmModel(pep_info.Chroms_pep,pep_info.scan_pep,w_mzi,w_mzf,w_rti,w_rtf,peak_dist,pep_info.pep_seq,settings,mz_peaks_pep);
                    [ pep_info.scan_partner mz_partner w_mzi w_mzf w_rti w_rtf]=peakDetection(pep_info.Chroms_partner,settings,peak_scan,mz_shift,distr_width,isot_distr_partner,peak_dist);
                    [pep_info.Chroms_partner gmm_partner tmp_partner s_partner pep_info.mz_partner NlogL_partner(j) noisy_gaussian_partner(j) ]=getGmmModel(pep_info.Chroms_partner,pep_info.scan_partner,w_mzi,w_mzf,w_rti,w_rtf,peak_dist,pep_info.partner_seq,settings,mz_peaks_partner);
                    %% Quantitation
                    [results.quants]=quantitation(results.quants,npeaks,pep_info,j,k,replicate);
                    %% adding current peptide to the list of already analyzed peptides
                    peptides{j}=[pepname num2str(pep_info.pep_charge)];
                    j=j+1;
                    pep_idx=pep_idx+1
                elseif ~flagPartner
                    pep_idx=pep_idx+1
                end
            else
                pep_idx=pep_idx+1
            end
        else
            pep_idx=pep_idx+1
        end
    catch ME
%                 rethrow(ME)
        pep_idx=pep_idx+1
        continue
    end
    save([savename 'in_progress']) % useful to start it from last iteration in case of any failure
    clear pep_info
end

save(savename)	
