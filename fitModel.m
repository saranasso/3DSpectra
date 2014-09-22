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
function [spectrum flag]=fitModel(peak_distrA,isot_distr,varargin)
% it fits the isotopic distribution model on the final spectrum using WLLS

spectrum(1,:,:)=peak_distrA';
w=isot_distr.prob_peaks;

I1A=lscov(isot_distr.i_peaks',peak_distrA',w); % weighted least squares
spectrum=isot_distr.i_peaks*I1A; % spectrum fitted by theoretical isotopic distribution

if ~isempty(varargin)
    pep_seq=varargin{1};
    isot_distr=distribution(pep_seq,2); % it takes into account the 2 main peaks
    isot_distr.prob_peaks=[isot_distr.prob_peaks, 0, 0];
end

[distr_sort idx_sort]=sort(isot_distr.i_peaks,'descend');
spectrum=spectrum.*isot_distr.prob_peaks; % each spectrum component is weighted by its probability
spectrum=spectrum(idx_sort);

if min(spectrum)>0
    flag=true;
else
    flag=false;
end

	
