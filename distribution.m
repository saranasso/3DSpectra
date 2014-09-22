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
function [isot_distr]=distribution(sequence,npeaks)
% it computes the theoretical isotopic distribution. Output has 2 fields:
% i_peaks and prob_peaks

path3=[pwd filesep 'smconfig.xml'];
path4=[pwd filesep 'smconfig.custom.xml'];
parser = javaObject('at.tugraz.genome.maspectras.parser.spectrummill.SpectrumMillAAConfigParser',path3,path4);
javaMethod('parse',parser);
intensities = parser.calculatePeptideIntensityDistributionFromSequence(sequence,'carbamidomethylation','1=Oxidized-Methionine,2=ICPL_heavy,3=ICPL_light', npeaks);
i=intensities.toString();
i=char(i);
format long e
isot_distr.i_peaks=str2num(i);
probabilities = parser.getProbabilities();
p=probabilities.toString();
p=char(p);
isot_distr.prob_peaks=str2num(p);



	
