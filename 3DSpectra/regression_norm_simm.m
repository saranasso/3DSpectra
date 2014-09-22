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
function [statistics]=regression_norm_simm(commons,ratio,all)
% it performs regression on all peptides and print it out on a file

[light cetera]= strtok(ratio,'to');
heavy=fliplr(strtok(fliplr(cetera),'ot'));

%% regression for maspectras
if all
    lightj=commons.light_J_all;
    heavyj=commons.heavy_J_all;
else
    lightj=commons.light_J;
    heavyj=commons.heavy_J;
end
max_dati=max([lightj;heavyj]);
lightj=lightj./max_dati;
heavyj=heavyj./max_dati;
N_data_2D=numel(lightj);

h=figure;
subplot 211
if ~all
    heading=[{['Ratio: ' light 'l' heavy 'h']};{'Maspectras: ASAPRatio'}; {'Common peptides'}];
else
    heading=[{['Ratio: ' light 'l' heavy 'h']};{'Maspectras: ASAPRatio'}; {'All peptides'}];
end
title(heading)
statistics.ASAPRatio=getRegression(heavyj,lightj);


%% regression for 3DSpectra
if all
    lights=commons.light_S_all;
    heavys=commons.heavy_S_all;
else
    lights=commons.light_S;
    heavys=commons.heavy_S;
end
max_dati=max([lights;heavys]);
lights=lights./max_dati;
heavys=heavys./max_dati;
N_data_3D=numel(lights);

subplot 212
statistics.Spectra3D=getRegression(heavys,lights);


if ~all
    [h1 z_value]=corr_diff(statistics.Spectra3D.R2,N_data_3D,statistics.ASAPRatio.R2,N_data_2D);
    print_file=[pwd filesep 'regressions_commons.ps'];
    if h1==1
        heading=[{['Ratio: ' light 'l' heavy 'h']}; '3DSpectra' ; {'Common peptides'};{'Statistically different at 5%'}];
    else
        heading=[{['Ratio: ' light 'l' heavy 'h']}; '3DSpectra' ; {'Common peptides'};{'Not statistically different at 5%'}];
    end
else
    heading=[{['Ratio: ' light 'l' heavy 'h']}; '3DSpectra' ; {'All peptides'}];
    print_file=[pwd filesep 'regressions_all.ps'];
end

title(heading)

%% printing to file
print(h,'-dpsc','-cmyk','-append',print_file)
pause(2)
close all	
