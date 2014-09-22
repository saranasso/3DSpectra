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
function [commons]=commons(file_xls_S,results,file_xls_J,Foglio_J,outliers_removal,hasReplicates)
% it finds common peptides between 3DSpectra and ASAPRatio (MASPECTRAS) quants

%% reading 3DSpectra results
Foglio_S=char(inputname(2));
idxs_ = findstr('_', Foglio_S);
Foglio_S(idxs_)=' ';
[nums,texts,raws]=xlsread([pwd '\' file_xls_S],Foglio_S);

for i=1:size(texts,1)
    pep_quant{i,1}=strtok(texts{i,1}, '0');
    replicate{i,1}=texts{i,8};
end
pep_quant_no_rep=unique(pep_quant);

if hasReplicates
    for i=1:length(pep_quant)
        pep_quant{i,1}=[pep_quant{i,1} replicate{i,1}];
    end
end
heavys_all=cell2mat(raws(:,6));
idxs_to_keep= intersect(find(~isnan(heavys_all)),find((heavys_all~=0)));
heavys_all=heavys_all(idxs_to_keep);
lights_all=cell2mat(raws(idxs_to_keep,5));
quants_all=cell2mat(raws(idxs_to_keep,7));

[pep_quant_PRE quants_all_PRE quants_light_PRE quants_heavy_PRE charges_PRE replicate_PRE label_PRE id_PRE]=read3DSpectraResults(results);
for i=1:size(pep_quant_PRE,1)
    pep_quant_PRE{i,1}=strtok(pep_quant_PRE{i,1}, '0');
end

if hasReplicates
    for i=1:length(pep_quant_PRE)
        pep_quant_PRE{i,1}=[pep_quant_PRE{i,1} replicate_PRE{i,1}];
    end
end

%% reading MASPECTRAS results
cd ..
[pep_quant_j_PRE quantj_all_PRE lightj_all heavyj_all]=readMASPECTRASresults(file_xls_J,Foglio_J,hasReplicates);

if outliers_removal
    ratio=fliplr(strtok(fliplr(file_xls_J),'\'));
    k = strfind(file_xls_J, 'outliers_included');
    file_xls_J_no_out=[file_xls_J(1:k-1) ratio];
    try
        [pep_quant_j quantj_all lightj_all heavyj_all]=readMASPECTRASresults(file_xls_J_no_out,Foglio_J,hasReplicates);
    catch
        disp('cd to the 3DSpectra folder before running this code!')
        return
    end
else
   pep_quant_j=pep_quant_j_PRE;
   quantj_all=quantj_all_PRE;
end

%% determine intersection pre and post outlier removal
[pep_commons_PRE i_s i_j]=intersect(pep_quant_PRE,pep_quant_j_PRE);

[pep_commons i_s i_j]=intersect(pep_quant,pep_quant_j);

idxs_s=find(ismember(pep_quant,pep_commons)==1);
idxs_j=find(ismember(pep_quant_j,pep_commons)==1);

quants=cell2mat(raws(idxs_s,7));
heavys=cell2mat(raws(idxs_s,6));
lights=cell2mat(raws(idxs_s,5));

quantj=quantj_all(idxs_j);
lightj=lightj_all(idxs_j);
heavyj=heavyj_all(idxs_j);

%% to evaluate % correction by 3D
removed=setdiff(pep_commons_PRE,pep_commons);
[in3D i_s i_r]=intersect(pep_quant,removed);
[in2D i_j i_r]=intersect(pep_quant_j,removed);

quants_PRE=cell2mat(raws(i_s,7));
quantj_PRE=quantj_all(i_j);

mean_s_PRE=mean(quants_PRE);
mean_j_PRE=mean(quantj_PRE);
std_s_PRE=std(quants_PRE);
std_j_PRE=std(quantj_PRE);

%%to evaluate additional # of quants by 3DSpectra compared to 2D
plus3D=setdiff(pep_quant,pep_commons);
plus2D=setdiff(pep_quant_j,pep_commons);

%% statistics on common peptides
mean_s=mean(quants);
sd_s=std(quants);
cv_s=std(quants)/mean(quants)*100;
N_s=numel(quants);
quantj=quantj(quantj~=0);
mean_j=mean(quantj);
sd_j=std(quantj);
cv_j=std(quantj)/mean(quantj)*100;
N_j=numel(quantj);
stat_S=[mean_s;sd_s;cv_s;N_s];
stat_J=[mean_j;sd_j;cv_j;N_j];
info=[{'Mean'};{'SD'};{'CV %'};{'# quants'}];
if outliers_removal
    cd('without outliers')
else
    cd('with outliers')
end
xlswrite(file_xls_S,info,Foglio_S,'J2');
xlswrite(file_xls_S,{'3DSpectra'},Foglio_S,'K1');
xlswrite(file_xls_S,stat_S,Foglio_S,'K2');
xlswrite(file_xls_S,{'ASAPRatio'},Foglio_S,'L1');
xlswrite(file_xls_S,stat_J,Foglio_S,'L2');
if hasReplicates
    for i=1:length(pep_commons)
        tmp=pep_commons{i,1};
        pep_commons{i,1}=[tmp(1:end-1)];
    end
    pep_commons=unique(pep_commons);
end
xlswrite(file_xls_S,pep_commons,Foglio_S,'N2');
xlswrite(file_xls_S,{'Common peptides:' num2str(numel(pep_commons))},Foglio_S,'N1');
xlswrite(file_xls_S,{'ASAPRAtio quantifications'},Foglio_S,'R1');
xlswrite(file_xls_S,pep_quant_j,Foglio_S,'R2');

xlswrite('correction',{'# HOLD'},file_xls_S,'A2');
xlswrite('correction',{'in3D'},file_xls_S,'B1');
xlswrite('correction',{'in2D'},file_xls_S,'C1');
xlswrite('correction',{'# XOR'},file_xls_S,'A3');
xlswrite('correction',numel(in3D),file_xls_S,'B2');
xlswrite('correction',numel(in2D),file_xls_S,'C2');
xlswrite('correction',numel(removed),file_xls_S,'B3');
xlswrite('correction',{'Correction %'},file_xls_S,'A4');
xlswrite('correction',(numel(in3D)/numel(removed))*100,file_xls_S,'B4');
xlswrite('correction',(numel(in2D)/numel(removed))*100,file_xls_S,'C4');
xlswrite('correction',{'% plus'},file_xls_S,'A6');
xlswrite('correction',(numel(plus3D)-numel(plus2D))/numel(pep_quant_j)*100,file_xls_S,'B6');
xlswrite('correction',[{'mean'};{'std'}],file_xls_S,'A7');
xlswrite('correction',[mean_s_PRE;std_s_PRE],file_xls_S,'B7');
xlswrite('correction',[mean_j_PRE;std_j_PRE],file_xls_S,'C7');
xlswrite('correction',{'# plus'},file_xls_S,'A5');
xlswrite('correction',numel(plus3D),file_xls_S,'B5');
xlswrite('correction',numel(plus2D),file_xls_S,'C5');

if outliers_removal
    stats=[mean(quantj_all);std(quantj_all);(std(quantj_all)/mean(quantj_all))*100;numel(quantj_all);numel(quantj_all)/numel(quantj_all_PRE)];
else
    stats=[mean(quantj_all);std(quantj_all);(std(quantj_all)/mean(quantj_all))*100;numel(quantj_all);100];
end

xlswrite(file_xls_S,{'ASAPRatio'},Foglio_S,'X1');
xlswrite(file_xls_S,stats,Foglio_S,'X2:X6');

%% output
commons.commons=pep_commons;
commons.pep_quant_S=pep_quant;
commons.quant_S=quants;
commons.light_S=lights;
commons.heavy_S=heavys;
commons.pep_quant_J=pep_quant_j;
commons.quant_J=quantj;
commons.quant_J_all=quantj_all;
commons.quant_S_all=quants_all;
commons.light_J=lightj;
commons.heavy_J=heavyj;
commons.light_J_all=lightj_all;
commons.heavy_J_all=heavyj_all;
commons.light_S_all=lights_all;
commons.heavy_S_all=heavys_all;
	
