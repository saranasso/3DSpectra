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
function [h z_value]=corr_diff(r1,n1,r2,n2)
% author: Sara Nasso
% it finds statistical significative difference (.05 level) between 2 correlation
% coefficients on two independent samples (check it through corrcoef)

%% Fisher's z-score transformation of Pearson's r. Fisher's transformation
% reduces skew and makes the sampling distribution more normal as sample size increases.
% Z = ln(abs((r+1)/r-1)))/2
%% Table of Z-score conversions for Pearson's r (r_to_z_score)
%    r             z'
r_to_z_score=[
0.0000        0.0000
0.0100        0.0100
0.0200        0.0200
0.0300        0.0300
0.0400        0.0400
0.0500        0.0500
0.0600        0.0601
0.0700        0.0701
0.0800        0.0802
0.0900        0.0902
0.1000        0.1003
0.1100        0.1104
0.1200        0.1206
0.1300        0.1307
0.1400        0.1409
0.1500        0.1511
0.1600        0.1614
0.1700        0.1717
0.1800        0.1820
0.1900        0.1923
0.2000        0.2027
0.2100        0.2132
0.2200        0.2237
0.2300        0.2342
0.2400        0.2448
0.2500        0.2554
0.2600        0.2661
0.2700        0.2769
0.2800        0.2877
0.2900        0.2986
0.3000        0.3095
0.3100        0.3205
0.3200        0.3316
0.3300        0.3428
0.3400        0.3541
0.3500        0.3654
0.3600        0.3769
0.3700        0.3884
0.3800        0.4001
0.3900        0.4118
0.4000        0.4236
0.4100        0.4356
0.4200        0.4477
0.4300        0.4599
0.4400        0.4722
0.4500        0.4847
0.4600        0.4973
0.4700        0.5101
0.4800        0.5230
0.4900        0.5361
0.5000        0.5493
0.5100        0.5627
0.5200        0.5763
0.5300        0.5901
0.5400        0.6042
0.5500        0.6184
0.5600        0.6328
0.5700        0.6475
0.5800        0.6625
0.5900        0.6777
0.6000        0.6931
0.6100        0.7089
0.6200        0.7250
0.6300        0.7414
0.6400        0.7582
0.6500        0.7753
0.6600        0.7928
0.6700        0.8107
0.6800        0.8291
0.6900        0.8480
0.7000        0.8673
0.7100        0.8872
0.7200        0.9076
0.7300        0.9287
0.7400        0.9505
0.7500        0.9730
0.7600        0.9962
0.7700        1.0203
0.7800        1.0454
0.7900        1.0714
0.8000        1.0986
0.8100        1.1270
0.8200        1.1568
0.8300        1.1881
0.8400        1.2212
0.8500        1.2562
0.8600        1.2933
0.8700        1.3331
0.8800        1.3758
0.8900        1.4219
0.9000        1.4722
0.9100        1.5275
0.9200        1.5890
0.9300        1.6584
0.9400        1.7380
0.9500        1.8318
0.9600        1.9459
0.9700        2.0923
0.9800        2.2976
0.9900        2.6467];

%% Significance of the difference between two correlations from two independent samples
% To compute the significance of the difference between two correlations from independent samples,
% such as a correlation for males vs. a correlation for females, follow these steps:

%% Use the table of z-score conversions or convert the two correlations to z-scores, as outlined above or the formula above. 
% Note that if the correlation is negative (we are considering only positive correlations!), the z value should be negative.
% z1=r_to_z_score(find(r_to_z_score(:,1)==r1),2);
z1 = log(abs((r1+1)/(r1-1)))/2;
% z2=r_to_z_score(find(r_to_z_score(:,1)==r2),2);
z2 = log(abs((r2+1)/(r2-1)))/2;
diff_z1_z2=abs(z1-z2);
%% Estimate the standard error of difference between the two correlations as:
% SE = SQRT[(1/(n1 - 3) + 1/(n2 - 3)]
% where n1 and n2 are the sample sizes of the two independent samples
SE = sqrt(1/(n1 - 3) + 1/(n2 - 3));
 
%% Divide the difference between the two z-scores by the standard error.
z_value=diff_z1_z2/SE;

%% If the z value for the difference computed in step 3 is 1.96 or higher, the difference in the correlations is significant at the .05 level. Use a 2.58 cutoff for significance at the .01 level.

if z_value >= 1.96 % soglia significatività al 5%   
    h=1;   
else   
    h=0;
end
% Example. Let a sample of 15 males have a correlation of income and education of .60, and let a sample of 20 females
% have a correlation of .50. We wish to test if this is a significant difference.
% The z-score conversions of the two correlations are .6931 and .5493 respectively, for a difference of .1438. 
% The SE estimate is SQRT[(1/12)+(1/17)] = SQRT[.1422] = .3770. 
% The z value of the difference is therefore .1438/.3770 = .381, much smaller than 1.96 and thus not significant at the .05 level. 
% (Source, Hubert Blalock, Social Statistics, NY: McGraw-Hill, 1972: 406-407.)	
