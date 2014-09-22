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
function statistics=getRegression(heavy,light)
% it performs linear regression if there is a significant linear
% relationship at 5% and if the Pearson's coefficient is not below 0.7
% (minimum level for good enough linearity) 

axis([0,max(heavy),0,max(light)])
axis manual
axis(axis)
hold on
plot(heavy,light,'o')
xlabel('Heavy abundances')
ylabel('Light abundances')
x_SD=sqrt(heavy);
y_SD=sqrt(light);
N_dati=numel(light);
[r p] = corr(heavy,light,'tail','gt','rows','pairwise');
statistics.R2=r^2;
statistics.pValue_R2=p;

if p<0.05
    [m,m_SE,q,q_SE,min_quad_distXY,residual,covar_post,J,SE]=fitlineXY(heavy,x_SD,light,y_SD,0.1,0);
    hold on
    x=(0:0.1:max(heavy));
    plot(x,m*x+q,'r')
    legend({['Quantifications #=' num2str(N_dati) ' R2xy=' num2str(r^2,'%.2f') ' pValue=' num2str(p)],['Regression line: m=' num2str(m,'%.2f') ' q=' num2str(q,'%.2f')  ' SE(m)=' num2str(m_SE,'%.2f') ' SE(q)=' num2str(q_SE,'%.2f') ' RMSE=' num2str(SE,'%.2f')]})
    statistics.RMSE=SE;
    statistics.m_SE=m_SE;
    statistics.q_SE=q_SE;
elseif p>0.05
    legend({['Quantifications #=' num2str(N_dati) ' R2xy=' num2str(r^2,'%.2f') ' pValue=' num2str(p)],'Linear relationship is not significant'})
end
	
