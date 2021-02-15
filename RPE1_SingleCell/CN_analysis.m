addpath /czlab/chzhang/CodeBase/matlab/hg38
load hg38.bins.mat
load hg38.mat
load ChrArms.mat

addpath /singlecellcenter/RPE-1/CN_analysis/ 

ReadDepthNormalization;
AllelicDepthProcessing;

%load SampleInfo.mat
%load SI_binnedCov.mat
%load SI_BinnedAllelicRatio.mat 

CalculateAllelicCN;
if ~exist('GenomeWideCovPlots')
	mkdir('GenomeWideCovPlots');
end
cd('GenomeWideCovPlots');
COV=AllelicCN_1Mb;
for i=1:length(Samples) 
    try
		%Cov=COV{i};
		Cov=AllelicCN_1Mb{i};
		GenomeWideCovPlot(Cov);
		title(Samples{i},'FontSize',24,'Interpreter','None')
   		print(gcf,'-dpdf',[Samples{i},'_allelicCov_1Mb.pdf']);
    	close all;
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end
cd ..
