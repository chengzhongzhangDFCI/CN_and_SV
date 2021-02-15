%% Individual sample CN plot
addpath /singlecellcenter/RPE-1/CN_analysis/ 
addpath /czlab/chzhang/CodeBase/matlab/hg38

load hg38.mat
load hg38.bins.mat
load ../SI_binnedCov.mat
Bins_10kb=hg38_10kb_bins(1:find(strcmp(hg38_10kb_bins.Chr,'chrX'),1,'last'),:);
[Bins_10kb,sortIdx]=sortrows(Bins_10kb,{'Chr','Start'});

DepthSampleNames=load('../SI_SampleNames.mat');

%load('SNP_binning_10kb.mat','Bins_10kb');
load SI_allelic_SampleNames.mat

% 1Mb
plotBinSize=1e6;
load('SI_BinnedAllelicRatio.mat','AlleleARatio_1Mb','SNPCount_1Mb');

%% 30x samples still use depth of coverage from 5x samples for better normalization of recurrent bias

DepthSamples=[DepthSampleNames.SI_5x_SampleNames,DepthSampleNames.SI_10x_SampleNames,DepthSampleNames.SI_30x_SampleNames];
SampleCov=[SI_5x_1Mb,SI_10x_1Mb,SI_30x_1Mb];
SampleCov=SampleCov(sortIdx,:);
SampleAlleleARatio=AlleleARatio_1Mb;
SNPCount=SNPCount_1Mb;
Samples=SampleNames;
Samples(SI_60x_idx)={'MN_SI_140609_P9_E3','MN_SI_140609_P9_E4','MN_SI_140609_P9_E5','MN_SI_140609_P9_E7','MN_SI_140609_P9_E8'};
AllelicCN_1Mb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	DepthSampleIdx=DepthSampleIdx(1);
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_1Mb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end
save SI_allelicCN_1Mb.mat AllelicCN_1Mb AllelicCN_Chr Samples SI_*_idx

% 250kb
plotBinSize=250000;
load('SI_BinnedAllelicRatio.mat','AlleleARatio_250kb','SNPCount_250kb');

%% 30x samples still use depth of coverage from 5x samples for better normalization of recurrent bias

%DepthSamples=[DepthSampleNames.SI_5x_SampleNames,DepthSampleNames.SI_10x_SampleNames,DepthSampleNames.SI_30x_SampleNames];
SampleCov=[SI_5x_250kb,SI_10x_250kb,SI_30x_250kb];
SampleCov=SampleCov(sortIdx,:);
SampleAlleleARatio=AlleleARatio_250kb;
SNPCount=SNPCount_250kb;
%Samples=SampleNames;
%Samples(SI_60x_idx)={'MN_SI_140609_P9_E3','MN_SI_140609_P9_E4','MN_SI_140609_P9_E5','MN_SI_140609_P9_E7','MN_SI_140609_P9_E8'};
AllelicCN_250kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	DepthSampleIdx=DepthSampleIdx(1);
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_250kb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end
save SI_allelicCN_250kb.mat AllelicCN_250kb AllelicCN_Chr Samples SI_*_idx

% 50kb
plotBinSize=50000;
load('SI_BinnedAllelicRatio.mat','AlleleARatio_50kb','SNPCount_50kb');

%% 30x samples still use depth of coverage from 5x samples for better normalization of recurrent bias

%DepthSamples=[DepthSampleNames.SI_5x_SampleNames,DepthSampleNames.SI_10x_SampleNames,DepthSampleNames.SI_30x_SampleNames];
SampleCov=[SI_5x_50kb,SI_10x_50kb,SI_30x_50kb];
SampleCov=SampleCov(sortIdx,:);
SampleAlleleARatio=AlleleARatio_50kb;
SNPCount=SNPCount_50kb;
%Samples=SampleNames;
%Samples(SI_60x_idx)={'MN_SI_140609_P9_E3','MN_SI_140609_P9_E4','MN_SI_140609_P9_E5','MN_SI_140609_P9_E7','MN_SI_140609_P9_E8'};
AllelicCN_50kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	DepthSampleIdx=DepthSampleIdx(1);
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_50kb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end

POS=AllelicCN_50kb{1}(:,1:2);
for i=1:length(AllelicCN_50kb)
	AllelicCN_50kb{i}=AllelicCN_50kb{i}(:,3:4);
end
save SI_allelicCN_50kb.mat AllelicCN_50kb POS  AllelicCN_Chr Samples SI_*_idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5-10x depth samples

DepthSamples=[DepthSampleNames.SI_5x_SampleNames,DepthSampleNames.SI_10x_SampleNames];
SampleCov=[SI_5x_1Mb(sortIdx,:),SI_10x_1Mb(sortIdx,:)];
SampleAlleleARatio=AlleleARatio_1Mb(:,SI_5x_idx);
SNPCount=SNPCount_1Mb(:,SI_5x_idx);
Samples=SampleNames(SI_5x_idx);
AllelicCN_1Mb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_1Mb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end
save SI_5x_allelicCN_1Mb.mat AllelicCN_1Mb Samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Genome-wide plots
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 250kb
plotBinSize=250000;
load('SI_BinnedAllelicRatio.mat','AlleleARatio_250kb','SNPCount_250kb');

%% 30x depth samples

DepthSamples=DepthSampleNames.SI_30x_SampleNames;
SampleCov=SI_30x_250kb(sortIdx,:);
SampleAlleleARatio=AlleleARatio_250kb(:,SI_30x_idx);
SNPCount=SNPCount_250kb(:,SI_30x_idx);
Samples=SampleNames(SI_30x_idx);
AllelicCN_250kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_250kb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end
save SI_30x_allelicCN_250kb.mat AllelicCN_250kb Samples

%% 5-10x depth samples

DepthSamples=[DepthSampleNames.SI_5x_SampleNames,DepthSampleNames.SI_10x_SampleNames];
SampleCov=[SI_5x_250kb(sortIdx,:),SI_10x_250kb(sortIdx,:)];
SampleAlleleARatio=AlleleARatio_250kb(:,SI_5x_idx);
SNPCount=SNPCount_250kb(:,SI_5x_idx);
Samples=SampleNames(SI_5x_idx);
AllelicCN_250kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_250kb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end
save SI_5x_allelicCN_250kb.mat AllelicCN_250kb Samples

%% 50kb
plotBinSize=0.5e5;
load('SI_BinnedAllelicRatio.mat','AlleleARatio_50kb','SNPCount_50kb');

%% 30x depth samples

DepthSamples=DepthSampleNames.SI_30x_SampleNames;
SampleCov=SI_30x_50kb(sortIdx,:);
SampleAlleleARatio=AlleleARatio_50kb(:,SI_30x_idx);
SNPCount=SNPCount_50kb(:,SI_30x_idx);
Samples=SampleNames(SI_30x_idx);
AllelicCN_50kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	DepthSampleIdx=find(strcmp(DepthSamples,SampleName));
	Depth=SampleCov(:,DepthSampleIdx);;
	AllelicRatio=SampleAlleleARatio(:,i);
	phasedSNPCounts=SNPCount(:,i);
	Cov=[Bins_10kb,dataset(Depth*2.*AllelicRatio,Depth*2.*(1-AllelicRatio),phasedSNPCounts,'VarNames',{'alleleA','alleleB','snp_count'})];
	WholeGenomeAllelicCN=[];
	ChrAllelicCN=[];
	for chr=1:23
		if chr<23
			chrlabel=['chr' int2str(chr)];
		else
			chrlabel='chrX';
		end
		cov_chr=Cov(strcmp(Cov.Chr,chrlabel),:);
		maxLength=cov_chr.End(end);
		maxScale=log10(double(maxLength)/1e6);
		XScale=maxScale*500;
		idx=round(plotBinSize/10000*(0.5:1:(maxLength/plotBinSize)));	
		snp_count=cov_chr.snp_count(idx);
		plotPos=double(cov_chr.Start(idx));
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);	
		plotAlleleA(snp_count<10)=nan;
		plotAlleleB(snp_count<10)=nan;
		ChrAllelicCN=[ChrAllelicCN;dataset({chrlabel},median(plotAlleleA,'omitnan'),median(plotAlleleB,'omitnan'),'VarNames',{'chr','alleleA','alleleB'})];
		WholeGenomeAllelicCN=[WholeGenomeAllelicCN;dataset(repmat({chrlabel},length(plotPos),1),plotPos,plotAlleleA,plotAlleleB,'VarNames',{'chr','pos','alleleA','alleleB'})];
	end	
	meanCN=median([ChrAllelicCN.alleleA;ChrAllelicCN.alleleB],'omitnan');
	WholeGenomeAllelicCN.alleleA=WholeGenomeAllelicCN.alleleA/meanCN;
	WholeGenomeAllelicCN.alleleB=WholeGenomeAllelicCN.alleleB/meanCN;
	ChrAllelicCN.alleleA=ChrAllelicCN.alleleA/meanCN;
	ChrAllelicCN.alleleB=ChrAllelicCN.alleleB/meanCN;
	AllelicCN_50kb{i}=WholeGenomeAllelicCN;
	AllelicCN_Chr{i}=ChrAllelicCN;
	fprintf(1,'done.\n');
catch
	fprintf(1,'Failed!\n');
end
end
save SI_30x_allelicCN_50kb.mat AllelicCN_50kb SampleNames

%% individual chromosome plot
for si=1:length(SampleNames)
	cov_chr=AllelicCN_1Mb{si};
	chrlist=unique(cov_chr.chr);
	for chri=1:length(chrlist)
		chrlabel=chrlist{chri};
		idx=strcmp(cov_chr.chr,chrlabel);
		plotPos=cov_chr.pos(idx);
		plotAlleleA=cov_chr.alleleA(idx);
		plotAlleleB=cov_chr.alleleB(idx);
		ChrCovPlotAllelic(chrlabel,plotPos,plotAlleleA,plotAlleleB);
		%title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
		print(gcf,'-dpdf',[SampleNames{si} '_' chrlabel '_1Mb.pdf']);
		close all
	end
end
