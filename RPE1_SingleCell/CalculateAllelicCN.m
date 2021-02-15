% 1Mb
plotBinSize=1e6;
SampleCov=SampleCov_1Mb;
Bins_10kb=hg38_10kb_bins(1:find(strcmp(hg38_10kb_bins.Chr,'chrX'),1,'last'),:);
[Bins_10kb,sortIdx]=sortrows(Bins_10kb,{'Chr','Start'});
SampleCov=SampleCov(sortIdx,:);
SampleAlleleARatio=AlleleARatio_1Mb;
SNPCount=SNPCount_1Mb;
Samples=SampleInfo.SampleNames;
AllelicCN_1Mb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	Depth=SampleCov(:,i);
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
save SI_allelicCN_1Mb.mat AllelicCN_1Mb AllelicCN_Chr Samples

% 250kb
plotBinSize=.25e6;
SampleCov=SampleCov_250kb;
Bins_10kb=hg38_10kb_bins(1:find(strcmp(hg38_10kb_bins.Chr,'chrX'),1,'last'),:);
[Bins_10kb,sortIdx]=sortrows(Bins_10kb,{'Chr','Start'});
SampleCov=SampleCov(sortIdx,:);
SampleAlleleARatio=AlleleARatio_250kb;
SNPCount=SNPCount_250kb;
Samples=SampleInfo.SampleNames;
AllelicCN_250kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	Depth=SampleCov(:,i);
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
        segA=SegmentCN(plotAlleleA,1,4,4);
        segB=SegmentCN(plotAlleleB,1,4,4);
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
save SI_allelicCN_250kb.mat AllelicCN_250kb AllelicCN_Chr Samples

% 50kb
plotBinSize=5e4;
SampleCov=SampleCov_50kb;
Bins_10kb=hg38_10kb_bins(1:find(strcmp(hg38_10kb_bins.Chr,'chrX'),1,'last'),:);
[Bins_10kb,sortIdx]=sortrows(Bins_10kb,{'Chr','Start'});
SampleCov=SampleCov(sortIdx,:);
SampleAlleleARatio=AlleleARatio_50kb;
SNPCount=SNPCount_50kb;
Samples=SampleInfo.SampleNames;
AllelicCN_50kb=cell(1,length(Samples));
AllelicCN_Chr=cell(1,length(Samples));
for i=1:length(Samples)
	SampleName=Samples{i};
	fprintf(1,'%s..',SampleName);
try
	Depth=SampleCov(:,i);
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
save SI_allelicCN_50kb.mat AllelicCN_50kb AllelicCN_Chr Samples

