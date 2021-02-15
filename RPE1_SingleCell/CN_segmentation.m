%addpath /czlab/chzhang/CodeBase/matlab/hg38
load hg38.mat
load ChrArms.mat
CHRs=ChrLabels(1:23);
load hg38.bins.mat
Bins_10kb=hg38_10kb_bins(1:find(strcmp(hg38_10kb_bins.Chr,'chrX'),1,'last'),:);

% Allelic coverage (single site and 10kb binned average)
load SI_allelicCN_250kb.mat

% Segmenting allelic coverage
ploidy=1;
MAbins=4; % 1Mb bins
minCNLen=6; % 1.5Mb minimum CN
AllelicSeg={};
AllelicSeg2={};
for i=1:length(Samples)	
	A_SEG=[];
	A_SEG2=[];
	B_SEG=[];
	B_SEG2=[];
	for ci=1:length(CHRs)
		chr=CHRs(ci);	
		AllelicDepth=AllelicCN_250kb{i}(strcmp(AllelicCN_250kb{i}.chr,chr),:);
		centromere_gap=Centromeres(strcmp(Centromeres.Chr,chr),:);
		if (centromere_gap.Start>0) % p-arm for non-acrocentric chromosomes
			binPos=AllelicDepth.pos(AllelicDepth.pos<=centromere_gap.Start)-120000;
			% alleleA
			meanDepth=AllelicDepth.alleleA(AllelicDepth.pos<=centromere_gap.Start,:);
			[A_seg,A_cn]=SegmentCN(meanDepth,ploidy,MAbins,minCNLen);
			A_seg.Start=binPos(A_seg.Start);A_seg.End=binPos(A_seg.End);
			A_seg.End(end)=centromere_gap.Start;
			A_SEG=[A_SEG;dataset(repmat(chr,size(A_seg.Start)),'VarNames',{'Chr'}),A_seg];
			[A_seg,A_cn]=SegmentCN(meanDepth,2*ploidy,MAbins,2*minCNLen);
			A_seg.Start=binPos(A_seg.Start);A_seg.End=binPos(A_seg.End);
			A_seg.End(end)=centromere_gap.Start;
			A_SEG2=[A_SEG2;dataset(repmat(chr,size(A_seg.Start)),'VarNames',{'Chr'}),A_seg];
			% alleleB
			meanDepth=AllelicDepth.alleleB(AllelicDepth.pos<=centromere_gap.Start,:);
			[B_seg,B_cn]=SegmentCN(meanDepth,ploidy,MAbins,minCNLen);
			B_seg.Start=binPos(B_seg.Start);B_seg.End=binPos(B_seg.End);
			B_seg.End(end)=centromere_gap.Start;		
			B_SEG=[B_SEG;dataset(repmat(chr,size(B_seg.Start)),'VarNames',{'Chr'}),B_seg];	
			[B_seg,B_cn]=SegmentCN(meanDepth,2*ploidy,MAbins,2*minCNLen);
			B_seg.Start=binPos(B_seg.Start);B_seg.End=binPos(B_seg.End);
			B_seg.End(end)=centromere_gap.Start;
			B_SEG2=[B_SEG2;dataset(repmat(chr,size(B_seg.Start)),'VarNames',{'Chr'}),B_seg];
		end
		% q-arm
		binPos=AllelicDepth.pos(AllelicDepth.pos>=centromere_gap.End)-120000;
		% alleleA
		meanDepth=AllelicDepth.alleleA(AllelicDepth.pos>=centromere_gap.End);
		[A_seg,A_cn]=SegmentCN(meanDepth,ploidy,MAbins,minCNLen);
		A_seg.Start=binPos(A_seg.Start);A_seg.End=binPos(A_seg.End);
		A_seg.Start(1)=centromere_gap.End;
		A_SEG=[A_SEG;dataset(repmat(chr,size(A_seg.Start)),'VarNames',{'Chr'}),A_seg];
		[A_seg,A_cn]=SegmentCN(meanDepth,2*ploidy,MAbins,2*minCNLen);
		A_seg.Start=binPos(A_seg.Start);A_seg.End=binPos(A_seg.End);
		A_seg.Start(1)=centromere_gap.End;
		A_SEG2=[A_SEG2;dataset(repmat(chr,size(A_seg.Start)),'VarNames',{'Chr'}),A_seg];
		% alleleB
		meanDepth=AllelicDepth.alleleB(AllelicDepth.pos>=centromere_gap.End,:);
		[B_seg,B_cn]=SegmentCN(meanDepth,ploidy,MAbins,minCNLen);
		B_seg.Start=binPos(B_seg.Start);B_seg.End=binPos(B_seg.End);
		B_seg.Start(1)=centromere_gap.End;		
		B_SEG=[B_SEG;dataset(repmat(chr,size(B_seg.Start)),'VarNames',{'Chr'}),B_seg];	
		[B_seg,B_cn]=SegmentCN(meanDepth,2*ploidy,MAbins,2*minCNLen);
		B_seg.Start=binPos(B_seg.Start);B_seg.End=binPos(B_seg.End);
		B_seg.Start(1)=centromere_gap.End;
		B_SEG2=[B_SEG2;dataset(repmat(chr,size(B_seg.Start)),'VarNames',{'Chr'}),B_seg];
	end
	SEG.A=A_SEG;
	SEG.B=B_SEG;
	AllelicSeg{i}=SEG;
	SEG.A=A_SEG2;
	SEG.B=B_SEG2;
	AllelicSeg2{i}=SEG;
end
save AllelicCNSeg.mat AllelicSeg AllelicSeg2;

return
