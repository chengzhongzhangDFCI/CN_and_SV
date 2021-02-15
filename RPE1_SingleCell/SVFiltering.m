BKPS=[INTER_BKPS;INTRA_BKPS];
bkps=BKPS(BKPS.TotalCount-BKPS.TCount<=1 & BKPS.maq>=30 & BKPS.NCount==0 & BKPS.TCount>=2,:);
bkps=sortrows(bkps,{'Sample','chr1','pos1'});

%save BKPS.mat SampleInfo bkps BKPS INTRA_SUPP_READS INTER_SUPP_READS %FB_SUPP_READS TD_SUPP_READS DL_SUPP_READS SRBKPS srbkps 
SampleNames=SampleInfo.SampleNames;
for i=1:length(SampleNames)
	Sample=SampleNames{i};	
    SV=bkps(strncmp(bkps.Sample,Sample,length(Sample)),:);
    sampleFields=find(strncmp(fieldnames(SV),Sample,length(Sample)))';	
	SV=SV(:,[1:11,13,14,sampleFields]);
	SV=SV(SV.TotalCount>1 | (SV.TotalCount>2 & ~strcmp(SV.chr1,SV.chr2)),:);	
	SV=SV(SV.(SV.Properties.VarNames{end})>=4,:);
	SampleSV_Depth4{i}=sortrows(SV,{'chr1','pos1'});
	SV=SV(SV.(SV.Properties.VarNames{end})>=6,:);
	SampleSV_Depth6{i}=sortrows(SV,{'chr1','pos1'});
	SV=SV(SV.(SV.Properties.VarNames{end})>=8,:);
	SampleSV_Depth8{i}=sortrows(SV,{'chr1','pos1'});
	SV=SV(SV.(SV.Properties.VarNames{end})>=10,:);
	SampleSV{i}=sortrows(SV,{'chr1','pos1'});
end
save SampleSV.mat SampleSV* SampleNames

