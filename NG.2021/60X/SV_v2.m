load SampleInfo.mat
addpath('/singlecellcenter/RPE-1/Analysis/scripts/');
loadDiscordantReads;
minCount=4;
ProcessingBySampleGroup_IntraBkps;
ProcessingBySampleGroup_InterBkps;
fprintf(1,'done\n');
SVCleaning;
SVlist=unique(AllSupp.SVidx);
SVs=[];SUPP=[];
for i=1:length(SVlist)
    supp=AllSupp(AllSupp.SVidx==SVlist(i),:);
    sv=[INTRA_BKPS(INTRA_BKPS.SVidx==SVlist(i),:);INTER_BKPS(INTER_BKPS.SVidx==SVlist(i),:)];
    SVs=[SVs;dataset(supp.SVGroup(1),'VarName','SVGroup'),sv];
    SUPP=[SUPP;supp];
end
SVs=sortrows(SVs,{'Group','Sample','SVGroup','SVidx'});

SampleSVs={};SampleSVSupp={};
for i=1:length(SampleInfo)
	Sample=SampleInfo.SampleNames(i);
	SampleSVs{i}=sortrows(SVs(strcmp(SVs.Sample,Sample),:),{'chr1','pos1'});
    SampleSVSupp{i}=SUPP(strcmp(SampleInfo.SampleNames(SUPP.SampleGroup),Sample),:);
end
save SampleSVs.mat SampleSVs SampleSVSupp;
