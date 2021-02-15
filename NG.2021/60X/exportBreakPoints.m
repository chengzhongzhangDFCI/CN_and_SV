load SampleInfo.mat
load SampleSVs.mat
for i=1:length(SampleInfo)
	Sample=SampleInfo.SampleNames{i};
	SV=SampleSVs{i};
	Supp=SampleSVSupp{i}(:,1:7);
	SV=SV(SV.suppSamples==1 & SV.NCount==0,1:13);
	SVsupp=[];
	for j=1:length(SV)
		SVsupp=[SVsupp;Supp(Supp.SVidx==SV.SVidx(j),:)];
	end
	export(SV,'File',[Sample '_bkpPairs.txt'],'delimiter','\t');
	export(SVsupp,'File',[Sample '_suppReads.txt'],'delimiter','\t');
end
