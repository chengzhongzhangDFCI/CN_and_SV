bam_folder='bamlinks';
bam_files=dir([bam_folder '/*.bam']);
SampleNames={};
SampleGroupNames={};
SampleInfo=[];
for i=1:length(bam_files)
	file=bam_files(i).name;
	SampleNames{i}=regexprep(file,'.bam','');
	info=split(SampleNames{i},'_');
	SampleInfo=[SampleInfo;dataset({[info{1} '_' info{2}]},info(3),info(4),info(5),'VarNames',{'Type','Date','Group','ID'})];
	SampleGroupNames{i}=[SampleInfo.Type{i} '_' SampleInfo.Date{i} '_' SampleInfo.Group{i}];
end
SampleInfo=[dataset(SampleNames',SampleGroupNames','VarNames',{'SampleNames','SampleGroupID'}),SampleInfo];
save SampleInfo.mat SampleInfo

addpath /singlecellcenter/RPE-1/Analysis/scripts
CN_analysis;
return
if ~exist('SampleSV.mat')
	SV_v2;
end
CN_segmentation;
addpath .
mkdir('IndividualChrPlots');
cd('IndividualChrPlots');
CN_plus_SV_plot;
cd('../CIRCOS');
exportCircos;
exit
