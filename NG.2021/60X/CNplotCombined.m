load SI_allelicCN_1Mb.mat ;
COV=AllelicCN_1Mb;
%load SampleSV.mat
%RA=SampleSV_Depth4;
load SampleInfo.mat;
SampleGroups=unique(SampleInfo.SampleGroupID);
for i=1:length(SampleGroups)
%	try
		Name=SampleGroups{i};
		id=find(strcmp(SampleInfo.SampleGroupID,Name));
		Cov=COV{id(1)};
		for j=2:length(id)
			cov1=COV{id(j)};
			Cov.alleleA=Cov.alleleA+cov1.alleleA;
			Cov.alleleB=Cov.alleleB+cov1.alleleB;
		end
		Cov.alleleA=Cov.alleleA/length(id);
		Cov.alleleB=Cov.alleleB/length(id);
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
%		SV=RA{i};
%		SV=SV(SV.TotalCount-SV.SplitCount>0 & (SV.SplitCount>0 | strcmp(SV.chr1,SV.chr2)),:);
%		lrbkps=SV(strcmp(SV.chr1,SV.chr2) & SV.pos2-SV.pos1>=1e6,:);
%		itbkps=SV(~strcmp(SV.chr1,SV.chr2),:);
%		lrbkps.pos1=double(lrbkps.pos1)/1e6;lrbkps.pos2=double(lrbkps.pos2)/1e6;
%		itbkps.pos1=double(itbkps.pos1)/1e6;itbkps.pos2=double(itbkps.pos2)/1e6;

		chrlist=unique(Cov.chr);
	%	chrlist={'chr4'};
		for chri=1:length(chrlist)
			chrlabel=chrlist{chri};
			idx=strcmp(Cov.chr,chrlabel);
			plotPos=Cov.pos(idx);
			plotAlleleA=Cov.alleleA(idx);
			plotAlleleB=Cov.alleleB(idx);
			ChrCovPlotAllelic(chrlabel,plotPos,plotAlleleA,plotAlleleB);
%			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1~=lrbkps.str2,:),4,1,'k');
%			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1==lrbkps.str2,:),4,-1,'k');
%			bkps1=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==-1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==-1)];
%			bkps2=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==1)];
%			for ki=1:length(bkps1)
%				plot(bkps1(ki)*[1,1],[4,5],'m');
%			end
%			for ki=1:length(bkps2)
%				plot(bkps2(ki)*[1,1],[4,3],'m');
%			end
			%title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
			print(gcf,'-dpdf',[Name '_' chrlabel '_1Mb.pdf']);
			close all
		end
%	catch
%		fprintf(1,'Failed at sample %s\n',Samples{i});
%	end
end

load SI_allelicCN_250kb.mat ;
COV=AllelicCN_250kb;
%load SampleSV.mat
%RA=SampleSV_Depth4;
Samples=SampleInfo.SampleNames;
for i=1:length(Samples)
	try
		Name=Samples{i};
		Cov=COV{i};
		Cov=[Cov,dataset(Cov.alleleA+Cov.alleleB,'VarNames','alleleCN')];
%		SV=RA{i};
%		SV=SV(SV.TotalCount-SV.SplitCount>0 & (SV.SplitCount>0 | strcmp(SV.chr1,SV.chr2)),:);
%		lrbkps=SV(strcmp(SV.chr1,SV.chr2) & SV.pos2-SV.pos1>=1e6,:);
%		itbkps=SV(~strcmp(SV.chr1,SV.chr2),:);
	%	export2circos([lrbkps;itbkps],Cov,[Samples{i}]);	
	%	export2circos([lrbkps;itbkps],[],[Samples{i}]);	
	%	lrbkps.pos1=double(lrbkps.pos1)/1e6;lrbkps.pos2=double(lrbkps.pos2)/1e6;
	%	itbkps.pos1=double(itbkps.pos1)/1e6;itbkps.pos2=double(itbkps.pos2)/1e6;

		chrlist=unique(Cov.chr);
		for chri=1:length(chrlist)
			chrlabel=chrlist{chri};
			idx=strcmp(Cov.chr,chrlabel);
			plotPos=Cov.pos(idx);
			plotAlleleA=Cov.alleleA(idx);
			plotAlleleB=Cov.alleleB(idx);
			ChrCovPlotAllelic(chrlabel,plotPos,plotAlleleA,plotAlleleB);
%			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1~=lrbkps.str2,:),4,1,'k');
%			plot_RA(lrbkps(strcmp(lrbkps.chr1,chrlabel) & lrbkps.str1==lrbkps.str2,:),4,-1,'k');
%			bkps1=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==-1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==-1)];
%			bkps2=[itbkps.pos1(strcmp(itbkps.chr1,chrlabel) & itbkps.str1==1);itbkps.pos2(strcmp(itbkps.chr2,chrlabel) & itbkps.str2==1)];
%			for ki=1:length(bkps1)
%				plot(bkps1(ki)*[1,1],[4,5],'m');
%			end
%			for ki=1:length(bkps2)
%				plot(bkps2(ki)*[1,1],[4,3],'m');
%			end
			%title([SampleNames{i} '_' chrlabel],'FontSize',18,'Interpreter','None');
			print(gcf,'-dpdf',[Samples{i} '_' chrlabel '_250kb.pdf']);
			close all
		end
	catch
		fprintf(1,'Failed at sample %s\n',Samples{i});
	end
end


