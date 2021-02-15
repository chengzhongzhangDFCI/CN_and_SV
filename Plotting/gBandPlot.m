function gBandPlot(chrlabel,loc,height)

load hg38.mat
gBands=GiemsaBands(strcmp(GiemsaBands.Chr,chrlabel),:);
gBands.Start=double(gBands.Start)/1e6;gBands.End=double(gBands.End)/1e6;
gBandColor=gBands.Stain;
gBandColor(strcmp(gBandColor,'acen'))={[0.9,0,0]};
gBandColor(strcmp(gBandColor,'gpos25'))={[0.75,0.75,0.75]};
gBandColor(strcmp(gBandColor,'gpos50'))={[0.5,0.5,0.5]};
gBandColor(strcmp(gBandColor,'gpos75'))={[0.25,0.25,0.25]};
gBandColor(strcmp(gBandColor,'gpos100'))={[0,0,0]};
gBandColor(strcmp(gBandColor,'stalk'))={[0,0.5,0]};
gBandColor(strcmp(gBandColor,'gvar'))={[0.5,0.25,0]};
gBandColor(strcmp(gBandColor,'gneg'))={[1,1,1]-0.01};
rectangle('Position',[gBands.Start(1),loc,gBands.End(end)-gBands.Start(1),height],'LineWidth',1.5,'FaceColor','none','EdgeColor',[0,0,0]);
for gi=1:length(gBands)
	rectangle('Position',[gBands.Start(gi),loc,gBands.End(gi)-gBands.Start(gi),height],'FaceColor',gBandColor{gi},'EdgeColor','none');
end
