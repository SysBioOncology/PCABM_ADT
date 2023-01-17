function visualizeSystem(TUprop,mySystem)

cellfun(@(f) evalin('caller',[f ' = TUprop.' f ';']), fieldnames(TUprop));
cellfun(@(f) evalin('caller',[f ' = mySystem.params.' f ';']), fieldnames(mySystem.params));
cellfun(@(f) evalin('caller',[f ' = mySystem.grid.' f ';']), fieldnames(mySystem.grid));
cellfun(@(f) evalin('caller',[f ' = mySystem.TU.' f ';']), fieldnames(mySystem.TU));
cellfun(@(f) evalin('caller',[f ' = mySystem.F.' f ';']), fieldnames(mySystem.F));
cellfun(@(f) evalin('caller',[f ' = mySystem.M1.' f ';']), fieldnames(mySystem.M1));
cellfun(@(f) evalin('caller',[f ' = mySystem.M2.' f ';']), fieldnames(mySystem.M2));

F = ones(N,M,3);

%PLOT NORMAL TUMOR CELLS
TUcellnormal = TUcells(~isRes);
F(TUcellnormal) = 1;
F(TUcellnormal+(N*M)) = 0;
F(TUcellnormal+(N*M)*2) = 0;

%PLOT RESISTANT TUMOR CELLS
TUcellres = TUcells(isRes);
F(TUcellres) = 0;
F(TUcellres+(N*M))=0;
F(TUcellres+(N*M)*2)=0;

%PLOT FIBROBLASTS
F(Fcells) = 0;
F(Fcells+(N*M)) = 204/255;
F(Fcells+(N*M)*2) = 0;

%PLOT M1
F(M1cells) = 0;
F(M1cells+(N*M)) = 128/255;
F(M1cells+(N*M)*2) = 255/255;

%PLOT M2
F(M2cells) = 51/255;
F(M2cells+(N*M)) = 255/255;
F(M2cells+(N*M)*2) = 255/255;

%PLOT FIGURE PROPERTIES
imshow(F,'InitialMagnification',600)
beginSc=100;
text(beginSc,beginSc+15,[num2str(round(mySystem.grid.StepsDone*4)),' hours'],...
    'Color','k','FontWeight','bold','FontSize',18,'VerticalAlignment','top')

end


