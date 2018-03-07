resultsFolder = 'C:\Users\smr2we\Box Sync\Michaela''s Files - Fibrosis Project\Fibrosis Project\';
ssFile = [resultsFolder 'Steady State Network\5x5_TGFB_noIL6_feedback.csv'];

steadystate = real(csvread(ssFile));
gridHeight = 5;
gridWidth = 5;

mcolI = steadystate(:,88);
mcolIII = steadystate(:,89);

MMP1 = steadystate(:,82);
MMP2 = steadystate(:,83);
MMP9 = steadystate(:,84);

    for x=1:length(mcolI)
        avgdep(x) = (mcolI(x) + mcolIII(x))/2;
        avgdeg(x) = (MMP1(x) + MMP2(x) + MMP9(x))/3;
    end
    
    %initialize grids
    depImSS = zeros(gridHeight,gridWidth);
    degImSS = zeros(gridHeight,gridWidth);
    
    position = 1;
    for y=1:gridHeight
        for x=1:gridWidth
            depImSS((gridHeight + 1) - x,y) = avgdep(position);
            degImSS((gridHeight + 1) - x,y) = avgdeg(position);
            position = position + 1;
        end
    end
    
    figure
    imagesc(depImSS);
    colormap gray;
    caxis([0 1]);
    colorbar;
    axis off;
    title('Deposition','FontSize',16);
    
    figure
    imagesc(degImSS);
    colormap gray;
    caxis([0 1]);
    colorbar;
    axis off;
    title('Degradation','FontSize',16);