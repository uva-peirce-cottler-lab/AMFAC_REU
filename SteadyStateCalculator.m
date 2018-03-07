%%Calculates steady state network states for all grid elements when there
%%is a horizontal linear IL1, IL6, TNFa gradient and a vertical linear TGFb
%%gradient. Shouold be in the Eclipse workspace so it has access to
%%network.met
%Tommy Athey Aug 2017
tic
gridWidth = 10;
gridHeight = 10;
netSize = 91;
load network.mat
tspan = [0 168]; %hopefully enough to reach steady state
ssnetwork = 'C:\Users\smr2we\Documents\Fibrosis Project\Steady State Network\10x10_noTGFB_IL6_feedback.csv';
states = zeros(gridWidth*gridHeight,netSize); %all networks will start at 0
for x=1:gridWidth %for every grid
    for y=1:gridHeight
        in = zeros(1,netSize);
        %all inputs will be 0.25 except for the inflammatory/fibrotic
        %cytokines, just like in the model
        params{1}(1,1:11) = [0.25 (y-1)/gridHeight 0.25 (x-1)/gridWidth (x-1)/gridWidth (x-1)/gridWidth 0.25 0.25 0.25 0.25 0.25];
        params{1}(1,12:13) = [0 0]; %turns off latent TGF-B feedback
        %params{1}(1,15:17) = [0 0 0]; %turns off IL-6 feedback
        [~,v] = ode23(@ODE,tspan,in,options,params,ODElist);
        states((x-1)*gridWidth+y,:) = v(end,:); %retrieve final network state
        disp((gridWidth*(x-1)+y)/(gridWidth*gridHeight)); %display progress
    end
end

%csvwrite(ssnetwork,states);
toc
