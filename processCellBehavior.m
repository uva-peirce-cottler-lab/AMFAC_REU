%%This script is used by the Repast model to process cell networks. It
%%should therefore be in the Eclipse workspace
%Tommy Athey Aug 2017

%runs the loop in parallel using parfor
%the network mat file should have been loaded already
%needs the variable "states" of the network states of the cells to process
parfor i=1:length(states(:,1)) %each row of the states variable is an incoming network state
    p = params;
    in = states(i,:); %set the initial values
    TGFB = TGFBweights(i); %TGFB weight for this fibroblast
    IL1 = IL1weights(i);
    IL6 = IL6weights(i);
    TNFa = TNFaweights(i);
    %set the weights (or steady state values) of the inputs to the initial
    %values as well
    %p{1}(1,1:13) = [in(1) in(20) in(32) in(39) in(42) in(44) in(12) in(26) in(7) in(28) in(14) 0 0];
    
    p{1}(1,1:11) = [0.25 TGFB 0.25 IL6 IL1 TNFa 0.25 0.25 0.25 0.25 0.25];
    p{1}(1,12:13) = [0 0]; %turns off latent TGF-B feedback
    %p{1}(1,15:17) = [0 0 0]; %turns off IL-6 feedback
    
    [t,y]=ode23(@ODE,tspan,in,options,p,ODElist); %simulate a time step of the network
    states(i,:) = y(length(y(:,1)),:); %change that network state to its new value
end
%the states variable can now be read elsewhere for the updated network
%states