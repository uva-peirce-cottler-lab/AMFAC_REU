for i=1:length(states(:,1)) %each row of the states variable is an incoming network state
    in = states(i,:); %set the initial values
    %set the weights (or steady state values) of the inputs to the initial
    %values as well
    
    params{1}(1,1:11) = [0.25 in(20) 0.25 in(39) in(42) in(44) 0.25 0.25 0.25 0.25 0.25];
    params{1}(1,12:13) = [0 0]; %turns off latent TGF-B feedback
    params{1}(1,15:17) = [0 0 0]; %turns off IL-6 feedback
end