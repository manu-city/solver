function [nonDimParams] = non_dim(dimParams, L, V)
%% Non-dimensionalisation of the domain, square domain is assumed. Units of inputs is in meters.

% L = domain length
% H = domain height
% o_x = obstacle x length
% o_y = obstacle y length
% clear = clearance between obstacles
% scale = length scale 

%% 
parameters = fieldnames(dimParams);

for i = 1:numel(parameters)
    name = parameters{i};
    name_ = append(name, '_');

    % Check for velocity data, requires velocity data to be names with
    % 'vel' in its name.
    if contains(name, 'vel') == 1
        nonDimParams.(name_) = dimParams.(name)/V;
    
    % Check for number data, does not need to be non-dimensionalised. Must
    % contain 'num' in its name
    elseif contains(name, 'num') == 1
        continue
      
    % All other data are lengths
    else
        nonDimParams.(name_) = dimParams.(name)/L;
    end
    
end

% Obtain Reynolds Number

nonDimParams.reynolds_ = (dimParams.bulkvel * dimParams.o_y)/dimParams.v;
