function narot = calc_narot(rubisco, chain_lead_vector)
% This function gets the rotation angle of a rubisco in a chain
%
% Inputs
% rubisco - a rubisco class object in a chain
% chain_lead_vector - the average vector of the chain containing rubisco
%
% Outputs
% narot - the rotation angle of the rubisco

    % if the rubisco orientation vector and the chain_lead_vector point in
    % the same direction
    if dot(chain_lead_vector, rubisco.vector) > 0
        narot = rubisco.narot; % use the rubisco's original narot
    else
        narot = -1 * rubisco.narot; % otherwise use the negative of the rubisco's narot
    end
end