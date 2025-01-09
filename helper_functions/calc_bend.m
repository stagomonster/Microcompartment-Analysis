function bend = calc_bend(rubisco1, rubisco2)
% This function calculates the amount of bend present between adjacent
% rubiscos in a chain.
%
% Inputs
% rubisco1 - a rubisco object paired with rubisco2
% rubisco2 - a rubisco object paired with rubisco1
%
% Outputs
% bend - the angle between a rubisco's orientation and the orientation
% of the next rubisco in the chain in degrees
%
% calc_bend.m © 2025 is licensed under CC BY-NC-SA 4.0

    % Normalize both rubisco orientation vectors
    rubisco_vec_1 = rubisco1.vector / norm(rubisco1.vector);
    rubisco_vec_2 = rubisco2.vector / norm(rubisco2.vector);

    tmpBend = dot(rubisco_vec_1, rubisco_vec_2); % calculate dot product of vectors
    if tmpBend > 0 % take dot product if they point in the same direction
        bend = tmpBend;
    else % else flip one vector and use that dot product
        bend = dot(rubisco_vec_1, -1 * rubisco_vec_2);
    end

    % Return bend angle in degrees between 0 and 90
    bend = acos(bend) * 180/pi;
end