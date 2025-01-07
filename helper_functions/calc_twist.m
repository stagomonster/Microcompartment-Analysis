function twist = calc_twist(rubisco1, rubisco2, chain_lead_vec)
% This function calculates the twist between two consecutive rubiscos in a
% chain.
%
% Inputs
% rubisco1 - a rubisco class object in a chain right before rubisco2
% rubisco2 - a rubisco class object in a chain right after rubisco1
% chain_lead_vec - the average vector of the chain containing rubisco1 and 
%                  rubisco2
%
% Outputs
% twist - the difference in the narot angles of rubisco1 and rubisco2, 
%         taking into account the symmetry of a rubisco. It is between -45 
%         and +45 degrees

    % get the narot angle of each rubisco
    narot1Raw = calc_narot(rubisco1, chain_lead_vec);
    narot2Raw = calc_narot(rubisco2, chain_lead_vec);

    % rubisco has c4 symmetry, so set limits on narot to between 0 and 90
    if narot1Raw < 0
        narot1Sym = mod((360 + mod(narot1Raw, -360)), 90);
    else
        narot1Sym = mod(mod(narot1Raw, 360), 90);
    end

    % rubisco has c4 symmetry, so set limits on narot to between 0 and 90
    if narot2Raw < 0
        narot2Sym = mod((360 + mod(narot2Raw, -360)), 90);
    else
        narot2Sym = mod(mod(narot2Raw, 360), 90);
    end
    
    % calculate the difference between the rubiscos' narot angles
    dAngle = narot2Sym - narot1Sym;

    % rubisco symmetry requires difference be between -45 and 45 degrees
    if dAngle > 45
        twist = 90 - dAngle;
    else
        if dAngle < -45
            twist = 90 + dAngle;
        else
            twist = dAngle;
        end
    end
end