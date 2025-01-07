% Important Constants to use throughout OHQ Scripts
function CONSTANTS = constants(pixel_size, diameter)
    persistent RUBISCO_DIAMETER_M PIXEL_SIZE AVOGADRO_NUMBER
    
    if isempty(RUBISCO_DIAMETER_M)
        RUBISCO_DIAMETER_M = 9.936e-9; % Rubisco Diameter, in meters. Change this according to particle.
        PIXEL_SIZE = 2.54e-10; % Pixel Size, in meters (always in bin2). Change this according to data set and particle.
        AVOGADRO_NUMBER = 6.02214076e+23;
    end

    if nargin > 0 && ~isempty(pixel_size) && ~isempty(diameter)
        PIXEL_SIZE = pixel_size;
        RUBISCO_DIAMETER_M = diameter;
    end
    
    CONSTANTS.RUBISCO_DIAMETER_M = RUBISCO_DIAMETER_M;
    CONSTANTS.PIXEL_SIZE = PIXEL_SIZE;
    CONSTANTS.AVOGADRO_NUMBER = AVOGADRO_NUMBER;
end