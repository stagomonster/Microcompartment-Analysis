function carb_indices = get_carboxysomes_indices_from_tomogram(data)
% This function takes in a read tomogram data table and returns all the
% unique carboxysome indexes found in column 21 of the data table.
    carb_indices = unique(data(:,21))';
end