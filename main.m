function carboxysome_data = main()
% This main function will execute the complete script pipeline that must
% be followed to obtain the chain and lattice output from an input table.
%
% Outputs:
% carboxysome_data - An array of carboxysome objects populated with data
%                    from every function in the main pipeline.
%
% main.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    % Run Dynamo to allow specialized functions
    run dynamo_activate.m

    % Add matlab paths to access scripts in this package
    addpath('helper_functions/');
    addpath('main_pipeline/');
    addpath('classes/');
    addpath('data/');

    % Obtain and process input data from file
    num_files = input("Is your input in a single file, or is it separated into even and odd files? Enter 1 for a single file or 2 for split files: ");
    if num_files == 1
        filename = input("Enter the name of your data file. Please ensure that your input data is in bin1 or 2: ", "s");
        outfile = input("What would you like to call your processed input data file? Please end the name with .tbl: ", "s");
        clean_value = input("What would you like the clean size to be? ");

        is_bin2 = input("Is your input file in bin2? Enter 1 if your data is in bin1 or 2 if your data is bin2: ");
        if is_bin2 == 2
            is_bin2 = true;
        elseif is_bin2 == 1
            is_bin2 = false;
        else
            error("Invalid input. Enter 1 or 2.");
        end

        preprocess_datasets(filename, outfile, clean_value, is_bin2)
    elseif num_files == 2
        odd_filename = input("Enter the name of the odd input data file: ", "s");
        even_filename = input("Enter the name of the even input data file: ", "s");
        outfile = input("What would you like to call your processed input data file? Please end the name with the extension .tbl: ", "s");
        clean_value = input("What would you like the clean to be? ");

        is_bin2 = input("Is your input file in bin2? Enter 1 if your data is in bin1 or 2 if your data is bin2: ");
        if is_bin2 == 2
            is_bin2 = true;
        elseif is_bin2 == 1
            is_bin2 = false;
        else
            error("Invalid input. Enter 1 or 2.");
        end

        preprocess_datasets(odd_filename, even_filename, outfile, clean_value, is_bin2);
    else
        error("Invalid input. Enter 1 or 2.");
    end

    % Set the filename to be worked on to be the name of the output
    % processed data file.
    filename = outfile;

    % Set the correct pixel size in the constants file
    pixel_size = input("Enter the pixel size, in meters, that your data uses. Please use scientific notation (ex. 2.54e-10) and keep in bin1 if you indicated that your input data was bin1: ");
    if ~is_bin2
        pixel_size = pixel_size * 2;
    end
    
    diameter = input("Enter the diameter, in meters, of your particle. Please use scientific notation again: ");
    CONSTANTS = constants(pixel_size, diameter);

    % SETUP READY. Begin Pipeline
    %% STEP 1: CONVEX HULL AND VOLUME
    fprintf('\nConvex Hull and Volume:\n');
    carb_volumes = convex_hull_and_volume(filename);

    %% STEP 2: LOCAL GLOBAL ALIGNMENT
    fprintf('\nLocal Global Alignment:\n');
    local_max_distance = input('Enter the maximum distance (in pixels) around a particle you want to search for potential oligimerization partners (to use default value of 2 * rubisco diameter enter "default"): ', 's');
    
    if strcmp(local_max_distance, 'default')
        carb_local = local_global_alignment(filename, carb_volumes);
    else
        carb_local = local_global_alignment(filename, carb_volumes, str2double(local_max_distance));
    end

    %% STEP 3: LINKAGES
    fprintf('\nLinkages:\n');
    linkages_choice = input('Would you like to use a custom set of parameters or a default set? Enter "custom" or "default": ', 's');

    if strcmp(linkages_choice, 'default')
        carb_linkages = linkages(filename, carb_local);
    else
        min_factor = input('Enter the minimum distance to search along the particle axis in terms of a multiple of the particle diameter: ');
        max_factor = input('Enter the maximum distance to search along the particle axis in terms of a multiple of the particle diameter: ');
        rad_factor = input('Enter the maximum distance to search radially from the particle axis in terms of a multiple of the particle diameter: ');
        max_angle = input('Enter the maximum angle between particle axes to search in (in degrees): ');
        mutual = input("Do you want to require that both particles in a given pair be in each other's search volumes? Enter 1 for yes, 0 for no: ");
        carb_linkages = linkages(filename, carb_local, min_factor, max_factor, rad_factor, max_angle, mutual);
    end

    %% STEP 4: CHAIN MAKER
    fprintf('\nChain Maker:\n');
    chain_min_length = input('Enter the minimum length to be considered a chain or "default" to use the default value of 3: ', 's');

    if strcmp(chain_min_length, 'default')
        carb_chains = chain_maker(filename, carb_linkages);
    else
        carb_chains = chain_maker(filename, carb_linkages, str2double(chain_min_length));
    end

    start_over = input('\nWould you like to terminate the program to use the new table (y/n)? ', 's');
    if strcmp(start_over, 'n')
        %% STEP 5: CHAIN LINKAGES
        fprintf('\nChain Linkages:\n');
        chain_links_choice = input('Would you like to use a custom set of parameters or the default set? Enter "custom" or "default": ', 's');
        
        if strcmp(chain_links_choice, 'default')
            carb_chain_links = chain_linkages(filename, carb_chains);
        else
            min_dist = input('Enter the minimum distance to allow between linked chains in nanometers: ');
            max_dist = input('Enter the maximum distance to allow between linked chains in nanometers: ');
            max_angle = input('Enter the maximum angle to allow between linked chains in degrees: ');
            carb_chain_links = chain_linkages(filename, carb_chains, min_dist, max_dist, max_angle);
        end    
    
        %% STEP 6: LATTICE GEN
        fprintf('\nLattice Gen:\n');
        lattice_threshold = input('Enter the minimum number of linked chains to constitute a lattice or "default" to use the default value of 5: ', 's');
    
        if strcmp(lattice_threshold, 'default')
            carboxysome_data = lattice_gen(filename, carb_chain_links);
        else
            carboxysome_data = lattice_gen(filename, carb_chain_links, str2double(lattice_threshold));
        end
    else
        carboxysome_data = carb_chains;
    end
end