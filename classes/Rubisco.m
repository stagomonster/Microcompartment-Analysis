classdef Rubisco < handle
% Model for storing data read from the .tbl file in a convinient way
%
% Rubisco.m © 2025 is licensed under CC BY-NC-SA 4.0
    
    properties
        tag uint32 = nan; % The rubisco's identifying number
        index uint16 = nan; % Rubisco's place in its carboxysome
        aligned = 0; % A boolean, 1 if the particle underwent alignment
        averaged = 0; % A boolean, 1 if the particle was included in an average
        dx = nan; % x, y, and z offsets from the center of its subtomogram
        dy = nan;
        dz = nan;
        tdrot = nan; % euler angles of the rubisco orientation
        tilt = nan;
        narot = nan;
        cc = nan; % cross correlation
        cc2 = nan;
        cpu = nan; % cpu used for alignment
        ftype = nan; % Fourier sampling
        ymintilt = nan; % y tilt range
        ymaxtilt = nan;
        xmintilt = nan; % x tilt range
        xmaxtilt = nan;
        fs1 = nan; % Fourier sampling free parameters
        fs2 = nan;
        tomo uint16 = nan; % tomogram it belongs to
        reg uint16 = nan; % carboxysome it belongs to
        class uint16 = nan; % particle class
        annotation uint16 = nan; % arbitrary label
        x = nan; % x, y, and z coordinates of subtomogram box center
        y = nan;
        z = nan;
        dshift = nan; % norm of shift vector
        daxis = nan; % difference in axis and narot orientations
        dnarot = nan;
        dcc = nan; % difference in cc
        otag = nan; % original tag from subboxing
        npar= nan; % number of added particles divided by subunit label
        ref= nan; % reference (for multireference projects)
        sref= nan; % subreference
        apix= nan; % angstroms per pixel
        def= nan; % defocus in microns
        eig1= nan; % eigencoefficients
        eig2= nan;
        ave_normal = []; % average of surrounding vectors normal to carboxysome surface
        vector= []; % orientation vector of rubisco
        inside = false; % inner or outer rubisco
        rubisco_above_me = nan; % the tag of the rubisco above this one
        rubisco_below_me = nan; % the tag of the rubisco below this one
        in_central_chain = false; % whether the rubisco is in the central chain of a lattice
    end

    methods
        function this_rubisco = Rubisco(tag1,aligned1,averaged1,dx1,dy1,dz1,tdrot1,...
                tilt1,narot1,cc1,cc21,cpu1,ftype1,ymintilt1,ymaxtilt1,xmintilt1,...
                xmaxtilt1,fs11,fs21,tomo1,reg1,class1,annotation1,x1,y1,z1,dshift1,...
                daxis1,dnarot1,dcc1,otag1,npar1,ref1,sref1,apix1,def1,eig11,eig21, ave_norm, vec)
            % Constructor: create an instance of container.
            if nargin > 0
                this_rubisco.tag = tag1;
                this_rubisco.aligned = aligned1;
                this_rubisco.averaged = averaged1;
                this_rubisco.dx = dx1;
                this_rubisco.dy = dy1;
                this_rubisco.dz = dz1;
                this_rubisco.tdrot = tdrot1;
                this_rubisco.tilt = tilt1;
                this_rubisco.narot = narot1;
                this_rubisco.cc = cc1;
                this_rubisco.cc2 = cc21;
                this_rubisco.cpu = cpu1;
                this_rubisco.ftype = ftype1;
                this_rubisco.ymintilt = ymintilt1;
                this_rubisco.ymaxtilt = ymaxtilt1;
                this_rubisco.xmintilt = xmintilt1;
                this_rubisco.xmaxtilt = xmaxtilt1;
                this_rubisco.fs1 = fs11;
                this_rubisco.fs2 = fs21;
                this_rubisco.tomo = tomo1;
                this_rubisco.reg = reg1;
                this_rubisco.class = class1;
                this_rubisco.annotation = annotation1;
                this_rubisco.x = x1;
                this_rubisco.y = y1;
                this_rubisco.z = z1;
                this_rubisco.dshift = dshift1;
                this_rubisco.daxis = daxis1;
                this_rubisco.dnarot = dnarot1;
                this_rubisco.dcc = dcc1;
                this_rubisco.otag = otag1;
                this_rubisco.npar = npar1;
                this_rubisco.ref = ref1;
                this_rubisco.sref = sref1;
                this_rubisco.apix = apix1;
                if nargin == 35
                    this_rubisco.def = NaN;
                    this_rubisco.eig1 = NaN;
                    this_rubisco.eig2 = NaN;
                    this_rubisco.ave_normal = [];
                    this_rubisco.vector = [];
                else
                    this_rubisco.def = def1;
                    this_rubisco.eig1 = eig11;
                    this_rubisco.eig2 = eig21;
                    this_rubisco.ave_normal = ave_norm;
                    this_rubisco.vector = vec;
                end
            end
        end

        function disp(self)
            % Print information about self to command window
            disp("Rubisco object with properties:")
            props = properties(self);
            for i = 1:length(props)
                name = props{i};
                value = self.(name);
                if isnan(value)
                    value = 'NaN';
                elseif ~isscalar(value)
                    value = mat2str(value);
                end
                disp("  " + name + " = " + value);
            end
        end

        function obj_copy = copy(obj)
            % Create a copy of a rubisco object
            obj_copy = Rubisco();
            props = properties(obj);
            for i = 1:length(props)
                obj_copy.(props{i}) = obj.(props{i});
            end
        end
    end
end

    