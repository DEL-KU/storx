%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing minimum feature size constraints             %
% based on convolution kernel used in shape and topology optimization       %
% (primarily used in level-set methods.)                                    %
%                                                                           %
% This Matlab code was written by:                                          %
% - Amir M. Mirzendehdel, Aerospace Engineering Department, KU              %
% - Krishnan Suresh, Mechanical Engineering Department, UW-Madison          %
%                                                                           %
% Please send your comments to: amirzend@ku.edu                             %
%                                                                           %
% The code is intended for educational purposes and theoretical details     %
% are discussed in the textbook:                                            %
% Introduction to Shape and Topology Optimization using MATLAB              %
%                                                                           %
% Disclaimer:                                                               %
% The authors reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event        %
% caused by the use of the program.                                         %
%                                                                           %
% License:                                                                  %
% This software is used, copied and distributed under the licensing         %
% agreement contained in the file LICENSE in the top directory of           %
% the distribution.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef  minimumFeatureSize_conv < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_filterKernel; % minimum feature size filter kernel
    end
    methods(Static)
        function K = makeFilterKernel(R)
            % R: integer radius in elements (R=1 => 3x3, R=2 => 5x5, ...)
            [x, y] = meshgrid(-R:R, -R:R);
            d = abs(x) + abs(y);        % Manhattan distance from center
            w = max(0, R+1 - d);        % triangular weight, zero outside "diamond"
            K = w / sum(w(:));          % normalize to sum = 1
        end

    end
    methods
        %% CONSTRUCTOR
        function obj = minimumFeatureSize_conv(solver,R)
            %MINIMUMFEATURESIZE_CONV Constructor for convolution-based minimum feature size filter.
            %
            %   obj = minimumFeatureSize_conv(solver, R)
            %
            %   This class enforces a minimum feature size in topology optimization by
            %   smoothing (filtering) the design sensitivities using a convolution
            %   kernel. The filter radius is specified in units of elements.
            %
            %   INPUTS:
            %     solver : Handle to the simulation / topology optimization solver.
            %              It is expected to provide:
            %                 - mesh information (element size, number of elements)
            %                 - existing element mask (m_existingElems) used to
            %                   avoid "bleeding" sensitivities outside the design domain.
            %
            %     R      : Integer filter radius in elements (R >= 1).
            %              The corresponding convolution kernel will have size
            %              (2*R+1)-by-(2*R+1). For example:
            %                  R = 1 -> 3x3 kernel
            %                  R = 2 -> 5x5 kernel
            %              Larger R produces stronger smoothing and enforces a larger
            %              minimum feature size.
            %
            %   OUTPUT:
            %     obj    : Instance of minimumFeatureSize_conv with an internally
            %              constructed, normalized filter kernel (stored in
            %              obj.m_filterKernel) that is later used to filter
            %              sensitivities (or densities) via 2D convolution.
            %
            %   NOTE:
            %     - The physical minimum length scale r_min (in length units) may be
            %       related to R through the element size, e.g.
            %           R ≈ ceil(r_min / max(hx, hy)).
            %     - The filter preserves the total "mass" locally by normalizing the
            %       kernel so that sum(m_filterKernel(:)) = 1.
            %

            % check if solver is valid
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid

            if nargin < 2, R = 1; end
            % constructor based on superclass
            obj = obj@mfgConstraints(solver);

            % feature size filter
            obj.m_filterKernel = obj.makeFilterKernel(R);
        end


        function [filteredDesign] = filterDesign(obj, design)
            % filter the design variables
            % input: obj, design variables
            % output: obj, filtered design variables
            filteredDesign = design .* obj.m_solver.m_existingElems;
        end

        function filteredSensitivity = filterSensitivity(obj, ~, sensField)
            % filter the sensitivity field
            % input:  obj, sensField
            % output: filteredSensitivity (same size as sensField)

            existingElems = obj.m_solver.m_existingElems;   % 1 = active, 0 = inactive

            % ---------------------------------------------------------------------
            % 1) Zero sensitivity outside the design domain
            % ---------------------------------------------------------------------
            sensField = sensField .* existingElems;

            % ---------------------------------------------------------------------
            % 2) Determine padding from kernel size (supports arbitrary odd-sized kernels)
            % ---------------------------------------------------------------------
            [kh, kw] = size(obj.m_filterKernel);
            padH = floor(kh/2);
            padW = floor(kw/2);

            % ---------------------------------------------------------------------
            % 3) Pad fields
            %    - sensitivities: replicate (to avoid artificial zeros at boundary)
            %    - mask: zeros (so neighbors outside domain are not counted)
            % ---------------------------------------------------------------------
            paddedSensField   = padarray(sensField,      [padH, padW], 'replicate');
            paddedExisting    = padarray(existingElems,  [padH, padW], 0);

            % ---------------------------------------------------------------------
            % 4) Convolution
            % ---------------------------------------------------------------------
            convSensField     = conv2(paddedSensField,  obj.m_filterKernel, 'valid');
            convExistingElems = conv2(paddedExisting,   obj.m_filterKernel, 'valid');

            % ---------------------------------------------------------------------
            % 5) Normalize ONLY where we actually have existing elements
            % ---------------------------------------------------------------------
            filteredSensitivity = zeros(size(sensField));

            % valid locations: inside design domain and with at least one neighbor
            validMask = (existingElems ~= 0) & (convExistingElems > 0);

            % avoid division by zero
            denom = max(convExistingElems(validMask), 1e-12);
            filteredSensitivity(validMask) = convSensField(validMask) ./ denom;
        end
    end
end
