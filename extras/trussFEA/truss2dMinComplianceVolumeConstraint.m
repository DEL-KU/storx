classdef truss2dMinComplianceVolumeConstraint < truss2d 
    % This class encapsulates methods to minimize the compliance a 2d truss system
    % assumes that the truss members are circular rods
    % The design variables are the cross-sectional areas
    properties(GetAccess = 'public', SetAccess = 'private')
        myInitialVolume;
        myInitialArea;
        myInitialCompliance;
        myInitialStress;
        myFinalArea;
        myFinalCompliance;
        myFinalVolume;
        myFinalStress;
        myYieldStress;
    end
    methods
        function obj = truss2dMinComplianceVolumeConstraint(xy,connectivity)
            obj = obj@truss2d(xy,connectivity);
            obj.myYieldStress(1:obj.myNumTrussBars) = 100e6; % default
        end
        function JRelative = complianceObjective(obj,x)
            Area = x.*obj.myInitialArea;
            obj = obj.assignA(Area);
            obj = obj.assemble();
            obj = obj.solve();
            J = obj.getCompliance();
            JRelative = J/obj.myInitialCompliance;
        end  
        function obj = assignYieldStress(obj,yieldStress,members)
            % assign sMax to one or more members
            % if members is not given, then assign to all 
            if (nargin == 2)
                members = 1:obj.myNumTrussBars;
            else
                assert(max(members) <= obj.myNumTrussBars);
                assert(min(members) >=  1);
            end
            obj.myYieldStress(members) = yieldStress;
        end 
        function obj = initialize(obj)
            obj.myInitialArea = obj.myArea;
            obj.myInitialVolume = sum(obj.myArea.*obj.myL);
            obj = obj.assemble();
            obj = obj.solve();
            obj.myInitialCompliance =  obj.getCompliance(); 
            obj.myInitialStress = obj.myStress;
        end
        function obj = optimize(obj)
            % initialize quantities of interest
            obj = obj.initialize();
            x0 = ones(1,obj.myNumTrussBars); % unitless quantities 
            LB = 1e-12*ones(1,obj.myNumTrussBars); % small non-zero values
            AinEq = (obj.myInitialArea.*obj.myL)/obj.myInitialVolume;
            BinEq = 1; 
%             opt = optimset('fmincon');
%             opt = optimset(opt,'GradObj','off');
            [xMin,~,~,Output]  = fmincon(@obj.complianceObjective,x0, ...
                    AinEq,BinEq,[],[],LB);
            obj = obj.assignA(xMin.*obj.myInitialArea);
            obj = obj.assemble();
            obj = obj.solve();
            obj.myFinalArea = obj.myArea;
            obj.myFinalCompliance =  obj.getCompliance(); 
            obj.myFinalVolume = sum(obj.myArea.*obj.myL);
            obj.myFinalStress = obj.myStress;
            disp(Output)
        end
    end
end

