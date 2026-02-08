classdef truss2dMinVolumeComplianceConstraint < truss2d 
    % compliance consrained volume minimization
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
        function obj = truss2dMinVolumeComplianceConstraint(xy,connectivity)
            obj = obj@truss2d(xy,connectivity);
            obj.myYieldStress(1:obj.myNumTrussBars) = 100e6; % default
        end
        function volRelative = volumeObjective(obj,x)
            Area = x.*obj.myInitialArea;
            obj = obj.assignA(Area);
            vol = sum(obj.myArea.*obj.myL);
            volRelative = vol/obj.myInitialVolume;
            
        end 
        function [cineq,ceq] = complianceConstraint(obj,x)
            Area = x.*obj.myInitialArea;
            obj = obj.assignA(Area);
            obj = obj.assemble();
            obj = obj.solve();
            J = obj.getCompliance();
            cineq = [];
            ceq = J/obj.myInitialCompliance - 1;
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
            obj = obj.initialize();
            x0 = ones(1,obj.myNumTrussBars); % unitless quantities      
            LB = 1e-8*ones(1,obj.myNumTrussBars); % small non-zero values
            opt = optimoptions('fmincon','PlotFcns',@optimplotfval, 'MaxFunEvals',25000);%
            [xMin,~,~,~]  = fmincon(@obj.volumeObjective,x0, ...
                   [],[],[],[],LB,[],@obj.complianceConstraint,opt);
            obj = obj.assignA(xMin.*obj.myInitialArea);
            obj = obj.assemble();
            obj = obj.solve();
            obj.myFinalArea = obj.myArea;
            obj.myFinalVolume = sum(obj.myArea.*obj.myL);
            obj.myFinalCompliance =  obj.getCompliance();
            obj.myFinalStress = obj.myStress;
        end
    end
end

