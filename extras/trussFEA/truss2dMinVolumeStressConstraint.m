classdef truss2dMinVolumeStressConstraint < truss2d 
    % stress-constrained volume minimization
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
        myLambda;
    end
    methods
        function obj = truss2dMinVolumeStressConstraint(xy,connectivity)
            obj = obj@truss2d(xy,connectivity);
            obj.myYieldStress(1:obj.myNumTrussBars) = 100e6; % default
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
        function volRelative = volumeObjective(obj,x)
            Area = x.*obj.myInitialArea;
            obj = obj.assignA(Area);
            vol = sum(obj.myArea.*obj.myL);
            volRelative = vol/obj.myInitialVolume;
        end 
        function [cineq,ceq] = stressConstraint(obj,x)
            Area = x.*obj.myInitialArea;
            obj = obj.assignA(Area);
            obj = obj.assemble();
            obj = obj.solve();
            nConstraints = 2*obj.myNumTrussBars; % two stress constraints per bar
            cineq = zeros(1,nConstraints);
            constraint = 1;
            for m = 1:obj.myNumTrussBars
                cineq(constraint) = obj.myStress(m)/obj.myYieldStress(m)-1;% tension 
                cineq(constraint+1) = -obj.myStress(m)/obj.myYieldStress(m) -1;%compression 
                constraint = constraint+2; % increment
            end 
            ceq = [];        
        end
        function obj = initialize(obj)
            obj.myInitialArea = obj.myArea;
            obj.myInitialVolume = sum(obj.myArea.*obj.myL);
            obj = obj.assemble();
            obj = obj.solve();
            obj.myInitialCompliance =  obj.getCompliance(); 
            obj.myInitialStress = obj.myStress;
        end
        function processLambda(obj)
            ineqnonlin = obj.myLambda.ineqnonlin;
            maxValue = max(abs(ineqnonlin));
            ineqnonlin = ineqnonlin/maxValue; % scaled
            for m = 1:obj.myNumTrussBars
                if (ineqnonlin(2*m-1) > 0.0001)
                    disp(['Bar ' num2str(m) ': Tensile stress active']);
                elseif (ineqnonlin(2*m) > 0.0001)
                    disp(['Bar ' num2str(m) ': Compressive stress active']);
                else
                    disp(['Bar ' num2str(m) ': Stress inactive']);
                end
            end
        end
        function obj = optimize(obj)
            obj = obj.initialize();
            x0 = ones(1,obj.myNumTrussBars); % unitless quantities 
            LB = 1e-12*ones(1,obj.myNumTrussBars); % small non-zero values
            [xMin,~,~,~,Lambda]  = fmincon(@obj.volumeObjective,x0, ...
                   [],[],[],[],LB,[],@obj.stressConstraint);
            obj = obj.assignA(xMin.*obj.myInitialArea);
            obj = obj.assemble();
            obj = obj.solve();
            obj.myFinalArea= obj.myArea;
            obj.myFinalVolume = sum(obj.myArea.*obj.myL);
            obj.myFinalCompliance =  obj.getCompliance(); 
            obj.myFinalStress = obj.myStress;
            obj.myLambda = Lambda;
        end
    end
end

