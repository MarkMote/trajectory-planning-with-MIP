classdef linOpt < handle
    %LINOPT translates high level MIP specifications into a format that can
    %be used with Gurobi
    
    properties (GetAccess = public, SetAccess = private)
        variables
        numVars   = 0;
        numStates = 0;
        numEqns   = 0;
        A;                  % size A = numEqns x numVars
        sense = '';              % size sense = numEqns
        b;                  % rhs
        vType = '';              % variable types
        epsilon = 10e-10    % Used to create inequality from strict inequality
        Mp = 10e6;
    end
    
    methods
        %% Base Functions
        function this = linOpt()
            % Stuff that happens when the object is created
            
        end
        
        function VarInds = define(this, name, elements, numAgents, type) % Defines variables
            %DEFINE defines a var and assigns each element a unique index
            if nargin == 3
                VarInds = this.numStates + (1:elements);  % return variable indicies
                
                this.numVars   = this.numVars + 1;        % One additional variable defined
                this.numStates = this.numStates+elements; % Add additional states
                this.variables{this.numVars,1} = name;
                this.variables{this.numVars,2} = elements;
                this.variables{this.numVars,3} = VarInds;
            else
                VarInds = zeros(numAgents,elements);
                for n = 1:numAgents
                    VarInds(n,:) = this.numStates + (1:elements);  % return variable indicies
                    
                    this.numVars   = this.numVars + 1;        % One additional variable defined
                    this.numStates = this.numStates+elements; % Add additional states
                    this.variables{this.numVars,1} = strcat(name,'(',num2str(n),',:)');
                    this.variables{this.numVars,2} = elements;
                    this.variables{this.numVars,3} = VarInds(n,:);
                end
            end
            
            if nargin == 5
                this.vType(VarInds) = type;
            else
                this.vType(VarInds) = 'C';
            end
            
        end
        
        
        function this = declareAtom( this, atom )
            lhs = atom{1};
            op  = atom{2};
            rhs = atom{3};
            this.numEqns = this.numEqns + 1 ;
            for i = 1:2:length(lhs)
                this.A(this.numEqns,lhs(i+1)) = lhs(i) ;
            end
            this.b(this.numEqns) = rhs;
            this.sense(this.numEqns) = op;
        end
        
        function this = atomImplies( this, a, b)
            if iscell(a) && ~iscell(b) % atom -> BDV
                if strcmp(a{2},'<')
                    this.declareAtom( { [a{1}, this.Mp, b] , '>' , a{3} } );
                elseif strcmp(a{2},'>')
                    this.declareAtom( { [a{1}, -this.Mp, b], '<' , a{3} } );
                else
                    fprintf('\nEquality not yet integrated')
                    pause(5);
                end
            elseif ~iscell(a) && iscell(b) % BDV -> atom
                if strcmp(b{2},'<')
                    this.declareAtom( { [b{1}, this.Mp, a] , '<' , b{3}+this.Mp } );
                elseif strcmp(b{2},'>')
                    this.declareAtom( { [b{1}, -this.Mp, a], '>' , b{3}-this.Mp} );
                    
                else
                    this.declareAtom( { [b{1}, this.Mp, a] , '<' , b{3}+this.Mp } );
                    this.declareAtom( { [b{1}, -this.Mp, a], '>' , b{3}-this.Mp} );
                end
            else % atom -> atom
                zeta = this.define('zetaSat',1,1,'B');
                if strcmp(a{2},'<')
                    this.declareAtom( { [a{1}, this.Mp, zeta] , '>' , a{3} } );
                    this.declareAtom( { [a{1}, this.Mp, zeta] , '<' , a{3}+this.Mp } ); % exp (zeta => a)
                elseif strcmp(a{2},'>')
                    this.declareAtom( { [a{1}, -this.Mp, zeta], '<' , a{3} } );
                    this.declareAtom( { [a{1}, -this.Mp, zeta], '>' , a{3}-this.Mp } ); % exp (zeta => a)
                else
                    fprintf('\nEquality not yet integrated')
                    pause(5);
                end
                if strcmp(b{2},'<')
                    this.declareAtom( { [b{1}, this.Mp, a] , '<' , b{3}+this.Mp } );
                elseif strcmp(b{2},'>')
                    this.declareAtom( { [b{1}, -this.Mp, a], '>' , b{3}-this.Mp} );
                else
                    this.declareAtom( { [b{1},  this.Mp, zeta] , '<' , b{3}+this.Mp} );
                    this.declareAtom( { [b{1},  -this.Mp, zeta], '>' , b{3}-this.Mp} );
                end
            end
        end
        
        
        function this = atomEquiv( this, a, b)
            if ~iscell(a) || ~iscell(b)
                this.atomImplies(a,b);
                this.atomImplies(b,a);
            else % For case of two atoms, use the same BDV
                zeta = this.define('zetaSat',1,1,'B');
                for loop = 1:2
                    if strcmp(a{2},'<')
                        this.declareAtom( { [a{1}, this.Mp, zeta] , '>' , a{3} } );
                        this.declareAtom( { [a{1}, this.Mp, zeta] , '<' , a{3}+this.Mp } ); % exp (zeta => a)
                    elseif strcmp(a{2},'>')
                        this.declareAtom( { [a{1}, -this.Mp, zeta], '<' , a{3} } );
                        this.declareAtom( { [a{1}, -this.Mp, zeta], '>' , a{3}-this.Mp } ); % exp (zeta => a)
                    else
                        fprintf('\nEquality not yet integrated')
                        pause(5);
                    end
                    if strcmp(b{2},'<')
                        this.declareAtom( { [b{1}, this.Mp, a] , '<' , b{3}+this.Mp } );
                    elseif strcmp(b{2},'>')
                        this.declareAtom( { [b{1}, -this.Mp, a], '>' , b{3}-this.Mp} );
                    else
                        this.declareAtom( { [b{1},  this.Mp, zeta] , '<' , b{3}+this.Mp} );
                        this.declareAtom( { [b{1},  -this.Mp, zeta], '>' , b{3}-this.Mp} );
                    end
                    c = a;
                    a = b;
                    b = c;
                end
            end
        end
        
        
        %% Helpful functions
        function this = normLeq(this, x, xmax, P)
            for p = 1:P
                this.declareAtom({ [sin(2*pi*p/P), x(1), cos(2*pi*p/P), x(2)] , '<', xmax});
            end
        end
        
        function this = avoidBox(this, x, box)
            psi1 = this.define('boxXmin',1,1,'B');
            psi2 = this.define('boxXmax',1,1,'B');
            psi3 = this.define('boxYmin',1,1,'B');
            psi4 = this.define('boxYmax',1,1,'B');
            this.atomImplies({ [1, x(1)], '>', box(1) }, psi1);
            this.atomImplies({ [1, x(1)], '<', box(2) }, psi2);
            this.atomImplies({ [1, x(2)], '>', box(3) }, psi3);
            this.atomImplies({ [1, x(2)], '<', box(4) }, psi4);
            this.declareAtom({ [1, psi1, 1, psi2, 1, psi3, 1, psi4], '<', 3});
        end
        
        function this = saturate(this, v, vc, vMax)
            vMin = - vMax;
            for i = 1:2
                zetaMax1 = this.define('zetaMax1',1,1,'B');
                zetaMin1 = this.define('zetaMin1',1,1,'B');
                if zetaMax1 == 105
                    1 
                end                
                % vc > vMax (<=>zetaMax1) => v=vMax
                this.atomEquiv({ [1, vc(i)], '>', vMax}, zetaMax1);
                this.atomImplies( zetaMax1 , { [1,v(i)],'=',vMax});

                
                % vc < vMin => v=vMin
                this.atomEquiv( zetaMin1 , { [1, vc(i)], '<', vMin});
                this.atomImplies( zetaMin1 , { [1,v(i)],'=',vMin});
                
                % (vc > vMin && vc < vMax) => (v=vc)
                zetaNom1 = this.define('zetaNom1',1,1,'B');
%                 fprintf('\n%d\n',zetaNom1);
                this.atomEquiv( {[1, zetaMin1, 1, zetaMax1], '<', 0.5}, zetaNom1);
                this.atomImplies( zetaNom1, {[1, v(i), -1, vc(i)], '=', 0});
            end
        end
        
        function this = expand_A(this)
            %EXPAND_A ensures that A has the right number of rows and cols
            [~, rowsA ] = size(this.A);
            if rowsA < this.numStates
                this.A(this.numEqns,this.numStates) = 0;
            end
        end
    end
    
end

