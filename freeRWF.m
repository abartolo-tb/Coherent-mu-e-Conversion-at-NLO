function [rValues, u1, u2] = freeRWF(m, Z, ...
                                          kappa, E, n, paramVec)
%freeRWF Calculates dirac wavefunction for free electron
%
%   Returns radial values, and the reduced spinor 
%   wavefunctions for the system.
%   Returned rValues are in units of [Fermi]^1
%   Returned u1, u2 are dimensionless
%
%   Inputs:
%   m = mass of the ejected particle in MeV (.5109989 for electron)
%   Z = atomic number
%   kappa = orbital eigenvalue of the Dirac Equation (+/- 1)
%   E = energy of the outgoing particle in MeV
%   n = integer specifying charge distribution used
%   paramVec = is a row vector of parameters characterizing 
%              the charge distribution of the nucleus
%
%   The following is a key for the integer n used to specifying
%   each charge distribution and the list of parameters expected
%
%   Key:
%
%   1 - Modified Harmonic Oscillator
%       paramVec = [alpha, a]
%       p(r) = norm*(1+alpha*(r/a)^2)*exp(-(r/a)^2)
%       a has dimensions of [Fermi]
%       alpha is dimensionless
%
%   2 - Two Parameter Fermi Distribution
%       paramVec = [c, y]
%       p(r) = norm/(1+exp((r-c)/y))
%       c has dimensions of [Fermi]
%       y has dimensions of [Fermi]
%
%   3 - Three Parameter Fermi Distribution
%       paramVec = [c, y, w]
%       p(r) = norm*(1+w*(r/c)^2)/(1+exp((r-c)/y))
%       c has dimensions of [Fermi]
%       y has dimensions of [Fermi]
%       w is dimensionless
%
%   4 - Two Parameter Gaussian Distribution
%       paramVec = [c, y]
%       p(r) = norm/(1+exp((r^2-c^2)/(y^2)))
%       c has dimensions of [Fermi]
%       y has dimensions of [Fermi]
%
%   5 - Three Parameter Gaussian Distribution
%       paramVec = [c, y, w]
%       p(r) = norm*(1+w*(r/c)^2)/(1+exp((r^2-c^2)/(z^2)))
%       c has dimensions of [Fermi]
%       y has dimensions of [Fermi]
%       w is dimensionless
%
%   6 - Sum of Gaussian Expansion
%       paramVec = [R_1, ..., R_n, Q_1, ..., Q_n, Gamma]
%       p(r) = Sum{i} A_i*[exp(-((r-R_i)/Gamma)^2) + ...
%                                           exp(-((r-R_i)/Gamma)^2)]
%       A_i = norm*Q_i/(2*pi^(3/2)*Gamma^3*(1+2(R_i/Gamma)^2))
%       R_i has dimensions of [Fermi]
%       Q_i is dimensionless
%       Gamma has dimensions of [Fermi]
%       We require Sum{i} Q_i = 1
%
%   7 - Fourier Bessel Expansion
%       paramVec = [R, a_1, a_2, ..., a_n]
%       p(r<=R) = Sum{k} a_k*j0(k*pi*r/R)
%       p(r>R) = 0
%       We use the spherical bessel function j0(z) = sin(z)/z
%       R has dimensions of [Fermi]
%       a_k has dimensions of [Fermi]^(-3)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ensure static inputs are physical

    % Mass of bound particle must be positive
    if(m<=0)
        error('particle mass must be positive');
    end

    % Nuclei need at least one proton
    if(Z<=0)
        error('nucleus must be charged');
    end
    
    % Highly unphysical nuclei will violate assumptions
    % made in choosing safe boundary conditions
    if(Z>=137)
        error('nucleus must be physical');
    end
    
    % Kappa must be a valid eigenvalue
    if (kappa~=-1)&&(kappa~=1)
        error('kappa must be +1 or -1');
    end
    
    % Energy of particle must be positive
    if (E<=0)
        error('energy must be positive');
    end


% Initialize relevant constants and the radial values
% Note: All calculations are done in dimensionless mass units

    % Fine Structure Constant
    alpha = 1/137.036;
    % Muon Mass (MeV)
    mMuon = 105.65837;
    % Free particle mass in units of muon mass
    m = m/(mMuon);
    % Energy in units of muon mass
    E = E/mMuon;
    % Muonic Bohr Radius
    rBohr = (1/(Z*alpha));
    % Extract distribution parameters from input
    params = paramExtract(mMuon, n, paramVec);

    % Find charge radius of the nucleus
    % Initialize temporary radial values
    numPoints = 5000;
    rMin = rBohr/10000;
    rMax = rBohr*10;
    logMin = log10(rMin/rBohr);
    logMax = log10(rMax/rBohr);
    scaledRValues = logspace(logMin, logMax, numPoints)';
    rValues = rBohr*scaledRValues;
    clear logMin logMax scaledRValues numPoints rMin rMax
    % Find charge distribution
    [~, chargeDist] = chargePotential(alpha, Z, n, params, rValues);
    % Find nuclear charge radius. This is currently defined as
    % the radius at which 99% of the charge is contained.
    fracCharge = .99;
    rNucleus = chargeRadius(chargeDist, rValues, fracCharge);
    clear rValues chargeDist fracCharge
    
    % Other values of interest (in dimensionless mass units)
    rMax = 10*max([rNucleus, rBohr]);
    rMin = rNucleus/10000;
    
    % Note, rBCi corresponds to the point where we may begin using
    % analytic solutions for the wavefunctions.
    rBCi = 10*max(rNucleus, Z*alpha/E);
    % Ensure rBCi is not greater than rMax
    if(rBCi>rMax)
        rBCi=rMax;
    end
    %clear rNucleus



% First, perform Grid-Search to find rough global minimum in 
% parameter space

    % Initialize range of radial wavefunction
    numPoints = 5000;
    logMin = log10(rMin/rBohr);
    logMax = log10(rMax/rBohr);
    scaledRValues = logspace(logMin, logMax, numPoints)';
    rValues = rBohr*scaledRValues;
    clear logMin logMax scaledRValues numPoints

    % Initialize potential over the range of rValues
    [potential, ~] = chargePotential(alpha, Z, n, params, rValues);
    
    % Find rMatch near rNucleus for Shoot-and-Match procedure
    temp = rValues(rValues >= rNucleus);
    rMatch = temp(1);
    clear temp
    
    % Find find rBC near rBCi to evaluate far-field
    % boundary conditions. If rBCi is less than rMatch,
    % set it equal to rMatch
    if(rBCi<rMatch)
        rBC = rMatch;
    else
        temp = rValues(rValues >= rBCi);
        if(isempty(temp))
            rBC = rValues(end);
        else
            rBC = temp(1);
        end
        clear temp
    end
    
    % Define convergence condition for purpose of determining
    % absolute and relative tolerance in ODEsolvers
    convError = .1;

    % Initialize search parameters
    scalingA = (0:.15:3)';
    scalingP = (0:.04:1)';
    % Create list of amplitudes and phases
    listA = 1*scalingA;
    listP = 2*pi*scalingP;
    clear scalingA scalingP
    % Find size of each list for control purposes
    numA = size(listA,1);
    numP = size(listP,1);
    
    % Prepare sliced variables to be used in parallel-for loop
    tempList = (0:1:(numA*numP-1))';
    tempSolutions = zeros(numA*numP, 1);
    % Use modular arithmetic to find appropriate indices for
    % traversing grid in parameter space
    indA = mod(tempList, numA) + 1;
    indP = floor(tempList/numA) + 1;
    % Generate lists of parameter pairs to be used
    tempA = listA(indA);
    tempP = listP(indP);
    clear tempList indA indP
    
    % Parallelize
    matlabpool open
    
    % Begin iteration
    parfor i = 1:(numA*numP)
        
        % Initialize solution space parameters
        a0 = tempA(i);
        phase = tempP(i);

        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, ...
                                    a0, E, kappa, m, convError);
        
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                                    E, phase, kappa, m, convError);
        
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecFree);
        
        % For purposes of the grid search, the mismatch is made
        % a scalar by adding each component in quadrature.
        delta = delU1^2 + delU1P^2 + delU2^2 + delU2P^2;
        
        % Store mismatch error and prepare for next loop
        tempSolutions(i) = delta;

    end

    % End parallelization
    matlabpool close
    
%     % Create 3d plot of mismatch over the parameter space
%     solutions = zeros( numA , numP );
%     for i=(1:numA*numP)
%         temp = tempSolutions(i);
%         indA = mod(i-1, numA)+1;
%         indP = floor((i-1)/numA)+1;
%         solutions(indA, indP) = temp;
%     end
%     newSol = solutions;
%     % Truncate solutions with a mismatch of over 30 for ease
%     % of viewing the plot
%     ind = (newSol>30);
%     newSol(ind)=30;
%     figure;
%     surf(listP, listA, newSol);
%     hold on
%     xlabel('Phase [Dimensionless]');
%     TeXString = texlabel('a_{0} [Dimensionless]');
%     ylabel(TeXString);
%     zlabel('Solution Mismatch in Quadrature');
%     title('Grid Search of Wavefunction Parameter Space');
%     hold off
%     clear newSol ind solutions temp indA indP TeXString
    
    
    % Return parameters that give the minimum from grid search
    [~, ind] = min(tempSolutions);
    minVec = [tempA(ind), tempP(ind)];
    
    clear ind rMatch convError i tempA tempP delta
    clear numA numP rValues tempSolutions listA listP



% Now, with a rough idea of the minimum in parameter space, find the
% best solution for a low number of radial points and low
% convergence criterion

    % Initialize range of radial wavefunction
    numPoints = 10000;
    logMin = log10(rMin/rBohr);
    logMax = log10(rMax/rBohr);
    scaledRValues = logspace(logMin, logMax, numPoints)';
    rValues = rBohr*scaledRValues;
    clear logMin logMax scaledRValues numPoints
    
    % Initialize potential over range of rValues
    [potential, ~] = chargePotential(alpha, Z, n, params, rValues);

    % Find rMatch near rNucleus for Shoot-and-Match procedure
    temp = rValues(rValues >= rNucleus);
    rMatch = temp(1);
    clear temp
    
    % Find find rBC near rBCi to evaluate far-field
    % boundary conditions. If rBCi is less than rMatch,
    % set it equal to rMatch
    if(rBCi<rMatch)
        rBC = rMatch;
    else
        temp = rValues(rValues >= rBCi);
        if(isempty(temp))
            rBC = rValues(end);
        else
            rBC = temp(1);
        end
        clear temp
    end
    
    % Phase and A0 from grid search minimum
    phase = minVec(2);
    a0 = minVec(1);
    clear minVec
        
    % Define convergence condition to be satisfied
    convError = .01;
    conv = false;
        
    % Begin iteration of wavefuctions. Continued until either
    % convergence conditions are satisfied or loop limit is reached.
    
    % Set maximum number of loops and loop counter
    maxLoop = 25;
    loopCount = 0;
    
    % Begin iteration
    while (conv == false)&&(loopCount < maxLoop)
        
        % Increment number of loops
        loopCount = loopCount + 1;
        
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, ...
                                a0, E, kappa, m, convError);
        
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                                E, phase, kappa, m, convError);
        
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecFree);
        
        % Test for convergence
        conv = convTest(delU1, delU1P, delU2, delU2P, convError);
        
        % If solutions have converged, exit loop
        if (conv==true)
            clear uVecNucl uVecFree
            clear delU1 delU1P delU2 delU2P
            break
        end
        
        % If solutions failed to converge, find better values for
        % the parameters a0 and phase
        
        % Pass current parameters to function which calculates the
        % N-R matrix
        M = nrMatrix(a0, phase, rValues, rBC, delU1, delU2, ...
                    potential, rMatch, E, kappa, m, convError);
        % From N-R matrix, find new values for a0 and phase
        newVec = [a0 ; phase] - (M\[delU1 ; delU2]);
        a0New = newVec(1);
        phaseNew = newVec(2);
        
        % Check if the solution's amplitude is approaching zero
        if(a0New<=.01)
            a0Min = .01;
            if(a0Min==a0)
                % In this situation, a0 is already at the 
                % minimum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the a0
                % component
                a0New = a0;
            else
                % Rescale gradient step so values are physical
                rescale = abs((a0Min - a0)/(a0New-a0));
                a0New = a0Min;
                phaseNew = phase + rescale*(phaseNew - phase);
            end
            clear rescale a0Min
        end
        % Check if the solution's amplitude is growing without
        % bound
        if(a0New>=3)
            a0Max = 3;
            if(a0Max==a0)
                % In this situation, a0 is already at the 
                % maximum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the a0
                % component
                a0New = a0;
            else
                % Rescale gradient step so values are physical
                rescale = abs((a0Max - a0)/(a0New-a0));
                a0New = a0Max;
                phaseNew = phase + rescale*(phaseNew - phase);
            end
            clear rescale a0Min
        end
        % Update a0
        a0 = a0New;
        
        % Update phase and check that it is in 
        % the interval [0, 2*pi)
        phase = mod(phaseNew, 2*pi);
        
        clear delU1 delU2 M newVec
        clear delU1P delU2P
        % Repeat loop
        
    end

    % Clear old temporary variables
    clear maxLoop loopCount potential rMatch



% At this point, the rough solutions to the radial wave equation
% have provided all the information they can. Now, with rough 
% values for the phase and scaling factor, we can solve the system
% in greater detail

    % Initialize range of radial wavefunction
    numPoints = 100000;
    logMin = log10(rMin/rBohr);
    logMax = log10(rMax/rBohr);
    scaledRValues = logspace(logMin, logMax, numPoints)';
    rValues = rBohr*scaledRValues;
    clear logMin logMax scaledRValues rMax rMin numPoints
    
    % Initialize potential over range of rValues
    [potential, ~] = chargePotential(alpha, Z, n, params, rValues);

    % Find rMatch near rNucleus for Shoot-and-Match procedure
    temp = rValues(rValues >= rNucleus);
    rMatch = temp(1);
    clear temp
    
    % Find find rBC near rBCi to evaluate far-field
    % boundary conditions. If rBCi is less than rMatch,
    % set it equal to rMatch
    if(rBCi<rMatch)
        rBC = rMatch;
    else
        temp = rValues(rValues >= rBCi);
        if(isempty(temp))
            rBC = rValues(end);
        else
            rBC = temp(1);
        end
        clear temp
    end
    
    % Define convergence condition to be satisfied
    convError = .001;
    conv = false;
        
    % Begin iteration of wavefuctions. Continued until either
    % convergence conditions are satisfied or loop limit is reached.
    
    % Set maximum number of loops and loop counter
    maxLoop = 25;
    loopCount = 0;
    
    % Begin iteration
    while (conv == false)&&(loopCount < maxLoop)
        
        % Increment number of loops
        loopCount = loopCount + 1;
        
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, ...
                                a0, E, kappa, m, convError);
        
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                                E, phase, kappa, m, convError);
        
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecFree);
        
        % Test for convergence
        conv = convTest(delU1, delU1P, delU2, delU2P, convError);
        
        % If solutions have converged, exit loop
        if (conv==true)
            % Unpack results from solvers
            u1Nucl = uVecNucl(:,1);
            u2Nucl = uVecNucl(:,2);
            u1Free = uVecFree(:,1);
            u2Free = uVecFree(:,2);
            % Merge results into complete wavefunctions
            u1 = [u1Nucl; u1Free(2:end)];
            u2 = [u2Nucl; u2Free(2:end)];
            clear uVecNucl uVecFree
            clear delU1 delU1P delU2 delU2P
            clear u1Nucl u1Free u2Nucl u2Free
            break
        end
        
        % If solutions failed to converge, find better values for
        % the parameters a0 and phase
        
        % Pass current parameters to function which calculates the
        % N-R matrix
        M = nrMatrix(a0, phase, rValues, rBC, delU1, delU2, ...
                    potential, rMatch, E, kappa, m, convError);
        % From N-R matrix, find new values for a0 and phase
        newVec = [a0 ; phase] - (M\[delU1 ; delU2]);
        a0New = newVec(1);
        phaseNew = newVec(2);
        
        % Check if the solution's amplitude is approaching zero
        if(a0New<=.01)
            a0Min = .01;
            if(a0Min==a0)
                % In this situation, a0 is already at the 
                % minimum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the a0
                % component
                a0New = a0;
            else
                % Rescale gradient step so values are physical
                rescale = abs((a0Min - a0)/(a0New-a0));
                a0New = a0Min;
                phaseNew = phase + rescale*(phaseNew - phase);
            end
            clear rescale a0Min
        end
        % Check if the solution's amplitude is growing without
        % bound
        if(a0New>=3)
            a0Max = 3;
            if(a0Max==a0)
                % In this situation, a0 is already at the 
                % maximum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the a0
                % component
                a0New = a0;
            else
                % Rescale gradient step so values are physical
                rescale = abs((a0Max - a0)/(a0New-a0));
                a0New = a0Max;
                phaseNew = phase + rescale*(phaseNew - phase);
            end
            clear rescale a0Min
        end
        % Update a0
        a0 = a0New;
        
        % Update phase and check that it is in 
        % the interval [0, 2*pi)
        phase = mod(phaseNew, 2*pi);

        clear delU1 delU2 M newVec
        clear delU1P delU2P
        % Repeat loop
        
    end



% At this point, the solution to the radial wave equation
% has either converged to its final value or the method
% has failed
    
    % Check to see if convergence failed
    if(conv==false)
        % If the solution failed to converge, throw an error
        error('freeElectronWaveFunction:failureToConverge', ...
                'After %d iterations, the %s %s %s %d .', ...
                loopCount, ...
                'radial wavefunctions and their slopes', ...
                ' failed to converge to less than a relative', ...
                ' difference of', convError);
    end
    
    % In this case, the solution convereged and the variables u1
    % and u2 are defined
    
    % Convert rValues to units of Fermi.
    rValues = rValues*((197.327)/mMuon);
    % Convert reduced wavefunction to units of (MeV*Fermi)^(-1/2)
    u1 = u1/sqrt(197.327);
    u2 = u2/sqrt(197.327);
    
%     % Display resulting wavefunctions
%     figure; 
%     plot(rValues, u1, '-r', rValues, u2, '-b');
%     hold on
%     title('Free Reduced Wavefunctions');
%     xlabel('Radius [Fermi]');
%     ylabel('Wavefunction [Fermi^{(-1/2)}*MeV^{(-1/2)}');
%     TeXString1 = texlabel('U_{1}');
%     TeXString2 = texlabel('U_{2}');
%     legend(TeXString1, TeXString2);
%     hold off
%     clear TeXString1 TeXString2
    
end

% Solves for nuclear part of the wavefunction from rMin to rMatch
function uVecNucl = nuclearSolver(rValues, potential, ...
                        rMatch, a0, E, kappa, Me, convError)
    % Find initial conditions for U1 and U2
    vMin = potEval(potential, rValues, rValues(1));
    B = E - vMin;
    if(kappa==-1)
        % Initial conditions for kappa = -1
        U2i = -a0*sqrt(B^2-Me^2)*rValues(1)*...
              ( ( sin(sqrt(B^2-Me^2)*rValues(1))/...
                    ((sqrt(B^2-Me^2)*rValues(1))^2) ) - ...
                ( cos(sqrt(B^2-Me^2)*rValues(1))/...
                    (sqrt(B^2-Me^2)*rValues(1)) )  );
        U1i = a0*sqrt((B+Me)/(B-Me))*...
                    sqrt(B^2-Me^2)*rValues(1)*...
                  ( sin(sqrt(B^2-Me^2)*rValues(1))/...
                    (sqrt(B^2-Me^2)*rValues(1)) );        
    else
        % Initial conditions for kappa = 1
        U2i = -a0*sqrt(B^2-Me^2)*rValues(1)*...
              ( sin(sqrt(B^2-Me^2)*rValues(1))/...
                (sqrt(B^2-Me^2)*rValues(1)) );
        U1i = -a0*sqrt((B+Me)/(B-Me))*...
               sqrt(B^2-Me^2)*rValues(1)*...
              ( ( sin(sqrt(B^2-Me^2)*rValues(1))/...
                    ((sqrt(B^2-Me^2)*rValues(1))^2) ) - ...
                ( cos(sqrt(B^2-Me^2)*rValues(1))/...
                    (sqrt(B^2-Me^2)*rValues(1)) )  );
    end
    clear vMin X B

    % Create ODE matrix
    odeMatrix = @(r,uVec) [ -kappa/r , ...
            (E - potEval(potential, rValues, r) + Me); ...
            -(E - potEval(potential, rValues, r) - Me), ...
               kappa/r]*uVec;
    
    % Set convergence options
    options = odeset('RelTol', convError/100, 'AbsTol', 10^(-12));
    
    % Try to solve system
    rValuesNucl = rValues(rValues<=rMatch);
    [~, uVecNucl] = ode45(odeMatrix , rValuesNucl, ...
                                            [U1i; U2i], options);
    clear U1i U2i odeMatrix rValuesNucl options
end

% Solves for free part of the wavefunction from rMax to rMatch
function uVecFree = freeSolver(rValues, potential, rMatch, ...
                        rBC, E, phase, kappa, Me, convError)

    % Separate out the regions where analytic solutions are viable
    % and where numerical solutions are needed
    rValuesAnalytic = rValues(rValues>=rBC);
    rValuesNum = rValues((rBC>rValues)&(rValues>=rMatch));
    
    % Find analytic solutions on the appropriate region
    if( kappa == -1)
        % Final conditions for kappa = -1
        U2analytic = -sqrt(2/(1+Me/E))*sqrt(E^2-Me^2)*...
                    rValuesAnalytic.*...
              ( ( sin(sqrt(E^2-Me^2)*rValuesAnalytic+phase)./...
                    ((sqrt(E^2-Me^2)*rValuesAnalytic).^2) ) - ...
                ( cos(sqrt(E^2-Me^2)*rValuesAnalytic+phase)./...
                    (sqrt(E^2-Me^2)*rValuesAnalytic) )  );
        U1analytic = sqrt(2/(1+Me/E))*sqrt((E+Me)/(E-Me))*...
                    sqrt(E^2-Me^2)*rValuesAnalytic.*...
                  ( sin(sqrt(E^2-Me^2)*rValuesAnalytic+phase)./...
                    (sqrt(E^2-Me^2)*rValuesAnalytic) );
    else
        % Final conditions for kappa = 1
        U2analytic = -sqrt(2/(1+Me/E))*...
                sqrt(E^2-Me^2)*rValuesAnalytic.*...
              ( sin(sqrt(E^2-Me^2)*rValuesAnalytic+phase)./...
                (sqrt(E^2-Me^2)*rValuesAnalytic) );
        U1analytic = -sqrt(2/(1+Me/E))*sqrt((E+Me)/(E-Me))*...
               sqrt(E^2-Me^2)*rValuesAnalytic.*...
              ( ( sin(sqrt(E^2-Me^2)*rValuesAnalytic+phase)./...
                    ((sqrt(E^2-Me^2)*rValuesAnalytic).^2) ) - ...
                ( cos(sqrt(E^2-Me^2)*rValuesAnalytic+phase)./...
                    (sqrt(E^2-Me^2)*rValuesAnalytic) )  );
    end

    % Check if the numerical region is actually null. If this is
    % the case, we only need the analytical solutions.
    if(isempty(rValuesNum))
            uVecFree = [U1analytic, U2analytic];
            return
    end
    
    % Remove the innermost value of the analytic solution and
    % use it as a boundary condition for the numerical solution
    U1f = U1analytic(1);
    U2f = U2analytic(1);
    U1analytic = U1analytic(2:end);
    U2analytic = U2analytic(2:end);
    rValuesNum = [rValuesNum; rValuesAnalytic(1)];
    
    % Unlike the Nuclear part of the wavefunction, the Free
    % part must be solved by integrating down from an initial
    % value to lower values of r. As such, we modify our
    % expressions so the ODE solver can still be used.
    
    % Create ODE matrix which is a function of rBar = rMax - r
    rMax = rValuesNum(end);
    odeMatrix = @(r,uVec) -[ -kappa/(rMax - r) , ...
       (E-potEval(potential, rValues,(rMax-r))+Me); ...
       -(E-potEval(potential, rValues, (rMax - r))-Me), ...
               kappa/(rMax - r)]*uVec;
    
    % Find radial values in terms of rBar
    rBarValues = sort(rMax - rValuesNum);
    clear rValuesCoul
    
    % Set convergence options
    options = odeset('RelTol', convError/100, 'AbsTol', 10^(-12));
    
    % Solve the system in terms of rBar
    [~, uVecRBar] = ode45(odeMatrix , rBarValues, ...
                                            [U1f; U2f], options);
    clear odeMatrix U1f U2f options
    
    % Convert back to regular radial values
    newRValues = rMax-rBarValues;
    clear rMax rBarValues
    
    % These values need to be sorted from rMatch to rMax and the
    % corresponding values of uVecRBar also must be changed
    uVecRBar1 = uVecRBar(:,1);
    uVecRBar2 = uVecRBar(:,2);
    [~, ind] = sort(newRValues);
    uVec1Num = uVecRBar1(ind);
    uVec2Num = uVecRBar2(ind);
    clear ind uVecRBar1 uVecRBar2 newRValues uVecRBar
    
    % Combine the numerical and analytic wavefunctions
    uVec1 = [uVec1Num; U1analytic];
    uVec2 = [uVec2Num; U2analytic];
    
    % Return results
    uVecFree = [uVec1, uVec2];
    clear uVec1 uVec2
end

% Calculates the mismatch of Nuclear and Free wavefunctions 
% and their slopes at rMatch
function [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecFree)
    % Unpack results from solvers
    u1Nucl = uVecNucl(:,1);
    u2Nucl = uVecNucl(:,2);
    u1Coul = uVecFree(:,1);
    u2Coul = uVecFree(:,2);
    rValuesNucl = rValues(rValues<=rMatch);
    rValuesCoul = rValues(rValues>=rMatch);
    clear uVecNucl uVecCoul
        
    % Calculate mismatch at rMatch
    u1diff = u1Coul(1) - u1Nucl(end);
    u1ave = (u1Coul(1)+u1Nucl(end))/2;
    u2diff = u2Coul(1) - u2Nucl(end);
    u2ave = (u2Coul(1)+u2Nucl(end))/2;
    delU1 = u1diff / u1ave;
    delU2 = u2diff / u2ave;
    clear u1diff u1ave u2diff u2ave
    
    % Calculate first derivatives at half step from rMatch
    u1PrimeNucl = (u1Nucl(end) - u1Nucl(end-1))/...
                    (rValuesNucl(end) - rValuesNucl(end-1));
    u2PrimeNucl = (u2Nucl(end) - u2Nucl(end-1))/...
                    (rValuesNucl(end) - rValuesNucl(end-1));
    u1PrimeCoul = (u1Coul(2) - u1Coul(1))/...
                    (rValuesCoul(2) - rValuesCoul(1));
    u2PrimeCoul = (u2Coul(2) - u2Coul(1))/...
                    (rValuesCoul(2) - rValuesCoul(1));

    % Calculate second derivatives at full step from rMatch
    u1DPrimeNucl = (u1Nucl(end)-2*u1Nucl(end-1)+u1Nucl(end-2))/...
                    (rValuesNucl(end)-rValuesNucl(end-1))^2;
    u2DPrimeNucl = (u2Nucl(end)-2*u2Nucl(end-1)+u2Nucl(end-2))/...
                    (rValuesNucl(end)-rValuesNucl(end-1))^2;
    u1DPrimeCoul = (u1Coul(3)-2*u1Coul(2)+u1Coul(1))/...
                    (rValuesCoul(2)-rValuesCoul(1))^2;
    u2DPrimeCoul = (u2Coul(3)-2*u2Coul(2)+u2Coul(1))/...
                    (rValuesCoul(2)-rValuesCoul(1))^2;
                
    % Calculate first derivatives at rMatch using corrections
    % from second derivatives
    u1PNucl = u1PrimeNucl + ...
        u1DPrimeNucl*(rValuesNucl(end) - rValuesNucl(end-1))/2;
    u2PNucl = u2PrimeNucl + ...
        u2DPrimeNucl*(rValuesNucl(end) - rValuesNucl(end-1))/2;
    u1PCoul = u1PrimeCoul + ...
        u1DPrimeCoul*(rValuesCoul(1) - rValuesCoul(2))/2;
    u2PCoul = u2PrimeCoul + ...
        u2DPrimeCoul*(rValuesCoul(1) - rValuesCoul(2))/2;
    clear u1PrimeNucl u1DPrimeNucl u2PrimeNucl u2DPrimeNucl
    clear u1PrimeCoul u1DPrimeCoul u2PrimeCoul u2DPrimeCoul
    clear rValuesNucl rValuesCoul
    
    % Calculate differences and averages of first derivative
    u1Pdiff = u1PNucl-u1PCoul;
    u1Pave = (u1PNucl+u1PCoul)/2;
    u2Pdiff = (u2PNucl-u2PCoul);
    u2Pave = (u2PNucl+u2PCoul)/2;
    clear u1PNucl u1PCoul p2PNucl p2PCoul
    
    % Find mismatches relative to average
    delU1P = u1Pdiff/u1Pave;
    delU2P = u2Pdiff/u2Pave;
    clear u1Pdiff u1Pave u2Pdiff u2Pave
end

% Tests if radial wave functions and their slopes have converged
function conv = convTest(delU1, delU1P, delU2, delU2P, ...
                                                    convError)
    differences = abs([delU1; delU1P; delU2; delU2P]);
    % See which differences are less than the maximum error
    logicals = (differences<=convError);
    % Test to see if all of the logicals are satisfied
    conv = all(logicals);
    clear differences logicals
end

% Calculates the matrix needed for the Newton-Raphson method
function M = nrMatrix(a0, phase, rValues, rBC, delU1, delU2, ...
                potential, rMatch, E, kappa, Me, convError)
    % The elements of this matrix are derivatives of delU1 or DelU2
    % with respect to changes in the arguements a0 or bGuess. In
    % order to approximate the derivatives, we will vary the
    % parameters by +/- .01%. This will result in two slopes for
    % each parameter, which we will average to estimate the local
    % slope.
    
    shift = .0001;
    % Find values of parameters to test
    a0min = a0*(1-shift);
    a0max = a0*(1+shift);
    phaseMin = phase*(1-shift);
    phaseMax = phase*(1+shift);
    
    % Find results for a0min
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, ...
                                a0min, E, kappa, Me, convError);
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                        E, phase, kappa, Me, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1a0min, ~, delU2a0min, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecFree);
    
    % Find results for a0max
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, ...
                                a0max, E, kappa, Me, convError);
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                        E, phase, kappa, Me, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1a0max, ~, delU2a0max, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecFree);
                                        
    % Find results for phaseMin
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0, ...
                                            E, kappa, Me, convError);
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                        E, phaseMin, kappa, Me, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1phaseMin, ~, delU2phaseMin, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecFree);
                                        
    % Find results for phaseMax
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0, ...
                                        E, kappa, Me, convError);
        % Find free wavefunction from solver
        uVecFree=freeSolver(rValues, potential, rMatch, rBC, ...
                        E, phaseMax, kappa, Me, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1phaseMax, ~, delU2phaseMax, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecFree);
          
                                        
    % Find slopes for each change in parameters
    dDelU1da0max = (delU1a0max - delU1)/(a0max - a0);
    dDelU2da0max = (delU2a0max - delU2)/(a0max - a0);
    dDelU1da0min = (delU1 - delU1a0min)/(a0 - a0min);
    dDelU2da0min = (delU2 - delU2a0min)/(a0 - a0min);
    dDelU1dphaseMax = (delU1phaseMax - delU1)/(phaseMax - phase);
    dDelU2dphaseMax = (delU2phaseMax - delU2)/(phaseMax - phase);
    dDelU1dphaseMin = (delU1 - delU1phaseMin)/(phase - phaseMin);
    dDelU2dphaseMin = (delU2 - delU2phaseMin)/(phase - phaseMin);
    
    % Average slopes with respect to the same parameters to get 
    % the local derivative
    dDelU1da0 = (dDelU1da0max + dDelU1da0min)/2;
    dDelU2da0 = (dDelU2da0max + dDelU2da0min)/2;
    dDelU1dphase = (dDelU1dphaseMax + dDelU1dphaseMin)/2;
    dDelU2dphase = (dDelU2dphaseMax + dDelU2dphaseMin)/2;
    
    % Return the N-R matrix
    M = [dDelU1da0, dDelU1dphase; dDelU2da0, dDelU2dphase];
end

% Extracts and scales input parameters from paramVec
function params = paramExtract(m, n, argVec)
    % Extract input parameters for charge distribution
    
    if (n==1)
        % Harmonic Oscillator model
        nArgs = size(argVec, 2);
        if(nArgs~=2)
            % Harmonic oscillator must have two parameters
            error('incorrect number of parameters for HO');
        end
        al = argVec(1);
        if(al<0)
            % Alpha parameter cannot be negative
            error('Alpha parameter for HO cannot be negative');
        end
        a = argVec(2);
        if(a<=0)
            % A parameter must be positive
            error('A parameter for HO must be positive');
        end
        % Convert "a" from fermi to dimensionless mass units
        a = a*(m/(197.327));
        params = [al;a];
        clear nArgs al a argCell
        
    elseif (n==2)||(n==3)||(n==4)||(n==5)
        
        % 2PF, 3PF, 2PG, or 3PG distribution
        nArgs = size(argVec, 2);
        if ((n==2)||(n==4))&&(nArgs~=2)
            % 2PF and 2PG need two parameters
            error('incorrect number of parameters for 2PF or 2PG');
        end
        if ((n==3)||(n==5))&&(nArgs~=3)
            % 3PF and 3PG need three parameters
            error('incorrect number of parameters for 3PF or 3PG');
        end
        params = zeros(nArgs, 1);
        for i = 1:nArgs
            params(i) = argVec(i);
        end
        if (params(1)<0)||(params(2)<=0)
            % C must be non-negative and Y must be positive
            error('must have C>=0 and Y>0 for NPF or NPG');
        end
        % Convert C, Y from fermi to dimensionless mass units
        params(1:2) = params(1:2)*(m/(197.327));
        clear nArgs i argCell
        
    elseif (n==6)
        
        % Sum of Gaussian distribution
        nArgs = size(argVec, 2);
        nPairs = floor((nArgs-1)/2);
        if(nPairs<=0)
            % At least one Gaussian must be parameterized
            error('SOG has no distribution parameter pairs');
        end
        if(nArgs~=(2*nPairs+1))
            % All R and Q parameters must come in pairs
            error('SOG parameters not paired correctly');
        end
        
        % Extract R-parameters
        rParams = zeros(nPairs, 1);
        for i=1:nPairs
            rParams(i) = argVec(i);
        end
        if any(rParams<0)
            % R parameters cannot be negative
            error('SOG R parameters must be non-negative');
        end
        
        % Extract Q-parameters
        qParams = zeros(nPairs, 1);
        for i=1:nPairs
            qParams(i) = argVec(i+nPairs);
        end
        if any(qParams<0)
            % Q parameters cannot be negative
            error('SOG Q parameters must be non-negative');
        end
        if (sum(qParams)<.9999)||(sum(qParams)>1.0001)
            % Q parameters should sum to one (up to rounding errors)
            error('SOG Q parameters do not sum to one');
        end
        
        % Extract Gamma
        gamma = argVec(end);
        if gamma<=0
            % RMS radius must be positive
            error('SOG gamma parameter must be positive');
        end
        
        % Return params
        params = [rParams; qParams; gamma];
        clear nArgs nPairs rParams qParams gamma i argCell
        
    elseif (n==7)
        
        % Fourier-Bessel Series
        nArgs = size(argVec, 2);
        if(nArgs<=1)
            error('too few parameters for Fourier-Bessel');
        end
        params = zeros(nArgs, 1);
        for i = 1:nArgs
            params(i) = argVec(i);
        end
        if params(1)<=0
            % Cutoff radius must be positive
            error('FB cutoff radius must be positive');
        end
        % Convert parameters from fermi to dimensionless mass units
        params(1) = params(1)*(m/(197.327));
        params(2:end) = params(2:end)*(197.327/m)^3;
        clear numArgs i cellArgs
        
    else
        
        % No valid charge distribution was specified
        error('invalid charge distribution specified');
        
    end
end

% Finds the electric potential and charge distribution over rValues
function [potential, chargeDist] = chargePotential(alpha, ...
                                    Z, n, params, rValues)
   % Extend rValues to the origin to minimize error from
   % integrating charge distribution to find potential.
   % Choose to use spacing between the first two points
   dif = rValues(2) - rValues(1);
   numPoints = round(rValues(1)/dif);
   extension = (0:(numPoints-1))';
   rValues = [dif*extension; rValues];
   clear dif extension

   % Determine the charge distribution
   if (n==1)
       % Harmonic Oscillator Model
       % params = [al ; a]
       al = params(1);
       a = params(2);
       chargeDist = (1+al*(rValues/a).^2).*exp(-(rValues/a).^2);
       clear al a
       
   elseif (n==2)
       % Two Parameter Fermi Distribution
       % params = [c ; y]
       c = params(1);
       y = params(2);
       chargeDist = 1./(1+exp((rValues-c)/y));
       clear c y
       
   elseif (n==3)
       % Three Parameter Fermi Distribution
       % params = [c ; y ; w]
       c = params(1);
       y = params(2);
       w = params(3);
       chargeDist = (1+w*(rValues/c).^2)./(1+exp((rValues-c)/y));
       clear c y w
       
   elseif (n==4)
       % Two Parameter Gaussian Distribution
       % params = [c ; y]
       c = params(1);
       y = params(2);
       chargeDist = 1./(1+exp((rValues.^2-c^2)/(y^2)));
       clear c y
       
   elseif (n==5)
       % Three Parameter Gaussian Distribution
       % params = [c ; y ; w]
       c = params(1);
       y = params(2);
       w = params(3);
       chargeDist = (1+w*(rValues/c).^2)./ ...
                                   (1+exp((rValues.^2-c^2)/(y^2)));
       clear c y w
       
   elseif (n==6)
       % Sum of Gaussian Distribution
       % params = [R1 ; R2 ; ... Q1 ; Q2 ; ... Gamma]
       
       nParams = size(params, 1);
       gamma = params(end);
       nComp = (nParams-1)/2;
       chargeDistComps = zeros( size(rValues, 1), nComp );
       
       % Find each component of the charge distribution
       for i=1:nComp
           Ri = params(i);
           Qi = params(i+nComp);
           Ai = Qi/((2*(pi^(3/2)))*(gamma^3)*(1+2*(Ri/gamma)^2));
           chargeDistComps(:,i) = ...
                       Ai*(exp(-((rValues-Ri)/gamma).^2) + ...
                                   exp(-((rValues+Ri)/gamma).^2));
       end
       
       % Sum the components
       chargeDist = sum(chargeDistComps, 2);
       
       clear nParams gamma nComp chargeDistComps Ri Qi Ai i
       
   elseif (n==7)
       % Fourier-Bessel Series
       % params = [R ; a1 ; a2 ; ...]
       
       R = params(1);
       nComp = size(params,1) - 1;
       chargeDistComps = zeros( size(rValues, 1), nComp );
       
       % Find each component of the charge distribution
       for i=1:nComp
           Ai = params(1+i);
           chargeDistComps(:,i) = ...
                       Ai*sin(i*pi*rValues/R)./(i*pi*rValues/R);
           % Remove Singularity at Origin
           chargeDistComps(1,i) = Ai;
       end
       
       % Sum the components
       chargeDist = sum(chargeDistComps, 2);
       
       % Impose cut-off radius R
       ind = (rValues>R);
       chargeDist(ind) = 0;
       
       clear R nComp chargeDistComps Ai i ind
       
   else
       
       % No valid charge distribution was specified
       error('invalid charge distribution specified');
       
   end   
   
   
   % Normalize the charge distribution to unity   
   norm = trapz(rValues, chargeDist.*(4*pi*rValues.^2));
   chargeDist = (1/norm)*chargeDist;
   clear n params norm
   
   % Find electric field everywhere by integrating charge
   % distribution
   elecField = (rValues.^(-2)).*cumtrapz(rValues, ...
                                        chargeDist.*(rValues.^2));
   % Account for artificial singularity at origin
   elecField(1) = 0;

   % Find potential everywhere by integrating electric field
   potential = -4*pi*Z*alpha*...
                   (trapz(rValues, elecField)+...
                    elecField(end)*rValues(end)-...
                    cumtrapz(rValues, elecField));
   clear elecField
   
   % Now remove extra points added near origin
   chargeDist = chargeDist((numPoints+1):end);
   potential = potential((numPoints+1):end);
   clear numPoints
end

% Finds the nuclear charge radius
function rCharge = chargeRadius(chargeDist, rValues, frac)
    % Find total charge as a function of radius
    chargeContent = cumtrapz(rValues, 4*pi*(rValues.^2).*chargeDist);
    % Determine indices of values which contain at least the
    % specified fraction of the charge
    ind = (chargeContent >= frac);
    temp = rValues(ind);
    % Now select the smallest rValue satisfying this condition
    rCharge = temp(1);
end

% Acts as short-cut for evaluating the potential
function potentialR = potEval(potential, rValues, r)
% This function is designed so that the ODEsolver can call this
% instead of re-evaluating all the integrals necessary to calculate
% the potential at each step.

    % Determine if called r is in rValues, and if so find the index
    ind = find(rValues == r , 1);
    if(isempty(ind))
        % If this value of r is not in rValues, need to 
        % estimate potential at r
        % Find points in rValues just above and below r
        indLess = (rValues < r);
        indGreater = (rValues > r);
        indLess = find(indLess, 1, 'last');
        indGreater = find(indGreater, 1, 'first');
        rLower = rValues(indLess);
        rGreater = rValues(indGreater);
        % Find potential at these points
        pLower = potential(indLess);
        pGreater = potential(indGreater);
        % Linearly interpolate between the two points
        slope = (pGreater - pLower)./(rGreater - rLower);
        potentialR = pLower + slope.*(r - rLower);        
    else
        % If this value of r is in rValues, just look up the
        % potential
        potentialR = potential(ind);
    end
end
