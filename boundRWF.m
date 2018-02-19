function [rValues, u1, u2, bindingE] = ...
                                    boundRWF(m, Z, n, paramVec)
%boundRWF Calculates dirac wavefunction for bound muon
%
%   Returns binding energy, radial values, and the 
%   reduced spinor wavefunctions for the system.
%   Returned rValues are in units of [Fermi]^1
%   Returned u1, u2 are in units of [Fermi]^(-1/2)
%   Returned binding energy is in units of MeV
%
%   Inputs:
%   m = mass of bound particle in MeV (105.65837 for muon)
%   Z = atomic number
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


% Initialize relevant constants, guesses, and the radial values
% Note: All calculations are done in dimensionless mass units

    % Fine Structure Constant
    alpha = 1/137.036;
    % Kappa from Dirac equation
    kappa = -1;
    % Bohr Radius
    rBohr = (1/(Z*alpha));
    % Extract and rescale distribution parameters from input
    params = paramExtract(m, n, paramVec);

    % Find charge radius of the nucleus
    % Initialize temporary radial values
    numPoints = 5000;
    rMin = rBohr/10000;
    rMax = rBohr*100;
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
    
    % If the nuclear radius is greater than the muonic Bohr
    % radius, we cannot trust this analysis
%    if(rNucleus>rBohr)
%        error('analysis is invalid');
%    end
    
    % Other values of interest (in dimensionless mass units)
    rMax = 10*max([rNucleus, rBohr]);
    rMin = rNucleus/10000;
    clear rNucleus



% FIRST, find the best solution for a low number of radial points
% and low convergence criterion
    
    % Initialize range of radial wavefunction
    numPoints = 10000;
    logMin = log10(rMin/rBohr);
    logMax = log10(rMax/rBohr);
    scaledRValues = logspace(logMin, logMax, numPoints)';
    rValues = rBohr*scaledRValues;
    clear logMin logMax scaledRValues numPoints

    % Find electric potential
    [potential, ~] = chargePotential(alpha, Z, n, params, rValues);

    % Find rMatch near rBohr for Shoot-and-Match procedure
    temp = rValues(rValues >= rBohr);
    rMatch = temp(1);
    clear temp
    
    % Define convergence condition to be satisfied
    convError = .01;
    conv = false;
    
    % Initialize guess for bound state energy. Use value for
    % an idealized hydrogen-like atom unless this binding energy
    % is unphysical.
    bHydro = (1/2)*(Z*alpha)^2;
    vMin = abs(potential(1));
    if (bHydro>vMin)
        bHydro = vMin*.9;
        bGuess = bHydro;
    else
        % Scaling factor has been added for faster convergence
        % at large Z.
        bGuess = bHydro*exp(-Z*alpha);
    end
    clear vMin
    
    % Set initial value for a0 in nuclear wavefunction
    % This value was chosen to roughly match the final value
    % of a0 for Aluminum
    a0 = .003*exp(-Z*alpha);

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
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0, ...
                                bGuess, kappa, convError);

        % Find Coulombic wavefunction from solver
        uVecCoul=coulombicSolver(rValues, potential, rMatch, ...
                                bGuess, kappa, convError);

        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecCoul);

        % Test for convergence
        conv = convTest(delU1, delU1P, delU2, delU2P, convError);

        % If solutions have converged, exit loop
        if (conv==true)
            clear conv delU1 delU1P delU2 delU2P uVecCoul uVecNucl
            break
        end

        % If solutions failed to converge, find better values for
        % the parameters a0 and bGuess

        % Pass current parameters to function which calculates the
        % Newton-Raphson matrix
        M = nrMatrix(a0, bGuess, rValues, delU1, delU2, ...
                        potential, rMatch, kappa, convError);
        % From N-R matrix, find new values for a0 and bGuess
        newVec = [a0 ; bGuess] - M\[delU1 ; delU2];
        a0New = newVec(1);
        bGuessNew = newVec(2);
        
        % Check that these new values correspond to physical
        % solutions. If they do not, rescale the descent vector
        % in a0-bGuess space such that the values are physical
        vMin = abs(potential(1));
        % Check if binding energy is greater than the potential
        if(bGuessNew > vMin)
            bGuessMax = .9*vMin;
            if(bGuess == bGuessMax)
                % In this situation, bGuess is already at the 
                % maximum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the bGuess
                % component
                bGuessNew = bGuessMax;
            else
                % Rescale gradient step so values are physical
                rescale = ...
                        abs((bGuessMax-bGuess)/(bGuessNew-bGuess));
                bGuessNew = bGuessMax;
                a0New = a0 + rescale*(a0New - a0);
            end
            clear rescale bGuessMax
        end
        % Check if binding energy is too close to the first
        % excited state
        if(bGuessNew <= bHydro*(1/3)*exp(-Z*alpha))
            bGuessMin = bHydro*(1/3)*exp(-Z*alpha);
            if(bGuessMin==bGuess)
                % In this situation, bGuess is already at the 
                % minimum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the bGuess
                % component
                bGuessNew = bGuessMin;
            else
            % Rescale gradient step so values are physical
                rescale = ...
                    abs((bGuessMin-bGuess)/(bGuessNew-bGuess));
                bGuessNew = bGuessMin;
                a0New = a0 + rescale*(a0New - a0);
            end
            clear rescale bGuessMin
        end
        % Check if the solution's amplitude is approaching zero
        if(a0New<=.00001)
            a0Min = .00001;
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
                bGuessNew = bGuess + rescale*(bGuessNew - bGuess);
            end
            clear rescale a0Min
        end

        % Update parameters
        a0 = a0New;
        bGuess = bGuessNew;
        clear a0New bGuessNew

        % Clear all temporary variables
        clear vMin delU1 delU2 M newVec uVecNucl uVecCoul
        clear delU1P delU2P
        
        % Repeat loop
    end

    


% At this point, the rough solutions to the radial wave equation
% have provided all the information they can. Now, with rough 
% values for the binding energy and scaling factor, we can 
% solve the system in greater detail
    
    % Clear old temporary variables
    clear maxLoop loopCount potential rMatch
    
    % Initialize range of radial wavefunction
    numPoints = 100000;
    logMin = log10(rMin/rBohr);
    logMax = log10(rMax/rBohr);
    scaledRValues = logspace(logMin, logMax, numPoints)';
    rValues = rBohr*scaledRValues;
    clear logMin logMax scaledRValues rMax rMin numPoints
    
    % Find electric potential everywhere
    [potential, ~] = chargePotential(alpha, Z, n, params, rValues);
    
    % Find rMatch near rBohr for Shoot-and-Match procedure
    temp = rValues(rValues >= rBohr);
    rMatch = temp(1);
    clear temp
    
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
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0, ...
                                bGuess, kappa, convError);
                            
        % Find Coulombic wavefunction from solver
        uVecCoul=coulombicSolver(rValues, potential, rMatch, ...
                                bGuess, kappa, convError);
        
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecCoul);
        
        % Test for convergence
        conv = convTest(delU1, delU1P, delU2, delU2P, convError);
        
        % If solutions have converged, exit loop
        if (conv==true)
            % Unpack results from solvers
            u1Nucl = uVecNucl(:,1);
            u2Nucl = uVecNucl(:,2);
            u1Coul = uVecCoul(:,1);
            u2Coul = uVecCoul(:,2);
            % Merge results
            u1 = [u1Nucl; u1Coul(2:end)];
            u2 = [u2Nucl; u2Coul(2:end)];
            bindingE = bGuess;
            clear uVecNucl uVecCoul
            clear delU1 delU1P delU2 delU2P
            clear u1Nucl u1Coul u2Nucl u2Coul
            clear loopCount maxLoop
            break
        end
        
        % If solutions failed to converge, find better values for
        % the parameters a0 and bGuess
        
        % Pass current parameters to function which calculates the
        % Newton-Raphson matrix
        M = nrMatrix(a0, bGuess, rValues, delU1, delU2, ...
                        potential, rMatch, kappa, convError);
        % From N-R matrix, find new values for a0 and B       
        newVec = [a0 ; bGuess] - M\[delU1 ; delU2];
        a0New = newVec(1);
        bGuessNew = newVec(2);
        
        % Check that these new values correspond to physical
        % solutions. If they do not, rescale the descent vector
        % in a0-bGuess space such that the values are physical
        vMin = abs(potential(1));
        % Check if binding energy is greater than the potential
        if(bGuessNew > vMin)
            bGuessMax = .9*vMin;
            if(bGuess == bGuessMax)
                % In this situation, bGuess is already at the 
                % maximum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the bGuess
                % component
                bGuessNew = bGuessMax;
            else
                % Rescale gradient step so values are physical
                rescale = ...
                        abs((bGuessMax-bGuess)/(bGuessNew-bGuess));
                bGuessNew = bGuessMax;
                a0New = a0 + rescale*(a0New - a0);
            end
            clear rescale bGuessMax
        end
        % Check if binding energy is too close to the first
        % excited state
        if(bGuessNew <= bHydro*(1/3)*exp(-Z*alpha))
            bGuessMin = bHydro*(1/3)*exp(-Z*alpha);
            if(bGuessMin==bGuess)
                % In this situation, bGuess is already at the 
                % minimum value. To avoid rescaling the descent
                % vector to zero, we just eliminate the bGuess
                % component
                bGuessNew = bGuessMin;
            else
            % Rescale gradient step so values are physical
                rescale = ...
                    abs((bGuessMin-bGuess)/(bGuessNew-bGuess));
                bGuessNew = bGuessMin;
                a0New = a0 + rescale*(a0New - a0);
            end
            clear rescale bGuessMin
        end
        % Check if the solution's amplitude is approaching zero
        if(a0New<=.00001)
            a0Min = .00001;
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
                bGuessNew = bGuess + rescale*(bGuessNew - bGuess);
            end
            clear rescale a0Min
        end
        
        % Update Parameters
        a0 = a0New;
        bGuess = bGuessNew;
        clear a0New bGuessNew

        % Clear Temporary Variables
        clear vMin delU1 delU2 M newVec uVecNucl uVecCoul
        clear delU1P delU2P
        
        % Repeat loop
    end



% At this point, the solution to the radial wave equation
% has either converged to its final value or the method
% has failed
    
    % Check to see if convergence failed
    if(conv==false)
        % If the solution failed to converge, throw an error
        error('After %d iterations, the %s %s %s %d .', ...
                loopCount, ...
                'radial wavefunctions and their slopes', ...
                ' failed to converge to a relative', ...
                ' difference of less than', convError);
    end
    
    % Otherwise, the solution convereged and the variables u1,
    % u2, and bindingE have been defined
    
    % Convert binding energy to units of MeV
    bindingE = m*bindingE;
    
    % Convert rValues to units of Fermi
    rValues = rValues*((197.327)/m);
    
    % Normalize the wavefunctions in units of Fermi^(-1/2) and 
    % return the results
    uVec = normalizer(rValues, u1, u2);
    u1 = uVec(:,1);
    u2 = uVec(:,2);
    clear uVec
    
%      % Display resulting wavefunctions
%      figure; 
%      plot(rValues, u1, '-b', rValues, u2, '-g');
%      hold on
%      title('Bound Reduced Wavefunctions');
%      xlabel('Radius [Fermi]');
%      TeXString0 = texlabel('Wavefunction [Fermi^(-1/2)]');
%      ylabel(TeXString0);
%      TeXString1 = texlabel('U_{1}');
%      TeXString2 = texlabel('U_{2}');
%      legend(TeXString1, TeXString2);
%      hold off
%      clear TeXString0 TeXString1 TeXString2
    
end

% Solves for nuclear part of the wavefunction from rMin to rMatch
function uVecNucl = nuclearSolver(rValues, potential, rMatch, ... 
                                a0, bGuess, kappa, convError)
    % Find Initial conditions for U1 and U2
    U1i = a0*(rValues(1));
    vMin = potEval(potential, rValues, rValues(1));
    U2i = ((bGuess + vMin)/(2 - kappa))*a0*(rValues(1))^2;
    clear vMin
    
    % Create ODE matrix
    odeMatrix = @(r,uVec) [ -kappa/r , ...
            (2-potEval(potential, rValues, r)-bGuess); ...
            (bGuess + potEval(potential, rValues, r)), ...
               kappa/r]*uVec;
    
    % Set convergence options
    options = odeset('RelTol', convError/100, 'AbsTol', 10^(-12));
           
    %Try to solve system
    rValuesNucl = rValues(rValues<=rMatch);
    [~, uVecNucl] = ode45(odeMatrix , rValuesNucl, ...
                                [U1i; U2i], options);
    clear U1i U2i odeMatrix rValuesNucl options
end

% Solves for Coulombic part of the wavefunction from rMax to rMatch
function uVecCoul = coulombicSolver(rValues, potential, rMatch, ...
                                bGuess, kappa, convError)

    % Define rMax where the boundary conditions will be evaluted
    rMax = rValues(end);
    
    % Restrict rValues to the region of interest for the far field
    % wavefunctions
    rValuesCoul = rValues(rValues>=rMatch);
    
    % Find the boundary conditions
    U1f = exp(-sqrt((2-bGuess)*bGuess)*rMax);
    U2f = -sqrt(bGuess/(2-bGuess))*...
                    exp(-sqrt((2-bGuess)*bGuess)*rMax);
    
    % Unlike the Nuclear part of the wavefunction, the Coulombic
    % part must be solved by integrating down from an initial
    % value to lower values of r. As such, we modify our
    % expressions so the ODE solver can still be used.
    
    % Create ODE matrix which is a function of rBar = rMax - r
    odeMatrix = @(r,uVec) -[ -kappa/(rMax - r) , ...
       (2-potEval(potential, rValues,(rMax-r))-bGuess); ...
       (bGuess+potEval(potential, rValues, (rMax - r))), ...
               kappa/(rMax - r)]*uVec;
    
    % Set convergence options
    options = odeset('RelTol', convError/100, 'AbsTol', 10^(-12));
           
    % Find radial values in terms of rBar
    rBarValues = sort(rMax - rValuesCoul);
    clear rValuesCoul
    
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
    uVec1 = uVecRBar1(ind);
    uVec2 = uVecRBar2(ind);
    clear ind uVecRBar1 uVecRBar2 newRValues uVecRBar
    
    % Return results
    uVecCoul = [uVec1, uVec2];
    clear uVec1 uVec2
end

% Calculates the mismatch of Nuclear and Coulombic wavefunctions 
% and their slopes at rMatch
function [delU1, delU1P, delU2, delU2P] = mismatch(rValues, ...
                                    rMatch, uVecNucl, uVecCoul)
    % Unpack results from solvers
    u1Nucl = uVecNucl(:,1);
    u2Nucl = uVecNucl(:,2);
    u1Coul = uVecCoul(:,1);
    u2Coul = uVecCoul(:,2);
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
function conv = convTest(delU1, delU1P, ...
                                    delU2, delU2P, convError)
    differences = abs([delU1; delU1P; delU2; delU2P]);
    % See which differences are less than the maximum error
    logicals = (differences<=convError);
    % Test to see if all of the logicals are satisfied
    conv = all(logicals);
    clear differences logicals
end

% Calculates the matrix needed for the Newton-Raphson method
function M = nrMatrix(a0, bGuess, rValues, delU1, delU2, ...
                            potential, rMatch, kappa, convError)
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
    bGmin = bGuess*(1-shift);
    bGmax = bGuess*(1+shift);
    
    
    % Find results for a0min
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0min,...
                                bGuess, kappa, convError);
        % Find Coulombic wavefunction from solver
        uVecCoul=coulombicSolver(rValues, potential, rMatch, ...
                                bGuess, kappa, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1a0min, ~, delU2a0min, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecCoul);
    
    % Find results for a0max
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0max,...
                                bGuess, kappa, convError);
        % Find Coulombic wavefunction from solver
        uVecCoul=coulombicSolver(rValues, potential, rMatch, ...
                                bGuess, kappa, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1a0max, ~, delU2a0max, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecCoul);
                                        
    % Find results for bGmin
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0, ...
                                bGmin, kappa, convError);
        % Find Coulombic wavefunction from solver
        uVecCoul=coulombicSolver(rValues, potential, rMatch, ...
                                bGmin, kappa, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1bGmin, ~, delU2bGmin, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecCoul);
                                        
    % Find results for bGmax
    
        % Find nuclear wavefunction from solver
        uVecNucl=nuclearSolver(rValues, potential, rMatch, a0, ...
                                bGmax, kappa, convError);
        % Find Coulombic wavefunction from solver
        uVecCoul=coulombicSolver(rValues, potential, rMatch, ...
                                bGmax, kappa, convError);
        % Find mismatch of wavefunctions and their slopes at
        % rMatch
        [delU1bGmax, ~, delU2bGmax, ~] = mismatch(rValues, ...
                                        rMatch, uVecNucl, uVecCoul);
          
                                        
    % Find slopes for each change in parameters
    dDelU1da0max = (delU1a0max - delU1)/(a0max - a0);
    dDelU2da0max = (delU2a0max - delU2)/(a0max - a0);
    dDelU1da0min = (delU1 - delU1a0min)/(a0 - a0min);
    dDelU2da0min = (delU2 - delU2a0min)/(a0 - a0min);
    dDelU1dbGmax = (delU1bGmax - delU1)/(bGmax - bGuess);
    dDelU2dbGmax = (delU2bGmax - delU2)/(bGmax - bGuess);
    dDelU1dbGmin = (delU1 - delU1bGmin)/(bGuess - bGmin);
    dDelU2dbGmin = (delU2 - delU2bGmin)/(bGuess - bGmin);
    
    % Average slopes with respect to the same parameters to get 
    % the local derivative
    dDelU1da0 = (dDelU1da0max + dDelU1da0min)/2;
    dDelU2da0 = (dDelU2da0max + dDelU2da0min)/2;
    dDelU1dbG = (dDelU1dbGmax + dDelU1dbGmin)/2;
    dDelU2dbG = (dDelU2dbGmax + dDelU2dbGmin)/2;
    
    % Return the N-R matrix
    M = [dDelU1da0, dDelU1dbG; dDelU2da0, dDelU2dbG];
end

% Normalizes radial wavefunctions
function uVec = normalizer(rValues, u1, u2)
    % Calculate normalization constant such that
    % Integral(U1^2 + U2^2, x=0, x=Inf) = 1
    norm = trapz(rValues, (u1.^2 +u2.^2) );
    u1Normed = u1/sqrt(norm);
    u2Normed = u2/sqrt(norm);
    uVec = [u1Normed, u2Normed];
    clear norm u1Normed u2Normed
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
