function overlapIntegralsNeutron()
%overlapIntegralsNeutron Calculates overlap integrals for neutron
% distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fIDStat = fopen('status1.txt', 'w');

    % Load previous results
    load('preMonteCarlo.mat')
    
    % Neutron distribution parameters and 1-sigma uncertainties
    cParam = 3.18;
    cSD = .19;
    zParam = .535;
    
    % Number of samples to be taken
    nSample = 10000;
    
    % Generate random samples for cParam
    cRands = normrnd(cParam, cSD, nSample, 1);
    
    % Prepare vectors to store overlap integral results
    Isn1vec = zeros(nSample,1);
    Isn2vec = zeros(nSample,1);
    Isn3vec = zeros(nSample,1);
    Isn4vec = zeros(nSample,1);
    tIsn1vec = zeros(nSample,1);
    tIsn2vec = zeros(nSample,1);
    tIsn3vec = zeros(nSample,1);
    tIsn4vec = zeros(nSample,1);
    Ivn1vec = zeros(nSample,1);
    Ivn2vec = zeros(nSample,1);
    Ivn3vec = zeros(nSample,1);
    Ivn4vec = zeros(nSample,1);
    
    save('workspace1.mat');
    
    % Loop over all selected parameter values
    for ind=1:nSample
        
        % Find neutron wavefunction
        pVecNeutron = [cRands(ind),zParam];
        PsiNeutron = wavefunction(A-Z, 2, pVecNeutron, rValues);
        
        % Fourier transform wavefunction
        PsiFNeutron = fourierTrans(PsiNeutron, rValues, kValues);
        
        % Square to find fourier transformed density distribution
        fDensityN = PsiFNeutron.^2;
        
        % Find overlap integrals
        Isn1vec(ind) = nucOverlap(kValues, fDensityN, Zs1);
        Isn2vec(ind) = nucOverlap(kValues, fDensityN, Zs2);
        Isn3vec(ind) = nucOverlap(kValues, fDensityN, Zs3);
        Isn4vec(ind) = nucOverlap(kValues, fDensityN, Zs4);
        tIsn1vec(ind) = nucOverlapMomentum(kValues, fDensityN, ...
                                                        Zs1, Mpi);
        tIsn2vec(ind) = nucOverlapMomentum(kValues, fDensityN, ...
                                                        Zs2, Mpi);
        tIsn3vec(ind) = nucOverlapMomentum(kValues, fDensityN, ...
                                                        Zs3, Mpi);
        tIsn4vec(ind) = nucOverlapMomentum(kValues, fDensityN, ...
                                                        Zs4, Mpi);
        Ivn1vec(ind) = nucOverlap(kValues, fDensityN, Zv1);
        Ivn2vec(ind) = nucOverlap(kValues, fDensityN, Zv2);
        Ivn3vec(ind) = nucOverlap(kValues, fDensityN, Zv3);
        Ivn4vec(ind) = nucOverlap(kValues, fDensityN, Zv4); 
        
        % Update saved values and post status
        save('workspace1.mat', 'Isn1vec', 'Isn2vec', 'Isn3vec', ...
                'Isn4vec', 'tIsn1vec', 'tIsn2vec', 'tIsn3vec', ...
                'tIsn4vec', 'Ivn1vec', 'Ivn2vec', 'Ivn3vec', ...
                'Ivn4vec', '-append')
        fprintf(fIDStat, 'Current Loop: %d \n', ind);
    end
    
    fclose(fIDStat);
    
end

function psi1 = wavefunction(Nnuc, n, params, rValues)

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
       
       nParams = size(params, 2);
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
       nComp = size(params, 2) - 1;
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
   
   
   % Normalize the charge distribution to Nnuc
   norm = trapz(rValues, chargeDist.*(4*pi*rValues.^2));
   chargeDist = (Nnuc/norm)*chargeDist;
   clear n params norm
   
   % Define wavefunction
   psi1 = sqrt(chargeDist);

end

function fPsi = fourierTrans(rPsi, rValues, kValues)
    % Convert rValues from Fermi to MeV^(-1)
    rValuesM = rValues/(197.327);
    % Convert Psi from Fermi^(-3/2) to MeV^(3/2)
    rPsiM = rPsi*(197.327)^(3/2);
    % Initialize fPsi
    fPsi = zeros(size(kValues,1),1);
    for ind=1:size(kValues,1)
        k = kValues(ind);
        integrand = 4*pi*rValuesM.^2.*sinc(rValuesM*k/pi).*rPsiM;
        fPsi(ind) = trapz(rValuesM,integrand);
    end
end

function Iw = nucOverlap(kValues, fDensity, Zw)
    % Create qtValues and qaValues = kValues
    qtValues = kValues;
    qaValues = kValues;
    % Create meshgrid values
    [qtMesh,qaMesh] = meshgrid(qtValues, qaValues);
    % Replicate Zw to correspond with meshgrid
    zMesh = repmat(Zw', size(kValues,1),1);
    % Interpolate fDensity for points on meshgrid
    fDint = interp1(qtValues, fDensity, ...
                        (1/2)*sqrt(qtMesh.^2+qaMesh.^2), 'spline');
    % Find integrand at each point
    integrand = qtMesh.^2.*qaMesh.^2.*zMesh.*fDint;
    % Integrate over all points
    Iw = trapz(kValues,trapz(kValues,integrand,2),1);
end

function Iw = nucOverlapMomentum(kValues, fDensity, Zw, Mpi)
    % Create qtValues and qaValues = kValues
    qtValues = kValues;
    qaValues = kValues;
    % Create meshgrid values
    [qtMesh,qaMesh] = meshgrid(qtValues, qaValues);
    % Replicate Zw to correspond with meshgrid
    zMesh = repmat(Zw', size(kValues,1),1);
    % Interpolate fDensity for points on meshgrid
    fDint = interp1(qtValues, fDensity, ...
                        (1/2)*sqrt(qtMesh.^2+qaMesh.^2), 'spline');
    % Find momentum dependent part on meshgrid points
    xt = (qtMesh/Mpi).^2;
    MD = ((2+xt)./sqrt(xt)).*acot(2./sqrt(xt));
    % Find integrand at each point
    integrand = qtMesh.^2.*qaMesh.^2.*zMesh.*fDint.*MD;
    % Integrate over all points
    Iw = trapz(kValues,trapz(kValues,integrand,2),1);
end
