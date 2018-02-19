function overlapIntegralsProton()
%overlapIntegralsProton Calculates overlap integrals for proton
% distribution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    fID = fopen('output.txt', 'w');
    
    % Electron Mass (MeV)
    Me = .5109989;
    % Muon Mass (MeV)
    Mmu = 105.65837;
    % Pion Mass (MeV)
    Mpi = 138.039;
    % Atomic Number
    Z = 13;
    % Atomic Weight
    A = 27;
    
    % Proton density distribution parameters
    pVecProton = [7.0, .43418*10^-1, .60298*10^-1, .28950*10^-2, ...
                    -.23522*10^-1, -.79791*10^-2, .23010*10^-2, ...
                    .10794*10^-2, .12574*10^-3, -.13021*10^-3, ...
                    .56563*10^-4, -.18011*10^-4, .42869*10^-5];
    
    % Extract reduced radial wavefunctions for electron and muon
    [rValues, u1mu, u2mu, bE] = boundRWF(Mmu, Z, 7, pVecProton);
    [~, u1em, u2em] = freeRWF(Me, Z, -1, Mmu - bE, 7, pVecProton);
    [~, u1ep, u2ep] = freeRWF(Me, Z, 1, Mmu - bE, 7, pVecProton);
    clear bE
    
    % Switch from reduced to normal radial wavefunctions
    gmu = u1mu./rValues;
    fmu = u2mu./rValues;
    gem = u1em./rValues;
    fem = u2em./rValues;
    gep = u1ep./rValues;
    fep = u2ep./rValues;
    
    % Define range of momentum values of interest (MeV)
    kValues = (.5:.5:400)';
    
    % Create the eight leptonic overlap functions
    Zs1 = lepOver(1,1,rValues,kValues,gem,fem,gmu,fmu);
    Zs2 = lepOver(1,2,rValues,kValues,gem,fem,gmu,fmu);
    Zs3 = lepOver(1,3,rValues,kValues,gep,fep,gmu,fmu);
    Zs4 = lepOver(1,4,rValues,kValues,gep,fep,gmu,fmu);
    Zv1 = lepOver(2,1,rValues,kValues,gem,fem,gmu,fmu);
    Zv2 = lepOver(2,2,rValues,kValues,gem,fem,gmu,fmu);
    Zv3 = lepOver(2,3,rValues,kValues,gep,fep,gmu,fmu);
    Zv4 = lepOver(2,4,rValues,kValues,gep,fep,gmu,fmu);
    
    t1 = toc;
    fprintf(fID, 'Check Point #1: %d \n', t1);
    
    % Clear unneeded variables
    clear gmu fmu gem fem gep fep t1 u1em u2em ...
                                u1ep u2ep u1mu u2mu bE
    
%%%%% Consider proton distributions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    % Find proton wavefunction
    PsiProton = wavefunction(Z, 7, pVecProton, rValues);
    
    % Fourier transform wavefunction
    PsiFProton = fourierTrans(PsiProton, rValues, kValues);
    
    % Square to find fourier transformed density distribution
    fDensityP = PsiFProton.^2;
    
    % Find overlap integrals
    Isp1 = nucOverlap(kValues, fDensityP, Zs1);
    Isp2 = nucOverlap(kValues, fDensityP, Zs2);
    Isp3 = nucOverlap(kValues, fDensityP, Zs3);
    Isp4 = nucOverlap(kValues, fDensityP, Zs4);
    tIsp1 = nucOverlapMomentum(kValues, fDensityP, Zs1, Mpi);
    tIsp2 = nucOverlapMomentum(kValues, fDensityP, Zs2, Mpi);
    tIsp3 = nucOverlapMomentum(kValues, fDensityP, Zs3, Mpi);
    tIsp4 = nucOverlapMomentum(kValues, fDensityP, Zs4, Mpi);
    Ivp1 = nucOverlap(kValues, fDensityP, Zv1);
    Ivp2 = nucOverlap(kValues, fDensityP, Zv2);
    Ivp3 = nucOverlap(kValues, fDensityP, Zv3);
    Ivp4 = nucOverlap(kValues, fDensityP, Zv4);
    
    % Print values of overlap integrals
    fprintf(fID, 'Isp1 = %d \n', Isp1);
    fprintf(fID, 'Isp2 = %d \n', Isp2);
    fprintf(fID, 'Isp3 = %d \n', Isp3);
    fprintf(fID, 'Isp4 = %d \n', Isp4);
    fprintf(fID, '~Isp1 = %d \n', tIsp1);
    fprintf(fID, '~Isp2 = %d \n', tIsp2);
    fprintf(fID, '~Isp3 = %d \n', tIsp3);
    fprintf(fID, '~Isp4 = %d \n', tIsp4);
    fprintf(fID, 'Ivp1 = %d \n', Ivp1);
    fprintf(fID, 'Ivp2 = %d \n', Ivp2);
    fprintf(fID, 'Ivp3 = %d \n', Ivp3);
    fprintf(fID, 'Ivp4 = %d \n', Ivp4);
    
    t2 = toc;
    fprintf(fID, 'Check Point #2: %d \n', t2);
    
    % Clear unneeded variables
    clear PsiProton pVecProton PsiFProton fDensityP ...
            Isp1 Isp2 Isp3 Isp4 tIsp1 tIsp2 tIsp3 tIsp4 ...
            Ivp1 Ivp2 Ivp3 Ivp4 t2
 
    fclose(fID);
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

function zFunc = lepOver(n1,n2,rValues,kValues,ge,fe,gmu,fmu)
    % n1 is used to distinguish the scalar case (n1=1) from
    % the vector case (n1=2).
    % n2 takes on values {1,2,3,4} for the four possible
    % overlap integrals.
    
    % Convert all functions from units of fm to units of MeV
    rValuesM = rValues/(197.327);
    geM = ge*(197.327)^(3/2);
    feM = fe*(197.327)^(3/2);
    gmuM = gmu*(197.327)^(3/2);
    fmuM = fmu*(197.327)^(3/2);
    % Initialize zFunc
    zFunc = zeros(size(kValues,1),1);
    if(n1==1)
        if(n2==1)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(2*pi^2))*rValuesM.^2 ...
                            .*sin(rValuesM*k)./(rValuesM*k)...
                            .*(geM.*gmuM+feM.*fmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        elseif(n2==2)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(8*pi))*rValuesM.^2 ...
                            .*(sin(rValuesM*k)./(rValuesM*k).^2 ...
                            - cos(rValuesM*k)./(rValuesM*k))...
                            .*(geM.*fmuM-feM.*gmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        elseif(n2==3)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(8*pi))*rValuesM.^2 ...
                            .*(sin(rValuesM*k)./(rValuesM*k).^2 ...
                            - cos(rValuesM*k)./(rValuesM*k))...
                            .*(geM.*gmuM+feM.*fmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        elseif(n2==4)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(2*pi^2))*rValuesM.^2 ...
                            .*sin(rValuesM*k)./(rValuesM*k)...
                            .*(feM.*gmuM-geM.*fmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        else
            error('Invalid value for n2.')
        end
    elseif(n1==2)
        if(n2==1)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(2*pi^2))*rValuesM.^2 ...
                            .*sin(rValuesM*k)./(rValuesM*k)...
                            .*(geM.*gmuM-feM.*fmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        elseif(n2==2)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(8*pi))*rValuesM.^2 ...
                            .*(sin(rValuesM*k)./(rValuesM*k).^2 ...
                            - cos(rValuesM*k)./(rValuesM*k))...
                            .*(geM.*fmuM+feM.*gmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        elseif(n2==3)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(8*pi))*rValuesM.^2 ...
                            .*(sin(rValuesM*k)./(rValuesM*k).^2 ...
                            - cos(rValuesM*k)./(rValuesM*k))...
                            .*(geM.*gmuM-feM.*fmuM);
                zFunc(ind) = trapz(rValuesM,integrand);
            end
        elseif(n2==4)
            for ind=1:size(kValues,1)
                k = kValues(ind);
                integrand = (1/(2*pi^2))*rValuesM.^2 ...
                            .*sin(rValuesM*k)./(rValuesM*k)...
                            .*(feM.*gmuM+geM.*fmuM);
                zFunc(ind) = -trapz(rValuesM,integrand);
            end
        else
            error('Invalid value for n2.')
        end        
    else
       error('Invalid value for n1.') 
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
