function xBar = reflectance_ih_v11(mdata,D,H,angle,pol)
%REFLECTANCE_IH calculates the spectral specular reflectance of a thin-film
%system with linear thickness inhomogeneity in layers. A semi-infinite substrate
%is assumed.
%
%Syntax:
%   R = reflectance_ih(mdata,D,H,angle,pol);
%
%Input:
%   'mdata' is a layer material data matrix with columns as follows:
%           [wl,ns,ks,n1,k1,...]. Here wl is wavelength [nm], n is index of
%           refraction and k is extinction coefficient. Suffix 's' stands
%           for substrate and numbers represent different layers counting
%           upwards from the substrate.
%
%   'D'     is a layer thickness vector as follows: [d1,...].
%
%   'H'     is a layer inhomogeneity vector as follows: [h1,...]. The
%           inhomogeneity value represents the relative change in
%           refractive index throughout the layer thickness. Negative h
%           value means that the refractive index is getting higher towards
%           the substrate.
%
%   'angle' is the angle of incidence [0,90).
%
%   'pol'   is the polarization state (either 's', 'p' or 'unpolarized').
%
%Output:
%   'R'     is a vector of reflectances corresponding to the wavelengths
%           given in wl.
%
% See also: reflectance, reflectance_d, reflectance_ft, reflectance_nu_soa,
% transmittance, transmittance_c, mldm



%/--- Checking input parameters ---/
assert(isreal(angle), 'Invalid angle.')
assert(all(angle >= 0) && all(angle < 90), ...
    'The angle must be non-negative and less than 90.');
assert(isreal(D), 'Thickness values must be real.')
%assert(all(D >= 0),'Layer thickness must be non-negative.')
assert(isreal(H), 'Inhomogeneity values must be real.')
assert(ismatrix(mdata) && isreal(mdata), ...
    'Invalid material data matrix.')
assert(((length(mdata(1,:))-3)/2) == numel(D), ...
    'The material data matrix has wrong size.');
% mdata(1,:))-3 mean we substract wl column and 2 substrate columns with n
% and k
assert(ischar(pol), 'Invalid polarization.')
pol = lower(pol);
assert(strcmp(pol,'s')||strcmp(pol,'p')||strcmp(pol,'unpolarized'), ...
    'pol must be either ''s'' or ''p''or ''unpolarized''.')


%/--- Data ---/
try
    % Substrate data
    wl = mdata(:,1);
    nsub = mdata(:,2);
    ksub = mdata(:,3);
    
    % Layer data
    layers = (length(mdata(1,:)) - 3)/2;
    Ln = mdata(:, 4 + ((1:layers) - 1) * 2);
    Lk = mdata(:, 5 + ((1:layers) - 1) * 2);
    
    % Medium (air)
    nm = 1.00027;
    km = 0;
    N0 = nm - 1i*km; 
    
    % Wave admittance (air)
    %permittivity = 8.854187817620e-12;
    %c0 = 299792458;
    %Y = permittivity*c0;
    
catch ERR
    fprintf(2, 'Error in reading input data.\n\nCaused by:\n');
    rethrow(ERR)
end

       


%/--- Calculate reflectance with matrix method ---/

    function Rs_pol = s_polarization(wl,nsub,ksub,Ln,Lk,N0)

            % Number of wavelengths
        nowl = length(wl);
        
        % Some initial values
        theta0 = ones(nowl,1) * (angle/180) * pi; %deg to rad
        sin_theta0 = sin(theta0);
        cos_theta0 = cos(theta0);
        eta0s = N0 .* cos_theta0; %(N0 * Y) .* cos(theta0)
        lambda = wl; %* 1e-9;
        sin_theta = sin_theta0;
        %As = [1 0;0 1];
        As11 = ones(nowl,1);
        As12 = zeros(nowl,1);
        As21 = As12;
        As22 = As11;
        Nprev = ones(nowl,1) * N0;
        N = Nprev; %seems like refractive indices if air
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Form and multiply layer matrices
        for indx = numel(D):-1:1; %starting from D matrix size to 1
            d = D(indx); %* 1e-9;
            h = H(indx);
            if d > 1e-2                
                J = floor(d/2) * (h ~= 0) + 1;
                d = d / J;
                for jj = 1:J;
                    hc = 1 - h*((jj-1)/max(J-1,1) - 0.5);
                    sin_thetaprev = sin_theta;
                    N = hc*((Ln(:,indx)) - 1i * Lk(:,indx));
                    sin_theta = (Nprev ./ N) .* sin_thetaprev;
                    cos_theta = sqrt(1 - sin_theta.^2);
                    etas = N .* cos_theta; %(N * Y) .* cos(theta);
                    sig = 2*pi * (d ./ lambda) .* N .* cos_theta;
                    
                    %Ms = [cos(sig), (1i*sin(sig))/etas; 1i*etas*sin(sig), cos(sig)];
                    Ms11 = cos(sig);
                    Ms12 = (1i*sin(sig)) ./ etas;
                    Ms21 = 1i*etas .* sin(sig);
                    Ms22 = Ms11;
                    
                    %As = As * Ms;
                    temp11 = As11.*Ms11 + As12.*Ms21;
                    temp12 = As11.*Ms12 + As12.*Ms22;
                    temp21 = As21.*Ms11 + As22.*Ms21;
                    temp22 = As21.*Ms12 + As22.*Ms22;
                    As11 = temp11;
                    As12 = temp12;
                    As21 = temp21;
                    As22 = temp22;
                    
                    Nprev = N;
                end
            end
        end
    
        
        % Take substrate into account
        Nsub = nsub - 1i*ksub;
        sin_thetam = (N ./ Nsub) .* sin_theta;
        cos_thetam = sqrt(1 - sin_thetam.^2);
        etams = Nsub .* cos_thetam; %(Nsub * Y) .* cos(thetam);
        % As = As * [1;etams] = [Bs;Cs]
        Bs = As11 + As12.*etams;
        Cs = As21 + As22.*etams;
        
        % Solve reflectance
        rhos = (eta0s.*Bs - Cs) ./ (eta0s.*Bs + Cs);
        Rs_pol = rhos .* conj(rhos);
    end

    
    function Rp_pol = p_polarization(wl,nsub,ksub,Ln,Lk,N0)

       % Number of wavelengths
        nowl = length(wl);

        % Some initial values
        theta0 = ones(nowl,1) * (angle/180) * pi;
        sin_theta0 = sin(theta0);
        cos_theta0 = cos(theta0);
        eta0p = N0 ./ cos_theta0; %(N0 * Y) ./ cos(theta0);
        lambda = wl; %* 1e-9;
        sin_theta = sin_theta0;
        %Ap = [1 0;0 1];
        Ap11 = ones(nowl,1);
        Ap12 = zeros(nowl,1);
        Ap21 = Ap12;
        Ap22 = Ap11;
        Nprev = ones(nowl,1) * N0;
        N = Nprev;
        
        % Form and multiply layer matrices
        for indx = numel(D):-1:1;
            d = D(indx); %* 1e-9;
            h = H(indx);
            if d > 1e-2
                J = floor(d/2) * (h ~= 0) + 1;
                d = d / J;
                for jj = 1:J;
                    hc = 1 - h*((jj-1)/max(J-1,1) - 0.5);
                    sin_thetaprev = sin_theta;
                    N = hc*(Ln(:,indx)) - 1i * Lk(:,indx);
                    sin_theta = (Nprev ./ N) .* sin_thetaprev;
                    cos_theta = sqrt(1 - sin_theta.^2);
                    etap = N ./ cos_theta; %(N * Y) ./ cos(theta);
                    sig = 2*pi * (d ./ lambda) .* N .* cos_theta;
                    
                    %Mp = [cos(sig), (1i*sin(sig))/etap; 1i*etap*sin(sig), cos(sig)];
                    Mp11 = cos(sig);
                    Mp12 = (1i*sin(sig)) ./ etap;
                    Mp21 = 1i*etap .* sin(sig);
                    Mp22 = Mp11;
                    
                    %Ap = Ap * Mp;
                    temp11 = Ap11.*Mp11 + Ap12.*Mp21;
                    temp12 = Ap11.*Mp12 + Ap12.*Mp22;
                    temp21 = Ap21.*Mp11 + Ap22.*Mp21;
                    temp22 = Ap21.*Mp12 + Ap22.*Mp22;
                    Ap11 = temp11;
                    Ap12 = temp12;
                    Ap21 = temp21;
                    Ap22 = temp22;
                    
                    Nprev = N;
                end
            end
        end
        
        % Take substrate into account
        Nsub = nsub - 1i*ksub;
        sin_thetam = (N./Nsub) .* sin_theta;
        cos_thetam = sqrt(1 - sin_thetam.^2);
        etamp = Nsub ./ cos_thetam; %(Nsub * Y) ./ cos(thetam);
        % Ap = Ap * [1;etamp]
        Bp = Ap11 + Ap12.*etamp;
        Cp = Ap21 + Ap22.*etamp;
        
        % Solve reflectance
        rhop = (eta0p.*Bp - Cp) ./ (eta0p.*Bp + Cp);
        Rp_pol = rhop .* conj(rhop);
    end

    

    % s-polarization
if pol == 's'
    Rs_pol = s_polarization(wl,nsub,ksub,Ln,Lk,N0);
   
% Result
xBar = Rs_pol;
end


    % p-polarization
if pol == 'p'
    Rp_pol = p_polarization(wl,nsub,ksub,Ln,Lk,N0);
   
% Result
xBar = Rp_pol;

 end


    % unpolarized-polarization
if pol == 'unpolarized'
    Rfun = (s_polarization(wl,nsub,ksub,Ln,Lk,N0)+p_polarization(wl,nsub,ksub,Ln,Lk,N0))/2;

% Result
xBar = Rfun;

end

end
