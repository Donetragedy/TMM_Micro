function R = refl_mod_wl_v11(Ms,NK,D,H,ang,pol)
%REFL_MOD_WL calculates the reflectance of a thin film stack as a function
%of the wavelength. Rough layer interfaces and refractive index
%inhomogeneity are taken into account in the model.
%
%Syntax:
%   R = refl_mod_wl(Ms,NK,D,Dr,H,ang,pol)
%
%Input:
%   'Ms'    is a substrate refractive index data matrix with columns as
%           follows: [wl,ns,ks]. Here wl is wavelength [nm], ns is
%           index of refraction and ks is extinction coefficient.
%
%   'NK'    is a layer refractive index data matrix with columns as
%           follows: [n1,k1,...]. Here n is index of refraction and k is
%           extinction coefficient. Number Suffixes stand for different
%           layers counting upwards from the substrate. The values should
%           correspond to the wavelengths in wl.
%
%   'D'     is a layer thickness vector as follows: [d1,...]. Number
%           Suffixes stand for different layers counting upwards from the
%           substrate.
%
%   'Dr'    is a interface roughness layer thickness vector as follows:
%           [dr1,...]. Number Suffixes stand for different interfaces
%           counting upwards from the substrate. There is always one
%           interface layer more than there are actual layers.
%
%   'H'     is a layer inhomogeneity vector as follows: [h1,...]. Number
%           Suffixes stand for different layers counting upwards from the
%           substrate. The inhomogeneity value represents the relative
%           change in refractive index throughout the layer thickness.
%           Negative h value means that the refractive index is getting
%           higher towards the substrate.
%
%   'angle' is the angle of incidence [0,90).
%
%   'pol'   is the polarization state (either 's' or 'p').
%
%Output:
%   'R'     is a vector of reflectances corresponding to the wavelengths
%           given in ang.
%
% See also: reflectance_ih, ilN, srN

% Version:          0.1.2
% Last Modified:    14-May-2014
% Author:           Henrik Mantynen


% answer = questdlg('Create interface and surface layers into the structure?', ...
% 	'', ...
% 	'Yes','No','-');
% 
% switch answer
% 
%     case 'Yes'
%                
%     % Number of layers
%     nol = length(D);
% 
%     % Refractive index matrix
%     NKsi = ilN(Ms(:,2:3),[NK(:,1),NK(:,2)]); %refr ind of interface between Si substrate and first layer above (1st columns n and 2nd k)
%     NKin1=ilN([NK(:,1),NK(:,2)],[NK(:,3),NK(:,4)]);
%     NKin2=ilN([NK(:,3),NK(:,4)],[NK(:,5),NK(:,6)]);
%     NKsr = srN([NK(:,end-1),NK(:,end)]); % refr ind of a surface roughness layer (calculated based on n and k of the top layer)
%     
%     M = [Ms,NKsi,[NK(:,1),NK(:,2)],NKin1,[NK(:,3),NK(:,4)],NKin2,[NK(:,5),NK(:,6)],NKsr];
%     
%     DD = zeros(1,3);
%     DD(:,1) = 2;
%     DD(:,2) = D(:,1);
%     DD(:,3) = 2;
%     DD(:,4) = D(:,2);
%     DD(:,5) = 2;
%     DD(:,6) = D(:,3);
%     DD(:,7) = 2;
%     
%     % Layer inhomogeneity vector
%     HH = zeros(1,2*nol+1);
%     HH(2:2:end) = H(:); 



%     case 'No'

    % Number of layers
    nol = length(D);

    % Refractive index matrix
    M = [Ms,NK];
    
    DD = zeros(1,numel(D));
    DD(:) = D(:);

    
    % Layer inhomogeneity vector
    HH = zeros(1,2*nol+1);
    HH(2:2:end) = H(:); 
    



% Calculate reflectance
R = reflectance_ih_v11(M,DD,HH,ang,pol);

end

