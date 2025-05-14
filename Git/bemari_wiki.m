function Neff = bemari(eps_base,eps_incl,vol_incl)
%BEMARI calculates the Bruggeman effective medium approximation for complex
%refractive index of a material with inclusions.
%
%Syntax:
%   Neff = bemari(N1,N2,f);
%
%Input:
%   'N1'    is the complex refractive index of the host material.
%
%   'N2'    is the complex refractive index of the inclusion material.
%
%   'f'     is the volume fraction (0 - 1) of the inclusion material.
%
%Output:
%   'Neff'  is the effective complex index of refraction.
%
%See also: bemaf, bemari3

% Version:          1.1.3
% Last Modified:    08-Jul-2013
% Author:           Henrik Mäntynen




% Changing complex indexes to dielectric constants

    factor_up = 2.*(1 - vol_incl).*eps_base + (1 + 2.* vol_incl).*eps_incl;
    factor_down = (2 + vol_incl).* eps_base + (1 - vol_incl).* eps_incl;
    Neff = eps_base.*factor_up./factor_down

end

