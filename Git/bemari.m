function Neff = bemari(N1,N2,f)
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


% Checking input
%assert((f <= 1) && (f>=0),'ERROR: Invalid input! f needs to be between 0 and 1.')
assert(all([isvector(N1),isvector(N2)]), ...
    'N1 and N2 must be vectors.')
assert(isequal(size(N1),size(N2)), ...
    'Vectors N1 and N2 must have the same size.')

% Changing complex indexes to dielectric constants
eps1 = N1.^2;
eps2 = N2.^2;

% Calculating the effective dielectric constant
p = sqrt(eps1./eps2);
b = ((3*f-1)*(1./p-p)+p)/4;
z = b + sqrt(b.^2+0.5);
epseff = z.*sqrt(eps1.*eps2);

% Changing back to complex index of refraction
Neff = sqrt(epseff);

end

