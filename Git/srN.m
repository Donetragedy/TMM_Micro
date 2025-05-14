function Nrough = srN(N)
%SRN calculates the effective complex refractive index of a surface
%roughness layer with Bruggeman effective medium approximation.
%
%Syntax:
%   Nrough = srN(N)
%
%Input:
%   'N'         is a matrix [n,k] representing the complex index of
%               refraction of the actual layer.
%
%Output:
%   'Nrough'    is a matrix [nr,kr] representing the complex index of
%               refraction of the surface roughness layer.
%
%See also: bemari, bemaf, ilN, incN

% Version:          1.0.6
% Last Modified:    08-Jul-2013
% Author:           Henrik Mäntynen


% Check input
assert(length(N(1,:))==2, 'N has to be a two column matrix.')

% Calculating
nair = zeros(length(N(:,1)),1) + 1.00027;
Nhost = N(:,1)-1i*N(:,2);
Neff = bemari(Nhost,nair,0.5);
Nrough = [real(Neff),-imag(Neff)];

end