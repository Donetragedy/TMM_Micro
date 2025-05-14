function Nmix = ilN(N1,N2,volfr)
%ILN calculates the effective complex refractive index of a interface
%mixture layer with Bruggeman effective medium approximation.
%
%Syntax:
%   Nmix = ilN(N1,N2,volfr)
%
%Input:
%   'N1'    is a matrix [n1,k1] representing the complex index of
%           refraction of the first layer.
%
%   'N2'    is a matrix [n2,k2] representing the complex index of
%           refraction of the second layer.
%
%   'volfr' shows how much volume fraction (0 - 1) of the inclusion material
%
%Output:
%   'Nmix'  is a matrix [nmix,kmix] representing the complex index of
%           refraction of the mixture layer.
%
%See also: bemari, bemaf, srN, incN

% Version:          1.0.4
% Last Modified:    08-Jul-2013
% Author:           Henrik Mäntynen


% Checking input
assert(length(N1(1,:))==2, ...
    'N1 has to be a two column matrix.')
assert(length(N1(1,:))==2, ...
    'N2 has to be a two column matrix.')
assert(all(size(N1)==size(N2)), ...
    'N1 and N2 must have the same size.')

% Calculating
Nhost = N1(:,1)-1i*N1(:,2);
Ninc = N2(:,1)-1i*N2(:,2);
Neff = bemari(Nhost,Ninc,volfr);
Nmix = [real(Neff),-imag(Neff)];

end