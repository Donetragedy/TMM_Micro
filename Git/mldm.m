function [ldmat,varargout] = mldm(wl,varargin)
%MLDM constructs a single layer data matrix from given material data
%matrices interpolating the values to the given wavelength range.
%
%Syntax:
%   ldmat = mdm(wl,mdata1,mdata2,...);
%   or
%   [ldmat,mdata1i,mdata2i,...] = mdm(wl,mdata1,mdata2,...);
%
%Input:
%   'wl'        is a column vector of the desired wavelength points [nm].
%
%   'mdataX'    are material data matrices with columns as follows:
%               [wavelength, n, k]. The wavelength ranges in these matrices
%               must be inside the range of the input parameter wl.
%
%Output:
%   'ldmat'     is the layer data matrix with columns as follows:
%               [wavelength,n1,k1,n2,k2,...]. Here all the given material
%               data matrices have been interpolated to the range of input
%               wl and concatenated to a single matrix.
%
%   'mdataXi'   are optional output matrices containing the interpolated
%               values of the corresponding input matrices 'mdataX'.
%
%See also: reflectance, transmittance

% Version:          1.0.3
% Last Modified:    08-Jul-2013
% Author:           Henrik Mäntynen


% Asserts
assert(iscolumn(wl), ...
    'The input wl must be a column vector.')
assert(isreal(wl) && all(wl > 0), ...
    'Wavelengths must be real and nonnegative.')
assert(size(varargin,2) > 0, ...
    'Too few input parameters.')

% Number of layers
nol = size(varargin,2); %?????????????????????????????????????????????????

% Initialize lmod and varargout
ldmat = [wl,zeros(length(wl),nol*2)];
varargout = cell(1,nol);

% Interpolate mdata and add it to lmod and varargout

for idx = 1:nol;
    mdata = varargin{idx};
    n = interp1(mdata(:,1),mdata(:,2),wl);
    k = interp1(mdata(:,1),mdata(:,3),wl);
    ldmat(:,(2 * idx)) = n;
    ldmat(:,(2 * idx + 1)) = k;

end

end

% This file is for interpolation the dispersion of material data to the
% desired spectrum of our measurement.


