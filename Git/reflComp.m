function reflComp(Rmeas,Rmod,varargin)
%REFLCOMP Plots a comparison of measured and modeled reflectances given as
%input.
%
%Syntax:
%   reflComp(Rmeas,Rmod);
%   or
%   reflComp(Rmeas,Rmod,'PropertyName',propertyvalue,...);
%
%Input:
%   'Rmeas' is a matrix of the measured reflectance as follows: [wl,R].
%           Here wl is a column vector of wavelengths [nm] and R is a
%           column vector of the corresponding reflectances.
%
%   'Rmod'  is a matrix of the modelled reflectance as follows: [wl,R].
%           Here wl is a column vector of wavelengths [nm] and R is a
%           column vector of the corresponding reflectances.
%
%Optional input:
%   'type'  is an optional input parameter that defines the type of the
%           input reflectance values: '%' for percentile [0-100] and 'abs'
%           for absolute [0-1] (The default is absolute).
%
%   'fnum'  is an optional input parameter that defines the number of the
%           figure where the plot is drawn.
%
%   'xlab'  is an optional input parameter that defines the x-axis label in
%           the plots. The nominal label is 'Wavelength [nm]'.
%
%See also: Rplot, nkplot

% Version:          1.2.5
% Last Modified:    09-May-2014
% Author:           Henrik Mäntynen


% Check input
assert(size(Rmeas,2) == 2, 'Invalid input Rmeas.')
assert(size(Rmod,2) == 2, 'Invalid input Rmod')

% Optional input
try
    mult = 100;
    newfig = true;
%     xlab = 'Wavelength [nm]';
    for n = 1:2:length(varargin);
        switch lower(varargin{n})
            case 'type'
                type = varargin{n+1};
                if strcmp(type,'%')
                    mult = 1;
                end
            case 'fnum'
                fn = varargin{n+1};
                newfig = false;
            case 'xlab'
                xlab = varargin{n+1};
        end
    end
    
catch ERR
    fprintf(2, ...
        'Failed to read optional input.\n\nCaused by:\n');
    rethrow(ERR)
end

% Calculate residuals and mean relative difference
try
    wlm = Rmeas(:,1);
    Rm = Rmeas(:,2)*mult;
    wl = Rmod(:,1);
    R = Rmod(:,2)*mult;
    resid = interp1(wl,R,wlm) - Rm;
    mrd = mean(abs(resid)./Rm)*100; %calculated as mean relative difference between sim and meas
    md = mean(abs(resid)); %calculated as mean difference between sim and meas
    
catch ERR
    fprintf(2, ...
        'Failed to calculate residuals and mean relative difference.\n\nCaused by:\n');
    rethrow(ERR)
end

% Create the plot
try
    if newfig
        h = figure();
    else
        h = figure(fn);
    end
    set(h,'Name','Comparison');
    
    subplot(2,1,1)
    plot(wl,R,'LineWidth',1), hold on
    plot(wlm,Rm,'.','MarkerSize',7)
    xlim([min(wl),max(wl)]);
    ylim([0,80]);
    grid on
    ax = gca;
    ax.FontSize = 15;
    legend('Modeled','Measured','FontSize', 15)
    xlabel('Wavelength [nm]','FontSize', 15)
    ylabel('R [%]','FontSize', 15)

%  Excel file saving
%     sim_data = table(wl,R);
%     writetable(sim_data,'4layer.xls');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      With_Si = [wl,R];
%      save('Cutted_input_data.mat', 'With_Si');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,1,2)
    plot(wlm,resid,'-x','MarkerSize',5)
    title(['\fontsize{15} Mean difference ',num2str(round(md,2)),'% and relative mean difference ',num2str(round(mrd,2)),'%'])
    xlim([min(wl),max(wl)]);
    grid on
    ax = gca;
    ax.FontSize = 15;
    xlabel('Wavelength [nm]','FontSize', 15)
    ylabel('Residual [% unit]','FontSize', 15)
    
catch ERR
    fprintf(2, ...
        'Failed to create the plot.\n\nCaused by:\n');
    rethrow(ERR)
end

end



