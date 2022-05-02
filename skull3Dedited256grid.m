
clearvars;

% =========================================================================
% DEFINE LITERALS
% =========================================================================
    

% medium parameters
c0              = 1500;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]
u               = 0.04;
cbrain          = 1560;     % speed of sound [m/s]
rhobrain        = 1040; 
cskullcortical  = 2800;     % speed of sound [m/s]
rhoskullcortical= 1850;
cskulltrabecular= 2300; % speed of sound [m/s]
rhoskulltrabecular= 1700; 
alphabrain = 1.2;
alphatrab = 32;
alphacort = 16;
alpha0 = 0;

% source parameters
source_f0       = 500000;   % source frequency [Hz]
source_roc      = 64e-3;    % bowl radius of curvature [m]
source_diameter = 64e-3;    % bowl aperture diameter [m]
source_mag      = rho0*c0*u;% source pressure [Pa]
offset          = source_roc - sqrt((source_roc^2)-(source_diameter/2)^2);

% grid parameters
axial_size      = 64e-3;    % total grid size in the axial dimension [m]
lateral_size    = 64e-3;   % total grid size in the lateral dimension [m]

% computational parameters
ppw             = 24;        % number of points per wavelength
t_end           = 120e-6;    % total compute time [s] (this must be long enough to reach steady state)
record_periods  = 3;        % number of periods to record
cfl             = 0.03;     % CFL number

% =========================================================================
% RUN SIMULATION
% =========================================================================

% --------------------
% GRID
% --------------------

% calculate the grid spacing based on the PPW and F0
%dx = c0 / (ppw * source_f0);   % [m]

%source_x_offset = round(offset/dx); %grid points for offset of source 

% compute the size of the grid
% Nx = roundEven(axial_size / dx) + source_x_offset;
% Ny = roundEven(lateral_size / dx);
% Nz = Ny
Nx = 256;
Ny = 256;
Nz = 256;

dx = axial_size/Nx;
dy = lateral_size/Ny;
dz = lateral_size/Nz;
% create the computational grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% compute points per temporal period
PPP = round(ppw / cfl);

% compute corresponding time spacing
dt = 1 / (PPP * source_f0);
%dt = (0.3*dx)/2800;
% create the time array using an integer number of points per period
Nt = round(t_end / dt);
kgrid.setTime(Nt, dt);

% calculate the actual CFL and PPW
disp(['PPW = ' num2str(c0 / (dx * source_f0))]);
disp(['CFL = ' num2str(c0 * dt / dx)]);


% --------------------
% SOURCE
% --------------------
% create time varying source
source_sig = createCWSignals(kgrid.t_array, source_f0, source_mag, 0);

bowl_pos  = [1, Ny/2, Nz/2];

focus_pos = [Nx, Ny/2, Nz/2];
grid_size = [Nx, Ny, Nz];

source_roc = 256;
source_diameter = 255;
% create bowl
bowl_binary = makeBowl(grid_size, bowl_pos, source_roc, source_diameter, focus_pos, 'Plot', true);

% assign source mask
source.p_mask = bowl_binary;

source.p = source_sig;
% --------------------
% MEDIUM
% --------------------
% import cranial geometry

load('skullmediumgeometry','mish')

imagesc(mish(:,:,128));
mediummm = rot90(mish);
mediummm = flip(mediummm,1);
medium.sound_speed = mediummm;
medium.sound_speed(medium.sound_speed == 0) = cbrain;
medium.sound_speed(medium.sound_speed == 3) = cskulltrabecular;
medium.sound_speed(medium.sound_speed == 2) = cskullcortical;
medium.sound_speed(medium.sound_speed == 1) = c0; % [m/s]

medium.density = mediummm;
medium.density(medium.density == 0) = rhobrain;
medium.density(medium.density == 3) = rhoskulltrabecular;
medium.density(medium.density == 2) = rhoskullcortical;
medium.density(medium.density == 1) = rho0; % [m/s]


medium.alpha_power = 2;
medium.alpha_coeff = mediummm;
medium.alpha_coeff(medium.alpha_coeff == 0) = alphabrain;
medium.alpha_coeff(medium.alpha_coeff == 3) = alphatrab;
medium.alpha_coeff(medium.alpha_coeff == 2) = alphacort;
medium.alpha_coeff(medium.alpha_coeff == 1) = alpha0; 



% --------------------
% SENSOR
% --------------------

% set sensor mask to record central plane, not including the source point
sensor.mask = zeros(Nx, Ny, Nz);
%sensor.mask(source_x_offset + 2:end, :, Nz/2 + 1) = 1;
sensor.mask(:, :,Nz/2 + 1) = 1;
% record the pressure
sensor.record = {'p'};

% average only the final few periods when the field is in steady state
sensor.record_start_index = kgrid.Nt - record_periods * PPP + 1;


% --------------------
% SIMULATION
% --------------------

% set input options
input_args = {...
    'PMLSize', 'auto', ...
    'PMLInside', false, ...
    'PlotPML', false, ...
    'DisplayMask', 'off'};

% run code

sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, ...
    input_args{:}, ...
    'DataCast', 'single', ...
    'PlotScale', [-1, 1] * source_mag);

% extract amplitude from the sensor data
amp = extractAmpPhase(sensor_data.p, 1/kgrid.dt, source_f0, ...
    'Dim', 2, 'Window', 'Rectangular', 'FFTPadding', 1);

% reshape data
%amp = reshape(amp, Nx - source_x_offset - 1, Ny);
amp = reshape(amp,Nx, Ny);
% extract pressure on axis
amp_on_axis = amp(:, Ny/2 + 1);

% define axis vectors for plotting
%x_vecc = kgrid.x_vec(source_x_offset + 2:end, :) - kgrid.x_vec(source_x_offset + 1);
x_vecc = kgrid.x_vec(2:end, :) - kgrid.x_vec(1);
y_vec = kgrid.y_vec;

% =========================================================================
% ANALYTICAL SOLUTION
% =========================================================================

% calculate the wavenumber
k = 2 * pi * source_f0 ./ c0;

% define radius and axis
x_max = (Nx * dx);
x_ref = 0:x_max/10000:x_max;

% calculate analytical solution
p_ref_axial = focusedBowlONeil(source_roc, source_diameter, source_mag / (c0 * rho0), source_f0, c0, rho0, x_ref, []);

% calculate analytical solution at exactly the same points as k-Wave
p_ref_axial_kw = focusedBowlONeil(source_roc, source_diameter, source_mag / (c0 * rho0), source_f0, c0, rho0, x_vecc, []);

% calculate error
%L2_error = 100 * sqrt( sum( (p_ref_axial_kw(:) - amp_on_axis(:)).^2 ) / sum( p_ref_axial_kw(:).^2 ) );
%Linf_error = 100 * max(abs(p_ref_axial_kw(:) - amp_on_axis(:))) / max(p_ref_axial_kw(:));
% =========================================================================
% Extracting Pressure and Time Data
% =========================================================================
x_vec = kgrid.x_vec(2:end, :) - kgrid.x_vec(1);
x_vec1 = kgrid.x_vec(1:end, :) - kgrid.x_vec(1);
timeperiod= 0:dt:t_end-dt;
periodofwave = 1/source_f0;
cyclelength = periodofwave/dt;
timeindex = Nt-(cyclelength-1);
wavetime = timeperiod(timeindex:Nt);
samplingfrequency = 1/dt;
freqv = (0:cyclelength/2-1)*samplingfrequency/cyclelength;
FFTarray = [];
FFTmag = [];
for i=(Ny/2)*Nx:((Ny/2)*Nx)+Nx-1
    pressuretime = sensor_data.p(i, 1:end);
    cycleextract = pressuretime(Nt-timeindex:end);
    %cycleextract = pressuretime;Nt-sensor.record_start_index
    FFT = fft(cycleextract);

    %symmetric so remove second half
    FFT = FFT(1:(round(cyclelength,0))/2);
    %magnitude of fft
    mFFT = (abs(FFT));
    scaledFFT = mFFT*(2/PPP);
    [FFTpeak, FFTindex] = max(scaledFFT);
    peakfreq = freqv(FFTindex);
    FFTarray = [FFTarray,peakfreq];
    FFTmag = [FFTmag,FFTpeak];
end

% =========================================================================
% VISUALISATION
% =========================================================================

% plot fft 
figure;
plot(x_vec1*10^3,FFTmag/10^3)
ylabel('Magnitude of Acoustic Pressure at Peak Frequency P(f) [kPa]')
xlabel('Axial Position (mm)')
title('Acoustic Pressure at Peak Frequency on Axial Positions after reaching steady-state with a bowl piston and skull geometry')


% =========================================================================
% VISUALISATION
% =========================================================================
% 
% % plot the pressure along the focal axis of the piston
% figure;
% plot(1e3 * x_ref, 1e-6 * p_ref_axial, 'k-');
% hold on;
% plot(1e3 * x_vec1, 1e-6 * amp_on_axis, 'b.');
% hold off;
% set(gca, 'XLim', [0, axial_size] * 1e3);
% xlabel('Axial Position [mm]');
% ylabel('Pressure [MPa]');
% legend('Exact', 'k-Wave', 'Location', 'Best');
% title('Axial Pressure');



% plot the pressure field 
figure;
imagesc(1e3 * x_vecc, 1e3 * y_vec, 1e-3*amp.');
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
axis image;
title('Acoustic Pressure Field with a heterogeneous medium of complex cranial geometry and focused bowl source generated using k-Wave');
a = colorbar;
a.Label.String = 'Acoustic Pressure (kPa)';
%% plane of skull
%sensor.mask(:, :,Nz/2 + 1) = 1;
xplane = squeeze(medium.sound_speed(:, :,Nz/2 + 1));

x_vector = [0;x_vec];
xplane = squeeze(medium.sound_speed(:, :,Nz/2 + 1));
xplane = flip(xplane);
x_vector = flip(x_vecc);
out = permute(xplane,[2,1]);

imagesc(x_vector*1e3,y_vec*1e3,out)
axis equal
title('Acoustic Pressure field medium')
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
colorbar


grid_weights = karray.getArrayGridWeights(kgrid);
xplanes = squeeze((grid_weights(:, :,Nz/2 + 1)));
xplanes = flip(xplanes);
outs = permute(xplanes,[2,1]);
imagesc(x_vector*1e3,y_vec*1e3,outs)
axis equal
title('Acoustic Pressure field source')
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
colorbar

visualisation = out + (outs.*1000);
imagesc(x_vector,y_vec,visualisation)
axis equal
title('Acoustic Pressure field source')
xlabel('Axial Position [mm]');
ylabel('Lateral Position [mm]');
colorbar



