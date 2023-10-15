f0=2.5e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/f0; % Wavelength
element_height=5/1000; % Height of element [m]
kerf=0.1/1000; % Kerf [m]
focus=[0 0 50]/1000; % Fixed focal point [m]
% Generate aperture for emission
emit_aperture = xdc_linear_array (256, lambda/2, element_height, kerf, 1, 1,focus);
% Set the impulse response and excitation of the emit aperture
impulse_response=sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*f0*(0:1/fs:2/f0));
xdc_excitation (emit_aperture, excitation);
% Generate aperture for reception
receive_aperture = xdc_linear_array (256, lambda/2, element_height, kerf, 1, 1,focus);
% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);
% Do phased array imaging
point_position=[0 0 70]/1000; % Position of the point to be imaged
no_lines=150; % Number of A-lines in image
sector=20 * pi/180; % Size of image sector
[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);
d_theta=sector/no_lines; % Increment in angle for 90 deg. image
% Pre-allocate some storage
image_data=zeros(1,no_lines);
theta= -sector/2;
dx=lambda/2;
for i=1:no_lines
% Set the focus for this direction
xdc_focus (emit_aperture, 0, [35*sin(theta) 0 70*cos(theta)]/1000);
xdc_focus (receive_aperture, 0, [35*sin(theta) 0 70*cos(theta)]/1000);
% Calculate the received response
[v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
% Store the result
image_data(1:max(size(v)),i)=v;
times(i) = t1;
% Steer in another angle
theta = theta + d_theta;
end
% Free space for apertures
xdc_free (emit_aperture)
xdc_free (receive_aperture)
% Adjust the data in time and display it as
% a gray scale image
min_sample=min(times)*fs;
for i=1:no_lines
rf_env=abs(hilbert([zeros(round(times(i)*fs-min_sample),1); image_data(:,i)]));
env(1:size(rf_env,1),i)=rf_env;
end
% make logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/fs*c/2;
x=((1:no_lines)-no_lines/2)*dx;
image(x*1000, depth*1000, env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('Image of cyst phantom (60 dB dynamic range)')
% Here the display of the data is inserted
%plot(image_data)
%% 
function [positions, amp] = cyst_phantom (N)
x_size = 40/1000; % Width of phantom [m] (x dim in the image)
y_size = 20/1000; % Transverse width of phantom [m] (thickness dimension of the phantom)
z_size = 40/1000; % Height of phantom [m]( depth dimension in the image)
z_start = 40/1000; % Start of phantom surface [m]( phantom starts at this depth)
% Create the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;
% Generate the amplitudes with a Gaussian distribution
amp=randn(N,1);
% Make the cyst and set the amplitudes to zero inside
r=10/1000; % Radius of cyst [m]
xc=0/1000; % Place of cyst [m]
zc=20/1000+z_start;
inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
amp = amp .* (1-inside);
% Place the point scatterers in the phantom
dz=z_size/10;
for i=N-9:N
x(i) = -15/1000;
y(i) = 0;
z(i) = z_start + (i-N+9)*dz;
amp(i) = 100;
end
% Return the variables
positions=[x y z];
end

