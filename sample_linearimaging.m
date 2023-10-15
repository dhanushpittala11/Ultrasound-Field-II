fc=2.5e6; % center frequency [Hz] of the transducer
f_sampl=100e6; % Sampling frequency [Hz]
c=1540; % Speed of sound [m/s]
lambda=c/fc; % Wavelength of the acoustic wave [m]

dim_x= lambda; % Width of element[m]
dim_y=5/1000; % Height of element [m]
w_x= dim_x/4; % Kerf gap between two physical elements[m]
focus=[0 0 50]/1000; % Fixed focal point [m]
N_elements=256; % Number of elements in the transducer
N_active=75; % Number of Active elements in the transducer

% Setting the sampling frequency
set_sampling(f_sampl);
% Generate aperture for emission
emit_aperture = xdc_linear_array (N_elements, dim_x, dim_y, w_x, 1, 3, focus);
waveform1=sin(2*pi*fc*(0:1/f_sampl:4/fc));
element_no=1;
ele_waveform (emit_aperture, element_no, waveform1);
[h,t] = calc_h (emit_aperture,[0 0 60]/1000);
plot((0:length(h)-1)/f_sampl+t,h)
ylabel('Response')
xlabel('Time [s]')
%% displaying the geometry of the emit aperture of the transducer
data = xdc_get(e_aperture,'rect');
[n,m]=size(data);

colormap(hot(128));
for i=1:m
x=[data(11,i), data(20,i); data(14,i), data(17,i)]*1000;
y=[data(12,i), data(21,i); data(15,i), data(18,i)]*1000;
z=[data(13,i), data(22,i); data(16,i), data(19,i)]*1000;
c= 1;
hold on
surf(x,y,z,c)
end
% Put som axis legends on
Hc = colorbar;
view(3)
xlabel('x [mm]')
ylabel('y [mm]')
zlabel('z [mm]')
grid
axis('image')
hold off
%% 
% Setting the impulse response and excitation of the emit aperture


impulse_response=sin(2*pi*fc*(0:1/f_sampl:2/fc));
impulse_response=impulse_response.*hanning(max(size(impulse_response))).';
xdc_impulse (emit_aperture, impulse_response);
excitation=sin(2*pi*fc*(0:1/f_sampl:2/fc));
xdc_excitation (emit_aperture, excitation);

% Generating aperture for reception
receive_aperture = xdc_linear_array (N_elements, dim_x, dim_y, w_x, 1, 3, focus);
% Set the impulse response for the receive aperture
xdc_impulse (receive_aperture, impulse_response);

% Loading the computer generated phantom
%[phantom_positions, phantom_amplitudes] = cyst_phantom(10000);
%we can also use the below command for scanning point phantom
%to know how the scan looks like
[phantom_positions, phantom_amplitudes] = point_phantom(30);
% performing linear array imaging

n_lines=N_elements-N_active+1; % Number of A-lines in image

dx=dim_x; % Increment for the image
z_focus=60/1000;

% Pre-allocate some storage
image_data=zeros(1,n_lines);
for i=1:n_lines
    x=(i-1-n_lines/2)*dx;
% Set the focus for this direction
    xdc_center_focus (emit_aperture, [x 0 0]);
    xdc_focus (emit_aperture, 0, [x 0 z_focus]);
    xdc_center_focus (receive_aperture, [x 0 0]);
    xdc_focus (receive_aperture, 0, [x 0 z_focus-20]);
% Set the active elements using the apodization
    apo=[zeros(1, i-1) hamming(N_active)' zeros(1, N_elements-N_active-i+1)];
    xdc_apodization (emit_aperture, 0, apo);
    xdc_apodization (receive_aperture, 0, apo);
% Calculate the received response
    [v, t1]=calc_scat(emit_aperture, receive_aperture, phantom_positions, phantom_amplitudes);
    % Store the result
    image_data(1:max(size(v)),i)=v;
    times(i) = t1;
 end
% Free space for apertures
xdc_free (emit_aperture)
xdc_free (receive_aperture)

% Adjust the data in time and display it as
% a gray scale image
min_sample=min(times)*f_sampl;
for i=1:n_lines
rf_env=abs(hilbert([zeros(round(times(i)*f_sampl-min_sample),1); image_data(:,i)]));
env(1:size(rf_env,1),i)=rf_env;
end
% making logarithmic compression to a 60 dB dynamic range
% with proper units on the axis
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=127*(env_dB+60)/60;
depth=((0:size(env,1)-1)+min_sample)/f_sampl*c/2;
x=((1:n_lines)-n_lines/2)*dx;
image(x*1000, depth*1000, env_gray)
xlabel('Lateral distance [mm]')
ylabel('Depth [mm]')
axis('image')
colormap(gray(128))
title('Image of point phantom (60 dB dynamic range)')

%% function for a computer generated phantom with a cyst 
function [positions, amp] = cyst_phantom (N)
x_size = 40/1000; % Width of phantom [m] (x dim in the image)
y_size = 20/1000; % Transverse width of phantom [m] (thickness dimension of the phantom)
z_size = 40/1000; % Height of phantom [m]( depth dimension in the image)
z_start = 40/1000; % Start of phantom surface [m]( phantom starts at this depth)

% Creating the general scatterers
x = (rand (N,1)-0.5)*x_size;
y = (rand (N,1)-0.5)*y_size;
z = rand (N,1)*z_size + z_start;

% Generating the amplitudes with a Gaussian distribution
amp=randn(N,1);

% Making the cyst and setting the amplitudes to zero inside.
r=10/1000; % Radius of cyst [m]
xc=0/1000; % Place of cyst [m]
zc=20/1000+z_start;
inside = ( ((x-xc).^2 + (z-zc).^2) < r^2);
amp = amp .* (1-inside);
% Placing the point scatterers in the phantom
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
%% function for generating the point phantom
function [positions, amp] = point_phantom(N)

dz=5/1000;        %  Distance between points [m]
z_start=10/1000;  %  Start of point [m]
       %  Number of point targets

%  Create the point scatterers positions

positions = [zeros(1,N); zeros(1,N); (1:N)*dz+z_start]';
amp=ones(N,1);

end
