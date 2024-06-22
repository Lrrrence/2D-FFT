cd('E:\OneDrive - University of Southampton\MATLAB Scripts\2D-FFT')

set(0,'DefaultLineLineWidth',1)
set(groot,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

close all
myfile = '1MHz_A1_fullwedge_y.csv';
M1 = readtable(myfile,'ReadVariableNames', false);
A1 = table2array(M1);
A1(1,:) = [];

A1 = A1(8:end,:);   %skip header
Fs = 1/((A1(2,1))-(A1(1,1)));   %sample rate
St = 1/Fs;     % Sampling period
A1 = A1(:,2:end); % remove time column

%% calc 2DFFT

%filterbw = [0.8e6 2e6]; %filter bandwidth
%A1 = bandpass(A1, filterbw, Fs,'ImpulseResponse','iir'); % filter 

A2 = padarray(A1(1:end,:),[100000 2000], 0,'post'); % zero pad

Nx = size(A2,2); % Number of samples collected along first dimension
Nt = size(A2,1); % Number of samples collected along second dimension
dx = 0.0008;  % Distance increment (i.e., Spacing between each column)
dt = St; % Time increment (i.e., Spacing between each row)
Nyq_k = 1/dx; % Nyquist of data in first dimension
Nyq_f = 1/dt; % Nyquist of data in second dimension
dk = 1/(Nx*dx);   % Wavenumber increment
df = 1/(Nt*dt);   % Frequency increment
k = -Nyq_k/2:dk:Nyq_k/2-dk; % wavenumber (m)
kr = k*(2*pi)/1000; % wavenumber (rad/mm)
f =-Nyq_f/2:df:Nyq_f/2-df; % frequency (Hz)

krhalf = kr(:,((size(kr,2))+1)/2:end); % take positive quadrant of kr
fhalf = f(:,((size(f,2))+1)/2:end); % take positive quadrant of f

fft2result = fftshift(fft2(A2))*dx*dt; %fft2
fft2resultB =fft2result.'; %transpose
fft2resultB = flip(fft2resultB,2); %flip

z = abs(fft2resultB); %absolute value

fft2data = z((size(z,1)+1)/2:end,(size(z,2)+1)/2:end); % take positive quadrant of z

znorm = rescale(fft2data,0,1); % normalise 0-1

%% plot
figure(1);

imagesc(fhalf/1e6,krhalf,znorm);
set(gca,'YDir','normal')
colorbar;
%colormap(brewermap([],'Spectral'))
colormap(flipud(colormap('viridis_white')))

xlim([0 2])
ylim([0 3])

ylabel('Wavenumber (rad/mm)')
xlabel('Frequency (MHz)')
hold on
%grid on
box on
ax = gca;
ax.LineWidth = 2;

multiplier = 2.5;

plot(A_Lamb.("A0 fd (MHz*mm)")/multiplier,A_Lamb.("A0 Wavenumber (rad/mm)"),'-','Color','k')
plot(S_Lamb.("S0 fd (MHz*mm)")/multiplier,S_Lamb.("S0 Wavenumber (rad/mm)"),'--','Color','k')
plot(A_Lamb.("A1 fd (MHz*mm)")/multiplier,A_Lamb.("A1 Wavenumber (rad/mm)"),'-','Color','k')
plot(S_Lamb.("S1 fd (MHz*mm)")/multiplier,S_Lamb.("S1 Wavenumber (rad/mm)"),'--','Color','k')

legend('A_{0}','S_{0}','A_{1}','S_{1}')

x0=0;
y0=0;
width=1200;
height=1200;
set(gcf,'position',[x0,y0,width,height])
set(gca,'fontsize', 20)
set(gca, 'FontName', 'HelveticaLTStd-Roman')
set(gcf, 'Color', 'w');

export_fig 2DFFT_A0_2.5mm.png -m2
