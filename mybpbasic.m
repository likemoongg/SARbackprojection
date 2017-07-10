function data = mybpbasic(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs a basic Backprojection operation. The         %
% following fields need to be populated:                               %
%                                                                      %
% data.Nfft: Size of the FFT to form the range profile                 %
% data.deltaF: Step size of frequency data (Hz)                        %
% data.minF: Vector containing the start frequency of each pulse (Hz)  %
% data.x_mat: The x-position of each pixel (m)                         %
% data.y_mat: The y-position of each pixel (m)                         %
% data.z_mat: The z-position of each pixel (m)                         %
% data.AntX: The x-position of the sensor at each pulse (m)            %
% data.AntY: The y-position of the sensor at each pulse (m)            %
% data.AntZ: The z-position of the sensor at each pulse (m)            %
% data.R0: The range to scene center (m)                               %
% data.phdata: Phase history data (frequency domain)                   %
% Fast time in rows, slow time in columns                              %
%                                                                      %
% The output is:                                                       %
% data.im_final: The complex image value at each pixel                 %
%                                                                      %
% Written by LeRoy Gorham, Air Force Research Laboratory, WPAFB, OH    %
% Email: leroy.gorham@wpafb.af.mil                                     %
% Date Released: 8 Apr 2010                                            %
%                                                                      %
% Gorham, L.A. and Moore, L.J., "SAR image formation toolbox for       %
% MATLAB," Algorithms for Synthetic Aperture Radar Imagery XVII        %
% 7669, SPIE (2010).                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define speed of light (m/s)
c = 299792458;

% Determine the size of the phase history data
data.K = size(data.phdata,1); % The number of frequency bins per pulse
data.Np = size(data.phdata,2); % The number of pulses

% Determine the azimuth angles of the image pulses (radians)
data.AntAz = unwrap(atan2(data.AntY,data.AntX));

% Determine the average azimuth angle step size (radians)
data.deltaAz = abs(mean(diff(data.AntAz)));

% Determine the total azimuth angle of the aperture (radians)
data.totalAz = max(data.AntAz) - min(data.AntAz);

% Determine the maximum scene size of the image (m)
data.maxWr = c/(2*data.deltaF);                            %c/2ΔF
data.maxWx = c/(2*data.deltaAz*mean(data.minF));           %c/（2Δθ*最小f均值）

% Determine the resolution of the image (m)
data.dr = c/(2*data.deltaF*data.K);                        %%c/2B
data.dx = c/(2*data.totalAz*mean(data.minF));              %%c/（2θ*最小f均值）

% Display maximum scene size and resolution
fprintf('Maximum Scene Size: %.2f m range, %.2f m cross range\n',data.maxWr,data.maxWx);
fprintf('Resolution: %.2fm range, %.2f m cross range\n',data.dr,data.dx);

% Calculate the range to every bin in the range profile (m)
data.r_vec = linspace(-data.Nfft/2,data.Nfft/2-1,data.Nfft)*data.maxWr/data.Nfft;

% Initialize the image with all zero values
data.im_final = zeros(size(data.x_mat));

% Set up a vector to keep execution times for each pulse (sec)
t = zeros(1,data.Np);

% Loop through every pulse
for ii = 1:data.Np

% Display status of the imaging process
if ii > 1
    t_sofar = sum(t(1:(ii-1)));
    t_est = (t_sofar*data.Np/(ii-1)-t_sofar)/60;%Np/ii-1=总时间/目前耗费时间，得出剩余时间
    fprintf('Pulse %d of %d, %.02f minutes remaining\n',ii,data.Np,t_est);
else
    fprintf('Pulse %d of %d\n',ii,data.Np);
end
tic

% Form the range profile with zero padding added
rc = fftshift(ifft(data.phdata(:,ii),data.Nfft));
%% 距离跃迁
qianyue(ii,:)=rc;
if ii==2
    figure
    plot(data.r_vec,20*log10(abs(rc')/max(abs(rc(:)))));
    axis tight;
    axis([-1 1 -50 0])
    title('After ifft');
    xlabel('距离向x（m）','fontsize',16)
    ylabel('归一化幅度（dB）','fontsize',16)
    figure
    guiyi=plot(data.r_vec,abs(rc')/max(abs(rc(:))));
    axis tight;
    xlim([-1,1])
    title('After ifft 方位角0.2°');
    xlabel('距离向x（m）','fontsize',16)
    ylabel('归一化幅度','fontsize',16)
    
    rc_temp=fftshift(ifft(data.phdata(:,ii)));
    vec_temp= linspace(-size(data.phdata,1)/2,size(data.phdata,1)/2-1,size(data.phdata,1))*data.maxWr/size(data.phdata,1);
    ivec_temp=linspace(-data.maxWr/2,data.maxWr/2,4096);
    irc_temp=interp1(vec_temp,rc_temp,ivec_temp,'spline');
    figure
    plot(ivec_temp,abs(irc_temp)/max(abs(irc_temp(:))));
    axis tight;
    xlim([-1,1])
    title('After ifft 方位角0.2°做了插值之后')
    
end
if ii==900
    figure
    guiyi=plot(data.r_vec,abs(rc')/max(abs(rc(:))));
    axis tight;
    xlim([-1,1])
    title('After ifft 方位角90°');
    xlabel('距离向x（m）','fontsize',16)
    ylabel('归一化幅度','fontsize',16)    
end


% Calculate differential range for each pixel in the image (m)
dR = sqrt((data.AntX(ii)-data.x_mat).^2 + ...
 (data.AntY(ii)-data.y_mat).^2 + ...
 (data.AntZ(ii)-data.z_mat).^2)-data.R0(ii);

% Calculate phase correction for image
phCorr = exp(1i*4*pi*data.minF(ii)/c*dR);

% Determine which pixels fall within the range swath
I = find(and(dR > min(data.r_vec), dR < max(data.r_vec)));

% Update the image using linear interpolation
data.im_final(I) = data.im_final(I) + interp1(data.r_vec,rc,dR(I),'linear') .* phCorr(I);

if ii==20
    figure
    ttemp=zeros(size(data.x_mat));
    ttemp(I)=ttemp(I)+interp1(data.r_vec,rc,dR(I),'linear') .* phCorr(I);
    surf(data.xtemp,data.ytemp,abs(ttemp));
    shading flat;
    xlabel('距离向')
    ylabel('方位向')
    title('在第20个脉冲组进行距离补偿之后的单组回波成像结果')
    figure
    surf(data.xtemp,data.ytemp,abs(data.im_final));
    shading flat;
    xlabel('距离向')
    ylabel('方位向')
    title('方位累加20个效果')
end

if ii==450
        figure
    ttemp=zeros(size(data.x_mat));
    ttemp(I)=ttemp(I)+interp1(data.r_vec,rc,dR(I),'linear') .* phCorr(I);
    surf(data.xtemp,data.ytemp,abs(ttemp));
    shading flat;
    xlabel('距离向')
    ylabel('方位向')
    title('在第450个脉冲组进行距离补偿之后的单组回波成像结果')
    figure
    surf(data.xtemp,data.ytemp,abs(data.im_final));
    shading flat;
    xlabel('距离向')
    ylabel('方位向')
    title('方位累加450个(共45°)效果(最终成像结果)') 
end
% Determine the execution time for this pulse
t(ii) = toc;
end
%% 距离跃迁成像
    figure
    surf(data.r_vec,(1:1:data.Np),abs(qianyue));
    xlim([-1,1])
    shading flat;
    title('距离跃迁')
    xlabel('距离向')
    ylabel('慢时间采样点数')
return
