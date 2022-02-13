clear;  
clc;
close all;
% Generates a hartmannograms for further processing in the WFS7

%%  SHWFS parameters
Lambda=0.65e-6;                 % wavelength [m]
k=2*pi/Lambda;                  % wavenumber

Pixels_X=500;                   % image pixels on the X axis
Pixels_Y=500;                   % image pixels on the Y axis
PixelSize=5.5e-6;               % pixel size [m]
Pitch=137.5e-6;                 % microlens size [m]
ML_focal=3.2e-3;                % focal dist [m]
Bits=8;                         % camera bit rate
ApertDiam=Pixels_X*PixelSize;   % diameter of the full aperture [m]
ApertRad=ApertDiam/2;           % radius of the full aperture [m]

L_X=Pixels_X*PixelSize;                    % sensor area on the X axis [m]
L_Y=Pixels_Y*PixelSize;                    % sensor area on the Y axis [m]
ML_Pixels=round(Pitch/PixelSize);          % pixels in CCD to each microlens

NumberLens_X=floor(Pixels_X/ML_Pixels);    % number of lenses in the array 
NumberLens_Y=floor(Pixels_Y/ML_Pixels);    % on the X and Y axis

% Set through the distance from the source to the sensor:
L_def=0.4;                                          %distance between the source and the SHWFS [m]
Defocus=(L_def-sqrt(L_def^2-(ApertDiam/2)^2))/2;    %in [m]  
disp(Defocus*10^6);

%%  Reference
IntXY_Ref=zeros(ML_Pixels);
H_ref={};
% Creating integration limits
ab = 1:ML_Pixels; cd = ab;
A = -Pitch/2+PixelSize*(ab-1);
B = -Pitch/2+PixelSize*ab;
C = -Pitch/2+PixelSize*(cd-1);
D = -Pitch/2+PixelSize*cd;
% Diffraction spot in one subaperture
for i=1:ML_Pixels
    Intensity=@(x,y) (sinc(Pitch*x/Lambda/ML_focal)).^2.*(sinc(Pitch*y/Lambda/ML_focal)).^2;
    for k=1:ML_Pixels
        IntXY_Ref(i,k)=integral2(Intensity,A(i),B(i),C(k),D(k));
    end
end
% Сonversion to digitized intensity in the camera (by bit)
MaxIntens=max(max(IntXY_Ref(:,:),[],1));
MinIntens=MaxIntens/(2^Bits-1);
IntensityPixelised=fix(IntXY_Ref(:,:)./MinIntens).*MinIntens/max(max(IntXY_Ref)); 
% Сombining into a single matrix (reference hartmannogram)
for N=1:NumberLens_X
    for M=1:NumberLens_Y
        H_ref{N,M}=IntensityPixelised;
    end
end
H_ref=cell2mat(H_ref);
% Output of the reference hartmannogram
figure('Name','Reference hartmannogram');
imshow(H_ref,[]);
% Save to the "Hartmannograms" folder 
imwrite(H_ref,'Hartmannograms\Hartmannogram_REF_500.bmp'); 

%%  Introduction slopes
% Сoordinates of subaperture centers 
ML_Center = Pitch/2;
k=1;
for i=0:NumberLens_X-1
    MLA_Centerx(k) = -L_X/2+i*Pitch+ML_Center;
    k=k+1;
end
MLA_Centery = MLA_Centerx;
% Matrices of centers on axes X and Y:
[MLA_CenterX,MLA_CenterY]=meshgrid(MLA_Centerx,MLA_Centery);

% Wavefront calculation (using the "deriv.m" program
[WF_derivX,WF_derivY] = deriv(MLA_CenterX,MLA_CenterY,L_def);

% Shift calculation
ShiftX = WF_derivX * ML_focal;
ShiftY = WF_derivY * ML_focal;

% For setting shifts individually:
% ShiftX = ones(20)*PixelSize;      % the offset is equal to the pixel size 
% ShiftY = zeros(20);               % the offset is equal to zero

%%  Hartmannogram of the distorted wavefront
tic
IntXY=zeros(ML_Pixels); 
H={};
for N=1:NumberLens_X
    if (mod(N,2)==0)
        wb=waitbar(N/NumberLens_X);
    end
    for M=1:NumberLens_Y
        Intensity=@(x,y) (sinc(Pitch*(x-ShiftX(N,M))/Lambda/ML_focal)).^2.*...
            (sinc(Pitch*(y-ShiftY(N,M))/Lambda/ML_focal)).^2;
        for i=1:ML_Pixels
            for k=1:ML_Pixels
                IntXY(i,k)=integral2(Intensity,A(i),B(i),C(k),D(k));
            end
        end
    % Сonversion to digitized intensity in the camera (by bit)
    MaxIntens=max(max(IntXY(:,:),[],1));
    MinIntens=MaxIntens/(2^Bits-1);
    IntensityPixelised=fix(IntXY(:,:)./MinIntens).*MinIntens/max(max(IntXY)); 
    % Recording an intensity pattern for one subaperture (N,M)
    H{N,M}=IntensityPixelised;
    end
end
delete(wb);
toc
% Сombining into a single matrix (hartmannogram)
H=cell2mat(H);
% Output of the reference hartmannogram
figure('Name','Processed hartmannogram');
imshow(H,[]);
% Save to the "Hartmannograms" folder 
imwrite(H,'Hartmannograms\Hartmannogram_DEF_500.bmp');

%%  Writing reference and shifts to a file 
% Data for reference
k=1;
for i=0:NumberLens_X-1
    MLA_Centerx_px(k) = (0+i*Pitch+ML_Center)/PixelSize;
    k=k+1;
end
MLA_Centery_px = MLA_Centerx_px;
[MLA_CenterX_px,MLA_CenterY_px]=meshgrid(MLA_Centerx_px,MLA_Centery_px);

MLA_CenterXr = reshape(MLA_CenterX_px',400,1);
MLA_CenterYr = reshape(MLA_CenterY_px',400,1);
LeftTopX = MLA_CenterXr - Pitch/(2*PixelSize);
LeftTopY = MLA_CenterYr - Pitch/(2*PixelSize);
RightBottomX = MLA_CenterXr + Pitch/(2*PixelSize);
RightBottomY = MLA_CenterYr + Pitch/(2*PixelSize);
Data = zeros(2400,1);
Data(1:6:end,:) = MLA_CenterXr;
Data(2:6:end,:) = MLA_CenterYr;
Data(3:6:end,:) = LeftTopX;
Data(4:6:end,:) = LeftTopY;
Data(5:6:end,:) = RightBottomX;
Data(6:6:end,:) = RightBottomY;

% Data for shifts
Shift_X = reshape((ShiftX/PixelSize)',400,1);
Shift_Y = reshape((ShiftY/PixelSize)',400,1);
New_CenterX = MLA_CenterXr + Shift_X;
New_CenterY = MLA_CenterYr + Shift_Y;
Shift_Centers = zeros(800,1);
Shift_Centers(1:2:end,:) = New_CenterX;
Shift_Centers(2:2:end,:) = New_CenterY;

% Creating a reference file
fig1 = fopen('Files_txt\Reference.txt', 'w'); 
fprintf(fig1, '%d\n',Pixels_X);             % image width in pixels
fprintf(fig1, '%d\n',Pixels_Y);             % image height in pixels
fprintf(fig1, '%5.6f\n',Lambda*10^6);       % wavelength [um]
fprintf(fig1, '%5.6f\n',PixelSize*10^6);    % pixel size [um]
fprintf(fig1, '%5.6f\n',ML_focal*10^3);     % focal length [mm]
fprintf(fig1, '%5.6f\n',Pitch*10^3);        % pitch size [mm]
fprintf(fig1, '%1.10f\n',Data); 
fclose(fig1);

% Creating a processed file
fig2 = fopen('Files_txt\Shifted_centers.txt', 'w'); 
fprintf(fig2, '%5.10f\n', Shift_Centers);   % shifting the centers [m]
fclose(fig2);