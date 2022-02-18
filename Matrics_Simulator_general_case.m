close all; clear; clc;
% Generates a hartmannograms for further processing in the WFS7

%%  SHWFS parameters
Lambda=0.65e-6;                 % wavelength [m]
k=2*pi/Lambda;                  % wavenumber

Pixels_X=500;                   % image pixels on the X axis
Pixels_Y=500;                   % image pixels on the Y axis
PixelSize= 5.5e-6;               % pixel size [m]
Pitch=136e-6;                 % microlens size [m]
ML_focal=3.2e-3;                % focal dist [m]
Bits=8;                         % camera bit rate
ApertDiam=Pixels_X*PixelSize;   % diameter of the full aperture [m]
ApertRad=ApertDiam/2;           % radius of the full aperture [m]

L_X=Pixels_X*PixelSize;                    % sensor area on the X axis [m]
L_Y=Pixels_Y*PixelSize;                    % sensor area on the Y axis [m]
ML_Pixels=round(Pitch/PixelSize);          % pixels in CCD to each microlens

NumberLens_X=ceil(Pixels_X/ML_Pixels);    % number of lenses in the array 
NumberLens_Y=ceil(Pixels_Y/ML_Pixels);    % on the X and Y axis

Cut_flag = 1;                             % if = 1 - cut pixels outside Pixels_X(Y)

%%  Setting the zernike coefficients
% Set through the distance from the source to the sensor:
L_def=1.65;                                         %distance between the source and the SHWFS [m]
Defocus=(L_def-sqrt(L_def^2-(ApertDiam/2)^2))/2;    %in [m]  

% Set the parameters individually:
X_Tilt = 0;
Defocus = 1e-6;

% Aberration coefficients
z11=X_Tilt;         % X-Tilt
zi11=0;             % Y-Tilt
z02=Defocus;        % Defocus 
z22=0;              % Vertical astigmatism
zi22=0;             % Oblique astigmatism
z13=0;              % Horizontal coma
zi13=0;             % Vertical coma
z04=0;              % Primary spherical

%%  Reference
iRows = ceil(Pixels_X/(2*ML_Pixels)); 
iCols = iRows;
iRefRows = 2*iRows + 1;
iRefCols = iRefRows;

for n = 0 : 2*iRows
    dx(n+1) = Pitch/PixelSize*(iRows - n);
end
x_0 = Pixels_X/2 - dx;
y_0 = x_0;
x_0_rep = repelem(x_0,ML_Pixels);
y_0_rep = repelem(y_0,ML_Pixels);
[X_0,Y_0] = meshgrid(x_0_rep,y_0_rep);
p0_x = round(y_0 - ML_Pixels/2);
p0_y = round(x_0 - ML_Pixels/2);
P0_x = repelem(p0_x,ML_Pixels);
P0_y = repelem(p0_y,ML_Pixels);
[P0_X,P0_Y] = meshgrid(P0_x,P0_y);
i = 1:ML_Pixels;
k = 1:ML_Pixels;
[I,K] = meshgrid(i,k);
i_rep = repmat(I,NumberLens_X+1,NumberLens_X+1);
k_rep = repmat(K,NumberLens_X+1,NumberLens_X+1);

IX =  P0_X + k_rep - 1;
IY =  P0_Y + i_rep - 1;

X = IX - X_0;
Y = IY - Y_0;

VAL = ((sinc(Pitch*PixelSize*X/Lambda/ML_focal)).^2.*(sinc(Pitch*PixelSize*Y/Lambda/ML_focal)).^2);
MaxIntens=max(max(VAL(:,:),[],1));
MinIntens=MaxIntens/(2^Bits-1);
IntensityPixelised_VAL=fix(VAL(:,:)./MinIntens).*MinIntens/max(max(VAL));

if Cut_flag == 1
    All_pixels = size(IntensityPixelised_VAL,1);
    Cut_pixels = All_pixels - Pixels_X-1;
    IntensityPixelised_VAL((All_pixels - Cut_pixels/2):All_pixels,:) = [];
    IntensityPixelised_VAL(:, (All_pixels - Cut_pixels/2):All_pixels) = [];
    IntensityPixelised_VAL(1:Cut_pixels/2,:) = [];
    IntensityPixelised_VAL(:,1:Cut_pixels/2) = [];
end

figure('Name','Reference hartmannogram');
imshow(IntensityPixelised_VAL,[]);
% Save to the "Hartmannograms" folder 
imwrite(IntensityPixelised_VAL,'Hartmannograms\REF_2048_136_val.bmp'); 
                
%%  Introduction slopes
% Сoordinates of subaperture centers 
ML_Center = Pitch/2;
k=1;
for i=0:NumberLens_X
    MLA_Centerx(k) = -L_X/2+i*Pitch+ML_Center;
    k=k+1;
end
MLA_Centery = MLA_Centerx;
% Matrices of centers on axes X and Y:
[MLA_CenterX,MLA_CenterY]=meshgrid(MLA_Centerx,MLA_Centery);

% Wavefront calculation (using the "zernike_8_deriv.m" program
[WF_derivY,WF_derivX]=zernike_8_deriv(MLA_Centerx/ApertRad,...
    MLA_Centery/ApertRad,z11,zi11,z02,z22,zi22,z13,zi13,z04);

% Shift calculation
ShiftX=WF_derivX * ML_focal/ApertRad/PixelSize;
ShiftY=WF_derivY * ML_focal/ApertRad/PixelSize;

%%  Hartmannogram of the distorted wavefront
ShiftX_rep = repelem(ShiftX,ML_Pixels);
ShiftY_rep = repelem(ShiftY,ML_Pixels);
[SHIFT_X,SHIFT_Y] = meshgrid(ShiftX_rep,ShiftY_rep);

VAL_shift = (sinc(Pitch*PixelSize*(X+SHIFT_Y)/Lambda/ML_focal)).^2.*...
    (sinc(Pitch*PixelSize*(Y+SHIFT_X)/Lambda/ML_focal)).^2;
MaxIntens=max(max(VAL_shift(:,:),[],1));
MinIntens=MaxIntens/(2^Bits-1);
IntensityPixelised_VAL_shift=fix(VAL_shift(:,:)./MinIntens).*MinIntens/max(max(VAL_shift));

if Cut_flag == 1
    All_pixels = size(IntensityPixelised_VAL_shift,1);
    Cut_pixels = floor(All_pixels - Pixels_X);
    IntensityPixelised_VAL_shift((All_pixels - Cut_pixels/2):All_pixels,:) = [];
    IntensityPixelised_VAL_shift(:, (All_pixels - Cut_pixels/2):All_pixels) = [];
    IntensityPixelised_VAL_shift(1:Cut_pixels/2,:) = [];
    IntensityPixelised_VAL_shift(:,1:Cut_pixels/2) = [];
end

figure('Name','Processed hartmannogram');
imshow(IntensityPixelised_VAL_shift,[]);
% Save to the "Hartmannograms" folder 
imwrite(IntensityPixelised_VAL_shift,'Hartmannograms\DEF_2048_136_val_shift.bmp'); 

%%  Writing reference and shifts to a file 
% % Data for reference
% k=1;
% for i=0:NumberLens_X-1
%     MLA_Centerx_px(k) = (0+i*Pitch+ML_Center)/PixelSize;
%     k=k+1;
% end
% MLA_Centery_px = MLA_Centerx_px;
% [MLA_CenterX_px,MLA_CenterY_px]=meshgrid(MLA_Centerx_px,MLA_Centery_px);
% 
% MLA_CenterXr = reshape(MLA_CenterX_px',400,1);
% MLA_CenterYr = reshape(MLA_CenterY_px',400,1);
% LeftTopX = MLA_CenterXr - Pitch/(2*PixelSize);
% LeftTopY = MLA_CenterYr - Pitch/(2*PixelSize);
% RightBottomX = MLA_CenterXr + Pitch/(2*PixelSize);
% RightBottomY = MLA_CenterYr + Pitch/(2*PixelSize);
% Data = zeros(2400,1);
% Data(1:6:end,:) = MLA_CenterXr;
% Data(2:6:end,:) = MLA_CenterYr;
% Data(3:6:end,:) = LeftTopX;
% Data(4:6:end,:) = LeftTopY;
% Data(5:6:end,:) = RightBottomX;
% Data(6:6:end,:) = RightBottomY;
% 
% % Data for shifts
% Shift_X = reshape((ShiftX/PixelSize)',400,1);
% Shift_Y = reshape((ShiftY/PixelSize)',400,1);
% New_CenterX = MLA_CenterXr + Shift_X;
% New_CenterY = MLA_CenterYr + Shift_Y;
% Shift_Centers = zeros(800,1);
% Shift_Centers(1:2:end,:) = New_CenterX;
% Shift_Centers(2:2:end,:) = New_CenterY;
% 
% % Creating a reference file
% fig1 = fopen('Files_txt\Reference.txt', 'w'); 
% fprintf(fig1, '%d\n',Pixels_X);             % image width in pixels
% fprintf(fig1, '%d\n',Pixels_Y);             % image height in pixels
% fprintf(fig1, '%5.6f\n',Lambda*10^6);       % wavelength [um]
% fprintf(fig1, '%5.6f\n',PixelSize*10^6);    % pixel size [um]
% fprintf(fig1, '%5.6f\n',ML_focal*10^3);     % focal length [mm]
% fprintf(fig1, '%5.6f\n',Pitch*10^3);        % pitch size [mm]
% fprintf(fig1, '%1.10f\n',Data); 
% fclose(fig1);
% 
% % Creating a processed file
% fig2 = fopen('Files_txt\Shifted_centers.txt', 'w'); 
% fprintf(fig2, '%5.10f\n', Shift_Centers);   % shifting the centers [m]
% fclose(fig2);