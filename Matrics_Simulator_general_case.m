close all; clear; clc;
% Generates a hartmannograms for further processing in the WFS7

%%  SHWFS parameters
Lambda=0.65e-6;                 % wavelength [m]
k=2*pi/Lambda;                  % wavenumber

Pixels_X=500;                   % image pixels on the X axis
Pixels_Y=500;                   % image pixels on the Y axis
PixelSize= 5.5e-6;              % pixel size [m]
Pitch=136e-6;                   % microlens size [m]
ML_focal=3.2e-3;                % focal dist [m]
Bits=8;                         % camera bit rate
% ApertDiam=2.75e-3;
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

IX =  P0_X + i_rep;
IY =  P0_Y + k_rep;

Y = IX - X_0;
X = IY - Y_0;

VAL = ((sinc(Pitch*PixelSize*X/Lambda/ML_focal)).^2.*(sinc(Pitch*PixelSize*Y/Lambda/ML_focal)).^2);
MaxIntens=max(max(VAL(:,:),[],1));
MinIntens=MaxIntens/(2^Bits-1);
IntensityPixelised_VAL=fix(VAL(:,:)./MinIntens).*MinIntens/max(max(VAL));

% Cropping pixels
if Cut_flag == 1
    IntensityPixelised_VAL(1:floor(ML_Pixels/2),:) = [];
    IntensityPixelised_VAL(:,1:floor(ML_Pixels/2)) = [];
    All_pixels = size(IntensityPixelised_VAL,1);
    IntensityPixelised_VAL(Pixels_X+1:All_pixels,:) = [];
    IntensityPixelised_VAL(:, Pixels_X+1:All_pixels) = [];
end

figure('Name','Reference hartmannogram');
imshow(IntensityPixelised_VAL,[]);
% Save to the "Hartmannograms" folder 
imwrite(IntensityPixelised_VAL,'Hartmannograms\REF_500_136_cut.bmp'); 
                
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

% Cropping pixels
if Cut_flag == 1
    IntensityPixelised_VAL_shift(1:floor(ML_Pixels/2),:) = [];
    IntensityPixelised_VAL_shift(:,1:floor(ML_Pixels/2)) = [];
    All_pixels = size(IntensityPixelised_VAL_shift,1);
    IntensityPixelised_VAL_shift(Pixels_X+1:All_pixels,:) = [];
    IntensityPixelised_VAL_shift(:, Pixels_X+1:All_pixels) = [];
end

figure('Name','Processed hartmannogram');
imshow(IntensityPixelised_VAL_shift,[]);
% Save to the "Hartmannograms" folder 
imwrite(IntensityPixelised_VAL_shift,'Hartmannograms\DEF_500_136_cut.bmp'); 
