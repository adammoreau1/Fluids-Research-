%Clearing old code
clear;
clc;
clf;

%Reading in the Data
uvelfileID = fopen('c_orient','rb');
uvel = fread(uvelfileID,128*128*128,'double');
fclose(uvelfileID);
uve2 = reshape(uvel, [128 128 128]);

%Pre-sizing arrays to speed up assigning the data
uve = zeros(128);
uve3 = zeros(128);
uve4 = zeros(128);

%looping through indicies to switch i & j
for i = 1:128
  for j=1:128
    for k=1:128
        %Switching i & j
        uve(i,j,k) = uve2(j,i,k);
        %Switching i & j and then taking the vertical segments and mapping
        %them to the XY plane in order to run surf command
        uve3(i,j,k) = uve2(j,k,i);
        uve4(i,j,k) = uve2(k,i,j);
    end
  end
end

%Center Slice at index Z = 64
h = surf(64 * ones(size(uve2,1),size(uve2,2)), 'CData', uve(:,:,64), 'FaceColor', 'TextureMap');
set(h,'edgecolor','none');
hold on;

%Top Slice at index Z = 1
a = surf(ones(size(uve2,1),size(uve2,2)), 'CData', uve(:,:,1), 'FaceColor', 'TextureMap');
%Bottom Slice at index Z = 128
b = surf(128*ones(size(uve2,1),size(uve2,2)), 'CData', uve(:,:,128), 'FaceColor', 'TextureMap');
set(a,'edgecolor','none');
set(b,'edgecolor','none');

%Side Slice at index Y = 1
c = surf(ones(size(uve2,1),size(uve2,2)), 'CData', uve3(:,:,1), 'FaceColor', 'TextureMap');
%Side Slice at index Y = 128
d = surf(128*ones(size(uve2,1),size(uve2,2)), 'CData', uve3(:,:,128), 'FaceColor', 'TextureMap');
set(c,'edgecolor','none');
set(d,'edgecolor','none');

%Side Slice at index X = 1
e = surf(ones(size(uve2,1),size(uve2,2)), 'CData', uve4(:,:,1), 'FaceColor', 'TextureMap');
%Side Slice at index X = 128
f = surf(128*ones(size(uve2,1),size(uve2,2)), 'CData', uve4(:,:,128), 'FaceColor', 'TextureMap');
set(e,'edgecolor','none');
set(f,'edgecolor','none');

%Rotation Matricies
rotcd = [1,0,0];
rotef = [0,1,0];

%Rotating Vertical Slices around the associated rotation matrix
rotate(c,rotcd,90);
rotate(d,rotcd,90);
rotate(e,rotef,90);
rotate(f,rotef,90);

%Allowing for some transperentcy
alphas = 0.9;
alpha(a,alphas);
alpha(b,alphas);
alpha(c,alphas);
alpha(d,alphas);
alpha(e,alphas);
alpha(f,alphas);

%Creating a correctly sized matrix
oness = 1 * ones(size(uve2,1),size(uve2,2));
oneteight = 128 * oness;

%Resizing each array so the sides fit together well
set(a,'ZData',oness);
set(b,'ZData',oneteight);
set(c,'YData',oness);
set(d,'YData',oneteight);
set(e,'XData',oness);
set(f,'XData',oneteight);

%Setting up the display
colorbar;
lighting phong;
axis square;
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
title('Orientation Test')
