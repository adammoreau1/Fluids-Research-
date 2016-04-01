%**************************************************************************
%********************* Setup of Data & Display ****************************
%**************************************************************************

%Clearing old Code
clear;
clc;
clf;
cla;

%Reading in the Data
uvelfileID = fopen('corr2pt3dq.cbin','rb');
uvel = fread(uvelfileID,128*128*128,'double');
fclose(uvelfileID);
uve = reshape(uvel, [128 128 128]);

%//////////////////////////////////////////////////////
%/////// From here down i & j should be correct ///////
%//////////////////////////////////////////////////////

%Setting the Searchvalue (User Defined)
%range: ~0.048min (not very good...) -- ~0.7max (also pretty bad)
searchvalue = 0.1

%Setting the Alpha Value [0-1] (Transperency)
alphas = 0.4;

%Setting up the display
lighting phong;
axis vis3d;
xlabel('X Axis')
ylabel('Y Axis')
zlabel('Z Axis')
title('Elipsoid Fitting w/ Axis Remapped')
ax = gca;
ax.XTick = [-pi/2 0 pi/2];
ax.YTick = [-pi/2 0 pi/2];
ax.ZTick = [-pi/2 0 pi/2];
ax.XLim = [-pi/2 pi/2];
ax.YLim = [-pi/2 pi/2];
ax.ZLim = [-pi/2 pi/2];
ax.XTickLabel = {'-\pi/2','0','\pi/2'};
ax.YTickLabel = {'-\pi/2','0','\pi/2'};
ax.ZTickLabel = {'-\pi/2','0','\pi/2'};
grid on;
hold on;

%**************************************************************************
%*******************Tri-liner Interpolation Routine************************
%**************************************************************************

%The Interpolation works by reading in the isosurface values and then
%finding when two cross the searchvalue. The loops are run twice, once to
%count the number of elements (it is dependent upon the searchvalue
%selected above) and to then pre-size the arrays when are filled on the
%second loop through

%Initializing all values to zero
lastvalue = 0;
lastjvalue = 0;
lastivalue = 0;
newvalue = 0;
count = 0;
place = 0;
knew=0;
jnew=0;
inew =0;
maxi = 0;
maxj = 0;
maxk = 0;

%Loops
%Loop to dertermine number of elements or fill arrays
for n = 1:2
    % 3 loops through i,j,k axies
    for i = 1:128
        for j=1:128
            for k=1:128
                %setting the newvalue
                newvalue = uve(i,j,k);
             
                %Comparing newvalue to the search value (Lower to Higher)
                if (newvalue <= searchvalue) && (lastvalue >= searchvalue)
                    %If a crossing is found in the first loop it is counted
                    if n == 1
                        count = count + 1;
                    end
                    %If it is the second loop it will be interpolated
                    if n == 2
                        %place marks the index of the X,Y,Z arrays and
                        %increments each time
                        place = place + 1;
                        
                        %Interpolation in i,j,k
                        knew= k  - 1 +(searchvalue-lastvalue)/(newvalue-lastvalue);
                        %jnew= j  - 1 +(searchvalue-lastjvalue)/(newvalue-lastjvalue);
                        %inew= i  - 1 +(searchvalue-lastivalue)/(newvalue-lastivalue);
                        
                        %Putting the interpolated values into X,Y,Z arrays
                        X(place) = 1.0*i;
                        Y(place) = 1.0*j;
                        Z(place) = 1.0*knew;
                    end
                end
                %Comparing newvalue to the search value (Higher to Lower)
                if (newvalue >= searchvalue) && (lastvalue <= searchvalue)
                    %If a crossing is found in the first loop it is counted
                    if n == 1
                        count = count + 1;
                    end
                    %If it is the second loop it will be interpolated
                    if n == 2
                        %place marks the index of the X,Y,Z arrays and
                        %increments each time
                        place = place + 1;
                        
                        %Interpolation in i,j,k
                        knew= k - 1 + (searchvalue-lastvalue)/(newvalue-lastvalue);
                                             
                        %Putting the interpolated values into X,Y,Z arrays
                        X(place) = 1.0*i;
                        Y(place) = 1.0*j;
                        Z(place) = 1.0*knew;
                        
                    end
                end
                %assigning the current value looped through to the lastvalue
                lastvalue = newvalue;
            end
        end
    end
    %Loops through 
    for k = 1:128
        for i=1:128
            for j=1:128
                %setting the newvalue
                newvalue = uve(i,j,k);
             
                %Comparing newvalue to the search value (Lower to Higher)
                if (newvalue <= searchvalue) && (lastvalue >= searchvalue)
                    %If a crossing is found in the first loop it is counted
                    if n == 1
                        count = count + 1;
                    end
                    %If it is the second loop it will be interpolated
                    if n == 2
                        %place marks the index of the X,Y,Z arrays and
                        %increments each time
                        place = place + 1;
                        
                        %Interpolation in i,j,k
                        jnew= j  - 1 +(searchvalue-lastvalue)/(newvalue-lastvalue);
                                                
                        %Putting the interpolated values into X,Y,Z arrays
                        X(place) = 1.0*i;
                        Y(place) = 1.0*jnew;
                        Z(place) = 1.0*k;
                    end
                end
                %Comparing newvalue to the search value (Higher to Lower)
                if (newvalue >= searchvalue) && (lastvalue <= searchvalue)
                    %If a crossing is found in the first loop it is counted
                    if n == 1
                        count = count + 1;
                    end
                    %If it is the second loop it will be interpolated
                    if n == 2
                        %place marks the index of the X,Y,Z arrays and
                        %increments each time
                        place = place + 1;
                        
                        %Interpolation in i,j,k
                        jnew= j - 1 + (searchvalue-lastvalue)/(newvalue-lastvalue);
                                             
                        %Putting the interpolated values into X,Y,Z arrays
                        X(place) = 1.0*i;
                        Y(place) = 1.0*jnew;
                        Z(place) = 1.0*k;
                        
                    end
                end
                %assigning the current value looped through to the lastvalue
                lastvalue = newvalue;
            end
        end
    end
    for i = 1:128
        for j=1:128
            for k=1:128
                %setting the newvalue
                newvalue = uve(i,j,k);
             
                %Comparing newvalue to the search value (Lower to Higher)
                if (newvalue <= searchvalue) && (lastvalue >= searchvalue)
                    %If a crossing is found in the first loop it is counted
                    if n == 1
                        count = count + 1;
                    end
                    %If it is the second loop it will be interpolated
                    if n == 2
                        %place marks the index of the X,Y,Z arrays and
                        %increments each time
                        place = place + 1;
                        
                        %Interpolation in i,j,k
                        inew= i  - 1 +(searchvalue-lastvalue)/(newvalue-lastvalue);
                        
                        %Putting the interpolated values into X,Y,Z arrays
                        X(place) = 1.0*inew;
                        Y(place) = 1.0*j;
                        Z(place) = 1.0*k;
                    end
                end
                %Comparing newvalue to the search value (Higher to Lower)
                if (newvalue >= searchvalue) && (lastvalue <= searchvalue)
                    %If a crossing is found in the first loop it is counted
                    if n == 1
                        count = count + 1;
                    end
                    %If it is the second loop it will be interpolated
                    if n == 2
                        %place marks the index of the X,Y,Z arrays and
                        %increments each time
                        place = place + 1;
                        
                        %Interpolation in i,j,k
                        inew= i - 1 + (searchvalue-lastvalue)/(newvalue-lastvalue);
                                             
                        %Putting the interpolated values into X,Y,Z arrays
                        X(place) = 1.0*inew;
                        Y(place) = 1.0*j;
                        Z(place) = 1.0*k;
                        
                    end
                end
                %assigning the current value looped through to the lastvalue
                lastvalue = newvalue;
            end
        end
    end
    
    %After looping through all i,j,k elements this will size X,Y,Z to the
    %total count
    if n == 1
        X = [1:count];
        Y = [1:count];
        Z = [1:count];
    end
end

%**************************************************************************
%*********************Data Display and analysis****************************
%**************************************************************************

%print the total number of elements in the matrix
display(count);

%Formating Data for placement into ellipsoid_fit
out = [X; Y; Z];
out = out.';
xpi = pi/64 * X - pi;
ypi = pi/64 * Y - pi;
zpi = pi/64 * Z - pi;
out2 = [xpi; ypi; zpi];
out3 = out2.';

%Display a scatter plot of Points
scatter3( xpi, ypi, zpi, '.');

%ellipsoid_fit
[ center, radii, evecs, v ] = ellipsoid_fit( out3 );

%Print out the information about the ellipsoid
fprintf ('The center of the Ellipsoid is given by ( %.3g, %.3g, %.3g )', center );
fprintf( '\n\nEllipsoid radii : a = %.3g b = %.3g c = %.3g', radii );
fprintf( '\n\nEigenvectors :\n %.3g %.3g %.3g\n%.3g %.3g %.3g\n%.3g %.3g %.3g\n', ...
    evecs(1), evecs(2), evecs(3), evecs(4), evecs(5), evecs(6), evecs(7), evecs(8), evecs(9) );
fprintf( '\n' );
fn = '\nThe function of the Ellipsoid is given as: \n\n 1 = %f X^2 + %f Y^2 + %f Z^2 + %f *2*XY + %f *2*XZ + %f *2*YZ + %f *2*X + %f *2*Y +%f *2*Z';
fprintf(fn, v);
fprintf('\n\n');

%Using a meshgrid to visualize the data
radii2 = [radii(2), radii(1), radii(3)];
maxd = max( radii2 );
step = maxd / 50;
[ x, y, z ] = meshgrid( -maxd:step:maxd , -maxd:step:maxd , -maxd:step:maxd );

Ellipsoid = v(1) *x.*x +   v(2) * y.*y + v(3) * z.*z + ...
          2*v(4) *x.*y + 2*v(5)*x.*z + 2*v(6) * y.*z + ...
          2*v(7) *x    + 2*v(8)*y    + 2*v(9) * z;
p = patch( isosurface( x, y, z, Ellipsoid, 1 ) );
set( p, 'FaceColor', 'g', 'EdgeColor', 'none' );
alpha(p,alphas)

%**************************************************************************
%***************Eigenvector Interpretation and Display*********************
%**************************************************************************

%Rotating the Center to be horizantal
%Each of these values come from the output of ellipsoid fit, and thus will
%be slightly greater than 0
center1 = [center(1),center(2),center(3)];

% %Multiplying the colums of evecs the rows of radii
evax = evecs(1) * radii(1);
evay = evecs(2) * radii(1);
evaz = evecs(3) * radii(1);

evbx = evecs(4) * radii(2);
evby = evecs(5) * radii(2);
evbz = evecs(6) * radii(2);

evcx = evecs(7) * radii(3);
evcy = evecs(8) * radii(3);
evcz = evecs(9) * radii(3);

% %Adding the Center to each of the dirrection eigenvectors
eva = [evax + center(1), evay + center(2), evaz + center(3)];
evb = [evbx + center(1), evby + center(2), evbz + center(3)];
evc = [evcx + center(1), evcy + center(2), evcz + center(3)];

% %Creating a line from the center to the point at the end of the eigenvector
evlinea = [center1;eva]; %get the angle inv tan of z/x (~0) & y/x (y/x is most important)
evlineb = [center1;evb];
evlinec = [center1;evc];

%Solving for the angles
ang = atand(eva(2)/eva(1))
angz = atand(eva(3)/eva(1))

%Plotting the three axis
plot3(evlinea(:,1),evlinea(:,2),evlinea(:,3),'black');
plot3(evlineb(:,1),evlineb(:,2),evlineb(:,3),'blue');
plot3(evlinec(:,1),evlinec(:,2),evlinec(:,3),'red');
