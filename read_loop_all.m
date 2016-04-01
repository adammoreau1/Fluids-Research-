%read_loop_all: This version will compare all cases of Richardson Numbers
%
%**************************************************************************
%********************* Setup of Data & Display ****************************
%**************************************************************************

%Clearing old Code
clear;
clc;
clf;
cla;
close all;
figure;
figure;
%large loops
for loop = 1:7
    
    %Reading in the Data
    if loop == 1
        uvelfileID = fopen('corr2pt3dq0010.cbin','rb'); wm = 0;
        %uvelfileID = fopen('corr2pt3dwm0010.cbin','rb'); wm = 1;
        ri = 0.1;
    end
    if loop == 2
        uvelfileID = fopen('corr2pt3dq0020.cbin','rb');
        %uvelfileID = fopen('corr2pt3dwm0020.cbin','rb');
        ri = 0.2;
    end
    if loop == 3
        uvelfileID = fopen('corr2pt3dq0050.cbin','rb');
        %uvelfileID = fopen('corr2pt3dwm0050.cbin','rb');
        ri = 0.5;
    end
    if loop == 4
        uvelfileID = fopen('corr2pt3dq0100.cbin','rb');
        %uvelfileID = fopen('corr2pt3dwm0100.cbin','rb');
        ri = 1;
    end
    if loop == 5
        uvelfileID = fopen('corr2pt3dq0200.cbin','rb');
        %uvelfileID = fopen('corr2pt3dwm0200.cbin','rb');
        ri = 2;
    end
    if loop == 6
        uvelfileID = fopen('corr2pt3dq0500.cbin','rb');
        %uvelfileID = fopen('corr2pt3dwm0500.cbin','rb');
        ri = 5;
    end
    if loop == 7
        uvelfileID = fopen('corr2pt3dq1000.cbin','rb');
        %uvelfileID = fopen('corr2pt3dwm1000.cbin','rb');
        ri = 10;
    end
    uvel = fread(uvelfileID,128*128*128,'double');
    fclose(uvelfileID);
    uve = reshape(uvel, [128 128 128]);
    
    %//////////////////////////////////////////////////////
    %/////// From here down i & j should be correct ///////
    %//////////////////////////////////////////////////////
    
    %Looping the searchvalues
    %Range: 0.01-0.89
    loopnum = 500;
    ang = [1:loopnum];
    angz = [1:loopnum];
    svs = [1:loopnum];
    cnt = [1:loopnum];
    centerx = [1:loopnum];
    centery = [1:loopnum];
    centerz = [1:loopnum];
    radiix = [1:loopnum];
    radiiy = [1:loopnum];
    radiiz = [1:loopnum];
    evx = [1:loopnum];
    evy = [1:loopnum];
    evz = [1:loopnum];
    c=0;
    for m = 1:loopnum
        searchvalue = m/733.82-0.00863727;
        c = c +1
        svs(m) = searchvalue;
        %Setting the Alpha Value [0-1] (Transperency)
        alphas = 0.4;
        
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
        cnt(m) = count;
        %**************************************************************************
        %*********************Data Display and analysis****************************
        %**************************************************************************
        
        %print the total number of elements in the matrix
        %display(count);
        
        %Formating Data for placement into ellipsoid_fit
        out = [X; Y; Z];
        out = out.';
        xpi = pi/64 * X - pi;
        ypi = pi/64 * Y - pi;
        zpi = pi/64 * Z - pi;
        out2 = [xpi; ypi; zpi];
        out3 = out2.';
        
        %Display a scatter plot of Points
        % scatter3( xpi, ypi, zpi, '.');
        
        %ellipsoid_fit
        [ center, radii, evecs, v ] = ellipsoid_fit( out3 );
        
        %**************************************************************************
        %***************Eigenvector Interpretation and Display*********************
        %**************************************************************************
        
        %Rotating the Center to be horizantal
        %Each of these values come from the output of ellipsoid fit, and thus will
        %be slightly greater than 0
        center1 = [center(1),center(2),center(3)];
        
        centerx(m) = center(1);
        centery(m) = center(2);
        centerz(m) = center(3);
        radiix(m) = radii(1);
        radiiy(m) = radii(2);
        radiiz(m) = radii(3);
        evx(m) = evecs(1);
        evy(m) = evecs(2);
        evz(m) = evecs(3);
        
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
        ang(m) = atand(eva(2)/eva(1));
        angz(m) = atand(eva(3)/eva(1));
    end;
    fprintf('END');
    figure(1);scatter(svs,ang,10,'filled');
    hold on;
    %colorbar;
    title('Angle of y/x')
    xlabel('Search Value')
    ylabel('Angle in Degrees')
    figure(2);scatter(svs,angz,10,'filled');
    hold on;
    title('Angle of z/x')
    %colorbar;
    xlabel('Search Value')
    ylabel('Angle in Degrees')
    
    if wm == 1
        fileID = fopen(sprintf('wm_%2.0f_all.txt',loopnum),'a');
    end
    if wm == 0
        fileID = fopen(sprintf('q_%2.0f_all.txt',loopnum),'a');
    end
    
    fprintf(fileID,'%4s %10s %11s %11s %7s %10s %10s %10s %9s %9s %9s %16s %16s %16s \n','Ri','Isovalue','Y/X Angle','Z/X Angle','Count','Center X','Center Y','Center Z','Radii a','Radii b','Radii c','Eigenvectors A','Eigenvectors B','Eigenvectors C');
    for p=1:loopnum
        fprintf(fileID,'%4.2f %10.4f %11.4f %11.4f %7.0f %10.4f %10.4f %10.4f %9.4f %9.4f %9.4f %16.4f %16.4f %16.4f \n',[ri svs(p) ang(p) angz(p) cnt(p) centerx(p) centery(p) centerz(p) radiix(p) radiiy(p) radiiz(p) evx(p) evy(p) evz(p)]);
    end
    fclose(fileID);
end
