clc
clear
close all
try
    load('vascular network 50-50.mat'); %This must be in the same directory as the Main function
    main(nx,ny,inlet,outlet,z);
catch
    disp("WARNING: The 'vascular network 50-50.mat' must be in the same directory as the Main function")
end    


%The CAPS are 2D matrices
%The lowercase are 1D matrices
function main(nx,ny,kin,kout,Map)

ORIGINAL_MAP = Map; %Backup
d = 15e-3; %[mm] diameter of the vecells
mi = 3e-5; %[mmHg-s] viscosity
l = 1.12; %[mm] length of vessel between two ones

nk = (nx - 1) * nx + ny;
G_value=0; A=0; PRES=0; p=0; G=0; q=0; SP=0; vf=0; V_F=0;

%Case A
type = "Case A";
ZERO();
PRESSURES(FILLTHEG());
SPEEDS();
PRINT();

%Case B
type = "Case B";
ZERO();
NEWMAP();
PRESSURES(FILLTHEG());
SPEEDS();
PRINT();

%Case C
type = "Case C";
mi = 50e-5; %[mmHg-s] viscosity
Map = ORIGINAL_MAP;
ZERO();
PRESSURES(FILLTHEG());
SPEEDS();
PRINT();


%Correct
function ZERO()
    G_value = (pi * d^4) / (128 * mi * l);
    A = (pi * d^2) / 4; %[mm^2] CrossSectionArea
    PRES = zeros(nx,ny); %Pressure matrix
    p = zeros(nk);
    G = zeros(nk,nk); %G matrix
    q = zeros(nk);
    SP = zeros(nx,ny); %Speed of flow
    vf = zeros(nk); %Volume_Flow
    V_F = zeros(nx,ny);%Volume_Flow
end

%Correct
function NEWMAP()
    %Removes a chunck of the vessels from the midle
    for i=22:40
        for j=18:32
            Map(i,j) = 0;
        end
    end
end

%Correct
function Ge = FILLTHEG()
    %Checking the matrix left bounderies
    for i=2:(nx - 1)
        j=1;
        k = (i - 1) * nx + j;
            if Map(i,j) == 1 %The Point in Question
                if Map(i,j+1) == 1 %Right Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k+1) =  - 1;
                end
                if Map(i-1,j) == 1 %Top Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k-nx) = - 1;
                end
                if Map(i+1,j) == 1 %Bottom Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k+nx) =  - 1;
                end
            end
    end

    %Checking the matrix right bounderies
    for i=2:(nx - 1)
        j=ny;
        k = (i - 1) * nx + j;
            if Map(i,j) == 1 %The Point in Question
                if Map(i,j-1) == 1 %Left Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k-1) =  - 1;
                end
                if Map(i-1,j) == 1 %Top Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k-nx) = - 1;
                end
                if Map(i+1,j) == 1 %Bottom Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k+nx) =  - 1;
                end
            end
    end

    %Checking the matrix top bounderies
    for j=2:(ny - 1)
        i=1;
        k = (i - 1) * nx + j;
            if Map(i,j) == 1 %The Point in Question
                if Map(i,j+1) == 1 %Right Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k+1) =  - 1;
                end
                if Map(i,j-1) == 1 %Left Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k-1) =  - 1;
                end
                if Map(i+1,j) == 1 %Bottom Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k+nx) = - 1;
                end
            end
    end

    %Checking the matrix bottom bounderies
    for j=2:(ny - 1)
        i=nx;
        k = (i - 1) * nx + j;
            if Map(i,j) == 1 %The Point in Question
                if Map(i,j+1) == 1 %Right Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k+1) =  - 1;
                end
                if Map(i,j-1) == 1 %Left Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k-1) = - 1;
                end
                if Map(i-1,j) == 1 %Top Point
                    G(k,k) = G(k,k) + 1;
                    G(k,k-nx) =  - 1;
                end
            end
    end

    %Checking the Top Left Point Only
    i = 1; j = 1;
    k = (i - 1) * nx + j;
    if Map(i,j) == 1 %The Point in Question
        if Map(i,j+1) == 1 %Right Point
            G(k,k) = G(k,k) + 1;
            G(k,k+1) =  - 1;
        end
        if Map(i+1,j) == 1 %Bottom Point
            G(k,k) = G(k,k) + 1;
            G(k,k+nx) = - 1;
        end
    end

    %Checking the Top Right Point Only
    i = 1;
    j = ny;
    k = (i - 1) * nx + j;
    if Map(i,j) == 1 %The Point in Question
        if Map(i,j-1) == 1 %Left Point
            G(k,k) = G(k,k) + 1;
            G(k,k-1) = - 1;
        end
        if Map(i+1,j) == 1 %Bottom Point
            G(k,k) = G(k,k) + 1;
            G(k,k+nx) = - 1;
        end
    end

    %Checking the Bottom Left Point Only
    i = nx; 
    j = 1;
    k = (i - 1) * nx + j;
    if Map(i,j) == 1 %The Point in Question
        if Map(i,j+1) == 1 %Right Point
            G(k,k) = G(k,k) + 1;
            G(k,k+1) = - 1;
        end
        if Map(i-1,j) == 1 %Top Point
            G(k,k) = G(k,k) + 1;
            G(k,k-nx) = - 1;
        end
    end

    %Checking the Bottom Right Point Only
    i = nx; 
    j = ny;
    k = (i - 1) * nx + j;
    if Map(i,j) == 1 %The Point in Question
        if Map(i,j-1) == 1 %Left Point
            G(k,k) = G(k,k) + 1;
            G(k,k-1) =  - 1;
        end
        if Map(i-1,j) == 1 %Top Point
            G(k,k) = G(k,k) + 1;
            G(k,k-nx) = - 1;
        end
    end
    
    %Checking the midle of the matrix skipping the bounderies
    for i=2:(nx - 1)
        for j=2:(ny - 1)
            k = (i - 1) * nx + j;
                if Map(i,j) == 1 %The Point in Question
                    if Map(i,j+1) == 1 %Right Point
                        G(k,k) = G(k,k) + 1;
                        G(k,k+1) = - 1;
                    end
                    if Map(i,j-1) == 1 %Left Point
                        G(k,k) = G(k,k) + 1;
                        G(k,k-1) = - 1;
                   end
                    if Map(i-1,j) == 1 %Top Point
                        G(k,k) = G(k,k) + 1;
                        G(k,k-nx) =  - 1;
                    end
                    if Map(i+1,j) == 1 %Bottom Point
                        G(k,k) = G(k,k) + 1;
                        G(k,k+nx) =  - 1;
                   end
                end
        end
    end

    Ge = G;
end

%Correct
function PRESSURES(Ge)
    a = Ge;
    Ge = sparse(G_value*Ge); 
    Ge(kout,:) = 0;
    Ge(kout,kout) = 1;
    Ge(kin,:) = 0;
    Ge(kin,kin) = 1;
    
    q(kin) = 25; %[mmHg]
    q(kout) = 5; %[mmHg]    

    Ge = sparse(Ge);
    p = Ge\q;
    p(isnan(p))=0;

    
    %Make the 1D vectors into 2x2 matrix
    for i=1:nx
        for j=1:ny
           k = (i - 1) * nx + j;
           PRES(i,j) = p(k);          
        end
    end
    
end

%Correct___Same with the FILLTHEG find the dP, the flow, and then the speed
function SPEEDS()
    %Checking the matrix left bounderies
    for i=2:(nx - 1)
        j=1;
        k = (i - 1) * nx + j;
            if PRES(i,j) > 0 %The Point in Question
                if PRES(i,j+1) > 0 %Right Point
                    dP = PRES(i,j) - PRES(i,j+1); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i-1,j) > 0 %Top Point
                    dP = PRES(i,j) - PRES(i-1,j); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i+1,j) > 0 %Bottom Point
                    dP = PRES(i,j) - PRES(i+1,j); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
            else
                vf(k) = 0;
                sp(k) = 0;
            end
    end

    %Checking the matrix right bounderies
    for i=2:(nx - 1)
        j=ny;
        k = (i - 1) * nx + j;
            if PRES(i,j) > 0 %The Point in Question
                if PRES(i,j-1) > 0 %Left Point
                    dP = PRES(i,j) - PRES(i,j-1); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i-1,j) > 0 %Top Point
                    dP = PRES(i,j) - PRES(i-1,j); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i+1,j) > 0 %Bottom Point
                    dP = PRES(i,j) - PRES(i+1,j); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
            else
                vf(k) = 0;
                sp(k) = 0;
            end
    end

    %Checking the matrix top bounderies
    for j=2:(ny - 1)
        i=1;
        k = (i - 1) * nx + j;
            if PRES(i,j) > 0 %The Point in Question
                if PRES(i,j+1) > 0 %Right Point
                    dP = PRES(i,j) - PRES(i,j+1); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i,j-1) > 0 %Left Point
                    dP = PRES(i,j) - PRES(i,j-1); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i+1,j) > 0 %Bottom Point
                    dP = PRES(i,j) - PRES(i+1,j); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
            else
                vf(k) = 0;
                sp(k) = 0;
            end
    end

    %Checking the matrix bottom bounderies
    for j=2:(ny - 1)
        i=nx;
        k = (i - 1) * nx + j;
            if PRES(i,j) > 0 %The Point in Question
                if PRES(i,j+1) > 0 %Right Point
                    dP = PRES(i,j) - PRES(i,j+1); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i,j-1) > 0 %Left Point
                    dP = PRES(i,j) - PRES(i,j-1); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
                if PRES(i-1,j) > 0 %Top Point
                    dP = PRES(i,j) - PRES(i-1,j); 
                    vf(k) = vf(k) + G_value * abs(dP);
                    sp(k) = vf(k) / A;
                end
            else
                vf(k) = 0;
                sp(k) = 0;
            end
    end

    %Checking the Top Left Point Only
    i = 1; j = 1;
    k = (i - 1) * nx + j;
    if PRES(i,j) > 0 %The Point in Question
        if PRES(i,j+1) > 0 %Right Point
            dP = PRES(i,j) - PRES(i,j+1); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
        if PRES(i+1,j) > 0 %Bottom Point
            dP = PRES(i,j) - PRES(i+1,j); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
    else
        vf(k) = 0;
        sp(k) = 0;
    end

    %Checking the Top Right Point Only
    i = 1;
    j = ny;
    k = (i - 1) * nx + j;
    if PRES(i,j) > 0 %The Point in Question
        if PRES(i,j-1) > 0 %Left Point
            dP = PRES(i,j) - PRES(i,j-1); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
        if PRES(i+1,j) > 0 %Bottom Point
            dP = PRES(i,j) - PRES(i+1,j); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
    else
        vf(k) = 0;
        sp(k) = 0;
    end

    %Checking the Bottom Left Point Only
    i = nx; 
    j = 1;
    k = (i - 1) * nx + j;
    if PRES(i,j) > 0 %The Point in Question
        if PRES(i,j+1) > 0 %Right Point
            dP = PRES(i,j) - PRES(i,j+1); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
        if PRES(i-1,j) > 0 %Top Point
            dP = PRES(i,j) - PRES(i-1,j); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
    else
        vf(k) = 0;
        sp(k) = 0;
    end

    %Checking the Bottom Right Point Only
    i = nx; 
    j = ny;
    k = (i - 1) * nx + j;
    if PRES(i,j) > 0 %The Point in Question
        if PRES(i,j-1) > 0 %Left Point
            dP = PRES(i,j) - PRES(i,j-1); 
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
        if PRES(i-1,j) > 0 %Top Point
            dP = PRES(i,j) - PRES(i-1,j);
            vf(k) = vf(k) + G_value * abs(dP);
            sp(k) = vf(k) / A;
        end
    else
        vf(k) = 0;
        sp(k) = 0;
    end
    
    %Checking the midle of the matrix skipping the bounderies
    for i=2:(nx - 1)
        for j=2:(ny - 1)
            k = (i - 1) * nx + j;
                if PRES(i,j) > 0 %The Point in Question
                    if PRES(i,j+1) > 0 %Right Point
                        dP = PRES(i,j) - PRES(i,j+1); 
                        vf(k) = vf(k) + G_value * abs(dP);
                        sp(k) = vf(k) / A;
                    end
                    if PRES(i,j-1) > 0 %Left Point
                        dP = PRES(i,j) - PRES(i,j-1); 
                        vf(k) = vf(k) + G_value * abs(dP);
                        sp(k) = vf(k) / A;
                   end
                    if PRES(i-1,j) > 0 %Top Point
                        dP = PRES(i,j) - PRES(i-1,j); 
                        vf(k) = vf(k) + G_value * abs(dP);
                        sp(k) = vf(k) / A;
                    end
                    if PRES(i+1,j) > 0 %Bottom Point
                        dP = PRES(i,j) - PRES(i+1,j);
                        vf(k) = vf(k) + G_value * abs(dP);
                        sp(k) = vf(k) / A;
                    end
                else
                    vf(k) = 0;
                    sp(k) = 0;
                end
        end
    end   
   
    %Make the 1D vectors into 2x2 matrix
    for i=1:nx
        for j=1:ny
           k = (i - 1) * nx + j;
           V_F(i,j) = vf(k); 
           SP(i,j) = sp(k);
        end
    end
end

%Correct
function PRINT()
    figure("Name", "G Matrix " + type)
    spy(G)
    title("G Matrix " + type);
    
    figure("Name", "Flow Map " + type)
    subplot(2,2,1)
    surf(0:nx-1,0:ny-1,V_F,'FaceColor','interp','EdgeColor','r','FaceLighting','phong')
    title("FLow Map " + type);
    xlabel("i");
    ylabel("j");
    zlabel("Volumetric Flow [mm3/s]");
    
    %figure("Name", "Speed Map " + type)
    subplot(2,2,2)
    surf(0:nx-1,0:ny-1,SP,'FaceColor','interp','EdgeColor','r','FaceLighting','phong')
    title("Speed Map " + type);
    xlabel("i");
    ylabel("j");
    zlabel("Speed [mm/s]");
    
    %figure("Name","Vessel Map " + type)
    subplot(2,2,3)
    surf(0:nx-1,0:ny-1,Map,'FaceColor','interp','EdgeColor','r','FaceLighting','phong')
    title("Vessel Map " + type);
    xlabel("i");
    ylabel("j");
    zlabel("Vessel");

    %figure("Name", "Pressure Map " + type)
    subplot(2,2,4)
    surf(0:nx-1,0:ny-1,PRES,'FaceColor','interp','EdgeColor','r','FaceLighting','phong')
    title("Pressure Map " + type);
    xlabel("i");
    ylabel("j");
    zlabel("Pressure [mmHg]");
end

end
