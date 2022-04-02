%This code is designed to take the desired semi-major axis, the desired
%inclination, and the target latitude, and display the various observation
%times that various eccentricities would yield
clc;clf;clear

%Inputs
    a=0.48;
    Inclination = 80;
    TargetLat=65;
    LightnessFactor=0;%In case you want to simulate the effect of pointing the sail at the sun. if not, just set to 0
    plotOptions=[1,1];
    ConeVerticies=100;
    ShowCone=1;
    DisplayTable=0;
    ecc=0:0.1:0.9;
    

%Givens & Conversions
    rSun=0.00465047; %AU
    uSun=(1-LightnessFactor)*3.9640*10^-14; %AU^3/s^2
    Inclination = Inclination*pi/180; %Convert to radians
    TargetLat=TargetLat*pi/180;%Convert to radians

%Calculate Keplerian Parameters
    Period=2*pi*sqrt(a^3/uSun)/24/3600;
    P=a.*(1-ecc.^2);
    Theta=(0:0.01:2*pi)'; %vector of thetas
    r=(P./(1+ecc.*cos(Theta)))'; %Matrix of r values, where each row is a different eccentricity
    Theta=Theta'; %Turns Theta for multiplication reasons
    x=r.*cos(Theta); %2-D X positions of each orbit
    y=r.*sin(Theta); %2-D Y positions of each orbit
    myIndices=find(Theta>=pi); %index at highest radius
    myrP=r(:,1); %Perigee Values
    myrAP=r(:,myIndices(1)); %Apogee Values
    
%Rotate 2-D plane about Y axis
    A=cos(Inclination); %First element of rotation matrix
    B=sin(Inclination); %Other element of rotation matrix
    my3DX=x.*A;
    my3DY=y;
    my3DZ=-x.*B;

%Find points within the cone
    % I noticed this is exactly the same for all orbits if they are in the
    % same general orbit, but if you dont believe me uncomment 68 and 82
    
    OrbitNum=1; BinaryMatrix=0*my3DX; %Binary blank matrix, shows where intersection points are
    for OrbitNum = 1:length(ecc)
        ThisX=my3DX(OrbitNum,:); %this orbits X values
        ThisY=my3DY(OrbitNum,:); %This orbits Y values
        ThisZ=my3DZ(OrbitNum,:); %This orbits Z values
        ThisAxisR=sqrt(ThisX.^2 + ThisY.^2); %Radi of the points from the z-axis
        myRatio=tan(TargetLat); %ratio for the cone: z/x
        for myPoint=1:length(ThisX) %loops through every point in the orbit
           thisConeR=ThisZ(myPoint)/myRatio; %Radius of the cone from the current zaxis
           if (thisConeR>ThisAxisR(myPoint)) %if the point is within the cone
               BinaryMatrix(OrbitNum,myPoint) = 1;
           end
        end
    end
    
if ~sum(sum(BinaryMatrix))
    error("Does not intersect the observation cone")
else
    %Isolate Vectors and rotate them back to 2-D
        BinaryVec=BinaryMatrix(1,:); %due to geometry, all of these vectors are the same, so only need first
        VecIndicies=find(BinaryVec==1); %indices where the craft is within the observation period
        ObservationVecs=[my3DX(:,VecIndicies(1)) my3DY(:,VecIndicies(1)) my3DZ(:,VecIndicies(1))]; %3-D Vectors where observation can begin for each orbit
        ToFlatRotationVec=inv([A 0 B; 0 1 0; -B 0 A]); %Rotation vector to convert 3-D back to 2-D
        myNewXY=ToFlatRotationVec*ObservationVecs'; %2-D vectors of starting observation point
        newX=myNewXY(1,:);
        newY=myNewXY(2,:);

    %Identify Observation Time
        newTheta=pi+atan(newY./newX); %theta at starting observation point
        TargetBigE=2*atan( sqrt((1-ecc)./(1+ecc)) .*tan(newTheta/2)); %Eccentric anomaly
        TargetMe=TargetBigE-ecc.*sin(TargetBigE); %Mean Anomaly
        tDays=(Period - 2*( abs( TargetMe*Period/(2*pi) ) ) ); %Observation time in days

    %plot options, from here on out its irrelevant
        if plotOptions(1)
            figure (1) %Plot 2-D orbits
            plot(x(1,:),y(1,:))    
            for i=2:length(ecc)
                hold on;plot(x(i,:),y(i,:))
            end
            hold on;
            plot(newX,newY,'x')
            grid on;
        end

        if plotOptions(2)
            figure (sum(plotOptions(1:2))) %Plot 3-D orbits
            plot3(my3DX(1,:),my3DY(1,:),my3DZ(1,:))    
            for i=2:length(ecc)
                hold on; plot3(my3DX(i,:),my3DY(i,:),my3DZ(i,:))
            end
            %Plot Sun
            hold on; [SunX SunY SunZ]=sphere(10); SunX=SunX.*rSun; SunY=SunY.*rSun; SunZ=SunZ.*rSun; 
            surf(SunX, SunY, SunZ); 

            if ShowCone 
                %Plot Target Cone:
                rConeVec = linspace(0,max(myrAP)/2,ConeVerticies) ;
                th = linspace(0,2*pi,ConeVerticies);
                [RCone,T] = meshgrid(rConeVec,th) ;
                XCone = RCone.*cos(T) ;
                YCone = RCone.*sin(T) ;
                ZCone = sqrt(XCone.^2 + YCone.^2)/cos(TargetLat);
                hold on; surf(XCone,YCone,ZCone)
            end

            daspect([1 1 1]);
            pbaspect([1 1 1]); grid on;
            xlabel("X"); ylabel("Y"); zlabel("Z")
        end
        if DisplayTable
            fprintf("At a semi major axis of %.3f AU and an inclination of %.1f degrees, an orbit will take %.0f days\n",a,Inclination*180/pi,Period)
            eccentricityValues=ecc'; ObservationPeriod=tDays'; ObservationPercentage=(tDays/Period*100)'; Perigee=myrP; Apogee = myrAP;
            disp(table(eccentricityValues, ObservationPeriod, ObservationPercentage, Perigee, Apogee))
        end
end