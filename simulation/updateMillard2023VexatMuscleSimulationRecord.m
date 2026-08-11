function benchRecord = updateMillard2023VexatMuscleSimulationRecord(...
                         mtInfo,...
                         benchConfig,...
                         benchRecord,...
                         timeIndex,...
                         simulationIndex,...
                         numberOfSimulations)

npts = benchConfig.npts;

if(isempty(benchRecord)==1)
  benchRecord.time                        = zeros(npts,numberOfSimulations);
  benchRecord.activation                  = zeros(npts,numberOfSimulations);
  benchRecord.activationDot               = zeros(npts,numberOfSimulations);  
  benchRecord.cpuTime                     = zeros(numberOfSimulations);
  benchRecord.normFiberForceAlongTendon   = zeros(npts,numberOfSimulations);
  benchRecord.normFiberForce              = zeros(npts,numberOfSimulations);
  benchRecord.activeFiberForceAlongTendon = zeros(npts,numberOfSimulations);
  benchRecord.activeFiberForce            = zeros(npts,numberOfSimulations);   
  benchRecord.normActiveFiberForce        = zeros(npts,numberOfSimulations); 
  benchRecord.normPassiveFiberForce        = zeros(npts,numberOfSimulations);   
  benchRecord.fiberForceAlongTendon       = zeros(npts,numberOfSimulations);
  benchRecord.tendonForce                 = zeros(npts,numberOfSimulations);
  benchRecord.fiberForce                  = zeros(npts,numberOfSimulations);  
  benchRecord.normFiberLength             = zeros(npts,numberOfSimulations);    
  benchRecord.pennationAngle              = zeros(npts,numberOfSimulations);
  benchRecord.normFiberVelocity           = zeros(npts,numberOfSimulations);
  benchRecord.pennationAngularVelocity    = zeros(npts,numberOfSimulations);
  %benchRecord.fiberStiffnessAlongTendon   = zeros(npts,numberOfSimulations);
  %benchRecord.fiberDampingAlongTendon     = zeros(npts,numberOfSimulations);  
  benchRecord.tendonStiffness             = zeros(npts,numberOfSimulations);
  benchRecord.tendonDamping               = zeros(npts,numberOfSimulations);
  
  benchRecord.crossBridgeStiffness        = zeros(npts,numberOfSimulations);
  benchRecord.crossBridgeDamping          = zeros(npts,numberOfSimulations); 
  benchRecord.normCrossBridgeForce        = zeros(npts,numberOfSimulations); 
  benchRecord.normCrossBridgeForceAlongTendon = zeros(npts,numberOfSimulations); 



  benchRecord.normProximalTitinLength     = zeros(npts,numberOfSimulations);
  benchRecord.normProximalTitinForce      = zeros(npts,numberOfSimulations);
  benchRecord.normDistalTitinLength       = zeros(npts,numberOfSimulations);
  benchRecord.normDistalTitinForce        = zeros(npts,numberOfSimulations);

  benchRecord.pathLength                  = zeros(npts,numberOfSimulations);
  benchRecord.pathVelocity                = zeros(npts,numberOfSimulations);

  benchRecord.dSystemEnergyLessWork       = zeros(npts,numberOfSimulations);
  benchRecord.systemEnergyLessWork        = zeros(npts,numberOfSimulations);
  benchRecord.tendonPotentialEnergy       = zeros(npts,numberOfSimulations);
  benchRecord.fiberPotentialEnergy        = zeros(npts,numberOfSimulations);
  benchRecord.fiberActiveWork             = zeros(npts,numberOfSimulations);
  benchRecord.dampingWork                 = zeros(npts,numberOfSimulations);
  benchRecord.boundaryWork                = zeros(npts,numberOfSimulations);

  benchRecord.tendonPower                 = zeros(npts,numberOfSimulations);
  benchRecord.fiberParallelElementPower   = zeros(npts,numberOfSimulations);
  benchRecord.fiberActivePower            = zeros(npts,numberOfSimulations);
  benchRecord.dampingPower                = zeros(npts,numberOfSimulations);
  benchRecord.boundaryPower               = zeros(npts,numberOfSimulations);

  benchRecord.musculotendonStiffness  = zeros(npts,numberOfSimulations);
  benchRecord.musculotendonDamping    = zeros(npts,numberOfSimulations);
  benchRecord.fiberStiffness  = zeros(npts,numberOfSimulations);
  benchRecord.fiberDamping    = zeros(npts,numberOfSimulations);


  benchRecord.fiberActiveForceLengthMultiplier = zeros(npts,numberOfSimulations);
  benchRecord.fiberPassiveForceLengthMultiplier= zeros(npts,numberOfSimulations);

  benchRecord.eventTime             = zeros(npts,numberOfSimulations);
  benchRecord.eventTendonForce      = zeros(npts,numberOfSimulations);
  benchRecord.eventNormFiberForce   = zeros(npts,numberOfSimulations);
  benchRecord.eventNormActiveFiberForce   = zeros(npts,numberOfSimulations);  
  benchRecord.eventNormFiberForce   = zeros(npts,numberOfSimulations);
  benchRecord.eventNormFiberLength  = zeros(npts,numberOfSimulations);
  benchRecord.eventPathLength       = zeros(npts,numberOfSimulations);  
  benchRecord.eventNormFiberVelocity= zeros(npts,numberOfSimulations);
   

  benchRecord.extra =  zeros(  npts,...
                             numberOfSimulations,...
                             10 );
  benchRecord.extraLabels = {};
  

  benchRecord.state = zeros(  npts,...
                              numberOfSimulations,...
                              benchConfig.numberOfMuscleStates);
  benchRecord.dstate = zeros(  npts,...
                              numberOfSimulations,...
                              benchConfig.numberOfMuscleStates);                            
end

idxSim = simulationIndex;
j = timeIndex;


if(isempty(mtInfo.extra)==0)
    numberOfExtra = length(mtInfo.extra);
    benchRecord.extra(j,idxSim,1:numberOfExtra)=mtInfo.extra';
    benchRecord.extraLabels = mtInfo.extraLabels;
end

                         
benchRecord.state(j,idxSim,:)=mtInfo.state.value';
benchRecord.dstate(j,idxSim,:)=mtInfo.state.derivative';

tendonForce = mtInfo.muscleDynamicsInfo.tendonForce;
       
lceN      = mtInfo.muscleLengthInfo.normFiberLength;
alpha     = mtInfo.muscleLengthInfo.pennationAngle;
dlceN     = mtInfo.fiberVelocityInfo.normFiberVelocity;
dalpha    = mtInfo.fiberVelocityInfo.pennationAngularVelocity;
        
falN      = mtInfo.muscleLengthInfo.fiberActiveForceLengthMultiplier;
fpeN      = mtInfo.muscleDynamicsInfo.normPassiveFiberForce;

faeN      = mtInfo.muscleDynamicsInfo.normActiveFiberForce;
fceNAT    = mtInfo.muscleDynamicsInfo.normFiberForce*cos(alpha);
fceN      = mtInfo.muscleDynamicsInfo.normFiberForce;
fceAT     = mtInfo.muscleDynamicsInfo.fiberForce*cos(alpha);
fce       = mtInfo.muscleDynamicsInfo.fiberForce;    
kx    = mtInfo.muscleDynamicsInfo.crossBridgeStiffness;
dx    = mtInfo.muscleDynamicsInfo.crossBridgeDamping;    
kTendon   = mtInfo.muscleDynamicsInfo.tendonStiffness;
dTendon   = mtInfo.muscleDynamicsInfo.tendonDamping;
fxN = mtInfo.muscleDynamicsInfo.normCrossBridgeForce;
fx = mtInfo.muscleDynamicsInfo.crossBridgeForce;

km = mtInfo.muscleDynamicsInfo.musculotendonStiffness;
dm = mtInfo.muscleDynamicsInfo.musculotendonDamping;
kf = mtInfo.muscleDynamicsInfo.fiberStiffness;
df = mtInfo.muscleDynamicsInfo.fiberDamping;

%%
%force and kinematic information
%%
benchRecord.activation(j,idxSim)                  = mtInfo.muscleDynamicsInfo.activation;
benchRecord.activationDot(j,idxSim)               = mtInfo.muscleDynamicsInfo.activationDerivative;
benchRecord.tendonForce(j,idxSim)                 = tendonForce;
benchRecord.normFiberForceAlongTendon(j,idxSim)   = fceNAT;
benchRecord.normFiberForce(j,idxSim)              = fceN;
benchRecord.normActiveFiberForce(j,idxSim)        = faeN;
benchRecord.fiberForceAlongTendon(j,idxSim)       = fceAT;
benchRecord.fiberForce(j,idxSim)                  = fce;    
benchRecord.tendonForce(j,idxSim)                 = tendonForce;
benchRecord.normFiberLength(j,idxSim)             = lceN;    
benchRecord.pennationAngle(j,idxSim)              = alpha;
benchRecord.normFiberVelocity(j,idxSim)           = dlceN;
benchRecord.pennationAngularVelocity(j,idxSim)    = dalpha;

benchRecord.crossBridgeStiffness(j,idxSim)        = kx;    
benchRecord.crossBridgeDamping(j,idxSim)          = dx;   

if(contains(mtInfo.type,'VEXAT')==1)
    benchRecord.normProximalTitinLength(j,idxSim)     = ...
        mtInfo.muscleLengthInfo.normTitin1Length;

    benchRecord.normProximalTitinForce(j,idxSim)      = ...
        mtInfo.muscleDynamicsInfo.normTitin1Force;

    benchRecord.normDistalTitinLength(j,idxSim)       = ...
        mtInfo.muscleLengthInfo.normTitin2Length;

    benchRecord.normDistalTitinForce(j,idxSim)        = ...
        mtInfo.muscleDynamicsInfo.normTitin2Force;

end
benchRecord.tendonStiffness(j,idxSim)             = kTendon;
benchRecord.tendonDamping(j,idxSim)               = dTendon;

benchRecord.pathLength(j,idxSim)                  = ...
  mtInfo.muscleLengthInfo.pathLength;
benchRecord.pathVelocity(j,idxSim)                = ...
  mtInfo.muscleLengthInfo.pathLengthDerivative;

benchRecord.fiberActiveForceLengthMultiplier(j,idxSim) = falN;
benchRecord.fiberPassiveForceLengthMultiplier(j,idxSim) = fpeN;
benchRecord.normPassiveFiberForce(j,idxSim) = fpeN;


%%
%Energetics & power information
%%
fae = mtInfo.muscleDynamicsInfo.activeFiberForce;
benchRecord.activeFiberForce(j,idxSim)             = fae;
benchRecord.activeFiberForceAlongTendon(j,idxSim)  = fae*cos(alpha);

benchRecord.normCrossBridgeForce(j,idxSim) = fxN;
benchRecord.normCrossBridgeForceAlongTendon(j,idxSim) = fxN*cos(alpha);


benchRecord.musculotendonStiffness(j,idxSim)  = km;
benchRecord.musculotendonDamping(j,idxSim)    = dm;

benchRecord.fiberStiffness(j,idxSim)  = kf;
benchRecord.fiberDamping(j,idxSim)    = df;

benchRecord.tendonPower(j,idxSim)               = mtInfo.muscleDynamicsInfo.tendonPower;
benchRecord.fiberParallelElementPower(j,idxSim) = mtInfo.muscleDynamicsInfo.fiberParallelElementPower;
benchRecord.fiberActivePower(j,idxSim)          = mtInfo.muscleDynamicsInfo.fiberActivePower;
benchRecord.dampingPower(j,idxSim)              = mtInfo.muscleDynamicsInfo.dampingPower;
benchRecord.boundaryPower(j,idxSim)             = mtInfo.muscleDynamicsInfo.boundaryPower;

benchRecord.systemEnergyLessWork(j,idxSim)        = nan;
benchRecord.dSystemEnergyLessWork(j,idxSim)       = nan;
benchRecord.tendonPotentialEnergy(j,idxSim)       = ...
    mtInfo.musclePotentialEnergyInfo.tendonPotentialEnergy;
benchRecord.fiberPotentialEnergy(j,idxSim)        = ...
    mtInfo.musclePotentialEnergyInfo.fiberPotentialEnergy;
benchRecord.fiberActiveWork(j,idxSim)             = nan;
benchRecord.dampingWork(j,idxSim)                 = nan;
benchRecord.boundaryWork(j,idxSim)                = nan;
    