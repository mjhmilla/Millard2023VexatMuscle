function benchRecord = createEmptySimulationRecord(...
                         benchConfig,...
                         numberOfSimulations)

npts = benchConfig.npts;


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

    