      
    //START Clusterizer parameters info
      
        char clusteringParametersTag [200] = "ceneknew";
        //char clusteringParametersTag [200] = "narrow";
      
        // Whether clusterizer should save images
        // of individual segments  
        bool clusteringSaveImages = true; //true;
        //bool clusteringSaveImages = true;
        
        // Whether clusterizer should run in 'calibration mode'
        // in which raw signal from detector is saved. Used 
        // for determining calibration function and its parameters
        bool clusteringCalibrationRun = false;
        
        // Whether should clusterizer provide debug messages
        bool clusteringDebugMode = false; //true;
        
        char librariesPrefix[200] = "/home/apmnair18/alice/FOCAL/";
        
        // File with geometry information - path
        char detectorGeometryFile[200] = "/home/apmnair18/Documents/MSc_Project/Simulation/GeometryFiles/geometry.txt"; // directory of geometry.txt
        
        // File with parameters information - path
        char clusterizerParametersFile[200] = "parameters.txt"; // directory of parameters.txt
        
        // Which segments to use as 'seeds'
        // for creating clusters from sub-clusters
        Int_t  COARSE_SEEDS[3] = {2,4,0};
        Int_t  FINE_SEEDS[2] = {1,3};
        Int_t  ALL_SEEDS [5] = {2,4,0,1,3};
        
        // Which segments paramteres should be used for calibration
        Int_t COARSE_CALIB = 0;
        Int_t FINE_CALIB = 1;
        
        // Energy weights for individual segments
        Float_t  ENERGY_WEIGHTS[6] = {1,1,1,1,1,1};
            
        // Default parameters for clusterizer
        // LKegacy leftover; should be overridden by 'parameters.txt'
        Float_t distance = 4.;
        Float_t eThreshold = 100000;

        // Noise settings
        Float_t pixelNoiseProb = 0;
        Float_t padNoiseSigma = 0;

        Float_t logWeight = 4.5; // w_0 for log weights
        Float_t powWeight = 2.2; // p for power weights
        
    //END Clusterizer parameters info
