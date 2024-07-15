
    //START Data sample info
        
        // Will be included in the output file name
        // char dataSampleTag[200] = "config_pi0_p1-400";
        //char dataSampleTag[200] = "config_pi0_p50";
        char dataSampleTag[200] = "";
        
        // Whether we're analysing embedded events - works only with DigitLoader
        // Will change how event info is stored
        Bool_t embedded = false;
        
        // Tag of the digit files, this needs to be set according to the input files
        // The generic inputfile name format is: digits_[DIGITSTAG]_[NUMBER].root
        char digitsTag[200] = "";
        
        // Folder with simulated data. 
        // For hits it is expected that folders numbered 
        //   [startFolder-endFolder] will be found inside.
        // For digits it is expected that digit files will be found inside 
        //   (see generic name format above)
        // char simFolder[200] = "/home/focal/storage/sims/box/config9/pi0/config_pi0_p1-400";
        //char simFolder[200] = "/home/focal/storage/sims/box/config11/pi0/config_pi0_p50";
        char simFolder[200] = "/home/apmnair18/Documents/MSc-Project/Simulation_outputs/pdg111_nevents199";
        
        // Where output files should be stored
        char clusteringOutputFileDir[200] = "/home/apmnair18/Documents/MSc-Project/test_outputs";
        // char clusteringOutputFileDir[200] = "/home/apmnair18/Documents/MSc-Project/Analysis_outputs/pdg111_nevents199";
        
        // Where images should be stored
        char clusteringImagesDir[200] = "/home/apmnair18/Documents/MSc-Project/test_outputs";

//        // DONE NOW AS RUNCLUSTERIZER ARGUMENTS
//        // Which folders (or files in case of digits) and events should be analyzed
//        Int_t startFolder = 1;
//        Int_t endFolder = 2;
//        Int_t startEvent = 0;
//        Int_t endEvent = 1;
    
    // END Data sample information
  
