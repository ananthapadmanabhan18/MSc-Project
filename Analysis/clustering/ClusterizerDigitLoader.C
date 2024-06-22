        
        // START HITS LOADER 
        
        char filename[200];
        sprintf(filename,"%s/digits_%s_%i.root",simFolder,digitsTag,nfolder);
                
        //check galice exists
        Long_t id, size, flags, mt;
        if (gSystem->GetPathInfo(filename,&id,&size,&flags,&mt) == 1) {
            cout << "ERROR: Digit folder file does not exist: " << nfolder << endl;
            continue;        
        }
        
        inputFile = new TFile(filename,"READ");
        
        // Loop over events in the folder
        for (Int_t ievt = startEvent; ievt <= endEvent; ievt++) {
             
            /***************************************************************************
             *  LOAD EVENT
             ***************************************************************************/
            
            TDirectory * actualDir = fFOCALCluster->SetOutputForEvent(Form("Event%i",ievt));   

            // Should not be necessary            
//            if (treeD) {
//              delete treeD;
//              treeD = 0;
//            }
            
            // Necessary to keep memory clean (unloads all objects of inputFile from memory)
            inputFile->Delete("T*");
            
            inputFile->GetDirectory(Form("Event%i",ievt))->GetObject("DigitTree",treeD);
            
            cout << "\tEvent: " << ievt << endl;
            
            fFOCALCluster->SetInput(treeD, "DIGITS");
            fFOCALCluster->LoadDigits2Cell(); // Creates fCells list
            
            actualDir->cd();
            
            TObjString * sPar;
            TParameter<Int_t> * iPar;
            
            inputFile->GetDirectory(Form("Event%i",ievt))->GetObject("SignalFolder",sPar);
            sPar->Write("SimFolder");
                 
            inputFile->GetDirectory(Form("Event%i",ievt))->GetObject("SignalEvent",iPar);
            iPar->Write("EventNumber");
            
            if (embedded) {
              inputFile->GetDirectory(Form("Event%i",ievt))->GetObject("BackgroundFolder",sPar);
              sPar->Write("SimFolderBackground");
                   
              inputFile->GetDirectory(Form("Event%i",ievt))->GetObject("BackgroundEvent",iPar);
              iPar->Write("EventNumberBackground");
            }            
            // END HITS LOADER
