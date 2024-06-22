        
        // START HITS LOADER 
        
        char filename[200];
        sprintf(filename,"%s/%i/%s",simFolder,nfolder,"galice.root");
                
        //check galice exists (Does not work for ROOT6?)
        /*
        Long_t * id, *size, *flags, *mt;
        if (gSystem->GetPathInfo(filename,id,size,flags,mt) == 1) {
            cout << "ERROR: FOLDER: " << nfolder << endl;
            continue;        
        }
        */
     
        //Alice run loader
        if (fRunLoader) {
          fRunLoader->Delete();
        }
        fRunLoader = AliRunLoader::Open(filename);
        
        if (!fRunLoader) {
            cout << "ERROR: FOLDER: " << nfolder << endl;
            continue;   
        }
        
        if (!fRunLoader->GetAliRun()) fRunLoader->LoadgAlice();
        if (!fRunLoader->TreeE()) fRunLoader->LoadHeader();
        if (!fRunLoader->TreeK()) fRunLoader->LoadKinematics();

        gAlice = fRunLoader->GetAliRun();
        
        //Focal init
        AliFOCAL *fFOCAL  = (AliFOCAL*)gAlice->GetDetector("FOCAL");
        AliFOCALLoader *fFOCALLoader = dynamic_cast<AliFOCALLoader*>(fRunLoader->GetLoader("FOCALLoader"));
        fFOCALLoader->LoadHits("READ");
        
        // Loop over events in the folder
        Int_t maxEvent = TMath::Min(endEvent, fRunLoader->GetNumberOfEvents());
        for (Int_t ievt = startEvent; ievt <= maxEvent; ievt++) {
             
            /***************************************************************************
             *  LOAD EVENT
             ***************************************************************************/
            
            TDirectory * actualDir = fFOCALCluster->SetOutputForEvent(Form("Event%i",ievt));
            
            Int_t ie = fRunLoader->GetEvent(ievt);
            if (ie)
              continue;
            TTree* treeH = fFOCALLoader->TreeH();
            cout << "\tEvent: " << ievt << " with " << treeH->GetEntries() << " tracks" << endl;
            
            fFOCALCluster->SetInput(treeH, "HITS");
            fFOCALCluster->LoadHits2Cells(); // Creates fCells list
            
            actualDir->cd();
            
            TObjString * sPar;
            TParameter<Int_t> * iPar;
            
            sPar = new TObjString(Form("%s/%i",simFolder,nfolder));
            sPar->Write("SimFolder");
            
            sPar = new TObjString(simFolder);
            sPar->Write("GeneralSimFolder");
            
            iPar = new TParameter<Int_t>("FolderNumber",nfolder);
            iPar->Write("FolderNumber");
      
            iPar = new TParameter<Int_t>("EventNumber",ievt);
            iPar->Write("EventNumber");
            
            // END HITS LOADER
