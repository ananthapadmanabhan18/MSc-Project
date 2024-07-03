#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TF2.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <AliRunLoader.h>
#include <AliFOCAL.h>
#include <AliFOCALLoader.h>
#include <AliFOCALGeometry.h>
#include <AliFOCALDigitizer.h>
#include <AliFOCALdigit.h>
#include <AliRun.h>
#include <AliStack.h>
#include <THStack.h>
#include <TParticle.h>


using namespace std;


void mci(Int_t startFolder, Int_t endFolder, Int_t startEvent, Int_t endEvent)
{
	// Parameteres for getting weighting function info
	Float_t infoParameterSigma = 1.0;
	Float_t infoParameterFitRadius = 5;
	Float_t infoParameterHistRadius = 5;

	char dataSampleTag[200] = "";
	char simFolder[200] ="/home/apmnair18/Documents/MSc-Project/Simulation_outputs/pdg111_nevents199";
	char infoOutputFileDir[200] = "../Out";

	const Float_t profileRadius = 5; // [cm] distance for collecting digits for shower profile

	Long_t totalEvents = 0;

	AliFOCALGeometry *geometry = new AliFOCALGeometry();
	const Char_t *detectorGeometryFile = gSystem->ExpandPathName("/home/apmnair18/Documents/MSc-Project/Simulation/GeometryFiles/geometry.txt");
	geometry->Init(detectorGeometryFile);
	Int_t nSegments = geometry->GetVirtualNSegments();


    std::cout<<std::endl;
	std::cout<< "#################################"<<std::endl;
	std::cout<< "#####  Number of Segments="<<nSegments << " #####"<<std::endl;
	std::cout<< "#################################"<<std::endl;   
    std::cout<<std::endl;

	Float_t z_FOCAL_front = geometry->GetFOCALZ0() - geometry->GetFOCALSizeZ() / 2 - 0.1;
	// Create digitizer
	AliFOCALDigitizer *digitizer = new AliFOCALDigitizer();
	digitizer->SetGeometry(geometry);
	digitizer->Initialize("HS");


	// Event info
	Char_t folder[200];
	Int_t event;

	// Incident particle info
	Int_t pdgCode[4];
	Float_t e[4];
	Float_t pt[4];
	Float_t phi[4];
	Float_t theta[4];
	Float_t Vx[4], Vy[4], Vz[4];
	Float_t cVx[4], cVy[4], cVz[4];
	Int_t conv_flag[4];

	// Fit info
	Float_t *totalE = new Float_t[2 * nSegments];
	Float_t *maxE = new Float_t[2 * nSegments];
	Float_t *sigma = new Float_t[2 * nSegments];
	Float_t *ampl = new Float_t[2 * nSegments];
	Float_t *chi2 = new Float_t[2 * nSegments];
	Int_t *ndef = new Int_t[2 * nSegments];

	// profile histos
	TH2F **showProf2DEvent = new TH2F *[nSegments];
	TH1F **showProf1DEvent = new TH1F *[nSegments];

	for(Int_t i=0; i<nSegments; i++)
	{

		Int_t nBins =5;
		Float_t cellsize =1;
		if(i==1 || i==3)
		{
			nBins =50;
			cellsize =0.1;
		}
		Float_t lim = 5 +0.5*cellsize; // limit of the histogram


		std::cout<<"The number of bins for bin number "<<i<<" is "<<nBins<<std::endl;
		std::cout<<"The cell size for bin number "<<i<<" is "<<cellsize<<std::endl;
		showProf2DEvent[i] = new TH2F(Form("prof2DEvent_%d", i), Form("Shower profile seg %d;x (cm);y (cm)", i), 2*nBins+1 , -lim , lim, 2*nBins+1 , -lim , lim);
		showProf1DEvent[i] = new TH1F(Form("prof1DEvent_%d", i),Form("Shower profile seg %d;r (cm);dN/rdr", i),nBins+1,-0.5*cellsize,lim);

	}



    std::cout<<std::endl;
	std::cout<< "#######################################"<<std::endl;
	std::cout<< "#### Starting to loop over folders ####"<<std::endl;
	std::cout<< "#######################################"<<std::endl;   
    std::cout<<std::endl;




	for(Int_t nfolder = startFolder;nfolder<=endFolder;nfolder++)
	{
		std::cout<<"####################"<<std::endl;
		std::cout<<"# Folder number: "<<nfolder<<" #"<<std::endl;
		std::cout<<"####################"<<std::endl;



		// Open the galice.root file
		char filename[200];
		sprintf(filename,"%s/%i/%s",simFolder,nfolder,"galice.root");
		printf("Opening file %s\n",filename);
		Long_t id,size,flags,mt;
		if (gSystem->GetPathInfo(filename, &id, &size, &flags, &mt) == 1)
		{
			cout << "Error opening galice.root in folder number: " << nfolder << endl;
			continue;
		}


		// Open the FOCAL.Hits.root file
		char filenameHits[200];
		sprintf(filenameHits, "%s/%i/%s", simFolder, nfolder, "FOCAL.Hits.root");
		printf("Opening file %s\n", filenameHits);
		if (gSystem->GetPathInfo(filenameHits, &id, &size, &flags, &mt) == 1)
		{
			cout << "Error opening FOCAL.HIts.root in folder number:  " << nfolder << endl;
			continue;
		}


		// Open the Kinematics.root file
		char filenameKin[200];
		sprintf(filenameKin, "%s/%i/%s", simFolder, nfolder, "Kinematics.root");
		printf("Opening file %s\n", filenameKin);
		if (gSystem->GetPathInfo(filenameKin, &id, &size, &flags, &mt) == 1)
		{
			cout << "Error opening Kinematics.root in folder number:  " << nfolder << endl;
			continue;
		}

		//Alice run loader
		AliRunLoader *fRunLoader = AliRunLoader::Open(filename);

		if(!fRunLoader)
		{
			cout << "Error opening galice.root in folder number: " << nfolder << endl;
			continue;
		}
		cout << "Opened galice.root in folder number: " << nfolder << endl;


		if (!fRunLoader->GetAliRun()){fRunLoader->LoadgAlice();}
		if (!fRunLoader->TreeE()){fRunLoader->LoadHeader();}
		if (!fRunLoader->TreeK()){fRunLoader->LoadKinematics();}

		gAlice = fRunLoader->GetAliRun();

		AliFOCAL *fFOCAL = (AliFOCAL *)gAlice->GetDetector("FOCAL");
		AliFOCALLoader *fFOCALLoader = dynamic_cast<AliFOCALLoader *>(fRunLoader->GetLoader("FOCALLoader"));
		fFOCALLoader->LoadHits("READ");

		// creating the output folder
		cout << "Open output file " << endl;
		gSystem->mkdir(infoOutputFileDir, true);


		TFile *f = new TFile(Form("%s/MCInfo_%s_%i.root", infoOutputFileDir, dataSampleTag, nfolder), "RECREATE");
		f->cd();


		// Create the tree MCInfo and the branches
		TTree *tInfo = new TTree("MCInfo", "MonteCarlo and other information");
		tInfo->Branch("Folder", &folder, "Folder[200]/C");
		tInfo->Branch("Event", &event, "Event/I");
		tInfo->Branch("PdgCode", pdgCode, "PdgCode[4]/I");
		tInfo->Branch("Energy", e, "Energy[4]/F");
		tInfo->Branch("Pt", pt, "Pt[4]/F");
		tInfo->Branch("Phi", phi, "Phi[4]/F");
		tInfo->Branch("Theta", theta, "Theta[4]/F");
		tInfo->Branch("Vx", Vx, "Vx[4]/F");
		tInfo->Branch("Vy", Vy, "Vy[4]/F");
		tInfo->Branch("Vz", Vz, "Vz[4]/F");
		tInfo->Branch("conv_flag", conv_flag, "conv_flag[4]/I");
		tInfo->Branch("cVx", cVx, "cVx[4]/F");
		tInfo->Branch("cVy", cVy, "cVy[4]/F");
		tInfo->Branch("cVz", cVz, "cVz[4]/F");
		tInfo->Branch("TotalEnergy", totalE, "TotalEnergy[40]/F");
		tInfo->Branch("MaxEnergy", maxE, "MaxEnergy[40]/F");
		tInfo->Branch("Sigma", sigma, "Sigma[40]/F");
		tInfo->Branch("Amplitude", ampl, "Amplitude[40]/F");
		tInfo->Branch("Chi2", chi2, "Chi2[40]/F");
		tInfo->Branch("NDEF", ndef, "NDEF[40]/I");



		// Create the profile histograms
		TProfile2D **showProf2D = new TProfile2D *[nSegments];
		TProfile **showProf1D = new TProfile *[nSegments];
		TProfile **showProf1D_dEdr = new TProfile *[nSegments];
		for(Int_t i=0;i<nSegments;i++)
		{
			Int_t nBins =5;
			Float_t cellSize =1;
			if(i==1 || i==3)
			{
				nBins =50;
				cellSize =0.1;
			}
			Float_t lim = 5 +0.5*cellSize; // limit of the histogram
			showProf2D[i] = new TProfile2D(Form("prof2D_%d", i), Form("Shower profile seg %d;x (cm);y (cm)", i), 2 * nBins + 1, -lim, lim, 2 * nBins + 1, -lim, lim);
			showProf1D[i] = new TProfile(Form("prof_%d", i), Form("Shower profile seg %d;r (cm);dN/rdr", i), nBins + 1, -0.5 * cellSize, lim);
			showProf1D_dEdr[i] = new TProfile(Form("prof_dEdr_%d", i), Form("Shower profile seg %d;r (cm);dN/dr", i), nBins + 1, -0.5 * cellSize, lim);
		}



		// Loop over events
		for(Int_t ievt = startEvent; ievt<=endEvent;ievt++)
		{
			totalEvents++;
			std::cout<<"######################"<<std::endl;
			std::cout<<"#  Event number: "<<ievt<<" #"<<std::endl;
			std::cout<<"######################"<<std::endl;

			strcpy(folder, Form("%s/%i", simFolder, nfolder));
			event = ievt;

			// Load the event
			Int_t ie = fRunLoader->GetEvent(ievt); //if ie=0 then event is found
			if (ie != 0){cout << "Error reading event " << ievt << endl;continue;}
			
			// Get the event info
			TTree *treeH = fFOCALLoader->TreeH();
			if (!treeH){std::cout << "TreeH is corrupt\n";break;continue;}

			//Print the number of tracks in the event
			cout << "Event: " << ievt << " with " << treeH->GetEntries() << " tracks" << endl;

			// Get the incident particle info
			TTree *treeK = fRunLoader->TreeK();

			//Initialize the variables to zero
			for(Int_t i=0;i<4;i++)
			{
				pdgCode[i] = 0;
				e[i] = 0;
				pt[i] = 0;
				phi[i] = 0;
				theta[i] = 0;
				Vx[i] = 0;
				Vy[i] = 0;
				Vz[i] = 0;
				conv_flag[i] = 0;
				cVx[i] = 0;
				cVy[i] = 0;
				cVz[i] = 0;
			}


			// Get the Stack for particle info
			AliStack *stack = fRunLoader->Stack();
			if (stack->GetNprimary() != 1){cout << "More than one primary found; this will lead to unexpected results... " << endl;}


			// the 0th index is for the primary particle
			auto primary = stack->Particle(0);
			pdgCode[0] = primary->GetPdgCode();
			e[0] = primary->Energy();
			pt[0] = primary->Pt();
			phi[0] = primary->Phi();
			theta[0] = primary->Theta();
			Vx[0] = primary->Vx();
			Vy[0] = primary->Vy();
			Vz[0] = primary->Vz();

			if(primary->GetFirstDaughter()<=0 && primary->GetLastDaughter()<=0){cout << "No daughter tracks... " << endl;}
			else
			{
				Int_t iTrk = 1;



				for(Int_t i = primary->GetFirstDaughter(); i <= primary->GetLastDaughter(); i++)
				{

					TParticle* part = stack->Particle(i);
					if(iTrk<4)
					{
						pdgCode[iTrk] = part->GetPdgCode();
						e[iTrk] = part->Energy();
						pt[iTrk] = part->Pt();
						phi[iTrk] = part->Phi();
						theta[iTrk] = part->Theta();
						Vx[iTrk] = part->Vx();
						Vy[iTrk] = part->Vy();
						Vz[iTrk] = part->Vz();						
					}
					else{cout << "WARNING: too many daughters" << endl;}

					if(part->GetPdgCode() == 22 && part->GetFirstDaughter()>=0)
					{

					}

					iTrk++;

				}//loop over all the daughters

			}// end of if statement that the primary has any daughters



		}//end of loop over events


		f->cd();
		f->Write();
		f->Close();
	}//end of loop over folders


	delete digitizer;
	delete[] maxE;
	delete[] sigma;
	delete[] ampl;
	delete[] chi2;
	delete[] ndef;
}    
