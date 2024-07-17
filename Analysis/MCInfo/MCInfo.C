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

void MCInfo(Int_t startFolder, Int_t endFolder, Int_t startEvent, Int_t endEvent)
{

	// Parameteres for getting weighting function info
	Float_t infoParameterSigma = 1.0;
	Float_t infoParameterFitRadius = 5;
	Float_t infoParameterHistRadius = 5;

	char dataSampleTag[200] = "";
	char simFolder[200] ="..//..//Simulation_outputs/pdg111_nevents199";
	char infoOutputFileDir[200] = "../../test_outputs/";

	const Float_t profileRadius = 5; // [cm] distance for collecting digits for shower profile

	Long_t totalEvents = 0;

	// Get geometry;
	AliFOCALGeometry *geometry = new AliFOCALGeometry();
	const Char_t *detectorGeometryFile = gSystem->ExpandPathName("..//..//Simulation/GeometryFiles/geometry.txt");
	geometry->Init(detectorGeometryFile);
	Int_t nSegments = geometry->GetVirtualNSegments();


	std::cout<< "#################################"<<std::endl;
	std::cout<< "Number of Segments"<<nSegments << std::endl;
	std::cout<< "#################################"<<std::endl;


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




	// Profile histos
	TH2F **showProf2DEvent = new TH2F *[nSegments];
	TH1F **showProf1DEvent = new TH1F *[nSegments];
	for (Int_t i = 0; i < nSegments; i++)
	{
		Int_t nBins = 5;
		Float_t cellSize = 1;
		if (i == 1 || i == 3)
		{
			nBins = 50;
			cellSize = 0.1;
		}
		Float_t lim = 5 + 0.5 * cellSize;
		showProf2DEvent[i] = new TH2F(Form("prof2DEvent_%d", i), Form("Shower profile seg %d;x (cm);y (cm)", i), 2 * nBins + 1, -lim, lim, 2 * nBins + 1, -lim, lim);
		showProf1DEvent[i] = new TH1F(Form("profEvent_%d", i), Form("Shower profile seg %d;r (cm);dN/rdr", i), nBins + 1, -0.5 * cellSize, lim);
	}





	// Loop over folders
	for (Int_t nfolder = startFolder; nfolder <= endFolder; nfolder++)
	{



		cout<<"FOLDER: "<<folder << endl;

		Long_t id, size, flags, mt;



		char filename[200];
		sprintf(filename, "%s/%i/%s", simFolder, nfolder, "galice.root");
		if (gSystem->GetPathInfo(filename, &id, &size, &flags, &mt) == 1){cout << "Error in galice.root in folder : " << nfolder << endl;continue;}




		char filenameHits[200];
		sprintf(filenameHits, "%s/%i/%s", simFolder, nfolder, "FOCAL.Hits.root");
		if (gSystem->GetPathInfo(filenameHits, &id, &size, &flags, &mt) == 1){cout << "Error in FOCAL.Hits.root in folder : " << nfolder << endl;continue;}



		char filenameKin[200];
		sprintf(filenameKin, "%s/%i/%s", simFolder, nfolder, "Kinematics.root");
		if (gSystem->GetPathInfo(filenameKin, &id, &size, &flags, &mt) == 1){cout << "Error in Kinematics.root in folder: " << nfolder << endl;continue;}




		// Alice run loader
		AliRunLoader *fRunLoader = AliRunLoader::Open(filename);

		if (!fRunLoader)
		{
			cout << "ERROR: FOLDER: " << nfolder << endl;
			continue;
		}
		cout << "Got runloader" << endl;

		if (!fRunLoader->GetAliRun()){fRunLoader->LoadgAlice();}
		if (!fRunLoader->TreeE()){fRunLoader->LoadHeader();}
		if (!fRunLoader->TreeK()){fRunLoader->LoadKinematics();}

		gAlice = fRunLoader->GetAliRun();

		// Focal init
		AliFOCAL *fFOCAL = (AliFOCAL *)gAlice->GetDetector("FOCAL");
		AliFOCALLoader *fFOCALLoader = dynamic_cast<AliFOCALLoader *>(fRunLoader->GetLoader("FOCALLoader"));
		fFOCALLoader->LoadHits("READ");

		// Creating the output File
		cout << "Open output file " << endl;
		gSystem->mkdir(infoOutputFileDir, true);




		/***************************************
		Creating the output TTree for the MCInfo
		****************************************/

		TFile *f = new TFile(Form("%s/MCInfo_%s_%i.root", infoOutputFileDir, dataSampleTag, nfolder), "RECREATE");
		f->cd();
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






		/******************************
		Creating the Profile Histograms
		*******************************/
		TProfile2D **showProf2D = new TProfile2D *[nSegments];
		TProfile **showProf1D = new TProfile *[nSegments];
		TProfile **showProf1D_dEdr = new TProfile *[nSegments];
		for (Int_t i = 0; i < nSegments; i++)
		{
			Int_t nBins = 5;
			Float_t cellSize = 1;
			if (i == 1 || i == 3)
			{
				nBins = 50;
				cellSize = 0.1;
			}
			Float_t lim = 5 + 0.5 * cellSize;
			showProf2D[i] = new TProfile2D(Form("prof2D_%d", i), Form("Shower profile seg %d;x (cm);y (cm)", i), 2 * nBins + 1, -lim, lim, 2 * nBins + 1, -lim, lim);
			showProf1D[i] = new TProfile(Form("prof_%d", i), Form("Shower profile seg %d;r (cm);dN/rdr", i), nBins + 1, -0.5 * cellSize, lim);
			showProf1D_dEdr[i] = new TProfile(Form("prof_dEdr_%d", i), Form("Shower profile seg %d;r (cm);dN/dr", i), nBins + 1, -0.5 * cellSize, lim);
		}






		/*******************************************************************
		Looping Over the Events and Filling the TTree and Profile Histograms
		********************************************************************/
		for (Int_t ievt = startEvent; ievt <= endEvent; ievt++)
		{

			totalEvents++;

			strcpy(folder, Form("%s/%i", simFolder, nfolder));
			event = ievt;

			Int_t ie = fRunLoader->GetEvent(ievt); // Get the event at index ievt
			if (ie != 0){cout << "Error reading event " << ievt << endl;continue;}
			

			TTree *treeH = fFOCALLoader->TreeH();
			if (!treeH){std::cout << "TreeH is corrupt\n";break;continue;}
			cout << "Event: " << ievt << " with " << treeH->GetEntries() << " tracks" << endl;


			TTree *treeK = fRunLoader->TreeK();
			if(!treeK){std::cout << "TreeK is corrupt\n";break;continue;}



			for (Int_t i = 0; i < 4; i++)	// Particle info
			{
				pdgCode[i] = 0; 			// PDG code
				e[i] = 0;					// Energy
				pt[i] = 0;					// Transverse momentum
				phi[i] = 0; 				// Azimuthal angle
				theta[i] = 0; 				// Polar angle
				Vx[i] = 0; 					// Vertex x
				Vy[i] = 0;					// Vertex y
				Vz[i] = 0; 					// Vertex z
				conv_flag[i] = 0;			// Conversion flag
				cVx[i] = 0;					// vertex if photon x
				cVy[i] = 0;					// vertex if photon y
				cVz[i] = 0;					// vertex if photon z			
			}



			AliStack *stack = fRunLoader->Stack();



			if (stack->GetNprimary() != 1){cout << "More than one primary found; this will lead to unexpected results... " << endl;}


			auto primary = stack->Particle(0);    // Primary particle
			pdgCode[0] = primary->GetPdgCode();
			e[0] = primary->Energy();
			pt[0] = primary->Pt();
			phi[0] = primary->Phi();
			theta[0] = primary->Theta();
			Vx[0] = primary->Vx();
			Vy[0] = primary->Vy();
			Vz[0] = primary->Vz();



			if (primary->GetFirstDaughter() <= 0 && primary->GetLastDaughter() <= 0)
			{
				cout << "No daughter tracks... " << endl;
			}
			else
			{
				Int_t iTrk = 1;
				for (Int_t i = primary->GetFirstDaughter(); i <= primary->GetLastDaughter(); i++)
				{
					TParticle *part = stack->Particle(i);
					if (iTrk < 4)
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
					else
					{
						cout << "WARNING: too many daughters" << endl;
					}

					if (part->GetPdgCode() == 22 && part->GetFirstDaughter() >= 0) // if the particle is a photon
					{
						auto cdaughter = stack->Particle(part->GetFirstDaughter());
						if (cdaughter->Vz() < z_FOCAL_front)
						{
							conv_flag[iTrk] = 1;
						}
						cVx[iTrk] = cdaughter->Vx();
						cVy[iTrk] = cdaughter->Vy();
						cVz[iTrk] = cdaughter->Vz();
					}
					iTrk++;

				}// end of loop over daughters

			}// End of if statement for daughters




			// Call the digitizer //
			digitizer->Hits2Digits(treeH->GetBranch("FOCAL"));

			
			// Get the digits//
			TClonesArray *digitsArray = digitizer->GetSDigits();


			// Creating the histograms //
			TH2F **histograms = new TH2F *[nSegments];
			TH1D **histograms1D = new TH1D *[nSegments];


			// Getting the Row, column and segment info //
			Int_t nCol, nRow;
			Float_t sizeX = geometry->GetFOCALSizeX();
			Float_t sizeY = geometry->GetFOCALSizeY();



			char name[20] = "Original";
			char name2[20] = "Projection";



			for (Int_t i = 0; i < nSegments; i++)
			{
				geometry->GetVirtualNColRow(i, nCol, nRow);
				Float_t x0, y0, z0, r;
				z0 = geometry->GetVirtualSegmentZ(i);
				r = TMath::Tan(theta[0]) * z0;
				x0 = TMath::Cos(phi[0]) * r;
				y0 = TMath::Sin(phi[0]) * r;
				histograms[i] = new TH2F(Form("%s_%i_%i_%i", name, nfolder, ievt, i), Form("%s_%i_%i_%i", name, nfolder, ievt, i), nCol * 2 * infoParameterHistRadius / sizeX, x0 - infoParameterHistRadius, x0 + infoParameterHistRadius, nRow * 2 * infoParameterHistRadius / sizeY, y0 - infoParameterHistRadius, y0 + infoParameterHistRadius);
				histograms[i]->SetStats(0);
			}

			cout << "Filling histograms" << endl;

			Int_t *rowSeed = new Int_t[nSegments];
			Int_t *colSeed = new Int_t[nSegments];
			Float_t *ampSeed = new Float_t[nSegments];
			for (Int_t iSeg = 0; iSeg < nSegments; iSeg++)
			{
				rowSeed[iSeg] = -1;
				colSeed[iSeg] = -1;
				ampSeed[iSeg] = -1;
			}


			TObjArray digitsForProfile;

			for (Int_t i = 0; i < digitsArray->GetEntries(); i++)
			{
				AliFOCALdigit *digit = (AliFOCALdigit *)digitsArray->UncheckedAt(i);
				Int_t col = digit->GetCol();
				Int_t row = digit->GetRow();
				Float_t x, y, z;
				Int_t segment = digit->GetSegment();

				Float_t x0, y0, z0, r;
				z0 = geometry->GetVirtualSegmentZ(segment);
				r = TMath::Tan(theta[0]) * z0; // Phi = 0 && x = 1 && y = 0;
				x0 = TMath::Cos(phi[0]) * r;
				y0 = TMath::Sin(phi[0]) * r;

				geometry->GetXYZFromColRowSeg(col, row, segment, x, y, z);
				Int_t energy = digit->GetAmp();
				histograms[segment]->Fill(x, y, (Float_t)energy);


				// prepare for profiles
				if (TMath::Abs(x - x0) < profileRadius && TMath::Abs(y - y0) < profileRadius)
				{
					digitsForProfile.AddLast(digit);
					if (energy > ampSeed[segment])
					{
						rowSeed[segment] = row;
						colSeed[segment] = col;
						ampSeed[segment] = energy;
					}
				}
			}// End of loop over digits




			// Fill profiles
			//
			// First fill per-event profiles,then average. This matter for the exmpy bins
			//  (alternative would be to use TProfile, but fill empty cells as a zero in the prof)

			cout << "Reset per event hists" << endl;
			for (Int_t i = 0; i < nSegments; i++)
			{
				showProf2DEvent[i]->Reset();
				showProf1DEvent[i]->Reset();
			}
			cout << "Fill digits" << endl;
			for (Int_t iDigit = 0; iDigit < digitsForProfile.GetEntriesFast(); iDigit++)
			{
				AliFOCALdigit *digit = (AliFOCALdigit *)digitsForProfile.UncheckedAt(iDigit);
				Int_t col = digit->GetCol();
				Int_t row = digit->GetRow();
				Float_t x, y, z;
				Int_t segment = digit->GetSegment();
				geometry->GetXYZFromColRowSeg(col, row, segment, x, y, z);
				Float_t x0, y0, z0;
				geometry->GetXYZFromColRowSeg(colSeed[segment], rowSeed[segment], segment, x0, y0, z0);
				Int_t energy = digit->GetAmp();

				showProf2DEvent[segment]->Fill(x - x0, y - y0, energy / ampSeed[segment]);
				Float_t r = TMath::Sqrt(pow(x - x0, 2) + pow(y - y0, 2));
				showProf1DEvent[segment]->Fill(r, energy / ampSeed[segment]);
			}

			cout << "Add hists" << endl;
			for (Int_t iSeg = 0; iSeg < nSegments; iSeg++)
			{
				Int_t xbins = showProf2DEvent[iSeg]->GetXaxis()->GetNbins();
				Int_t ybins = showProf2DEvent[iSeg]->GetYaxis()->GetNbins();
				// could use this to count # bins per r bin
				for (Int_t ibin = 0; ibin < xbins; ibin++)
				{
					for (Int_t jbin = 0; jbin < ybins; jbin++)
					{
						Float_t x = showProf2DEvent[iSeg]->GetXaxis()->GetBinCenter(ibin + 1);
						Float_t y = showProf2DEvent[iSeg]->GetYaxis()->GetBinCenter(jbin + 1);
						showProf2D[iSeg]->Fill(x, y, showProf2DEvent[iSeg]->GetBinContent(ibin + 1, jbin + 1));
					}
				}
				xbins = showProf1DEvent[iSeg]->GetNbinsX();
				for (Int_t ibin = 0; ibin < xbins; ibin++)
				{
					Float_t r = showProf1DEvent[iSeg]->GetXaxis()->GetBinCenter(ibin + 1);
					// Would be better to use actual (discretised) area instead of 1/r
					if (r != 0)
						showProf1D[iSeg]->Fill(r, showProf1DEvent[iSeg]->GetBinContent(ibin + 1) / r);
					showProf1D_dEdr[iSeg]->Fill(r, showProf1DEvent[iSeg]->GetBinContent(ibin + 1));
				}
			}

			cout << "Fitting histograms" << endl;

			for (Int_t i = 0; i < nSegments; i++)
			{
				// Compute initial parameters
				Float_t x0, y0, z0, r;
				z0 = geometry->GetVirtualSegmentZ(i);
				r = TMath::Tan(theta[0]) * z0; // Phi = 0 && x = 1 && y = 0;
				x0 = TMath::Cos(phi[0]) * r;
				y0 = TMath::Sin(phi[0]) * r;
				Float_t amp = histograms[i]->GetMaximum();
				// Define function
				TF2 *fc = new TF2("cauchy", "[3]/(1+((x-[0])*(x-[0]) + (y-[1])*(y-[1]))/([2]*[2]))", x0 - infoParameterFitRadius, x0 + infoParameterFitRadius, y0 - infoParameterFitRadius, y0 + infoParameterFitRadius);
				fc->SetParameters(x0, y0, infoParameterSigma, 1);
				fc->SetParNames("x0", "y0", "#sigma", "I");

				geometry->GetVirtualNColRow(i, nCol, nRow);
				Int_t yBin = nRow * infoParameterHistRadius / sizeY / 2;
				histograms1D[i] = histograms[i]->ProjectionX(Form("%s_%i", name2, i),yBin, yBin);
				histograms1D[i]->SetStats(0);
				TF1 *fc1D = new TF1("cauchy1D", "[2]/(1+TMath::Power((x-[0]),2)/([1]*[1]))", x0 - infoParameterFitRadius, x0 + infoParameterFitRadius);
				fc1D->SetParameters(x0, infoParameterSigma, 1);
				fc1D->SetParNames("x0", "#sigma", "I");

				cout << Form("Initial parameters: %.2e, %.2e, %.2e, %i", x0, y0, infoParameterSigma, 1) << endl;

				// Fit the histogram
				if (pdgCode[0] == 22)
				{
					histograms[i]->Fit("cauchy", "RN", "goff");
					histograms1D[i]->Fit("cauchy1D", "R", "goff");
				}

				// Get Results
				Double_t pars[4];
				fc->GetParameters(pars);
				totalE[i * 2] = histograms[i]->Integral();
				maxE[i * 2] = amp;
				sigma[i * 2] = TMath::Abs(pars[2]);
				ampl[i * 2] = pars[3];
				chi2[i * 2] = fc->GetChisquare();
				ndef[i * 2] = fc->GetNDF();
				cout << Form("Final parameters: %.2e, %.2e, %.2e, %.2e; Chi/n = %.2e/%i",
							 pars[0], pars[1], pars[2], pars[3], fc->GetChisquare(), fc->GetNDF())
					 << endl;
				fc1D->GetParameters(pars);
				totalE[i * 2 + 1] = 0;
				maxE[i * 2 + 1] = histograms1D[i]->GetMaximum();
				sigma[i * 2 + 1] = TMath::Abs(pars[1]);
				ampl[i * 2 + 1] = pars[2];
				chi2[i * 2 + 1] = fc->GetChisquare();
				ndef[i * 2 + 1] = fc->GetNDF();
				cout << Form("Final parameters1D: %.2e, %.2e, %.2e; Chi/n = %.2e/%i",
							 pars[0], pars[1], pars[2], fc1D->GetChisquare(), fc1D->GetNDF())
					 << endl;

				cout << "Saving..." << endl;

				fc->Delete();
				fc1D->Delete();
				histograms[i]->Delete();
				histograms1D[i]->Delete();
			}


			tInfo->Fill(); // Fill the tree

			delete[] rowSeed;
			delete[] colSeed;
			delete[] ampSeed;

			delete[] histograms;
			delete[] histograms1D;

		} // End of loop over Events

		f->cd();
		f->Write();
		f->Close();

		delete[] showProf2D;
		delete[] showProf1D;
		delete[] showProf1D_dEdr;
		fRunLoader->Delete();

	} // End of loop over folders

	delete digitizer;
	delete[] maxE;
	delete[] sigma;
	delete[] ampl;
	delete[] chi2;
	delete[] ndef;
}
