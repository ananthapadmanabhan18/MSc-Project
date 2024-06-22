#include "iostream"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TParameter.h"

#include "AliFOCALGeometry.h"
#include "AliFOCALCluster.h"
#include "AliFOCALRinger.h"
#include "ObjectMap.h"

// void Analysis(int startFolder, int endFolder, int startEvent, int endEvent, const Char_t *MCInfoDir, const Char_t *clustersFolder, const Char_t *dataSampleTag, const Char_t *outputdir)
void Analysis(int startFolder, int endFolder, int startEvent, int endEvent)
{




   const Char_t *MCInfoDir = "/home/apmnair18/Documents/MSc_Project/Analysis_outputs/pdg111_nevents199/";
   const Char_t *clustersFolder= "/home/apmnair18/Documents/MSc_Project/Analysis_outputs/pdg111_nevents199/";
   const Char_t *dataSampleTag="";
   const Char_t *outputdir="/home/apmnair18/Documents/MSc_Project/Analysis_outputs/pdg111_nevents199/";




  // Segment to analyze:
  //    >=0... clusters from specific segment
  //    -1 ... total clusters from whole detector
  //    -2 ... total clusters from coarse segments only
  //    -3 ... total clusters from fine segments only
  int segmentToAnalyze = -1;
  bool analysisDebugMode = true;
  int storeSegmentInfo = 0; // flag to store per-segment info
  float PI0_MASS = 0.1354;   // GeV
  float PI0_SIGMA = 0.01161; // GeV
  float PI0_TOLERANCE = 3;
  float MAX_PAIR_DISTANCE = 5;          // cm
  float MAX_REJECTION_DISTANCE = 0.135; // cm

  AliFOCALGeometry *geometry = new AliFOCALGeometry();
  const Char_t *detectorGeometryFile = gSystem->ExpandPathName("/home/apmnair18/Documents/MSc_Project/Simulation/GeometryFiles/geometry.txt");
  geometry->Init(detectorGeometryFile);

  Char_t folder[200];
  int event;

  int nClusters;

  int pdgCode[4];
  float e[4];
  float pt[4];
  float phi[4];
  float theta[4];
  float Vx[4], Vy[4], Vz[4];
  float cVx[4], cVy[4], cVz[4];
  int conv_flag[4];

  float clusterE[3];
  float clusterX[3];
  float clusterY[3];
  float clusterZ[3];
  float Width1[30];
  float Width2[30];
  float ClusterPhi[30];
  float E_long[30];
  float SeedE[30];
  int Ncells[30];
  float W1Cclust[3];   // overall width 1 from all coarse clusters -- these had been suggested by Marco
  float W2Cclust[3];   // overall width 2 from all coarse clusters
  float W1Fclust[3];   // overall width 1 from all fine clusters
  float W2Fclust[3];   // overall width 2 from all fine clusters
  float VW1Cclust[3];  // overall, vect-added width 1 from all coarse clusters -- using Marco's sigma1-sigma2 interpolation method
  float VW2Cclust[3];  // overall, vect-added width 2 from all coarse clusters
  float VW1Fclust[3];  // overall, vect-added width 1 from all fine clusters
  float VW2Fclust[3];  // overall, vect-added width 2 from all fine clusters
  float PhiRCclust[3]; // angular concentration of shower-proj ellipse in coarse segments
  float PhiRFclust[3]; // angular concentration of shower-proj ellipse in fine segments
  float PhiACclust[3]; // average angular orientation of shower ellipse in coarse segments
  float PhiAFclust[3]; // average angular orientation of shower ellipse in fine segments
  float ERclust[3];    // coarse energy ratio E5/(E2+E4)

  float HCALEnergy[3]; // Energy in HCAL
  float HCALEisoR2[3]; // isolation in HCAL R=0.2
  float HCALEisoR4[3]; // isolation in HCAL R=0.4

  bool allFit;
  float distMatch[3];
  bool primaryFits[3];
  int rejections;
  float invariantMass;

  TFile *inputFile = 0;
  TTree *tClusters = 0;
  TFile *fInfo = 0;



  if (analysisDebugMode)
    std::cout << "Creating cluster maps" << std::endl;



  int infoSegment = segmentToAnalyze;
  if (infoSegment < 0)
  {
    if ((infoSegment == -1) || (infoSegment == -3))
      infoSegment = 1;
    else
      infoSegment = 2;
  }





  if (analysisDebugMode)
    std::cout << Form("Using segment %i for geometry info.", infoSegment) << std::endl;
  int nCols = (int)(geometry->GetFOCALSizeX() / geometry->GetVirtualPadSize(infoSegment));
  int nRows = (int)(geometry->GetFOCALSizeY() / geometry->GetVirtualPadSize(infoSegment));
  ObjectMap *clusterMap = new ObjectMap(nCols, nRows);

  // Loop over folders
  for (int nfolder = startFolder; nfolder <= endFolder; nfolder++)
  {

    std::cout << "FOLDER: " << nfolder << std::endl;

    // START CLUSTERS LOADER

    char filename[200];
    sprintf(filename, "%s/clusters__%i.root", clustersFolder, nfolder);

    Long_t id, size, flags, mt;
    if (gSystem->GetPathInfo(filename, &id, &size, &flags, &mt) == 1)
    {
      std::cout << "ERROR: Cluster file does not exist: " << filename << std::endl;
      continue;
    }

    if (inputFile)
      inputFile->Close();

    inputFile = new TFile(filename, "READ");

    if (fInfo)
    {
      fInfo->Close();
      fInfo = 0;
    }
    fInfo = new TFile(Form("%s/MCInfo__%i.root", MCInfoDir, nfolder), "READ");

    if (!inputFile->IsOpen() || !fInfo->IsOpen())
    {
      std::cout << "Cluster file " << inputFile->GetName() << " or info file " << fInfo->GetName() << " not found" << std::endl;
      continue;
    }

    TTree *tInfo;
    fInfo->GetObject("MCInfo", tInfo);

    tInfo->SetBranchAddress("Folder", folder);
    tInfo->SetBranchAddress("Event", &event);

    tInfo->SetBranchAddress("PdgCode", pdgCode);
    tInfo->SetBranchAddress("Energy", e);
    tInfo->SetBranchAddress("Pt", pt);
    tInfo->SetBranchAddress("Phi", phi);
    tInfo->SetBranchAddress("Theta", theta);
    tInfo->SetBranchAddress("Vx", Vx);
    tInfo->SetBranchAddress("Vy", Vy);
    tInfo->SetBranchAddress("Vz", Vz);
    tInfo->SetBranchAddress("conv_flag", conv_flag);
    tInfo->SetBranchAddress("cVx", cVx);
    tInfo->SetBranchAddress("cVy", cVy);
    tInfo->SetBranchAddress("cVz", cVz);

    TFile *outputFile = new TFile(Form("%s/analysis__%i.root",outputdir, nfolder),"RECREATE");

    outputFile->cd();

    TTree *analysis = new TTree("analysis", "Clusterizer analysis results");

    analysis->Branch("Folder", folder, "Folder[200]/C");
    analysis->Branch("Event", &event, "Event/I");
    analysis->Branch("nClusters", &nClusters, "nClusters/I");

    // Particles, mother + daughters
    analysis->Branch("PdgCode", pdgCode, "PdgCode[4]/I");
    analysis->Branch("Energy", e, "Energy[4]/F");
    analysis->Branch("Pt", pt, "Pt[4]/F");
    analysis->Branch("Phi", phi, "Phi[4]/F");
    analysis->Branch("Theta", theta, "Theta[4]/F");
    analysis->Branch("conv_flag", conv_flag, "conv_flag[4]/I");
    analysis->Branch("Vx", Vx, "Vx[4]/F");
    analysis->Branch("Vy", Vy, "Vy[4]/F");
    analysis->Branch("Vz", Vz, "Vz[4]/F");
    analysis->Branch("cVx", cVx, "cVx[4]/F");
    analysis->Branch("cVy", cVy, "cVy[4]/F");
    analysis->Branch("cVz", cVz, "cVz[4]/F");

    analysis->Branch("ClusterE", clusterE, "ClusterE[3]/F");
    analysis->Branch("ClusterX", clusterX, "ClusterX[3]/F");
    analysis->Branch("ClusterY", clusterY, "ClusterY[3]/F");
    analysis->Branch("ClusterZ", clusterZ, "ClusterZ[3]/F");
    if (storeSegmentInfo)
    {
      analysis->Branch("Width1", Width1, "Width1[30]/F");
      analysis->Branch("Width2", Width2, "Width2[30]/F");
      analysis->Branch("ClusterPhi", ClusterPhi, "ClusterPhi[30]/F");
      analysis->Branch("E_long", E_long, "E_long[30]/F");
      analysis->Branch("SeedE", SeedE, "SeedE[30]/F");
      analysis->Branch("Ncells", Ncells, "Ncells[30]/F");
    }
    analysis->Branch("W1Cclust", W1Cclust, "W1Cclust[3]/F");
    analysis->Branch("W2Cclust", W2Cclust, "W2Cclust[3]/F");
    analysis->Branch("W1Fclust", W1Fclust, "W1Fclust[3]/F");
    analysis->Branch("W2Fclust", W2Fclust, "W2Fclust[3]/F");
    analysis->Branch("VW1Cclust", VW1Cclust, "VW1Cclust[3]/F");
    analysis->Branch("VW2Cclust", VW2Cclust, "VW2Cclust[3]/F");
    analysis->Branch("VW1Fclust", VW1Fclust, "VW1Fclust[3]/F");
    analysis->Branch("VW2Fclust", VW2Fclust, "VW2Fclust[3]/F");
    analysis->Branch("PhiRCclust", PhiRCclust, "PhiRCclust[3]/F");
    analysis->Branch("PhiRFclust", PhiRFclust, "PhiRCclust[3]/F");
    analysis->Branch("PhiACclust", PhiACclust, "PhiACclust[3]/F");
    analysis->Branch("PhiAFclust", PhiAFclust, "PhiACclust[3]/F");
    analysis->Branch("ERclust", ERclust, "ERclust[3]/F");

    analysis->Branch("HCALE", HCALEnergy, "HCALE[3]/F");
    analysis->Branch("HCALEisoR2", HCALEisoR2, "HCALEisoR2[3]/F");
    analysis->Branch("HCALEisoR4", HCALEisoR4, "HCALEisoR4[3]/F");

    analysis->Branch("AllFit", &allFit, "AllFit/O");
    analysis->Branch("DistMatch", distMatch, "DistMatch[3]/F");
    analysis->Branch("PrimaryFits", primaryFits, "PrimaryFits[3]/O");
    analysis->Branch("Rejections", &rejections, "Rejections/I");
    analysis->Branch("InvariantMass", &invariantMass, "InvariantMass/F");
    //
    // Loop over events in the folder
    for (int ievt = startEvent; ievt <= endEvent; ievt++)
    {

      /***************************************************************************
       *  LOAD EVENT
       ***************************************************************************/

      inputFile->GetDirectory(Form("Event%i", ievt))->GetObject("fTreeR", tClusters);

      TBranch *bClusters;
      if (segmentToAnalyze == -1)
        bClusters = tClusters->GetBranch("AliFOCALCluster");
      else
        bClusters = tClusters->GetBranch("AliFOCALClusterItr");

      TClonesArray *clustersArray = 0;
      bClusters->SetAddress(&clustersArray);
      bClusters->GetEvent(0);

      if (analysisDebugMode)
        std::cout << "\tEvent: " << ievt << std::endl;

      TObjString *sPar;
      TParameter<int> *iPar;

      inputFile->GetDirectory(Form("Event%i", ievt))->GetObject("SimFolder", sPar);
      inputFile->GetDirectory(Form("Event%i", ievt))->GetObject("EventNumber", iPar);

      tInfo->GetEntry(ievt);

      if ((strcmp(folder, sPar->GetString().Data()) != 0) || (iPar->GetVal() != event))
      {
        std::cout << Form("Event info doesn't fit for folder:%i, event: %i. Skipping event", nfolder, ievt) << std::endl;
        continue;
      }

      // END CLUSTERS LOADER

      for (int i = 0; i < 3; i++)
      {
        clusterE[i] = -1; // set this to -1 like Davide, perhaps
        clusterX[i] = 0;
        clusterY[i] = 0;
        clusterZ[i] = -1; // set this to -1 like Davide, perhaps
        distMatch[i] = 999;
        primaryFits[i] = false;
      }
      for (int i = 0; i < 30; i++)
      {
        Width1[i] = -1;
        Width2[i] = -1;
        ClusterPhi[i] = -6;
        E_long[i] = -1;
        SeedE[i] = -1;
        Ncells[i] = -1;
      }

      for (int i = 0; i < 3; i++)
      {
        W1Cclust[i] = -1; // suggested by Marco
        W2Cclust[i] = -1; // probably not needed, since calculated below based on Width1, Width2 initialized above
        W1Fclust[i] = -1;
        W2Fclust[i] = -1;
        VW1Cclust[i] = -1;
        VW2Cclust[i] = -1;
        VW1Fclust[i] = -1;
        VW2Fclust[i] = -1;
        PhiRCclust[i] = -1;
        PhiRFclust[i] = -1;
        PhiACclust[i] = -1;
        PhiAFclust[i] = -1;
        ERclust[i] = -1;
        HCALEnergy[i] = 0;
        HCALEisoR2[i] = 0;
        HCALEisoR4[i] = 0;
      }

      rejections = 0;
      invariantMass = 0;

      // Sort them
      float tmpF;
      int tmpI;
      bool modified = true;
      while (modified)
      {
        modified = false;
        for (int i = 3; i > 1; i--)
        {

          if (e[i] > e[i - 1])
          {

            modified = true;

            tmpI = pdgCode[i];
            pdgCode[i] = pdgCode[i - 1];
            pdgCode[i - 1] = tmpI;

            tmpF = e[i];
            e[i] = e[i - 1];
            e[i - 1] = tmpF;

            tmpF = pt[i];
            pt[i] = pt[i - 1];
            pt[i - 1] = tmpF;

            tmpF = phi[i];
            phi[i] = phi[i - 1];
            phi[i - 1] = tmpF;

            tmpF = theta[i];
            theta[i] = theta[i - 1];
            theta[i - 1] = tmpF;

            tmpF = Vx[i];
            Vx[i] = Vx[i - 1];
            Vx[i - 1] = tmpF;

            tmpF = Vy[i];
            Vy[i] = Vy[i - 1];
            Vy[i - 1] = tmpF;

            tmpF = Vz[i];
            Vz[i] = Vz[i - 1];
            Vz[i - 1] = tmpF;

            tmpI = conv_flag[i];
            conv_flag[i] = conv_flag[i - 1];
            conv_flag[i - 1] = tmpI;

            tmpF = cVx[i];
            cVx[i] = cVx[i - 1];
            cVx[i - 1] = tmpF;

            tmpF = cVy[i];
            cVy[i] = cVy[i - 1];
            cVy[i - 1] = tmpF;

            tmpF = cVz[i];
            cVz[i] = cVz[i - 1];
            cVz[i - 1] = tmpF;
          }
        }
      }

      // First create and initializer maps for fast-searching of clusters
      clusterMap->ResetMap();

      AliFOCALRinger ringer;
      ringer.Init(100);

      if (analysisDebugMode)
        std::cout << "Registering clusters into maps" << std::endl;

      // Now register all sub-clusters in the map and reset their flag
      nClusters = clustersArray->GetEntries();
      if (analysisDebugMode)
        std::cout << "Event has " << nClusters << " found clusters." << std::endl;
      for (int i = 0; i < nClusters; i++)
      {

        AliFOCALCluster *cluster = (AliFOCALCluster *)clustersArray->UncheckedAt(i);
        cluster->SetFlag(true);

        int segment = cluster->Segment();

        if (segment != segmentToAnalyze)
          continue;

        int colt = (int)((cluster->X() + geometry->GetFOCALSizeX() / 2) / geometry->GetVirtualPadSize(infoSegment));
        int rowt = (int)((cluster->Y() + geometry->GetFOCALSizeY() / 2) / geometry->GetVirtualPadSize(infoSegment));

        if (!clusterMap->Get(colt, rowt))
        {
          clusterMap->InsertAt(colt, rowt, cluster);
          if (analysisDebugMode)
            std::cout << Form("\tRegistering cluster at XY=[%.2f,%.2f] to ColRow=[%i,%i]", cluster->X(), cluster->Y(), colt, rowt) << std::endl;
        }
        else
        {
          bool found = false;
          ringer.SetRing(1);
          int xCol, yCol;
          while (ringer.GetNext(xCol, yCol))
          {
            if (!clusterMap->Get(colt + xCol, rowt + yCol))
            {
              clusterMap->InsertAt(colt + xCol, rowt + yCol, cluster);
              if (analysisDebugMode)
                std::cout << Form("\tRegistering cluster at XY=[%.2f,%.2f] to ColRow=[%i,%i]", cluster->X(), cluster->Y(), colt + xCol, rowt + yCol) << std::endl;
              found = true;
              break;
            }
          }
          if (!found && analysisDebugMode)
            std::cout << Form("WARNING cluster at [%.2f,%.2f] not registered in the map, no space left in vicinity!",cluster->X(), cluster->Y())<< std::endl;
        }
      }

      // Now, go through clusters near positions of first and second decay photon (resp. primary photon)
      // and see if their clusters were found

      if (analysisDebugMode)
        std::cout << "Looking for clusters matching single/embedded particles" << std::endl;

      bool singlePhotonEvent = (pdgCode[0] == 22);

      float generalZ = geometry->GetVirtualSegmentZ(2); // Use second segment Z; approx shower max (same used for cluster positions)
      if (segmentToAnalyze >= 0)
        generalZ = geometry->GetVirtualSegmentZ(segmentToAnalyze);

      float x[3] = {0, 0, 0};
      float y[3] = {0, 0, 0};

      // Get x,y info of original particles (photons, primary or decay)
      if (singlePhotonEvent)
      {
        x[0] = TMath::Cos(phi[0]) * TMath::Tan(theta[0]) * generalZ;
        y[0] = TMath::Sin(phi[0]) * TMath::Tan(theta[0]) * generalZ;
      }
      else
      {
        x[0] = TMath::Cos(phi[1]) * TMath::Tan(theta[1]) * generalZ;
        y[0] = TMath::Sin(phi[1]) * TMath::Tan(theta[1]) * generalZ;
        x[1] = TMath::Cos(phi[2]) * TMath::Tan(theta[2]) * generalZ;
        y[1] = TMath::Sin(phi[2]) * TMath::Tan(theta[2]) * generalZ;
        if (e[3] > 0)
        {
          x[2] = TMath::Cos(phi[3]) * TMath::Tan(theta[3]) * generalZ;
          y[2] = TMath::Sin(phi[3]) * TMath::Tan(theta[3]) * generalZ;
        }
      }

      // Recalculate x,y to col,row for map use
      int col[3] = {-1, -1, -1};
      int row[3] = {-1, -1, -1};
      int realCol[3] = {-1, -1, -1};
      int realRow[3] = {-1, -1, -1};

      // create holders for primary (resp. decay) particles
      AliFOCALCluster **primary = new AliFOCALCluster *[3];
      primary[0] = 0;
      primary[1] = 0;
      primary[2] = 0;

      int particlesToStudy = 0;
      if (singlePhotonEvent)
        particlesToStudy = 1;
      else if (e[3] == 0)
        particlesToStudy = 2;
      else
        particlesToStudy = 3;

      TLorentzVector lor[3];

      // Locate closest cluster for each particleToStudy
      for (int particle = 0; particle < particlesToStudy; particle++)
      {

        col[particle] = (int)((x[particle] + geometry->GetFOCALSizeX() / 2) / geometry->GetVirtualPadSize(infoSegment));
        row[particle] = (int)((y[particle] + geometry->GetFOCALSizeX() / 2) / geometry->GetVirtualPadSize(infoSegment));

        if (analysisDebugMode)
        {
          std::cout << Form("\tPrimary at XY=[%.2f,%.2f], looking around ColRow[%i,%i]:",x[particle], y[particle], col[particle], row[particle])<< std::endl;
        }

        int xCol, yRow;
        // Search for cluster in the map around col,row position
        for (int ring = 0; ring < 5; ring++)
        {
          ringer.SetRing(ring);
          while (ringer.GetNext(xCol, yRow))
          {
            //                std::cout << Form("\t...at ColRow[%i,%i]",col[particle]+xCol,row[particle]+yRow) << std::endl;
            AliFOCALCluster *c = (AliFOCALCluster *)clusterMap->Get(col[particle] + xCol, row[particle] + yRow);
            // Check that cluster has not been assigned to other particle
            bool isUsed = 0;
            for (int ipart = 0; ipart < particle; ipart++)
            {
              if (primary[ipart] == c)
                isUsed = 1;
            }
            if (isUsed)
              continue;
            if (!primary[particle])
              if (c)
              {
                if (!primary[particle])
                {
                  primary[particle] = c;
                  realCol[particle] = col[particle] + xCol;
                  realRow[particle] = row[particle] + yRow;
                  distMatch[particle] = TMath::Sqrt(TMath::Power(c->X() - x[particle], 2) +
                                                    TMath::Power(c->Y() - y[particle], 2));
                }
                else
                {
                  // Distance of previously found cluster from original particle
                  float oldDistance = distMatch[particle];
                  // Distance of actual found cluster from original particle
                  float newDistance = TMath::Sqrt(TMath::Power(c->X() - x[particle], 2) +
                                                  TMath::Power(c->Y() - y[particle], 2));
                  if (newDistance < oldDistance)
                  {
                    realCol[particle] = col[particle] + xCol;
                    realRow[particle] = row[particle] + yRow;
                    primary[particle] = c;
                    distMatch[particle] = newDistance;
                  }
                }
              }
          } // end while ringer
        }   // end for rings
        if (primary[particle])
        {

          float distance = distMatch[particle];
          if (distance <= MAX_REJECTION_DISTANCE)
            primaryFits[particle] = true;
          else
            primaryFits[particle] = false;

          if (analysisDebugMode)
          {
            std::cout << Form("\t\tFound XY=[%.2f,%.2f] at ColRow=[%i,%i]; distanceFits=%o",
                              primary[particle]->X(), primary[particle]->Y(),
                              realCol[particle], realRow[particle], primaryFits[particle])
                      << std::endl;
          }

          clusterX[particle] = primary[particle]->X();
          clusterY[particle] = primary[particle]->Y();
          clusterZ[particle] = primary[particle]->Z();
          clusterE[particle] = primary[particle]->E();

          HCALEnergy[particle] = primary[particle]->GetHCALEnergy();
          HCALEisoR2[particle] = primary[particle]->GetIsoEnergyR2();
          HCALEisoR4[particle] = primary[particle]->GetIsoEnergyR4();

          for (int iseg = particle * 10; iseg < (particle * 10) + 10; iseg++)
          { // these 6 lines just added from Davide
            Width1[iseg] = primary[particle]->GetWidth1(iseg % 10);
            Width2[iseg] = primary[particle]->GetWidth2(iseg % 10);
            ClusterPhi[iseg] = primary[particle]->GetPhi(iseg % 10);
            E_long[iseg] = primary[particle]->GetSegmentEnergy(iseg % 10);
            SeedE[iseg] = primary[particle]->GetSeedEnergy(iseg % 10);
          }

          // Cluster shapes
          W1Cclust[particle] = 0;
          W2Cclust[particle] = 0;
          W1Fclust[particle] = 0;
          W2Fclust[particle] = 0;
          //	 PhiRCclust[particle] = 0;
          //	 PhiRFclust[particle] = 0;

          ERclust[particle] = E_long[particle * 10 + 2] + E_long[particle * 10 + 4]; // put this variable back in and everywhere above and in the ntuple (done)
          if (ERclust[particle] > 0)
            ERclust[particle] = E_long[particle * 10 + 5] / ERclust[particle];
          else
            ERclust[particle] = 1e3;
          float totEC = 0, totEF = 0;
          float totPhiEC = 0, totPhiEF = 0;
          // intermediate values for calcuating the angular concentration measure (radius):
          float PhiRXF = 0, PhiRYF = 0, PhiRXC = 0, PhiRYC = 0;
          for (int iSeg = 0; iSeg < 6; iSeg++)
          {
            if (iSeg == 1 || iSeg == 3)
            {
              if (Width1[particle * 10 + iSeg] > 0 && Width2[particle * 10 + iSeg] > 0)
              { // Skip clusters where width could not be determined
                W1Fclust[particle] += E_long[particle * 10 + iSeg] * Width1[particle * 10 + iSeg];
                W2Fclust[particle] += E_long[particle * 10 + iSeg] * Width2[particle * 10 + iSeg];
                totEF += E_long[particle * 10 + iSeg];

                float CP = ClusterPhi[particle * 10 + iSeg];

                if (TMath::Abs(CP) < 3.15 / 2.0)
                {

                  PhiRXF += E_long[particle * 10 + iSeg] * TMath::Cos(2.0 * CP); // !!! adjust Phi
                  PhiRYF += E_long[particle * 10 + iSeg] * TMath::Sin(2.0 * CP); // !!! adjust Phi
                  totPhiEF += E_long[particle * 10 + iSeg];
                }
              }
            }
            else
            {
              if (Width1[particle * 10 + iSeg] > 0 && Width2[particle * 10 + iSeg] > 0)
              { // Skip clusters where width could not be determined
                W1Cclust[particle] += E_long[particle * 10 + iSeg] * Width1[particle * 10 + iSeg];
                W2Cclust[particle] += E_long[particle * 10 + iSeg] * Width2[particle * 10 + iSeg];
                totEC += E_long[particle * 10 + iSeg];

                float CP = ClusterPhi[particle * 10 + iSeg];

                if (TMath::Abs(CP) < 3.15 / 2.0)
                {

                  PhiRXC += E_long[particle * 10 + iSeg] * TMath::Cos(2.0 * CP); // !!! adjust Phi
                  PhiRYC += E_long[particle * 10 + iSeg] * TMath::Sin(2.0 * CP); // !!! adjust Phi
                  totPhiEC += E_long[particle * 10 + iSeg];
                }
              }
            }
          }
          if (totEC > 0)
          {
            W1Cclust[particle] /= totEC;
            W2Cclust[particle] /= totEC;
          }
          if (totEF > 0)
          {
            W1Fclust[particle] /= totEF;
            W2Fclust[particle] /= totEF;
          }
          if (totPhiEC > 0)
          {
            PhiRXC /= totPhiEC;
            PhiRYC /= totPhiEC;
          }
          if (totPhiEF > 0)
          {
            PhiRXF /= totPhiEF;
            PhiRYF /= totPhiEF;
          }

          PhiRCclust[particle] = TMath::Sqrt(PhiRXC * PhiRXC + PhiRYC * PhiRYC);
          PhiRFclust[particle] = TMath::Sqrt(PhiRXF * PhiRXF + PhiRYF * PhiRYF);

          PhiACclust[particle] = TMath::ATan2(PhiRYC, PhiRXC) * 0.5;
          PhiAFclust[particle] = TMath::ATan2(PhiRYF, PhiRXF) * 0.5;

          // Cluster shapes making use of angular information
          VW1Cclust[particle] = 0;
          VW2Cclust[particle] = 0;
          VW1Fclust[particle] = 0;
          VW2Fclust[particle] = 0;

          // intermediate values: for calculating the vectorially-added widths:
          float W1, W2, CP, CPA, CosDP, SinDP, W1P, W2P;

          if (analysisDebugMode)
            std::cout << "Calculate vectorially added width " << std::endl;
          for (int iSeg = 0; iSeg < 6; iSeg++)
          {
            if (iSeg == 1 || iSeg == 3)
            {

              CP = ClusterPhi[particle * 10 + iSeg];
              W1 = Width1[particle * 10 + iSeg];
              W2 = Width2[particle * 10 + iSeg];
              if (W1 > 0 && W2 > 0)
              {
                CPA = PhiAFclust[particle];

                W1P = 0.0;
                W2P = 0.0;

                if (TMath::Abs(CP) < 3.15 / 2.0)
                {
                  // operations with phi difference come here
                  CosDP = TMath::Cos(CP - CPA);
                  SinDP = TMath::Sin(CP - CPA);
                  W1P = 1.0 / TMath::Sqrt(TMath::Power(CosDP / W1, 2) + TMath::Power(SinDP / W2, 2));
                  W2P = 1.0 / TMath::Sqrt(TMath::Power(SinDP / W1, 2) + TMath::Power(CosDP / W2, 2));
                }

                VW1Fclust[particle] += E_long[particle * 10 + iSeg] * W1P;
                VW2Fclust[particle] += E_long[particle * 10 + iSeg] * W2P;
              }
            }
            else
            {

              CP = ClusterPhi[particle * 10 + iSeg];
              W1 = Width1[particle * 10 + iSeg];
              W2 = Width2[particle * 10 + iSeg];
              if (W1 > 0 && W2 > 0)
              {
                CPA = PhiACclust[particle];

                W1P = 0.0;
                W2P = 0.0;

                if (TMath::Abs(CP) < 3.15 / 2.0)
                {
                  // operations with phi difference come here
                  CosDP = TMath::Cos(CP - CPA);
                  SinDP = TMath::Sin(CP - CPA);
                  W1P = 1.0 / TMath::Sqrt(TMath::Power(CosDP / W1, 2) + TMath::Power(SinDP / W2, 2));
                  W2P = 1.0 / TMath::Sqrt(TMath::Power(SinDP / W1, 2) + TMath::Power(CosDP / W2, 2));
                }

                VW1Cclust[particle] += E_long[particle * 10 + iSeg] * W1P;
                VW2Cclust[particle] += E_long[particle * 10 + iSeg] * W2P;
              }
            }
          }
          if (totPhiEC > 0)
          {
            VW1Cclust[particle] /= totPhiEC;
            VW2Cclust[particle] /= totPhiEC;
          }
          if (totPhiEF > 0)
          {
            VW1Fclust[particle] /= totPhiEF;
            VW2Fclust[particle] /= totPhiEF;
          }

          // computations of vectorially-added sigma_1 and sigma_2 should come here (above)

          double primaryE = primary[particle]->E();
          // Get Particle vector
          // Use double precision
          double cx = primary[particle]->X();
          double cy = primary[particle]->Y();
          double cz = primary[particle]->Z();
          double primaryD = TMath::Sqrt(cx * cx + cy * cy + cz * cz);
          // Create lorentz vector for this particle
          lor[particle].SetPxPyPzE(primaryE * cx / primaryD,
                                   primaryE * cy / primaryD,
                                   primaryE * cz / primaryD, primaryE);
        }
      } // end for particlesToStudy
      if (analysisDebugMode)
        std::cout << "End loop over particles" << std::endl;

      // Check if all fit
      allFit = true;
      for (int particle = 0; particle < particlesToStudy; particle++)
      {

        if (!primaryFits[particle])
          allFit = false;
      }

      // Calculate invariant mass
      if (particlesToStudy > 1)
      {

        bool allFound = true;
        for (int particle = 0; particle < particlesToStudy; particle++)
        {

          if (!primary[particle])
            allFound = false;
        }
        // Continue only if all primary particles has been matched with a cluster
        if (allFound)
        {
          // std::cout << "Calculating invariant mass" << std::endl;
          TLorentzVector resultVector;
          resultVector.SetPxPyPzE(0, 0, 0, 0);
          for (int particle = 0; particle < particlesToStudy; particle++)
          {
            resultVector += lor[particle];
          }

          invariantMass = resultVector.M();
          if (invariantMass < 0)
          {
            std::cout << "Invariant mass negative: E " << resultVector.E() << " p " << resultVector.P() << std::endl;
            for (int particle = 0; particle < particlesToStudy; particle++)
            {
              lor[particle].Print();
            }
          }
        }
        else
        {
          // std::cout << "Skipping calculation of invariant mass, not all decay products were found." << std::endl;
          invariantMass = 0;
        }
      }
      else
      {
        // std::cout << "Skipping calculation of invariant mass - single-photon event." << std::endl;
        invariantMass = 0;
      }

      // Now find out, how many 'fake' rejections would be there (using only the "primary" clusters)
      // std::cout << "Looking for invariant matches from background" << std::endl;
      for (int particle = 0; particle < particlesToStudy; particle++)
      {

        if (!primary[particle])
          continue;

        int xCol, yRow;
        // Search for cluster in the map around col,row position
        int maxRing = (int)(MAX_PAIR_DISTANCE / geometry->GetVirtualPadSize(infoSegment));
        for (int ring = 1; ring <= maxRing; ring++)
        {
          ringer.SetRing(ring);
          while (ringer.GetNext(xCol, yRow))
          {
            AliFOCALCluster *c = (AliFOCALCluster *)clusterMap->Get(realCol[particle] + xCol, realRow[particle] + yRow);
            if (c)
            {

              double secondaryE = c->E();
              // Get Particle vector
              double secondaryD = TMath::Sqrt(c->X() * c->X() +
                                              c->Y() * c->Y() +
                                              c->Z() * c->Z());
              // Create lorentz vector for this particle
              TLorentzVector secondaryLor;
              secondaryLor.SetPxPyPzE(secondaryE * c->X() / secondaryD,
                                      secondaryE * c->Y() / secondaryD,
                                      secondaryE * c->Z() / secondaryD, secondaryE);
              TLorentzVector resultVector;
              resultVector.SetPxPyPzE(0, 0, 0, 0);
              resultVector = lor[particle] + secondaryLor;
              if (TMath::Abs(resultVector.M() - PI0_MASS) < PI0_TOLERANCE * PI0_SIGMA / 2)
                rejections++;
            }
          } // end while ringer
        }   // end for ringer
      }     // end for particlesToStudy

      analysis->Fill();

    } // end events

    outputFile->cd();
    // analysis->Write();
    outputFile->Write();
    outputFile->Close();

  } // end folders

  delete clusterMap;
  std::cout << "DONE" << std::endl;
}
