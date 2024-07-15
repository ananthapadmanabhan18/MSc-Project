#include "TROOT.h"
#include "TSystem.h"



void RunClusterizer(Int_t startRun, Int_t endRun, Int_t startEvent, Int_t endEvent) {
    #include "ClusterizerParameters.C"
    gROOT->ProcessLine(Form(".x %s/LoadFOCAL.C(\"%s\")",librariesPrefix,librariesPrefix));

    gROOT->ProcessLine(Form(".x Clusterizer.C(%d,%d,%d,%d)",startRun, endRun, startEvent, endEvent));
}
