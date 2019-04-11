#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <array>
#include <utility>

#include "Riostream.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TArrow.h"

#include "CHit2D.hxx"
#include "CHit3D.hxx"
#include "CCluster3D.hxx"
#include "CInputRead.hxx"
#include "CCreate3DHits.hxx"
#include "CCluster3DHits.hxx"

int CRecon(){
    
    if (!TClass::GetDict("CHit2D")) {
        gROOT->ProcessLine(".L CHit2D.cxx++");
    }
    
    if (!TClass::GetDict("CHit3D")) {
        gROOT->ProcessLine(".L CHit3D.cxx++");
    }
    
    if (!TClass::GetDict("CCluster3D")) {
        gROOT->ProcessLine(".L CCluster3D.cxx++");
    }
    
    TFile* hfile = new TFile("testEvent_3DST+emptyECAL_event_222_sampleT_v2.root","READ");
    TTree* tree = (TTree*)hfile->Get("EDepSimTree");
    
    int nprocessEvents = 1;
    int skip=0;
    nprocessEvents+=skip;
    CInputRead(tree,nprocessEvents,skip);
    
    hfile->Close();
    delete hfile;
    
    TFile* hfile2D = new TFile("FileWith2DHits.root","READ");
    TTree* tree2D = (TTree*)hfile2D->Get("treeWith2DHits");
    
    CCreate3DHits(tree2D);

    hfile2D->Close();
    delete hfile2D;
    
    TFile* hfile3D = new TFile("FileWith3DHits.root","READ");
    TTree* tree3D = (TTree*)hfile3D->Get("treeWith3DHitt");
    
    CCluster3DHits(tree3D);
    
   // hfile3D->Close();
   // delete hfile3D;
    
    
    return 0;
}
