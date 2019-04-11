#ifndef CInputRead_hxx_seen
#define CInputRead_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include "CHit2D.hxx"
#include <iostream>
#include <vector>

#include "TH2F.h"
#include "TPad.h"
#include "TH1.h"
#include "TVector3.h"

bool SortHits3DTrue(const TVector3& lhs,  const TVector3& rhs){
    
    if(lhs.X()!=rhs.X()){
        return lhs.X()<rhs.X();
    }else if(lhs.X()==rhs.X()){
        if(lhs.Z()!=rhs.Z()){
            return lhs.Z()<rhs.Z();
        }else if(lhs.Z()==rhs.Z()){
            return lhs.Y()<rhs.Y();
        }
    }
    return false;
    
};



void CInputRead(TTree *tree,int& nevts,int& skip){
    

    std::vector<CHit2D> hitsXY;
    std::vector<CHit2D> hitsXZ;
    std::vector<CHit2D> hitsYZ;


    
    int event;
    double hitLocation[1000][3],hitT[1000][3];
    int hitPE[1000][3];
    
    tree->SetBranchAddress("event",&event);
    tree->SetBranchAddress("hitLocation",&hitLocation);
    tree->SetBranchAddress("hitPE_measure",&hitPE);
    tree->SetBranchAddress("hitT",&hitT);
    
    TFile* hfile2 = new TFile("FileWith2DHits.root","RECREATE");
    TTree *newTree = new TTree("treeWith2DHits","tree with 2D hits");
    
    newTree->Branch("2DHitsXY","vector<CHit2D>",&hitsXY,8000,1);
    newTree->Branch("2DHitsXZ","vector<CHit2D>",&hitsXZ,8000,1);
    newTree->Branch("2DHitsYZ","vector<CHit2D>",&hitsYZ,8000,1);

    int nevents = tree->GetEntries();
    if(nevts>nevents){
        std::cout<<"process max possible number of events : "<<nevents<<std::endl;
        nevts=nevents;
    }
    
  
    for(int evt=skip;evt<nevts;++evt){
        tree->GetEntry(evt);
        

    TH2F* chargeXYplot = new TH2F("chargeXY","chargeXY",240,0,240,240,0,240);
        chargeXYplot->GetXaxis()->SetTitle("X,cm");
        chargeXYplot->GetYaxis()->SetTitle("Y,cm");
    TH2F* chargeXZplot = new TH2F("chargeXZ","chargeXZ",240,0,240,240,0,240);
        chargeXZplot->GetXaxis()->SetTitle("X,cm");
        chargeXZplot->GetYaxis()->SetTitle("Z,cm");
    TH2F* chargeYZplot = new TH2F("chargeYZ","chargeYZ",240,0,240,240,0,240);
        chargeYZplot->GetXaxis()->SetTitle("Y,cm");
        chargeYZplot->GetYaxis()->SetTitle("Z,cm");
    
    TH2F* timeXYplot = new TH2F("timeXY","timeXY",240,0,240,240,0,240);
        timeXYplot->GetXaxis()->SetTitle("X,cm");
        timeXYplot->GetYaxis()->SetTitle("Y,cm");
    TH2F* timeXZplot = new TH2F("timeXZ","timeXZ",240,0,240,240,0,240);
        timeXZplot->GetXaxis()->SetTitle("X,cm");
        timeXZplot->GetYaxis()->SetTitle("Z,cm");
    TH2F* timeYZplot = new TH2F("timeYZ","timeYZ",240,0,240,240,0,240);
        timeYZplot->GetXaxis()->SetTitle("Y,cm");
        timeYZplot->GetYaxis()->SetTitle("Z,cm");

    std::vector<std::pair<int,std::shared_ptr<std::vector<int>>>> trueInfoXY;
     std::vector<std::pair<int,std::shared_ptr<std::vector<int>>>> trueInfoXZ;
     std::vector<std::pair<int,std::shared_ptr<std::vector<int>>>> trueInfoYZ;
    
    for(std::size_t i=0;i<240*240+10;++i){
        std::shared_ptr<std::vector<int>> vect1(new std::vector<int>);
        trueInfoXY.push_back(std::make_pair(i,vect1));
        
       std::shared_ptr<std::vector<int>> vect2(new std::vector<int>);
        trueInfoXZ.push_back(std::make_pair(i,vect2));
        
        std::shared_ptr<std::vector<int>> vect3(new std::vector<int>);
        trueInfoYZ.push_back(std::make_pair(i,vect3));
    }
    
    for(int i=0;i<1000;i++){
        if(hitLocation[i][0]>0){

            TVector3 position(hitLocation[i][0],hitLocation[i][1],hitLocation[i][2]);
            

            double chargeXY = (double)hitPE[i][2];
            double chargeXZ = (double)hitPE[i][1];
            double chargeYZ = (double)hitPE[i][0];
            
            double timeXY = hitT[i][2];
            double timeXZ = hitT[i][1];
            double timeYZ = hitT[i][0];
            int ibinXY=chargeXYplot->FindBin(position.X()/10,position.Y()/10);
            chargeXYplot->AddBinContent(chargeXYplot->FindBin(position.X()/10,position.Y()/10),chargeXY);
            double pT1 = timeXYplot->GetBinContent(timeXYplot->FindBin(position.X()/10,position.Y()/10));
            if(pT1!=0 && pT1>timeXY)timeXYplot->SetBinContent(timeXYplot->FindBin(position.X()/10,position.Y()/10),timeXY);
            if(pT1==0)timeXYplot->SetBinContent(timeXYplot->FindBin(position.X()/10,position.Y()/10),timeXY);
            
            for(std::size_t j=0;j<trueInfoXY.size();++j){
                if(ibinXY==trueInfoXY[j].first){
                    trueInfoXY[j].second->push_back(i);
                    break;
                }
            }
 

            
            int ibinXZ=chargeXZplot->FindBin(position.X()/10,position.Z()/10);
            chargeXZplot->AddBinContent(chargeXZplot->FindBin(position.X()/10,position.Z()/10),chargeXZ);
            double pT2 = timeXZplot->GetBinContent(timeXZplot->FindBin(position.X()/10,position.Z()/10));
            if(pT2!=0 && pT2>timeXZ)timeXZplot->SetBinContent(timeXZplot->FindBin(position.X()/10,position.Z()/10),timeXZ);
            if(pT2==0)timeXZplot->SetBinContent(timeXZplot->FindBin(position.X()/10,position.Z()/10),timeXZ);
            
            for(std::size_t j=0;j<trueInfoXZ.size();++j){
                if(ibinXZ==trueInfoXZ[j].first){
                    trueInfoXZ[j].second->push_back(i);
                    break;
                }
            }
            
            int ibinYZ=chargeYZplot->FindBin(position.Y()/10,position.Z()/10);
            chargeYZplot->AddBinContent(chargeYZplot->FindBin(position.Y()/10,position.Z()/10),chargeYZ);
            double pT3 = timeYZplot->GetBinContent(timeYZplot->FindBin(position.Y()/10,position.Z()/10));
            if(pT3!=0 && pT3>timeYZ)timeYZplot->SetBinContent(timeYZplot->FindBin(position.Y()/10,position.Z()/10),timeYZ);
            if(pT3==0)timeYZplot->SetBinContent(timeYZplot->FindBin(position.Y()/10,position.Z()/10),timeYZ);
            
            for(std::size_t j=0;j<trueInfoYZ.size();++j){
                if(ibinYZ==trueInfoYZ[j].first){
                    trueInfoYZ[j].second->push_back(i);
                    break;
                }
            }
        }
    }
        int nXY=0;
        int nXZ=0;
        int nYZ=0;
    for(int i=1; i<240*240+2;++i){
        if(chargeXYplot->GetBinContent(trueInfoXY[i-1].first)>0){
           CHit2D hit;
        int nx  = chargeXYplot->GetXaxis()->GetNbins()+2;
        int ny  = chargeXYplot->GetYaxis()->GetNbins()+2;
        int binx = trueInfoXY[i-1].first%nx;
        int biny = ((trueInfoXY[i-1].first-binx)/nx)%ny;
            hit.SetRow(chargeXYplot->GetXaxis()->GetBinLowEdge(binx)+0.5);
            hit.SetColumn(chargeXYplot->GetYaxis()->GetBinLowEdge(biny)+0.5);
            hit.SetCharge(chargeXYplot->GetBinContent(trueInfoXY[i-1].first));
            hit.SetTime(timeXYplot->GetBinContent(trueInfoXY[i-1].first));
            hit.SetConstituents(*trueInfoXY[i-1].second);
            hit.SetPlane(2);
            hit.SetId(nXY);
            hitsXY.push_back(hit);
            nXY++;
        }
        
        if(chargeXZplot->GetBinContent(trueInfoXZ[i-1].first)>0){
            CHit2D hit;
            int nx  = chargeXZplot->GetXaxis()->GetNbins()+2;
            int ny  = chargeXZplot->GetYaxis()->GetNbins()+2;
            int binx = trueInfoXZ[i-1].first%nx;
            int biny = ((trueInfoXZ[i-1].first-binx)/nx)%ny;
            hit.SetRow(chargeXZplot->GetXaxis()->GetBinLowEdge(binx)+0.5);
            hit.SetColumn(chargeXZplot->GetYaxis()->GetBinLowEdge(biny)+0.5);
            hit.SetCharge(chargeXZplot->GetBinContent(trueInfoXZ[i-1].first));
            hit.SetTime(timeXZplot->GetBinContent(trueInfoXZ[i-1].first));
            hit.SetConstituents(*trueInfoXZ[i-1].second);
            hit.SetPlane(1);
            hit.SetId(nXZ);
            hitsXZ.push_back(hit);
            nXZ++;
        }
        
        if(chargeYZplot->GetBinContent(trueInfoYZ[i-1].first)>0){
            CHit2D hit;
            int nx  = chargeYZplot->GetXaxis()->GetNbins()+2;
            int ny  = chargeYZplot->GetYaxis()->GetNbins()+2;
            int binx = trueInfoYZ[i-1].first%nx;
            int biny = ((trueInfoYZ[i-1].first-binx)/nx)%ny;
            hit.SetRow(chargeYZplot->GetXaxis()->GetBinLowEdge(binx)+0.5);
            hit.SetColumn(chargeYZplot->GetYaxis()->GetBinLowEdge(biny)+0.5);
            hit.SetCharge(chargeYZplot->GetBinContent(trueInfoYZ[i-1].first));
            hit.SetTime(timeYZplot->GetBinContent(trueInfoYZ[i-1].first));
            hit.SetConstituents(*trueInfoYZ[i-1].second);
            hit.SetPlane(0);
            hit.SetId(nYZ);
            hitsYZ.push_back(hit);
            nYZ++;
        }
        
    }
        

        delete chargeXZplot;
        delete chargeXYplot;
        delete chargeYZplot;
        delete timeXYplot;
        delete timeXZplot;
        delete timeYZplot;
        
        newTree->Fill();
        hfile2->Write();
        
        hitsXY.clear();
        hitsXZ.clear();
        hitsYZ.clear();
        
        trueInfoXY.clear();
        trueInfoXZ.clear();
        trueInfoYZ.clear();
        
        
    }

    
    hfile2->Close();
    delete hfile2;

    
};


#endif
