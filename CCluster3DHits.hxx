#ifndef CCluster3DHits_hxx_seen
#define CCluster3DHits_hxx_seen

#include "TROOT.h"
#include "TObject.h"
#include "TMath.h"
#include "CHit3D.hxx"
#include "CCluster3D.hxx"
#include <iostream>
#include <vector>
#include "TVector3.h"

#include "TTmplDensityCluster.hxx"
#include "TTmplMinimalSpanningTree.hxx"

#define DebugPlot


template<typename T>
class CPositionMetric{
public:
    double operator()(T lhs,T rhs){
        double x1=lhs.GetPosition().X();
        double y1=lhs.GetPosition().Y();
        double z1=lhs.GetPosition().Z();
        
        double x2=rhs.GetPosition().X();
        double y2=rhs.GetPosition().Y();
        double z2=rhs.GetPosition().Z();
        
        return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
    }
    
};

void CCluster3DHits(TTree* tree3D){
    
    const unsigned int minHits=2;
    const double maxDist=3;
    
    std::vector<CHit3D>* hits3D=0;
    std::vector<CHit2D>* hits2D=0;

#ifdef DebugPlot
    TCanvas* c6 = new TCanvas("c6","A Simple Graph Example6",200,10,700,500);
    
    TH3F* plot3D = new TH3F("3DClustersRECON","3DClustersRECON",150,30,180,80,160,240,120,120,240);
    //TH3F* plot3D = new TH3F("3DHitsRECON","3DHitsRECON",240,0,240,240,0,240,240,0,240);
    plot3D->GetXaxis()->SetTitle("X,cm");
    plot3D->GetYaxis()->SetTitle("Y,cm");
    plot3D->GetZaxis()->SetTitle("Z,cm");
#endif
    
    tree3D->SetBranchAddress("3DHits",&hits3D);
    tree3D->SetBranchAddress("Unused2DHits",&hits2D);

    
    std::vector<CCluster3D> clusters3D;
    std::vector<CHit3D> unused3D;
    std::vector<CHit2D> unused2D;
    
    TFile* hfile = new TFile("FileWith3DClusters.root","RECREATE");
    TTree *newTree = new TTree("treeWith3DClusters","tree with 3D clusters");
    
    newTree->Branch("3DClusters","vector<CCluster3D>",&clusters3D,8000,1);
    newTree->Branch("Unused3DHits","vector<CHit3D>",&unused3D,8000,1);
    newTree->Branch("Unused2DHits","vector<CHit2D>",&unused2D,8000,1);
    
    int nentry = tree3D->GetEntries();
    for(int i=0;i<nentry;++i){
        
        unused2D=(*hits2D);
        
        tree3D->GetEntry(i);
        
        typedef TTmplDensityCluster<CHit3D,CPositionMetric<CHit3D>> AlgDBScan;
        typedef TTmplMinimalSpanningTree<CHit3D,CPositionMetric<CHit3D>,CHit3D> AlgMST;
        
        
        std::shared_ptr<AlgDBScan> clus(new AlgDBScan(minHits,maxDist));
        
        clus->Cluster(*hits3D);
        
        unsigned int nClusters=clus->GetClusterCount();
        
        std::cout<<"Number of Clusters Formed="<<nClusters<<std::endl;
        
        std::vector<CHit3D> used3DHits;
        
        for(std::size_t c=0;c<nClusters;++c){

            std::vector<CHit3D> clusteredHits=clus->GetPoints(c);

            std::shared_ptr<AlgMST> mstAlg(new AlgMST);
            std::cout<<"Cluster Size="<<clusteredHits.size()<<std::endl;
            mstAlg->AddVertices(clusteredHits.begin(),clusteredHits.end());
            mstAlg->MakeTree(*clusteredHits.begin());
            
            const AlgMST::MinimalSpanningTree& tree = mstAlg->GetTree();
            
            int placeInTree=-2;
            int topDepth=-2;
            
            for (int t = 0; t < tree.size(); ++t) {
                if(topDepth<tree[t].VertexDepth){
                    topDepth=tree[t].VertexDepth;
                    placeInTree=t;
                }
            }
            
            std::shared_ptr<AlgMST> mstAlg2(new AlgMST);
            mstAlg2->AddVertices(clusteredHits.begin(),clusteredHits.end());
            mstAlg2->MakeTree(*std::find(clusteredHits.begin(),clusteredHits.end(),tree[placeInTree].Object));
            const AlgMST::MinimalSpanningTree& tree2 = mstAlg2->GetTree();
            
            int placeInTree2=-2;
            int topDepth2=-2;
            
            for (int t = 0; t < tree2.size(); ++t) {
                if(topDepth2<tree2[t].VertexDepth){
                    topDepth2=tree2[t].VertexDepth;
                    placeInTree2=t;
                }
            }
#ifdef DebugPlot
            std::cout<<"ROOT:"<<tree2[0].VertexDepth<<" X="<<tree2[0].Object.GetPosition().X()<<" Y="<<tree2[0].Object.GetPosition().Y()<<" Z="<<tree2[0].Object.GetPosition().Z()<<std::endl;
            
                        std::cout<<"TOPH:"<<tree2[placeInTree2].VertexDepth<<" X="<<tree2[placeInTree2].Object.GetPosition().X()<<" Y="<<tree2[placeInTree2].Object.GetPosition().Y()<<" Z="<<tree2[placeInTree2].Object.GetPosition().Z()<<std::endl;
#endif
            for(std::size_t h=0;h<clusteredHits.size();++h){
                used3DHits.push_back(clusteredHits[h]);
            plot3D->Fill(clusteredHits[h].GetPosition().X(),clusteredHits[h].GetPosition().Y(),clusteredHits[h].GetPosition().Z(),(c+1)*50);
            }
        }
        
        
        for(std::vector<CHit3D>::iterator h=hits3D->begin();h!=hits3D->end();++h){
            if(std::find(used3DHits.begin(),used3DHits.end(),*h)!=used3DHits.end())continue;
            unused3D.push_back(*h);
        }
        
#ifdef DebugPlot
        c6->cd();
        plot3D->Draw("LEGO2Z");
        gPad->Print("Clusters3D_rec.gif");
        for(int i=0;i<360;i=i+50){
            std::cout<<i<<std::endl;
            gPad->SetPhi(30+i);
            gPad->Print("Clusters3D_rec.gif+");
        }
#endif
        
        newTree->Fill();
        
        unused3D.clear();
        unused2D.clear();
        clusters3D.clear();
    }
    hfile->Close();
    delete hfile;
};


#endif
