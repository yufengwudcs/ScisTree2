//
//  ScistPerfPhyImp.hpp
//  
//
//  Created by Yufeng Wu on 7/27/18.
//
//

#ifndef ScistPerfPhyImp_hpp
#define ScistPerfPhyImp_hpp

#include <set>
#include <vector>
#include <map>
#include "ScistPerfPhyUtils.hpp"
#include "ScistGenotype.hpp"
#include "TreeBuilder.h"
#include "UtilsNumerical.h"

class PhyloDistance;
class MarginalTree;

// *************************************************************************************
// Utiltiies

void OutputMutationTree( const char *filenameMT, const string &strMutTree, bool fLabel );


// *************************************************************************************
// Build phylogeny by tree search

class ScistPerfPhyProbOnTree;

class ScistPerfPhyMLE
{
public:
    ScistPerfPhyMLE(ScistGenGenotypeMat &genos, const std::string *pNJTree = NULL);
    double Infer( std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces = NULL, std::string *pstrTreeNW = NULL );
    double InferFastNNI( std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces = NULL, std::string *pstrTreeNW = NULL, map<pair<int,int>, set<int> > *pmapChangedPosHistory = NULL );
    double InferFastSPR( const std::string &treeInit, std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces = NULL, std::string *pstrTreeNW = NULL, map<pair<int,int>, set<int> > *pmapChangedPosHistory = NULL );
    double InferFastSPRInitTrees( const std::vector<std::string> &listInitTreesMLE, std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces = NULL, std::string *pstrTreeNW = NULL );
    double InferReroot( const string &strTreeInit, std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces = NULL, std::string *pstrTreeNW = NULL );
    void SetVerbose(bool f) { fVerbose = f; }
    void SetBrOpt(bool f) { fOptBrLen = f; }
    void SetOutput(bool f) { fOutput = f; }
    void SetPPOut(bool f) { fOutputPPWithEdgeLabels = f; }
    void SetPPOutLabel(bool f) { fOutputLabel = f; }
    void SetSPR(bool f) { fSPR = f; }
    //void SetHeuSPRNum(int n) { maxHeuSPRNum = n; }
    void SetExactSPR(bool f) { fExactSPR = f; }
    void SetReroot(bool f) { fReroot = f; }
    void SetSPRSrcFrac(double f, int t) { fracSPRSrc = f; thresSPRDropStop = t; }
    void SetIterateMode(bool f) { fIterate = f; }
    void SetNumThreads(int n) { numThreads = n; }
    void SetCellNames( const std::vector<std::string> &listCellNamesIn ) { listCellNames=listCellNamesIn; }
    void SetSiteNames( const std::vector<std::string> &listSiteNamesIn ) { listSiteNames=listSiteNamesIn; }
    void SetMutTreeFileName(const std::string &strMutTreeFileNameIn) { this->strMutTreeFileName = strMutTreeFileNameIn; }
    void SetInitTrees(const std::vector<std::string> &listTreesIn) { this->listInitTrees = listTreesIn; }
    static void GetNgbrTreesFrom(int numHaps, const std::string &strTree, std::set<std::string> &setNgbrTrees );
    static void GetNgbrTreesFromSPR(int numHaps, const std::string &strTree, std::set<std::string> &setNgbrTrees );
    static std::string RemapLeafLbls(int numHaps, const std::string &strTree, const std::map<int,int> &mapLabels );
    
private:
    void Init();
    std::string ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const;
    void FindChangedGenos( int site, const std::pair<ScistPerfPhyCluster,ScistPerfPhyCluster> &clusToAdd, std::set< std::pair<std::pair<int,int>, int> > & listChangedPlaces ) const;
    static std::string RemapLeafLbls(int numHaps, const std::string &strTree0Based, const vector<int> &listLblsOld );
    double ScoreTree(const string &strTree, std::vector<std::pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > &listChangedCluster) const;
    double ScoreTreeMT(const string &strTree, std::vector<std::pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > &listChangedCluster) const;
    double ScoreTreeWithSite(ScistPerfPhyProbOnTree &probTree, MarginalTree &tree, int site, ScistPerfPhyCluster &clusChanged1, ScistPerfPhyCluster &clusChanged2) const;
    double CalcMaxProbUpperBound() const;
    double OptBranchLens( const std::string &strTree, std::string &strTreeBrOpt );
    double CalcChangedGenosProb(const std::set< std::pair<std::pair<int,int>, int> > & listChangedPlaces) const;
    std::string ConvCellTreeStr(const std::string &strTree) const;
    std::string ConvMutTreeStr(const std::string &strTree) const;
    double FastFindBestNNITreeFrom(const string &strTreeFrom, string &strBestNNITree);
    double FastFindBestRerootTreeFrom(const string &strTreeFrom, string &strBestNNITree);
    void CollectUncertainGenosFrom( const std::map<std::pair<int,int>, std::set<int> > &mapChangedPlacesHist, std::set<std::pair<int,int> > &setUncertainPos ) const;
    
    ScistGenGenotypeMat &genosInput;
    std::vector<ScistPerfPhyCluster> listClusMutsInputHetero;
    std::vector<ScistPerfPhyCluster> listClusMutsInputHomo;
    std::vector<int> listInputColMulti;
    ScistPerfPhyGuideTree treeGuide;
    bool fVerbose;
    bool fOptBrLen;
    bool fOutput;
    bool fOutputPPWithEdgeLabels;
    bool fOutputLabel;
    bool fSPR;
    bool fExactSPR;
    bool fReroot;
    double fracSPRSrc;
    int thresSPRDropStop;
    bool fIterate;      // if on, would run it again during search
    //int maxHeuSPRNum;
    int numThreads;
    std::string strNJTree;
    std::vector<double> listSitePriorScore;
    std::vector<std::string> listCellNames;
    std::vector<std::string> listSiteNames;
    std::string strMutTreeFileName;
    std::vector<std::string> listInitTrees;
};

// *************************************************************************************
// Build phylogeny by tree search with branch length

class ScistFullPerfPhyMLE : public NumericalAlgoUtils
{
public:
    ScistFullPerfPhyMLE(ScistGenGenotypeMat &genos);
    void Infer();
    void SetVerbose(bool f) { fVerbose = f; }
    virtual double EvaluateAt(double pt, void *pParam);
    double OptBranchLens(MarginalTree &tree);
    
private:
    void Init();
    double CalcLikelihoodOf(MarginalTree &tree) const;
    double CalcLikelihoodOf(ScistPerfPhyProbOnTree &sppp, int site, MarginalTree &tree, double totEdgeLen, const std::vector<std::set<int> > &listClades) const;
    std::string ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const;
    
    ScistGenGenotypeMat &genosInput;
    //std::vector<ScistPerfPhyCluster> listClusMutsInput;
    std::vector<ScistPerfPhyCluster> listClusMutsInputHetero;
    std::vector<ScistPerfPhyCluster> listClusMutsInputHomo;
    std::vector<int> listInputColMulti;
    ScistPerfPhyGuideTree treeGuide;
    bool fVerbose;
    std::vector<std::map<std::set<int>, double> > cacheProbMutClades;
    MarginalTree *pMargTreeOptBrLen;
    int brOptIndex;
};

// *************************************************************************************
// Tree probability

class ScistPerfPhyProbOnTree
{
public:
    ScistPerfPhyProbOnTree( ScistGenGenotypeMat &genos, MarginalTree &mtreeIn );
    double CalcProbMaxForSite(int site, ScistPerfPhyCluster &clusChangedMut, ScistPerfPhyCluster &clusChangedHomoMut) const;
    void PrepareProbMaxQ(std::vector< std::vector<double> > &listProbsSite, std::vector<double> &probPre ) const;
    double CalcProbMaxForSiteForTree(int site, MarginalTree &treeEval);
    void PrepareProbMaxQ( );
    double CalcProbForSite(int site, double totEdgeLen, const std::vector<std::set<int> > &listClades) const;
    
private:
    void Init();
    double CalcProbMaxForSiteHap(int site, ScistPerfPhyCluster &clusChanged) const;
    double CalcProbMaxForSiteGeno(int site, ScistPerfPhyCluster &clusChangedHetero, ScistPerfPhyCluster &clusChangedHomo) const;
    double CalcProbForSiteHap(int site, double totEdgeLen, const std::vector<std::set<int> > &listClades) const;
    double CalcProbForSiteGeno(int site, double totEdgeLen, const std::vector<std::set<int> > &listClades) const;
    
    ScistGenGenotypeMat &genosInput;
    ScistHaplotypeMat genosInputHap;
    MarginalTree &mtree;
    std::vector<double> listSitePriorScore;
    std::vector< std::vector<double> > listProbsSiteRepetitive;
    std::vector<double> listSiteMaxProbLeaf;
};



#endif /* ScistPerfPhyImp_hpp */
