//
//  ScistPerfPhyImp.cpp
//  
//
//  Created by Yufeng Wu on 7/27/18.
//
//

#include "ScistPerfPhyImp.hpp"
#include "ScistGenotype.hpp"
#include "Utils3.h"
#include "Utils4.h"
#include "TreeBuilder.h"
#include "MarginalTree.h"
#include "RBT.h"
#include "PhylogenyTree.h"
#include <iomanip>
#include <cmath>
#include <thread>
#include "UtilsNumerical.h"
#include "ScistLocalSearch.hpp"

//const int MAX_SPR_OP = 1;


// *************************************************************************************
// Utiltiies

void OutputMutationTree( const char *filenameMT, const string &strMutTree, bool fLabel )
{
    PhylogenyTreeBasic treeMut;
    treeMut.ConsOnNewickEdgeLabelTree( strMutTree );
    if( fLabel )
    {
        treeMut.OutputGML(filenameMT);
    }
    else
    {
        treeMut.OutputGMLNoLabel(filenameMT);
    }
}

// *************************************************************************************
// Build phylogeny by tree search with branch length


ScistFullPerfPhyMLE :: ScistFullPerfPhyMLE(ScistGenGenotypeMat &genos) : genosInput(genos), fVerbose(false), pMargTreeOptBrLen(NULL), brOptIndex(-1)
{
    Init();
}

void ScistFullPerfPhyMLE :: Infer()
{
    set<ScistPerfPhyCluster> setClusAllGuide;
    this->treeGuide.GetAllClusters( setClusAllGuide );
    string strTreeOpt = ConsTreeFromSetClusters( setClusAllGuide );
    
    MarginalTree treeOpt;
    ReadinMarginalTreesNewickWLenString( strTreeOpt, this->genosInput.GetNumHaps(), treeOpt );
    treeOpt.InitUnitEdgelen();
    
    //double loglikeliOptInit = CalcLikelihoodOf(treeOpt);
    
    // optimize branch length
    double loglikeliOptBr = OptBranchLens(treeOpt);
    strTreeOpt = treeOpt.GetNewickSorted(true);
//cout << "Initial tree: "  << treeOpt.GetNewick() << ", log-likelihood: " << loglikeliOptBr << endl;
    
    set<string> setTreeSearchedBefore;
    setTreeSearchedBefore.insert( strTreeOpt );
    
    // now search for neighborhood of the current tree to optimize the tree
    while(true)
    {
        set<string> setNgbrTrees;
        //GetNgbrTreesFromSPR( this->genosInput.GetNumHaps(), strTreeOpt, setNgbrTrees );
        ScistPerfPhyMLE :: GetNgbrTreesFrom( this->genosInput.GetNumHaps(), strTreeOpt, setNgbrTrees );
        if( fVerbose)
        {
            cout << "Current best likelihood: " << loglikeliOptBr << ", current tree: " << treeOpt.GetNewickSorted(true) << ", tree neighborhood size: " << setNgbrTrees.size() << endl;
        }
        bool fCont = false;
        for(set<string> :: iterator it = setNgbrTrees.begin(); it != setNgbrTrees.end(); ++it)
        {
            if( setTreeSearchedBefore.find(*it) != setTreeSearchedBefore.end() )
            {
                continue;
            }
            setTreeSearchedBefore.insert(*it);
            
//cout << "Neighbor tree: " << *it << endl;
            MarginalTree treeStep;
            ReadinMarginalTreesNewickWLenString( *it, this->genosInput.GetNumHaps(), treeStep );
            //treeStep.InitUnitEdgelen();
//cout << "treeStep: " << treeStep.GetNewick() << endl;
            double loglikeliStep = OptBranchLens( treeStep );
            //double loglikeliStep = CalcLikelihoodOf( treeStep );
//cout << ", loglikeliStep (w/ branch length optimization): " << loglikeliStep << endl;
            if( loglikeliStep > loglikeliOptBr )
            {
                //cout << "BETTER.\n";
                loglikeliOptBr = loglikeliStep;
                strTreeOpt = *it;
                treeOpt = treeStep;
                fCont = true;
            }
        }
        if( fCont == false )
        {
            break;
        }
    }
    
    cout << "**** Optimal cost: " << loglikeliOptBr << endl;
    cout << "Constructed single cell phylogeny: " << treeOpt.GetNewickSorted(false) << endl;
    cout << "With branch length: " << treeOpt.GetNewickSorted(true) << endl;
}

void ScistFullPerfPhyMLE :: Init()
{
    //
    cacheProbMutClades.resize( genosInput.GetNumSites() );
    // get all clusters
    //listClusMutsInput.clear();
    //for(int s=0; s<genosInput.GetNumSites(); ++s)
    //{
    //    set<int> muts;
    //    genosInput.GetMutRowsHapAtSite(s, muts);
    //    ScistPerfPhyCluster clus(muts);
    //    listClusMutsInput.push_back(clus);
    //}
    listClusMutsInputHetero.clear();
    listClusMutsInputHomo.clear();
    for(int s=0; s<genosInput.GetNumSites(); ++s)
    {
        set<int> muts;
        genosInput.GetRowsWithGenoAtSite(s, 1, muts);
        ScistPerfPhyCluster clus(muts);
        listClusMutsInputHetero.push_back(clus);
        
        set<int> muts2;
        genosInput.GetRowsWithGenoAtSite(s, 2, muts2);
        ScistPerfPhyCluster clus2(muts2);
        listClusMutsInputHomo.push_back(clus2);
    }
    
    this->genosInput.GetColMultiplicityMap( listInputColMulti );
    
    // construct NJ tree as the initial tree
    string strNJ = this->genosInput.ConsNJTreeZeroRoot ();
    this->treeGuide.Init(strNJ);
}

double ScistFullPerfPhyMLE :: OptBranchLens(MarginalTree &tree)
{
    //
    this->pMargTreeOptBrLen = &tree;
    
    const double MIN_BR_LEN = 0.01;
    const double MAX_BR_LEN = 10.0;
    const double TOLNUM = 0.2;
    
    double loglikeliRes = -1.0*HAP_MAX_INT;
    
    // optimize branch of each once and only once
    for(int br = 0; br<tree.GetTotNodesNum(); ++br)
    {
        if( br == tree.GetRoot() )
        {
            continue;
        }
        this->brOptIndex = br;
        double brLen = tree.GetEdgeLen( br );
        double brNew = brLen;
        double likeliMax = -1.0*Func1DMinBrent(MIN_BR_LEN, brLen, MAX_BR_LEN, TOLNUM, &brNew );
        if( likeliMax > loglikeliRes )
        {
            loglikeliRes = likeliMax;
            tree.SetBranchLen(br, brNew);
        }
        else
        {
            tree.SetBranchLen(br, brLen);
        }
    }
    return loglikeliRes;
}

double ScistFullPerfPhyMLE :: EvaluateAt(double pt, void *pParam)
{
    //
    YW_ASSERT_INFO( pMargTreeOptBrLen != NULL, "Tree to opt branch: null" );
    YW_ASSERT_INFO( brOptIndex >= 0, "Branch opt not set");
    pMargTreeOptBrLen->SetBranchLen(brOptIndex, pt);
    return -1.0*CalcLikelihoodOf( *pMargTreeOptBrLen );
}

double ScistFullPerfPhyMLE :: CalcLikelihoodOf( MarginalTree &tree) const
{
    set<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > setClusDone;
    double res = 0.0;
    
    vector< set<int> > listClades;
    tree.ConsDecedentLeavesInfoLabels(listClades);
    for(int i=0; i<(int)listClades.size(); ++i)
    {
        DecAllNumInSet( listClades[i] );
//cout << "Tree clade: ";
//DumpIntSet(listClades[i]);
    }
    double totEdgeLen = tree.GetTotEdgeLen();
    ScistPerfPhyProbOnTree sppp(this->genosInput, tree);
    
    for(int site=0; site<genosInput.GetNumSites(); ++site)
    {
        pair<ScistPerfPhyCluster,ScistPerfPhyCluster> pp( listClusMutsInputHetero[site], listClusMutsInputHomo[site] );
        if( setClusDone.find(pp) != setClusDone.end() )
        {
            continue;
        }
        int multi = this->listInputColMulti[site];
        double loglikeliSite = CalcLikelihoodOf(sppp, site, tree, totEdgeLen, listClades);
        res += loglikeliSite*multi;
        setClusDone.insert( pp );
    }
    
    return res;
}

double ScistFullPerfPhyMLE :: CalcLikelihoodOf(ScistPerfPhyProbOnTree &sppp, int site, MarginalTree &tree, double totEdgeLen, const vector<set<int> > &listClades) const
{
    return sppp.CalcProbForSite( site, totEdgeLen, listClades );
}

std::string ScistFullPerfPhyMLE :: ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const
{
    //
    // now construct tree
    ScistInfPerfPhyUtils treeBuild;
    map<int, ScistPerfPhyCluster> mapPickedClus;
    int s = 0;
    for(set<ScistPerfPhyCluster>:: iterator it = setClusters.begin(); it != setClusters.end(); ++it)
    {
        mapPickedClus[s] = *it;
        ++s;
    }
    string strTree = treeBuild.ConsTreeWCombDistClus( this->genosInput, mapPickedClus );
    return strTree;
}

// *************************************************************************************
// Build phylogeny by tree search


ScistPerfPhyMLE :: ScistPerfPhyMLE(ScistGenGenotypeMat &genos, const std::string *pNJTree) : genosInput(genos), fVerbose(false), fOptBrLen(false), fOutput(true), fOutputPPWithEdgeLabels(false), fOutputLabel(true), fSPR(false), /*maxHeuSPRNum(0), */ fExactSPR(false), fReroot(false), fracSPRSrc(0.5), thresSPRDropStop(100), fIterate(false), numThreads(1)
{
    if( pNJTree != NULL )
    {
        strNJTree = *pNJTree;
    }
    
    Init();
}

double ScistPerfPhyMLE :: Infer( std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces, std::string *pstrTreeNW )
{
//cout << "ScistPerfPhyMLE :: Infer\n";

    // Jan 16, 2023: now switching to fast NNI search
    double loglikeliBest = MAX_NEG_DOUBLE_VAL;
    string strTreeBest;
    // keep track of changed position histories
    map< pair<int,int>, set<int> > mapChangedPosHistory;
    
    // YW: for now, just iteratively run once if iterative mode is on
    int numRounds = 0;
    //const int NUM_ROUNDS = 2;
    //while( ++numRounds <= NUM_ROUNDS )
    while(true)
    {
        ++numRounds;
        if( fVerbose && fIterate )
        {
            cout << "^^^^^^^ ITERATIVELY SEARCH ROUND " << numRounds << "^^^^^^^^^^^^\n";
            if( numRounds == 1 )
            {
                //cout << "Init tree: " << this->strNJTree << endl;
            }
        }
        // always perform NNI search first
        string strOptNNITree = strNJTree;
        std::set< std::pair<std::pair<int,int>, int> > listChangedPlacesStep;
        double loglikeliBestStep = MAX_NEG_DOUBLE_VAL;
        
        if( fSPR == false )
        {
            loglikeliBestStep = InferFastNNI(&listChangedPlacesStep, &strOptNNITree, &mapChangedPosHistory);
            if( fVerbose )
            {
                cout << "Maximum log-likelihood by NNI local search: " << loglikeliBestStep
/*<< ", maximum log-likelihood so far: " << loglikeliBest */ << endl;
                if( numRounds == 1 )
                {
                    cout << "Best NNI local search tree: " << strOptNNITree << endl;
                }
            }
        }
        
#if 0
        // always perform NNI search first
        string strOptNNITree;
        std::set< std::pair<std::pair<int,int>, int> > listChangedPlacesStep;
        double loglikeliBestStep = InferFastNNI(&listChangedPlacesStep, &strOptNNITree, &mapChangedPosHistory);
        if( fVerbose )
        {
            cout << "Maximum log-likelihood by NNI local search: " << loglikeliBestStep /*<< ", maximum log-likelihood so far: " << loglikeliBest */ << endl;
            if( numRounds == 1 )
            {
                cout << "Best NNI local search tree: " << strOptNNITree << endl;
            }
        }
        //if( pstrTreeNW != NULL)
        //{
        //    *pstrTreeNW = strOptNNITree;
        //}
#endif
        
        // if SPR is on, run SPR from the NNI tree
        if( fSPR )
        {
            string treeInit = strOptNNITree;
            //vector<string> listTrees;
            //listTrees.push_back(*pstrTreeNW);
            loglikeliBestStep = InferFastSPR(treeInit, &listChangedPlacesStep, &strOptNNITree, &mapChangedPosHistory);
            if( fVerbose )
            {
                //cout << "Best SPR local search tree: " << strOptNNITree << " with log likelihood:" << loglikeliBestStep << endl;
            }
        }
        
        // run reroot if turned on
        if( fReroot )
        {
            string treeInit = strOptNNITree;
            //vector<string> listTrees;
            //listTrees.push_back(*pstrTreeNW);
            loglikeliBestStep = InferReroot(treeInit, &listChangedPlacesStep, &strOptNNITree);
            if( fVerbose )
            {
                //cout << "Best tree afer re-rooting: " << strOptNNITree << " with log likelihood:" << loglikeliBestStep << endl;
            }
            //exit(1);
        }
        
        // if no better, stop
        const double MIN_INC_LOG_LIKELI = 0.01;
        if( loglikeliBestStep  <= loglikeliBest + MIN_INC_LOG_LIKELI )
        {
            break;
        }
        
        if( fVerbose )
        {
            cout << "+++ Likelihood improved to:" << loglikeliBestStep << endl;
        }
        
        // take this
        loglikeliBest = loglikeliBestStep;
        if( pstrTreeNW != NULL)
        {
            *pstrTreeNW = strOptNNITree;
        }
        if( plistChangedPlaces != NULL )
        {
            *plistChangedPlaces = listChangedPlacesStep;
        }
        //this->genosInput.ResetMaximalGenos();
        for(set<pair<pair<int,int>,int> > :: iterator it = listChangedPlacesStep.begin(); it != listChangedPlacesStep.end(); ++it )
        {
            this->genosInput.SetGenotypeAt(it->first.first, it->first.second, it->second);
        }
        
#if 0
        if( fVerbose )
        {
            if(  fOutput)
            {
                //cout << "Genotypes called by maximal single position probability\n";
                //const string strDesc = "Single-site maximal probability genotypes";
                //this->genosInput.OutputImput(&strDesc);
            }
#if 0
            cout << "List of corrected genotypes (site, cell, new genotype) in base-1 (from maximal single position probability genotypes): \n";
            for(set<pair<pair<int,int>,int> > :: iterator it = listChangedPlacesStep.begin(); it != listChangedPlacesStep.end(); ++it )
            {
                cout << "[ " << setw(6) << it->first.second+1 << " " << setw(6) << it->first.first+1 << " ]: " << it->second << endl;
            }
#endif
        }
        
#if 0
        // change genotype
        if( plistChangedPlaces != NULL )
        {
            for(set<pair<pair<int,int>,int> > :: iterator it = plistChangedPlaces->begin(); it != plistChangedPlaces->end(); ++it )
            {
                this->genosInput.SetGenotypeAt(it->first.first, it->first.second, it->second);
            }
            
            if( fVerbose )
            {
                if(  fOutput)
                {
                    //cout << "Genotypes called by maximal single position probability\n";
                    //const string strDesc = "Single-site maximal probability genotypes";
                    //this->genosInput.OutputImput(&strDesc);
                }
                
                cout << "List of corrected genotypes (site, cell, new genotype) in base-1 (from maximal single position probability genotypes): \n";
                for(set<pair<pair<int,int>,int> > :: iterator it = plistChangedPlaces->begin(); it != plistChangedPlaces->end(); ++it )
                {
                    cout << "[ " << setw(6) << it->first.second+1 << " " << setw(6) << it->first.first+1 << " ]: " << it->second << endl;
                }
            }
        }
#endif
#endif
        
        if( this->fIterate == false )
        {
            break;
        }
        
        // Iterative search for better trees
        this->genosInput.SetGenotypeCalled(true);       // use called genotypes to build NJ tree
        set<pair<int,int> > setUncertainGenoPos;
        CollectUncertainGenosFrom( mapChangedPosHistory, setUncertainGenoPos );
        if( fVerbose )
        {
            cout << "Set of uncertain genotype positions (cell, site): size:  " << setUncertainGenoPos.size() << "   " << endl;
            /*
            for(set<pair<int,int> > :: iterator itt4 = setUncertainGenoPos.begin(); itt4 != setUncertainGenoPos.end(); ++itt4)
            {
                cout << "[" << itt4->first << ", " << itt4->second << "]  ";
            }
            cout << endl; */
        }
        this->genosInput.SetUncertainGenoPos(setUncertainGenoPos);
#if 0
        // collect current clades
        set<set<int> > setCladesCurr;
        for(int sa=0; sa<this->genosInput.GetNumSites(); ++sa)
        {
            set<int> smuts;
            bool fNoChange = true;
            //this->genosInput.GetMutRowsHapAtSite(sa, smuts);
            for(int c=0; c<this->genosInput.GetNumGenos(); ++c)
            {
                int g = this->genosInput.GetGenotypeAt(c, sa);
                double p = this->genosInput.GetGenotypeProbAllele0At(c, sa);
                if( g != 0 )
                {
                    if( p >= 0.5 )
                    {
                        fNoChange = false;
                        break;
                    }
                    smuts.insert(c);
                }
                else if( p < 0.5 )
                {
                    fNoChange = false;
                    break;
                }
            }
            if( smuts.size() > 1 && fNoChange == true )
            {
                setCladesCurr.insert(smuts);
            }
        }

        if( fVerbose )
        {
            cout << "Sets of pre-chosen clades for initial tree re-construction number: " << setCladesCurr.size() << endl;
            for(set<set<int> > :: iterator itt3 = setCladesCurr.begin(); itt3 != setCladesCurr.end(); ++itt3)
            {
                DumpIntSet( *itt3 );
            }
        }
#endif
        
        //this->strNJTree = this->genosInput.ConsNJTreeZeroRoot( &setCladesCurr );
        this->strNJTree = this->genosInput.ConsNJTreeZeroRoot(  );
        Init();
        //if( plistChangedPlaces != NULL )
        //{
        //    plistChangedPlaces->clear();
        //}
    }
#if 0
    // change the genotype matrix
    if( plistChangedPlaces != NULL )
    {
        cout << "Number of corrected genotypes (site, cell, new genotype) in base-1 (from maximal single position probability genotypes): " << plistChangedPlaces->size() << "\n";
        //for(set<pair<pair<int,int>,int> > :: iterator it = plistChangedPlaces->begin(); it != plistChangedPlaces->end(); ++it )
        //{
        //    cout << "[ " << setw(6) << it->first.second+1 << " " << setw(6) << it->first.first+1 << " ]: " << it->second << endl;
        //}
        
        //this->genosInput.ResetMaximalGenos();
        this->genosInput.ChangeGenosAtPositions(*plistChangedPlaces);
        this->genosInput.Dump();
        bool fc = this->genosInput.AreAllGenosCompatible();
        YW_ASSERT_INFO( fc == true, "The genotypes are not compatible!" );
    }
#endif
    
    if( fOutput )
    {
        // output the matrix
        //YW_ASSERT_INFO(plistChangedPlaces != NULL, "Cannot be null");
        //YW_ASSERT_INFO(pstrTreeNW != NULL, "Cannot be null2");
        if( plistChangedPlaces != NULL && pstrTreeNW != NULL )
        {
            if( fVerbose )
            {
                cout << "Number of corrected genotypes (site, cell, new genotype) in base-1 (from maximal single position probability genotypes): " << plistChangedPlaces->size() << "\n";
                //for(set<pair<pair<int,int>,int> > :: iterator it = plistChangedPlaces->begin(); it != plistChangedPlaces->end(); ++it )
                //{
                //    cout << "[ " << setw(6) << it->first.second+1 << " " << setw(6) << it->first.first+1 << " ]: " << it->second << endl;
                //}
            }
            
            ScistGenGenotypeMat *pMatRes = this->genosInput.Copy();
            pMatRes->ChangeGenosAtPositions(*plistChangedPlaces);
            // output genotypes
            cout << "Called genotypes output to file: " << pMatRes->GetImpGenosFileName() << endl;
            pMatRes->OutputImput();
            
            if( fOutputPPWithEdgeLabels )
            {
                ScistHaplotypeMat *pMatResHap = dynamic_cast<ScistHaplotypeMat *>(pMatRes);
                if(pMatResHap == NULL)
                {
                    cout << "** Right now, only output perfect phylogeny for binary genotypes\n";
                }
                else
                {
                    string strTreeEdgeLabel = ConsRootedPerfectPhylogenyFromMat(pMatResHap->GetHapMat(), true, true);
                    cout << "** Perfect phylogeny (with sites labeled on edges) from the imputed genotypes: " << strTreeEdgeLabel << endl;
                    
                    string strMutTree = ConsEdgeLabeTree(strTreeEdgeLabel);
                    string strMutTreeConv = ConvMutTreeStr(strMutTree);
                    cout << "^^ Mutation tree: " << strMutTreeConv << endl;
                    
                    // output mutation tree file
                    OutputMutationTree( this->strMutTreeFileName.c_str(), strMutTreeConv, this->fOutputLabel );
                }
            }
            
            delete pMatRes;
            
            cout << "**** Maximum log-likelihood: " << loglikeliBest << ", number of changed genotypes: " << plistChangedPlaces->size() << endl;
//#if 0
            cout << "Computed log-lielihood from changed genotypes: " << CalcChangedGenosProb(*plistChangedPlaces) << endl;
//#endif
            //cout << "Minimum cost: " << CalcMaxProbUpperBound() - loglikeliBest << endl;
            
            string strTreeOptOut = ConvCellTreeStr(*pstrTreeNW);
            cout << "Constructed single cell phylogeny: " << *pstrTreeNW << endl;
        }
    }
    return loglikeliBest;
}

double ScistPerfPhyMLE :: InferFastNNI( std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces, std::string *pstrTreeNW, map< pair<int,int>, set<int> > *pmapChangedPosHistory )
{
    const double DEF_THRES_PROB_INC = log(1.00001);
    // fast NNI search
//cout << "ScistPerfPhyMLE :: InferFastNNI\n";
//cout << "Genotypes: \n";
//this->genosInput.Dump();
    //
    string strTreeOpt;
    if( strNJTree.length() > 0 )
    {
        strTreeOpt = strNJTree;
    }
    else
    {
        set<ScistPerfPhyCluster> setClusAllGuide;
        this->treeGuide.GetAllClusters( setClusAllGuide );
//cout << "Number of clusters: " << setClusAllGuide.size() << endl;
        strTreeOpt = ConsTreeFromSetClusters( setClusAllGuide );
    }
//cout << "Initial tree: " << strTreeOpt << endl;
    
    //set<ScistPerfPhyCluster> setClusAllGuideUse;
    //GetClustersFromTree(strTreeOpt, setClusAllGuideUse);
    std::vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersOpt;
    //double loglikeliBest = ScoreSetClusters( setClusAllGuideUse, listChangedClustersOpt );
    double loglikeliBest = ScoreTree( strTreeOpt, listChangedClustersOpt );
//cout << "Init likelihood: " << loglikeliBest << endl;
    set<string> setTreeSearchedBefore;
    setTreeSearchedBefore.insert( strTreeOpt );

//#if 0
    // now search for neighborhood of the current tree to optimize the tree
    while(true)
    {
        if( fVerbose)
        {
            //cout << "Current best likelihood: " << loglikeliBest << endl;
            //cout << "Current best likelihood: " << loglikeliBest << ", opt tree: " << strTreeOpt << endl;
        }
        string strBestNNITree;
        double loglikeliStep = FastFindBestNNITreeFrom( strTreeOpt, strBestNNITree);
//cout << "loglikeliStep : " << loglikeliStep << ", found tree: " << strBestNNITree << endl;
        
        // re-check NNI tree
        vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersStepTmp;
//double loglikeliStep3 = ScoreTree( strBestNNITree, listChangedClustersStepTmp );
//cout << " re-computed likelihood: " << loglikeliStep3 << endl;
//YW_ASSERT_INFO(std::fabs(loglikeliStep3-loglikeliStep) < 0.01, "WRONG");
        
        if( loglikeliStep   > loglikeliBest + DEF_THRES_PROB_INC )
        {
            loglikeliBest = loglikeliStep;
            strTreeOpt = strBestNNITree;
            if( fVerbose )
            {
                //cout << "BETTER: updated loglikeBest: " << loglikeliBest << endl;
            }
        }
        else
        {
//cout << "NOT BETTER\n";
            // has reached local optimum
            break;
        }
    }
//#endif
    
    // get the changed genotypes
    vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersStep;
    //double loglikeliStep = ScoreSetClusters( setClus, listChangedClustersStep);
//#if 0
    double loglikeliStep2 = ScoreTree( strTreeOpt, listChangedClustersStep );
//cout << "current opt tree: " << strTreeOpt << ", re-computed likelihood: " << loglikeliStep2 << ", loglikelibest:" << loglikeliBest << endl;
    YW_ASSERT_INFO(std::fabs(loglikeliStep2-loglikeliBest) < 0.01, "WRONG1156");
//#endif
    listChangedClustersOpt = listChangedClustersStep;
    
    // output the final tree
    //this->genosInput.ResetMaximalGenos();
    std::set< std::pair<std::pair<int,int>, int> > listChangedPlaces;
    for(int site = 0; site<this->genosInput.GetNumSites(); ++site)
    {
        FindChangedGenos(site, listChangedClustersOpt[site], listChangedPlaces);
    }
    if( plistChangedPlaces != NULL )
    {
        *plistChangedPlaces = listChangedPlaces;
    }
    if( pstrTreeNW != NULL)
    {
        *pstrTreeNW = strTreeOpt;
    }
    // update history of changed positions if needed
    if(pmapChangedPosHistory != NULL )
    {
        for(set<pair<pair<int,int>,int> > :: iterator itt2 = listChangedPlaces.begin(); itt2 != listChangedPlaces.end(); ++itt2)
        {
            (*pmapChangedPosHistory)[ itt2->first ].insert( itt2->second );
        }
    }
    
#if 0
    if( fVerbose )
    {
        if(fOutput)
        {
            //cout << "Genotypes called by maximal single position probability\n";
            //const string strDesc = "Single-site maximal probability genotypes";
            //this->genosInput.OutputImput(&strDesc);
        }
        
        cout << "List of corrected genotypes (site, cell, new genotype) in base-1 (from maximal single position probability genotypes): \n";
        for(set<pair<pair<int,int>,int> > :: iterator it = listChangedPlaces.begin(); it != listChangedPlaces.end(); ++it )
        {
            cout << "[ " << setw(6) << it->first.second+1 << " " << setw(6) << it->first.first+1 << " ]: " << it->second << endl;
        }
    }
    if( fOutput )
    {
        // output the matrix
        ScistGenGenotypeMat *pMatRes = this->genosInput.Copy();
        pMatRes->ChangeGenosAtPositions(listChangedPlaces);
        if( fVerbose )
        {
            cout << "Called genotypes\n";
            pMatRes->OutputImput();
        }
        if( fOutputPPWithEdgeLabels )
        {
            ScistHaplotypeMat *pMatResHap = dynamic_cast<ScistHaplotypeMat *>(pMatRes);
            if(pMatResHap == NULL)
            {
                cout << "** Right now, only output perfect phylogeny for binary genotypes\n";
            }
            else
            {
                string strTreeEdgeLabel = ConsRootedPerfectPhylogenyFromMat(pMatResHap->GetHapMat(), true, true);
                cout << "** Perfect phylogeny (with sites labeled on edges) from the imputed genotypes: " << strTreeEdgeLabel << endl;
                
                string strMutTree = ConsEdgeLabeTree(strTreeEdgeLabel);
                string strMutTreeConv = ConvMutTreeStr(strMutTree);
                cout << "^^ Mutation tree: " << strMutTreeConv << endl;
                
                // output mutation tree file
                OutputMutationTree( this->strMutTreeFileName.c_str(), strMutTreeConv, this->fOutputLabel );
            }
        }
        
        delete pMatRes;
    }
    
    // change genotype
    for(set<pair<pair<int,int>,int> > :: iterator it = listChangedPlaces.begin(); it != listChangedPlaces.end(); ++it )
    {
        this->genosInput.SetGenotypeAt(it->first.first, it->first.second, it->second);
    }
    
    double res = loglikeliBest;
    
    if( fOutput )
    {
        cout << "**** Maximum log-likelihood: " << loglikeliBest << ", number of changed genotypes: " << listChangedPlaces.size() << endl;
        cout << "Computed log-lielihood from changed genotypes: " << CalcChangedGenosProb(listChangedPlaces) << endl;
        //cout << "Minimum cost: " << CalcMaxProbUpperBound() - loglikeliBest << endl;
        
        string strTreeOptOut = ConvCellTreeStr(strTreeOpt);
        cout << "Constructed single cell phylogeny: " << strTreeOptOut << endl;
    }
#endif
    double res = loglikeliBest;
    if( fOptBrLen )
    {
        string strTreeBrOpt;
        double loglikeliBestBr = OptBranchLens( strTreeOpt, strTreeBrOpt );
        res = loglikeliBestBr;
        if( fOutput )
        {
            cout << "**** Maximum log-likelihood (with branch length optimization): " << loglikeliBestBr << endl;
            string strTreeBrOptOut = ConvCellTreeStr(strTreeBrOpt);
            cout << "Single cell phylogeny with branch length: " << strTreeBrOptOut << endl;
        }
    }
    
    return res;
}

double ScistPerfPhyMLE :: InferFastSPR( const string &strTreeInit, std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces, std::string *pstrTreeNW, map<pair<int,int>, set<int> > *pmapChangedPosHistory )
{
    // perform fast SPR local search
//cout << "ScistPerfPhyMLE :: InferFastSPR\n";
    // test code: YW
    //set<ScistPerfPhyCluster> setClusAllGuide;
    //this->treeGuide.GetAllClusters( setClusAllGuide );
//cout << "Number of clusters: " << setClusAllGuide.size() << endl;
    //string strTreeOpt2 = ConsTreeFromSetClusters( setClusAllGuide );
//cout << "Initial tree: " << strTreeOpt2 << endl;
    
#if 0
    //
    string strTreeInit;
    if( strNJTree.length() > 0 )
    {
        strTreeInit = strNJTree;
    }
    else
    {
        set<ScistPerfPhyCluster> setClusAllGuide;
        this->treeGuide.GetAllClusters( setClusAllGuide );
//cout << "Number of clusters: " << setClusAllGuide.size() << endl;
        strTreeInit = ConsTreeFromSetClusters( setClusAllGuide );
    }
#endif
    
    ScistFastSPRLocalSearchLoop optSPRLocalSearch(this->genosInput, strTreeInit);
    optSPRLocalSearch.SetNumThreads(this->numThreads);
    optSPRLocalSearch.SetVerbose(this->fVerbose);
    //optSPRLocalSearch.SetMaxHeuSPRsNum(this->maxHeuSPRNum);
    //const double FRAC_SPR = 0.5;
    optSPRLocalSearch.SetHeuristicMode(!fExactSPR, fracSPRSrc, thresSPRDropStop);
    string strTreeOpt;
    double probMax = optSPRLocalSearch.FindOpt(strTreeOpt);
    
    //if( plistChangedPlaces != NULL )
    //{
    //    cout << "Warning: matrix changes not implemented yet\n";
    //}
    
    // get the changed genotypes
    vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersStep;
    //double loglikeliStep = ScoreSetClusters( setClus, listChangedClustersStep);
//#if 0
    double loglikeliStep2 = ScoreTree( strTreeOpt, listChangedClustersStep );
//cout << "current opt tree: " << strTreeOpt << ", re-computed likelihood: " << loglikeliStep2 << ", probMax:" << probMax << endl;
    YW_ASSERT_INFO(std::fabs(loglikeliStep2-probMax) < 0.01, "WRONG1145");
//#endif
    //listChangedClustersOpt = listChangedClustersStep;
    
    // output the final tree
    //this->genosInput.ResetMaximalGenos();
    std::set< std::pair<std::pair<int,int>, int> > listChangedPlaces;
    for(int site = 0; site<this->genosInput.GetNumSites(); ++site)
    {
        FindChangedGenos(site, listChangedClustersStep[site], listChangedPlaces);
    }
    if( plistChangedPlaces != NULL )
    {
        *plistChangedPlaces = listChangedPlaces;
    }
    
    if( pstrTreeNW != NULL)
    {
        *pstrTreeNW = strTreeOpt;
    }
    
    // update history of changed positions if needed
    if(pmapChangedPosHistory != NULL )
    {
        for(set<pair<pair<int,int>,int> > :: iterator itt2 = listChangedPlaces.begin(); itt2 != listChangedPlaces.end(); ++itt2)
        {
            (*pmapChangedPosHistory)[ itt2->first ].insert( itt2->second );
        }
    }
    
    //if( fOutput )
    //{
    //    cout << "**** Maximum log-likelihood: " << probMax << ", number of changed genotypes: unknown "  << endl;
    //    cout << "Computed log-lielihood from changed genotypes: not computed." << endl;
    //    //cout << "Minimum cost: " << CalcMaxProbUpperBound() - loglikeliBest << endl;
    //    string strTreeOptOut = ConvCellTreeStr(strTreeOpt);
    //    cout << "Constructed single cell phylogeny: " << strTreeOptOut << endl;
    //}
    
    return probMax;
}

double ScistPerfPhyMLE :: InferFastSPRInitTrees( const std::vector<std::string> &listInitTreesMLE, std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces, std::string *pstrTreeNW )
{
    // perform fast SPR local search
    //cout << "ScistPerfPhyMLE :: InferFastSPR\n";
    
    double probMaxAll = MAX_NEG_DOUBLE_VAL;
    string strOptTreeAll;
    for(unsigned int i=0; i<listInitTreesMLE.size(); ++i)
    {
        // init guide tree with this string
        //this->treeGuide.InitDecAll(listInitTreesMLE[i]);
        string strTreeOpt;
        std::set< std::pair<std::pair<int,int>, int> > listChanges;
        double probMax = InferFastSPR( listInitTreesMLE[i],  &listChanges, &strTreeOpt);

        if(probMax > probMaxAll)
        {
            probMaxAll = probMax;
            strOptTreeAll = strTreeOpt;
            
            if( plistChangedPlaces != NULL )
            {
                *plistChangedPlaces = listChanges;
            }
        }
    }
    if( pstrTreeNW != NULL)
    {
        *pstrTreeNW = strOptTreeAll;
    }
    return probMaxAll;
}

double ScistPerfPhyMLE :: FastFindBestNNITreeFrom(const string &strTreeFrom, string &strBestNNITree)
{
    // fast NNI search of local trees
    //this->treeGuide.InitDecAll(strTreeFrom);
    //this->treeGuide.Init(strTreeFrom);
    
    ScistFastNNILocalSearch optNNI(this->genosInput, strTreeFrom, this->numThreads);
//    ScistFastNNILocalSearch2 optNNI(this->genosInput, strTreeFrom, this->numThreads);
    double res = optNNI.CalcCurrMaxProb(strBestNNITree);
//cout << "FastFindBestNNITreeFrom: curMaxProb: " << res << endl;
    return res;
    //YW_ASSERT_INFO(false, "Not implemented.");
    //return MAX_NEG_DOUBLE_VAL;
}

double ScistPerfPhyMLE :: FastFindBestRerootTreeFrom(const string &strTreeFrom, string &strBestRerootTree)
{
    //return MAX_NEG_DOUBLE_VAL;
//#if 0
    // has reach local opitma for NNI, now try SPR
    ScistFastRerootLocalSearch treeSearchCur(this->genosInput, strTreeFrom, this->numThreads );
    
    int nodeST = -1, nodeSTSub = -1;
    double probMaxStepReroot = treeSearchCur.GetOptSubtreeReroot(nodeSTSub, nodeST);
    
    // get the opt reroot tree
    MarginalTree mTreeCurr;
    treeSearchCur.GetCurrTree(mTreeCurr);
    mTreeCurr.RerootSubtree(nodeST, nodeSTSub);
    strBestRerootTree = mTreeCurr.GetNewickNoBrLen();
    
    return probMaxStepReroot;
//#endif
}


double ScistPerfPhyMLE :: InferReroot( const string &strTreeInit, std::set< std::pair<std::pair<int,int>, int> > * plistChangedPlaces, std::string *pstrTreeNW )
{
    const double DEF_THRES_PROB_INC = log(1.00001);
    // fast NNI search
//cout << "ScistPerfPhyMLE :: InferFastNNI\n";
//cout << "Genotypes: \n";
//this->genosInput.Dump();
    //
    string strTreeOpt = strTreeInit;
//cout << "Reroot optimization: Initial tree: " << strTreeOpt << endl;
    
    //set<ScistPerfPhyCluster> setClusAllGuideUse;
    //GetClustersFromTree(strTreeOpt, setClusAllGuideUse);
    std::vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersOpt;
    //double loglikeliBest = ScoreSetClusters( setClusAllGuideUse, listChangedClustersOpt );
    double loglikeliBest = ScoreTree( strTreeOpt, listChangedClustersOpt );
//cout << "Init likelihood: " << loglikeliBest << endl;

//#if 0
    // now search for neighborhood of the current tree to optimize the tree
    while(true)
    {
        if( fVerbose)
        {
            //cout << "Current best likelihood: " << loglikeliBest << endl;
            //cout << "Current best likelihood: " << loglikeliBest << ", opt tree: " << strTreeOpt << endl;
        }
        string strBestNNITree;
        double loglikeliStep = FastFindBestRerootTreeFrom( strTreeOpt, strBestNNITree);
//cout << "loglikeliStep : " << loglikeliStep << ", found tree: " << strBestNNITree << endl;
        
        // re-check NNI tree
        vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersStepTmp;
#if 0
double loglikeliStep3 = ScoreTree( strBestNNITree, listChangedClustersStepTmp );
cout << " re-computed likelihood: " << loglikeliStep3 << endl;
YW_ASSERT_INFO(std::fabs(loglikeliStep3-loglikeliStep) < 0.01, "WRONG");
#endif
        if( loglikeliStep   > loglikeliBest + DEF_THRES_PROB_INC )
        {
            loglikeliBest = loglikeliStep;
            strTreeOpt = strBestNNITree;
            if( fVerbose )
            {
                //cout << "BETTER: updated loglikeBest: " << loglikeliBest << endl;
            }
        }
        else
        {
//cout << "NOT BETTER\n";
            // has reached local optimum
            break;
        }
    }
//#endif
    
    // get the changed genotypes
    vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > listChangedClustersStep;
    //double loglikeliStep = ScoreSetClusters( setClus, listChangedClustersStep);
//#if 0
    double loglikeliStep2 = ScoreTree( strTreeOpt, listChangedClustersStep );
//cout << "current opt tree: " << strTreeOpt << ", re-computed likelihood: " << loglikeliStep2 << ", loglikelibest:" << loglikeliBest << endl;
    YW_ASSERT_INFO(std::fabs(loglikeliStep2-loglikeliBest) < 0.01, "WRONG1156");
//#endif
    listChangedClustersOpt = listChangedClustersStep;
    
    // output the final tree
    //this->genosInput.ResetMaximalGenos();
    std::set< std::pair<std::pair<int,int>, int> > listChangedPlaces;
    for(int site = 0; site<this->genosInput.GetNumSites(); ++site)
    {
        FindChangedGenos(site, listChangedClustersOpt[site], listChangedPlaces);
    }
    if( plistChangedPlaces != NULL )
    {
        *plistChangedPlaces = listChangedPlaces;
    }
    if( pstrTreeNW != NULL)
    {
        *pstrTreeNW = strTreeOpt;
    }
    
    return loglikeliBest;
}

double ScistPerfPhyMLE :: OptBranchLens( const std::string &strTree, std::string &strTreeBrOpt )
{
    //
    MarginalTree treeBrOpt;
    ReadinMarginalTreesNewickWLenString( strTree, this->genosInput.GetNumHaps(), treeBrOpt );
    ScistFullPerfPhyMLE sfpp(this->genosInput);
    double res = sfpp.OptBranchLens(treeBrOpt);
    strTreeBrOpt = treeBrOpt.GetNewickSorted(true);
    return res;
}

void ScistPerfPhyMLE :: Init()
{
//cout << "Init: ScistPerfPhyMLE\n ";
    //
    // get all clusters
    listClusMutsInputHetero.clear();
    listClusMutsInputHomo.clear();
    for(int s=0; s<genosInput.GetNumSites(); ++s)
    {
        set<int> muts;
        genosInput.GetRowsWithGenoAtSite(s, 1, muts);
        ScistPerfPhyCluster clus(muts);
        listClusMutsInputHetero.push_back(clus);
        
        set<int> muts2;
        genosInput.GetRowsWithGenoAtSite(s, 2, muts2);
        ScistPerfPhyCluster clus2(muts2);
        listClusMutsInputHomo.push_back(clus2);
    }
    
    // YW: don't get multiplicity. 01/13/23
    //this->genosInput.GetColMultiplicityMap( listInputColMulti );
//cout << "(2)Init: ScistPerfPhyMLE\n ";
    // construct NJ tree as the initial tree
    if( strNJTree.size() == 0 )
    {
        string strNJ = this->genosInput.ConsNJTreeZeroRoot();
        strNJTree = strNJ;
    }
//cout << "Guide tree: " << strNJ << endl;
    //string strNJ = this->genosInput.ConsNJTree();
//cout << "Zero-rooted initial tree: " << strNJ << endl;
//cout << "Genotype input: \n";
//this->genosInput.Dump();
    //
    this->treeGuide.Init(strNJTree);
//cout << "Init genotype prob...\n";
    // set the prior score to be zero
    listSitePriorScore.clear();
    for(int i=0; i<this->genosInput.GetNumSites(); ++i)
    {
        double logprobInit = 0.0;
        for(int h=0; h<this->genosInput.GetNumHaps(); ++h)
        {
            double p = this->genosInput.GetGenotypeProbAllele0At(h, i);
            logprobInit += log(p);
        }
        listSitePriorScore.push_back(logprobInit);
    }
//cout << "Done Init.\n";
}

std::string ScistPerfPhyMLE :: ConsTreeFromSetClusters( const std::set<ScistPerfPhyCluster> &setClusters ) const
{
//cout << "All the clusters: \n";
//for(set<ScistPerfPhyCluster> :: const_iterator it = setClusters.begin(); it != setClusters.end(); ++it)
//{
//it->Dump();
//}
    //
    // now construct tree
    ScistInfPerfPhyUtils treeBuild;
    map<int, ScistPerfPhyCluster> mapPickedClus;
    int s = 0;
    for(set<ScistPerfPhyCluster>:: iterator it = setClusters.begin(); it != setClusters.end(); ++it)
    {
        mapPickedClus[s] = *it;
        ++s;
    }
    string strTree = treeBuild.ConsTreeWCombDistClus( this->genosInput, mapPickedClus, false );
    return strTree;
}

void ScistPerfPhyMLE :: GetNgbrTreesFrom(int numHaps, const std::string &strTree, std::set<std::string> &setNgbrTrees )
{
//cout << "GetNgbrTreesFrom: numHaps: " << numHaps << ", tree: " << strTree << endl;
    //
    setNgbrTrees.clear();
    MarginalTree treeCurr;
    ReadinMarginalTreesNewickWLenString( strTree, numHaps, treeCurr );
    vector<MarginalTree> listNgbrTrees;
    FindOneNNIMTreesFrom( treeCurr, listNgbrTrees );
    for(int i=0; i<(int)listNgbrTrees.size(); ++i)
    {
        string strTree = listNgbrTrees[i].GetNewickSorted(false);
        setNgbrTrees.insert(strTree);
    }
}

void ScistPerfPhyMLE :: GetNgbrTreesFromSPR(int numHaps, const std::string &strTree, std::set<std::string> &setNgbrTrees )
{
    //
    setNgbrTrees.clear();
    MarginalTree treeCurr;
    ReadinMarginalTreesNewickWLenString( strTree, numHaps, treeCurr );
    string strSelf = treeCurr.GetNewickSorted(false);
//cout << "strTree: " << strTree << ", strSelf: " << strSelf << endl;
    
    // map to consecutive order as required by RBT
    vector<int> listLeafLblsOld;
    //treeCurr.MapLeafLblConsecutiveOrder( listLeafLblsOld );
    treeCurr.GetLabelList(listLeafLblsOld);
//cout << "Mapped leaves: ";
//DumpIntVec(listLeafLblsOld);
//cout << "Changed tree: " << treeCurr.GetNewick() << endl;
    
    // use RBT utility
    vector<int> listLbls;
    treeCurr.GetLabelList(listLbls);
//cout << "listLbss: ";
//DumpIntVec(listLbls);
    vector<int> parPosList;
    treeCurr.GetParPosInfo( parPosList );
//cout << "parPosList: ";
//DumpIntVec(parPosList);
    vector<double> listEdgeDistOut;
    treeCurr.GetTreeEdgeLen( listEdgeDistOut );
    RBT treeCurrRBT( numHaps, listLbls, parPosList, listEdgeDistOut);
    vector<RBT *> ngbrTrees;
    treeCurrRBT.FindSPRDistOneNgbrs( ngbrTrees );
    
//cout << "GetNgbrTreesFromSPR: init tree: " << strTree << endl;
    for(int i=0; i<(int)ngbrTrees.size(); ++i)
    {
        string strNW = ngbrTrees[i]->GetNewick();
        string strNWBack = RemapLeafLbls(numHaps, strNW, listLeafLblsOld );
        setNgbrTrees.insert(strNWBack);
//cout << "strNW: " << strNW  << ", SPR tree: " << strNWBack << endl;
    }
    // remove self
    setNgbrTrees.erase(strSelf);
    
    for(int i=0; i<(int)ngbrTrees.size(); ++i)
    {
        delete ngbrTrees[i];
    }
}

std::string ScistPerfPhyMLE :: RemapLeafLbls(int numHaps, const std::string &strTree0Based, const vector<int> &listLblsOld )
{
    //
    MarginalTree treeCurr;
    ReadinMarginalTreesNewickWLenString( strTree0Based, numHaps, treeCurr );
    map<int,int> mapLblsBack;
    for(int i=0;i<(int)listLblsOld.size(); ++i)
    {
        mapLblsBack[i] = listLblsOld[i];
    }
    treeCurr.RemapLeafLabels( mapLblsBack );
    return treeCurr.GetNewickSorted(false);
}

std::string ScistPerfPhyMLE :: RemapLeafLbls(int numHaps, const std::string &strTree, const std::map<int,int> &mapLabels )
{
    //
    MarginalTree treeCurr;
    ReadinMarginalTreesNewickWLenString( strTree, numHaps, treeCurr );
    treeCurr.RemapLeafLabels( mapLabels );
    return treeCurr.GetNewickSorted(false);
}

std::string ScistPerfPhyMLE :: ConvCellTreeStr(const std::string &strTree) const
{
    //
    if( this->listCellNames.size() == 0 )
    {
        // no conversion if no cell names specified
        return strTree;
    }
    
    TaxaMapper taxaMapper;
    for(int i=0; i<(int)listCellNames.size(); ++i)
    {
        taxaMapper.AddTaxaStringWithId( i+1, listCellNames[i] );
    }
    //
    return taxaMapper.ConvIdStringWithOrigTaxa( strTree );
}

std::string ScistPerfPhyMLE :: ConvMutTreeStr(const std::string &strTree) const
{
    //
    if( this->listSiteNames.size() == 0 )
    {
        // no conversion if no cell names specified
        return strTree;
    }
    
    TaxaMapper taxaMapper;
    for(int i=0; i<(int)listSiteNames.size(); ++i)
    {
        taxaMapper.AddTaxaStringWithId( i+1, listSiteNames[i] );
    }
    //
    return taxaMapper.ConvIdStringWithOrigTaxa( strTree );
}

void ScistPerfPhyMLE :: FindChangedGenos( int siteToAdd, const pair<ScistPerfPhyCluster,ScistPerfPhyCluster> &clusToAdd, set< pair<pair<int,int>, int> > & listChangedPlaces ) const
{
    // find list of positions where the genos are changed
    ScistPerfPhyCluster clusInt, clusThisOnly, clusRHSOnly;
    clusToAdd.first.IntersectWith( listClusMutsInputHetero[siteToAdd], clusInt, clusThisOnly, clusRHSOnly );
    ScistPerfPhyCluster clusInt2, clusThisOnly2, clusRHSOnly2;
    clusToAdd.second.IntersectWith( listClusMutsInputHomo[siteToAdd], clusInt2, clusThisOnly2, clusRHSOnly2 );
    // get changed 0
    set<int> setss;
    PopulateSetWithInterval(setss, 0, this->genosInput.GetNumHaps()-1);
    set<int> rows0Orig;
    this->genosInput.GetRowsWithGenoAtSite(siteToAdd, 0, rows0Orig);
    SubtractSets(setss, rows0Orig);
    ScistPerfPhyCluster clus0(setss);
    clus0.SubtractFrom( clusToAdd.first );
    clus0.SubtractFrom( clusToAdd.second );
    
    ScistPerfPhyClusterItor itor0(clus0);
    itor0.First();
    while(itor0.IsDone() == false )
    {
        int sc = itor0.GetCurrentSC();
        pair<int,int> pp(sc, siteToAdd);
        pair<pair<int,int>, int>  pp0(pp, 0);
        listChangedPlaces.insert(pp0);
        itor0.Next();
    }
    
    // This only: new mutants
    ScistPerfPhyClusterItor itor1(clusThisOnly);
    itor1.First();
    while(itor1.IsDone() == false )
    {
        int sc = itor1.GetCurrentSC();
        pair<int,int> pp(sc, siteToAdd);
        pair<pair<int,int>, int>  pp1(pp, 1);
        listChangedPlaces.insert(pp1);
        itor1.Next();
    }
    // RHS only: new wildtype
    ScistPerfPhyClusterItor itor2(clusThisOnly2);
    itor2.First();
    while(itor2.IsDone() == false )
    {
        int sc = itor2.GetCurrentSC();
        pair<int,int> pp(sc, siteToAdd);
        pair<pair<int,int>, int>  pp2(pp,2);
        listChangedPlaces.insert(pp2);
        itor2.Next();
    }
}

double ScistPerfPhyMLE :: ScoreTree(const string &strTree, std::vector<std::pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > &listChangedCluster) const
{
    // run multi-thread version if specified
    if( this->numThreads > 1 )
    {
        return ScoreTreeMT(strTree, listChangedCluster);
    }
    
//cout << "ScoreTree: tree: " << strTree << endl;
    // score the current tree
    MarginalTree treeToScore;
    ReadinMarginalTreesNewickWLenString(strTree, this->genosInput.GetNumHaps(), treeToScore);
//cout << "Score tree: " << treeToScore.GetNewick() << endl;
    set<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > setClusDone;
    //map<pair<ScistPerfPhyCluster,ScistPerfPhyCluster>, pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > mapChangedClus;
    double res = 0.0;
    ScistPerfPhyProbOnTree probTree( this->genosInput, treeToScore );
    
    // YW: do we really need the map as above? TBD
    for(int site=0; site<genosInput.GetNumSites(); ++site)
    {
//cout << "ScoreTree: site " << site << endl;
//cout << "Heterozygote clus: ";
//listClusMutsInputHetero[site].Dump();
//cout << "Homozygous clus: ";
//listClusMutsInputHomo[site].Dump();
        //pair<ScistPerfPhyCluster,ScistPerfPhyCluster> pp0( listClusMutsInputHetero[site], listClusMutsInputHomo[site] );
        //if( setClusDone.find(pp0) != setClusDone.end() )
        //{
        //    listChangedCluster.push_back( mapChangedClus[ pp0 ] );
        //    continue;
        //}
        //int multi = this->listInputColMulti[site];
        pair<ScistPerfPhyCluster,ScistPerfPhyCluster> clusChanged;
        double loglikeliSite = ScoreTreeWithSite(probTree, treeToScore, site, clusChanged.first, clusChanged.second);
        //mapChangedClus[ pp0 ] = clusChanged;
        listChangedCluster.push_back(clusChanged);
        res += loglikeliSite;
        //setClusDone.insert( pp0 );
//cout << "ScoreTree: site: " << site << ", loglikesite: " << loglikeliSite << endl;
//cout << "site prob: " << loglikeliSite << ": clusChanged: ";
//clusChanged.first.Dump();
//cout << "  and ";
//clusChanged.second.Dump();
    }
    
    return res;
}

// helper function for threading
// utility for multithreading
static void UtilFastTreeProbMT( ScistPerfPhyProbOnTree *pprobTree, int tindex, int siteBlkStart, int siteBlkEnd, vector<double> *plistBlkProbs, vector<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > *plistChangedClusterMT )
{
    double res = 0.0;
    for(int site=siteBlkStart; site<=siteBlkEnd; ++site)
    {
        pair<ScistPerfPhyCluster,ScistPerfPhyCluster> clusChanged;
        //double loglikeliSite = perfMLE.ScoreTreeWithSite(probTree, treeToScore, site, clusChanged.first, clusChanged.second);
        double loglikeliSite =  pprobTree->CalcProbMaxForSite( site, clusChanged.first, clusChanged.second );
        //mapChangedClus[ pp0 ] = clusChanged;
        // YW: hack: don't reeturn changed genos for testing
        (*plistChangedClusterMT).push_back(clusChanged);
        res += loglikeliSite;
    }
    (*plistBlkProbs)[tindex] = res;
}


double ScistPerfPhyMLE :: ScoreTreeMT(const string &strTree, std::vector<std::pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > &listChangedCluster) const
{
    listChangedCluster.clear();
//cout << "ScoreTreeMT: tree: " << strTree << endl;
    // score the current tree
    MarginalTree treeToScore;
    ReadinMarginalTreesNewickWLenString(strTree, this->genosInput.GetNumHaps(), treeToScore);
//cout << "Score tree: " << treeToScore.GetNewick() << endl;
    set<pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > setClusDone;
    //map<pair<ScistPerfPhyCluster,ScistPerfPhyCluster>, pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > mapChangedClus;
    double res = 0.0;
    ScistPerfPhyProbOnTree probTree( this->genosInput, treeToScore );
    
    // break sites into pieces
    int numThreadsUse = this->numThreads;
    if( genosInput.GetNumSites() < numThreadsUse )
    {
        numThreadsUse = genosInput.GetNumSites();
    }
    int szBlock = genosInput.GetNumSites()/numThreadsUse;
    int site = 0;
    vector<double> listResultsMT(numThreadsUse);
    vector<vector< pair<ScistPerfPhyCluster,ScistPerfPhyCluster> > > listClusChangedMT(numThreadsUse);
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        int siteBlkEnd = site + szBlock-1;
        if( siteBlkEnd >= genosInput.GetNumSites()  || t == numThreadsUse-1)
        {
            siteBlkEnd = genosInput.GetNumSites()-1;
        }
        YW_ASSERT_INFO(site <= siteBlkEnd, "WRONG111");

        // start thread
        thread *pthr = new thread(UtilFastTreeProbMT, &probTree, t, site, siteBlkEnd, &listResultsMT, &listClusChangedMT[t] );
        listPtrThreads.push_back(pthr);
        
        site += szBlock;
    }
    
    // wait to join
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      listPtrThreads[j]->join();
    }
    for(unsigned int j=0; j<listPtrThreads.size(); ++j)
    {
      delete listPtrThreads[j];
    }
    listPtrThreads.clear();
    
    // sum the prob
    res = 0.0;
    for(unsigned int i=0; i<listResultsMT.size(); ++i)
    {
        res += listResultsMT[i];
    }
    for(unsigned int i=0; i<listClusChangedMT.size(); ++i)
    {
        for(unsigned int j=0; j<listClusChangedMT[i].size(); ++j)
        {
            listChangedCluster.push_back( listClusChangedMT[i][j]);
        }
    }
//cout << "Done with ScoreTreeMT: prob: " << res << endl;
    return res;
}

double ScistPerfPhyMLE :: ScoreTreeWithSite(ScistPerfPhyProbOnTree &probTree, MarginalTree &tree, int site, ScistPerfPhyCluster &clusChanged1, ScistPerfPhyCluster &clusChanged2) const
{
//cout << "site: " << site << ", tree: " << tree.GetNewickSorted(false) << endl;
    return probTree.CalcProbMaxForSite( site, clusChanged1, clusChanged2 );
}

double ScistPerfPhyMLE :: CalcMaxProbUpperBound() const
{
    //
    double res = 0.0;
    for(int s=0; s<this->genosInput.GetNumSites(); ++s)
    {
        for(int h=0; h<this->genosInput.GetNumHaps(); ++h)
        {
            double p0 = this->genosInput.GetGenotypeProbAllele0At(h, s);
            double p1 = 1 - p0;
            if( p0 >= p1)
            {
                res += log(p0);
            }
            else
            {
                res += log(p1);
            }
        }
    }
    return res;
}

double ScistPerfPhyMLE :: CalcChangedGenosProb(const std::set< std::pair<std::pair<int,int>, int> > & listChangedPlaces) const
{
    //
    double res = 0.0;
    map<pair<int,int>, int> mapChangedPlaces;
    for( std::set< std::pair<std::pair<int,int>, int> > :: const_iterator it = listChangedPlaces.begin(); it != listChangedPlaces.end(); ++it)
    {
        mapChangedPlaces[it->first] = it->second;
    }
    
    for(int s=0; s<this->genosInput.GetNumSites(); ++s)
    {
        for(int h=0; h<this->genosInput.GetNumHaps(); ++h)
        {
            pair<int,int> pp(h, s);
            int allele = this->genosInput.GetGenotypeAt(h, s);
            std::map< std::pair<int,int>, int > :: const_iterator it = mapChangedPlaces.find(pp);
            if( it != mapChangedPlaces.end() )
            {
                int alleleAlt = it->second;
                YW_ASSERT_INFO(allele == alleleAlt, "Wrong");
                allele = alleleAlt;
            }
            
            double p0 = this->genosInput.GetGenotypeProbAllele0At(h, s);
            double p1 = 1 - p0;
            if( allele == 0 )
            {
                res += log(p0);
            }
            else
            {
                res += log(p1);
            }
        }
    }
    
    return res;
}

void ScistPerfPhyMLE :: CollectUncertainGenosFrom( const std::map<std::pair<int,int>, std::set<int> > &mapChangedPlacesHist, std::set<std::pair<int,int> > &setUncertainPos ) const
{
    //
    setUncertainPos.clear();
    for( std::map<std::pair<int,int>, std::set<int> > :: const_iterator it = mapChangedPlacesHist.begin(); it != mapChangedPlacesHist.end(); ++it )
    {
        if( it->second.size() >= 1 )
        {
            // this position has changed alleles, so mark it uncertain
            setUncertainPos.insert( it->first );
        }
    }
}

// *************************************************************************************
// Tree probability


ScistPerfPhyProbOnTree :: ScistPerfPhyProbOnTree( ScistGenGenotypeMat &genos, MarginalTree &mtreeIn ) : genosInput(genos), mtree(mtreeIn)
{
    // set the prior score to be zero
    listSitePriorScore.clear();
    for(int i=0; i<this->genosInput.GetNumSites(); ++i)
    {
        double logprobInit = 0.0;
        for(int h=0; h<this->genosInput.GetNumHaps(); ++h)
        {
            double p = this->genosInput.GetGenotypeProbAt(h, i, 0);
            logprobInit += log(p);
        }
        listSitePriorScore.push_back(logprobInit);
    }
    
    Init();
}

void ScistPerfPhyProbOnTree :: Init()
{
    //
    ScistTernaryMat *pGenoMat = dynamic_cast<ScistTernaryMat *>(&this->genosInput);
    if( pGenoMat == NULL )
    {
        return;     // only work with genotype data
    }
    this->genosInputHap.SetSize( this->genosInput.GetNumHaps(), this->genosInput.GetNumSites()*2 );
    for(int h=0; h<this->genosInput.GetNumHaps(); ++h)
    {
        for(int s=0; s<this->genosInput.GetNumSites(); ++s)
        {
            double p0 = pGenoMat->GetGenotypeProbAt( h, s, 0 );
            double p1 = pGenoMat->GetGenotypeProbAt( h, s, 1 );
            double p2 = pGenoMat->GetGenotypeProbAt( h, s, 2 );
            double p12 = p1+p2;
            double p01 = p0+p1;
            int allele0 = 0;
            if( p0 < p12 )
            {
                allele0 = 1;
            }
            this->genosInputHap.SetGenotypeAt(h, 2*s, allele0);
            this->genosInputHap.SetGenotypeProbAt(h, 2*s, p0);
            int allele1 = 0;
            if( p01 < p2 )
            {
                allele1 = 1;
            }
            this->genosInputHap.SetGenotypeAt(h, 2*s+1, allele1);
            this->genosInputHap.SetGenotypeProbAt(h, 2*s+1, p01);
        }
    }
}

double ScistPerfPhyProbOnTree :: CalcProbMaxForSite(int site, ScistPerfPhyCluster &clusChangedMut, ScistPerfPhyCluster &clusChangedHomoMut) const
{
    ScistHaplotypeMat *pHapMat = dynamic_cast<ScistHaplotypeMat *>(&this->genosInput);
    
    if( pHapMat != NULL )
    {
        clusChangedHomoMut.Clear();
        return CalcProbMaxForSiteHap(site, clusChangedMut);
    }
    else    // right now, must be of genotype
    {
        return CalcProbMaxForSiteGeno(site, clusChangedMut, clusChangedHomoMut);
    }
}

void ScistPerfPhyProbOnTree :: PrepareProbMaxQ( )
{
    // this is for repetitively evaluting different trees for the same matrix input
    // assume: leaf won't change during rep
    this->listProbsSiteRepetitive.clear();
    this->listProbsSiteRepetitive.resize( this->genosInput.GetNumSites() );
    this->listSiteMaxProbLeaf.clear();
    this->listSiteMaxProbLeaf.resize( this->genosInput.GetNumSites() );
    for(int site=0; site<(int)listProbsSiteRepetitive.size(); ++site)
    {
        this->listProbsSiteRepetitive[site].resize( this->mtree.GetTotNodesNum() );
        this->listSiteMaxProbLeaf[site] = MAX_NEG_DOUBLE_VAL;
        for(int node=0; node<this->genosInput.GetNumHaps(); ++node)
        {
            // a single leaf in the split
            int lvlbl = mtree.GetLabel(node)-1;
            //cout << "Leaf: " << lvlbl << endl;
            double p0 = this->genosInput.GetGenotypeProbAllele0At(lvlbl, site);
            double logpStep = log((1-p0)/p0);
            this->listProbsSiteRepetitive[site][node] = logpStep;
            //cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
            if( this->listSiteMaxProbLeaf[site] < logpStep )
            {
                this->listSiteMaxProbLeaf[site] = logpStep;
            }
        }
    }
}

double ScistPerfPhyProbOnTree :: CalcProbMaxForSiteForTree(int site, MarginalTree &treeEval)
{
    //cout << "ScoreTreeWithSite: tree: " << tree.GetNewick() << ", site: " << site << endl;
    // score the site wrt the tree (i.e. find the best split of the tree for this site)
    //double res = -1.0*HAP_MAX_INT;
    double res = listSiteMaxProbLeaf[site];
    // do a bottom up
    //vector<double> listNodeSplitProb;
    // init to be bad
    //for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    //{
    //    listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
    //}
    
//cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) << endl;
//mtree.Dump();
    
    for(int node=this->genosInput.GetNumHaps(); node<treeEval.GetTotNodesNum(); ++node)
    {
        double logpStep;

        // get the two children and add them up
        int childLeft = treeEval.GetLeftDescendant(node);
        int childRight = treeEval.GetRightDescendant(node);
//cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
        //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
        logpStep = this->listProbsSiteRepetitive[site][childLeft] + this->listProbsSiteRepetitive[site][childRight];
//cout << "log prob: " << logpStep << " for node: " << node << endl;
        //listNodeSplitProb[node] = logpStep;
        //listNodeSplitProb.push_back( logpStep );
        this->listProbsSiteRepetitive[site][node] = logpStep;
        if( logpStep > res )
        {
            //cout << "Better at node: " << node << endl;
            res = logpStep;
        }
    }
    
    // if nothing is good, just take all-0
    if( res < 0.0 )
    {
        //
        res = 0;
    }
//cout << "Max prob at this site: " << res << ", added prior:" << res + this->listSitePriorScore[site] << " at site " << nodeOpt << endl;
    //cout << "clust changed: ";
    //clusChanged.Dump();
    return res + this->listSitePriorScore[site];
}

double ScistPerfPhyProbOnTree :: CalcProbMaxForSiteHap(int site, ScistPerfPhyCluster &clusChanged) const
{
    //cout << "ScoreTreeWithSite: tree: " << tree.GetNewick() << ", site: " << site << endl;
    // score the site wrt the tree (i.e. find the best split of the tree for this site)
    double res = -1.0*HAP_MAX_INT;
    // do a bottom up
    vector<double> listNodeSplitProb;
    // init to be bad
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        listNodeSplitProb.push_back( -1.0*HAP_MAX_INT );
    }
    
//cout << "CalcProbMaxForSiteHap: mtree: " << mtree.GetNewickSorted(false) << endl;
//mtree.Dump();
    
    int nodeOpt = -1;
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        //cout << "node " << node << endl;
        //if( node == mtree.GetRoot() )
        //{
            //continue;
        //}
        double logpStep;
        if( mtree.IsLeaf(node) )
        {
            // a single leaf in the split
            int lvlbl = mtree.GetLabel(node)-1;
            //cout << "Leaf: " << lvlbl << endl;
            double p0 = this->genosInput.GetGenotypeProbAllele0At(lvlbl, site);
            if( p0 < YW_VERY_SMALL_FRACTION)
            {
                p0 = YW_VERY_SMALL_FRACTION;
            }
            else if( p0 > 1.0-YW_VERY_SMALL_FRACTION)
            {
                p0 = 1.0-YW_VERY_SMALL_FRACTION;
            }
            logpStep = log((1-p0)/p0);
//cout << "Set leaf " << node << " log prob to: " << logpStep << ", p0=" << p0 << endl;
        }
        else
        {
            // get the two children and add them up
            int childLeft = mtree.GetLeftDescendant(node);
            int childRight = mtree.GetRightDescendant(node);
//cout << "node: " << node << ", childLeft: " << childLeft << ", childRight: " << childRight << endl;
            //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
            
            YW_ASSERT_INFO( listNodeSplitProb[childLeft] > -1.0*HAP_MAX_INT, "Bad left" );
            YW_ASSERT_INFO( listNodeSplitProb[childRight] > -1.0*HAP_MAX_INT, "Bad right1" );
            logpStep = listNodeSplitProb[childLeft] + listNodeSplitProb[childRight];
        }
//cout << "log prob: " << logpStep << " for node: " << node << endl;
        listNodeSplitProb[node] = logpStep;
        if( logpStep > res )
        {
            //cout << "Better at node: " << node << endl;
            res = logpStep;
            nodeOpt = node;
        }
    }
    
    set<int> nodeOptSplitLbls;
    
    // if nothing is good, just take all-0
    if( res < 0.0 )
    {
        //
        res = 0;
        nodeOpt = -1;
    }
    else
    {
        YW_ASSERT_INFO(nodeOpt >= 0, "Node not found");
        set<int> nodeOptSplit;
        mtree.GetLeavesUnder(nodeOpt, nodeOptSplit);
        mtree.GetlabelsFor( nodeOptSplit, nodeOptSplitLbls );
        DecAllNumInSet(nodeOptSplitLbls);
    }
    ScistPerfPhyCluster clus(nodeOptSplitLbls);
    clusChanged = clus;
//cout << "Max prob at this site: " << res << ", added prior:" << res + this->listSitePriorScore[site] << " at site " << nodeOpt << endl;
    //cout << "clust changed: ";
    //clusChanged.Dump();
    return res + this->listSitePriorScore[site];
}

double ScistPerfPhyProbOnTree :: CalcProbMaxForSiteGeno(int site, ScistPerfPhyCluster &clusChangedHetero, ScistPerfPhyCluster &clusChangedHomo) const
{
    //
    set<int> setSC0, setSC1, setSC2;
    this->genosInput.GetRowsWithGenoAtSite(site, 0, setSC0);
    this->genosInput.GetRowsWithGenoAtSite(site, 1, setSC1);
    this->genosInput.GetRowsWithGenoAtSite(site, 2, setSC2);
    
    // first accumulate for each node, the sum of diff p1/p0
    vector<double> vecSumDiffP10, vecSumDiffP21;
    vector<double> vecMaxSumDiff21;
    vector<int> vecMaxSumDiff21Node;
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        double logpStep, logpStep2;
        if( mtree.IsLeaf(node) )
        {
            // a single leaf in the split
            int lvlbl = mtree.GetLabel(node)-1;
            //cout << "Leaf: " << lvlbl << endl;
            double p0 = this->genosInput.GetGenotypeProbAt(lvlbl, site, 0);
            double p1 = this->genosInput.GetGenotypeProbAt(lvlbl, site, 1);
            double p2 = this->genosInput.GetGenotypeProbAt(lvlbl, site, 2);
            logpStep = log(p1/p0);
            logpStep2 = log(p2/p1);
            vecMaxSumDiff21.push_back( logpStep2 );
            vecMaxSumDiff21Node.push_back(node);
        }
        else
        {
            // get the two children and add them up
            int childLeft = mtree.GetLeftDescendant(node);
            int childRight = mtree.GetRightDescendant(node);
            //cout << "childLeft: " << childLeft << ", right: " << childRight << endl;
            
            YW_ASSERT_INFO( vecSumDiffP10[childLeft] > -1.0*HAP_MAX_INT, "Bad left (geno)" );
            YW_ASSERT_INFO( vecSumDiffP10[childRight] > -1.0*HAP_MAX_INT, "Bad right2" );
            logpStep = vecSumDiffP10[childLeft] + vecSumDiffP10[childRight];
            logpStep2 = vecSumDiffP21[childLeft] + vecSumDiffP21[childRight];
            
            double maxSumLogp21 = logpStep2;
            int nodeMax = node;
            if( vecSumDiffP21[childLeft] > maxSumLogp21)
            {
                maxSumLogp21 = vecSumDiffP21[ childLeft ];
                nodeMax = vecMaxSumDiff21Node[childLeft];
            }
            if( vecSumDiffP21[childRight] > maxSumLogp21 )
            {
                maxSumLogp21 = vecSumDiffP21[ childRight ];
                nodeMax = vecMaxSumDiff21Node[ childRight ];
            }
            vecMaxSumDiff21.push_back( maxSumLogp21 );
            vecMaxSumDiff21Node.push_back( nodeMax );
        }
//cout << "log prob: " << logpStep << endl;
        vecSumDiffP10.push_back( logpStep );
        vecSumDiffP21.push_back( logpStep2 );
    }
    
    // do another scan to find the best
    double res = -1.0*HAP_MAX_INT;
    int node1 = -1, node2 = -1;
    for(int node=0; node<mtree.GetTotNodesNum(); ++node)
    {
        double p2Part = 0.0;
        double node2MaxUse = -1;
        if( vecMaxSumDiff21[node] > 0.0 )
        {
            p2Part = vecMaxSumDiff21[node];
            node2MaxUse = vecMaxSumDiff21Node[node];
        }
        if( vecSumDiffP10[node] + p2Part > res )
        {
            res = vecSumDiffP10[node] + p2Part;
            
            node1 = node;
            node2 = node2MaxUse;
        }
    }
    
    // figure out the genos
    set<int> dummy;
    ScistPerfPhyCluster clusDummy(dummy);
    if( res < 0.0 )
    {
        //
        clusChangedHetero = clusDummy;
        clusChangedHomo = clusDummy;
    }
    else
    {
        YW_ASSERT_INFO(node1 >= 0, "Wrong");
        set<int> nodeOptSplit, nodeOptSplitLbls;
        mtree.GetLeavesUnder(node1, nodeOptSplit);
        mtree.GetlabelsFor( nodeOptSplit, nodeOptSplitLbls );
        DecAllNumInSet(nodeOptSplitLbls);
        set<int> nodeOptSplitLbls2;
        if(node2 >= 0 )
        {
            set<int> nodeOptSplit2;
            mtree.GetLeavesUnder(node2, nodeOptSplit2);
            mtree.GetlabelsFor( nodeOptSplit2, nodeOptSplitLbls2 );
            DecAllNumInSet(nodeOptSplitLbls2);
        }
        SubtractSets(nodeOptSplitLbls, nodeOptSplitLbls2);
        
        ScistPerfPhyCluster clus1(nodeOptSplitLbls);
        clusChangedHetero = clus1;
        ScistPerfPhyCluster clus2(nodeOptSplitLbls2);
        clusChangedHomo = clus2;
    }
    
    return res + this->listSitePriorScore[site];
}

double ScistPerfPhyProbOnTree :: CalcProbForSite(int site, double totEdgeLen, const vector<set<int> > &listClades) const
{
    ScistHaplotypeMat *pHapMat = dynamic_cast<ScistHaplotypeMat *>(&this->genosInput);
    
    if( pHapMat != NULL )
    {
        return CalcProbForSiteHap(site, totEdgeLen, listClades);
    }
    else    // right now, must be of genotype
    {
        return CalcProbForSiteGeno(site, totEdgeLen, listClades);
    }
}

double ScistPerfPhyProbOnTree :: CalcProbForSiteHap(int site, double totEdgeLen, const vector<set<int> > &listClades) const
{
    vector<double> listCladeProb;
    for(int i=0; i<mtree.GetTotNodesNum(); ++i)
    {
        listCladeProb.push_back(-1.0*HAP_MAX_INT);
    }
    
    // get the sum of prob0
    double sumProb0 = 0.0;
    for(int h=0; h<this->genosInput.GetNumHaps(); ++h)
    {
        sumProb0 += log( this->genosInput.GetGenotypeProbAllele0At(h, site) );
    }
    
    double loglikeTot = -1.0*HAP_MAX_INT;
    for(int i=0; i<mtree.GetTotNodesNum(); ++i)
    {
        if( i == mtree.GetRoot() )
        {
            continue;
        }
        double brLen = mtree.GetEdgeLen(i);
        double probPrior = brLen/totEdgeLen;
        double probCladeOnly = 0.0;
        if( mtree.IsLeaf(i))
        {
            int lbl = *listClades[i].begin();
            double p0 = this->genosInput.GetGenotypeProbAllele0At(lbl, site);
            double p1 = 1-p0;
            double pr = log(p1/p0);
            probCladeOnly = pr;
        }
        else
        {
            int childLeft = mtree.GetLeftDescendant(i);
            int childRight = mtree.GetRightDescendant(i);
            probCladeOnly = listCladeProb[childLeft] + listCladeProb[childRight];
        }
//cout << "probPrior: " << probPrior << endl;
        listCladeProb[i] = probCladeOnly + log(probPrior);
        //double probClade = CalcProbMutClade(site, listClades[i] );
        // YW: need to check this
        //loglikeTot = GetLogSumOfTwo(loglikeTot, log(probPrior) + probCladeOnly);
        if( loglikeTot < listCladeProb[i])
        {
            loglikeTot = listCladeProb[i];
        }
    }
    //return loglikeTot + sumProb0;
    double res = loglikeTot + sumProb0;
//cout << "log prob at site: " << res << endl;
    return res;
}

double ScistPerfPhyProbOnTree :: CalcProbForSiteGeno(int site, double totEdgeLen, const vector<set<int> > &listClades) const
{
    ScistPerfPhyProbOnTree *pthis=const_cast<ScistPerfPhyProbOnTree *>(this);
    ScistPerfPhyProbOnTree spppt(pthis->genosInputHap, this->mtree);
    return spppt.CalcProbForSite( 2*site, totEdgeLen, listClades ) + spppt.CalcProbForSite(2*site+1, totEdgeLen, listClades);
}
