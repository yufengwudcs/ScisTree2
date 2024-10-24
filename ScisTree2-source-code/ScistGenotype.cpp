//
//  ScistGenotype.cpp
//  
//
//  Created by Yufeng Wu on 5/25/18.
//
//

#include "ScistGenotype.hpp"
#include "Utils3.h"
#include <cmath>
#include <queue>
#include <iomanip>
#include <thread>
#include "PhylogenyTree.h"
#include "TreeBuilder.h"
#include "MarginalTree.h"
#include "Utils4.h"
#include "UtilsNumerical.h"
#include "RerootTreeUtils.h"

// *************************************************************************************
// genotypes: integer matrix

const double DEF_DROPOUT_PROB = 0.2;

ScistGenGenotypeMat::ScistGenGenotypeMat() : thresSignifcant(0.0), numThreads(1), probDropout(DEF_DROPOUT_PROB), fGenoSet(false)
{
}

void ScistGenGenotypeMat :: TrimCliquesMaxDiff( std::set<std::set<int> > &listCliques, int maxToKeep ) const
{
//cout << "Entering trim, number of cliques: " << listCliques.size() << ", maxToKeep: " << maxToKeep << endl;
    // keep only the most different ones
    if( (int)listCliques.size() <= maxToKeep )
    {
        return;
    }
    // find the distance between two sets
    map< pair<  const set<int> *, const set<int> * >, int > mapPairCliqueDiff;
    for( set<set<int> > :: iterator it1 = listCliques.begin(); it1 != listCliques.end(); ++it1 )
    {
        set<set<int> > :: iterator it2 = it1;
        ++it2;
        for(; it2 != listCliques.end(); ++it2)
        {
            //
            set<int> sint;
            JoinSets(*it1, *it2, sint);
            pair<  const set<int> *, const set<int> * > pp1( &(*it1), &(*it2) ), pp2( &(*it2), &(*it1) );
            mapPairCliqueDiff[pp1] = it1->size()+it2->size()-2*sint.size();
            mapPairCliqueDiff[pp2] = mapPairCliqueDiff[pp1];
        }
    }
    
    // increamentally add the most different; first add the first clique
    set<const set<int> * > listCliquesNext;
    listCliquesNext.insert( &(*listCliques.begin()) );
    while( (int)listCliquesNext.size() < maxToKeep )
    {
        const set<int> * pcliqueToAdd = NULL;
        int diffMax = 0;
        for( set<set<int> > :: iterator it1 = listCliques.begin(); it1 != listCliques.end(); ++it1 )
        {
            //
            if( listCliquesNext.find(&(*it1)) != listCliquesNext.end() )
            {
                //
                continue;
            }
            
            //
            int diffCurr = 0;
            for( set<const set<int> * > :: iterator it2 = listCliquesNext.begin(); it2 != listCliquesNext.end(); ++it2 )
            {
                pair<  const set<int> *, const set<int> * > pp( *it2, &(*it1) );
                YW_ASSERT_INFO(mapPairCliqueDiff.find(pp) != mapPairCliqueDiff.end(), "Fail to find");
                diffCurr += mapPairCliqueDiff[pp];
            }
            if( diffCurr > diffMax)
            {
                diffMax = diffCurr;
                pcliqueToAdd = &(*it1);
            }
        }
        YW_ASSERT_INFO( pcliqueToAdd != NULL, "Cannot be null" );
        listCliquesNext.insert(pcliqueToAdd);
//cout << "In TrimCliquesMaxDiff: adding clique: ";
//DumpIntSet(*pcliqueToAdd);
    }
    set<set<int> > listCliquesNextUse;
    for(set<const set<int> *> :: iterator it = listCliquesNext.begin(); it != listCliquesNext.end(); ++it)
    {
        listCliquesNextUse.insert( *(*it) );
    }
    listCliques = listCliquesNextUse;
}

ScistGenGenotypeMat * ScistGenGenotypeMat :: SubMatrix( const std::set<int> &setRows, const std::set<int> &setSites ) const
{
    ScistGenGenotypeMat *pMatNew = CreateNewMat();
    pMatNew->SetSize( setRows.size(), setSites.size() );
    // set row name
    int rowCurr = 0;
    for( set<int> :: iterator it = setRows.begin(); it != setRows.end(); ++it )
    {
        int siteCurr = 0;
        for( set<int> :: iterator it2 = setSites.begin(); it2 != setSites.end(); ++it2 )
        {
            pMatNew->SetGenotypeAt( rowCurr, siteCurr, GetGenotypeAt(*it, *it2) );
            pMatNew->SetGenotypeProbAt( rowCurr, siteCurr, GetGenotypeProbAllele0At(*it, *it2) );
            ++siteCurr;
        }
        
        pMatNew->SetGenotypeName( rowCurr, GetGenotypeName(*it) );
        ++rowCurr;
    }
    return pMatNew;
}

std::string ScistGenGenotypeMat :: ConsNJTree() const
{
    //
    PhyloDistance dist;
    // setup pairwise hamming distance
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=i+1; j<GetNumHaps(); ++j)
        {
            //
            double d = CalcHammingDistBetwHaps(i, j);
            dist.SetDistance( i, j, d );
//cout << "Distance between (" << i << "," << j << "): " << d << endl;
        }
    }
    DistanceTreeBuilder dtb( dist );
    for(int i=0; i<GetNumHaps(); ++i)
    {
        int indexUse = i+1;
        string strIndexToUse = std::to_string(indexUse);
        dtb.SetTaxonName( i, strIndexToUse );
    }
    return dtb.NJ();
}

std::string ScistGenGenotypeMat :: ConsNJTreeZeroRoot( const std::set<std::set<int> > *psetCladesCurr ) const
{
//cout << "ScistGenGenotypeMat: numthreads " << this->numThreads << endl;
    //if( this->numThreads > 1 && GetNumHaps() > this->numThreads )
    if( this->numThreads > 1  )
    {
        return ConsNJTreeZeroRootMT( psetCladesCurr );
    }
    //
    PhyloDistance dist;
    // setup pairwise hamming distance
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=i+1; j<GetNumHaps(); ++j)
        {
            //
            double d = CalcHammingDistBetwHaps(i, j);
            dist.SetDistance( i, j, d );
            //cout << "Distance between (" << i << "," << j << "): " << d << endl;
        }
    }
    // add one more hap: all-0
    for(int i=0; i<GetNumHaps(); ++i)
    {
        //
        double d = 0.0;
        for( int s=0; s<GetNumSites(); ++s )
        {
            if(GetGenotypeAt(i, s) != 0 )
            {
                d += 1.0;
            }
        }
        d = d/GetNumSites();
        dist.SetDistance(i, GetNumHaps(), d);
    }
    
    // use vector distance instead for faster access
    int maxHapSize = 2*GetNumHaps()+1;
    dist.PrepareVecDists( maxHapSize );
    
    DistanceTreeBuilder dtb( dist );
    for(int i=0; i<=GetNumHaps(); ++i)
    {
        int indexUse = i+1;
        string strIndexToUse = std::to_string(indexUse);
        dtb.SetTaxonName( i, strIndexToUse );
    }
    string strNJWithRoot;
    if( psetCladesCurr == NULL )
    {
        strNJWithRoot = dtb.NJ();
    }
    else
    {
//cout << "Build NJ tree with constraints..\n";
        strNJWithRoot = dtb.ConstrainedNJ( *psetCladesCurr, GetNumHaps() );
    }
//cout << "strNJWithRoot: " << strNJWithRoot << endl;
    // reroot
    string strIdRoot = std::to_string(GetNumHaps()+1);
    char strNJWithRootBuf[102400];
    strcpy( strNJWithRootBuf, strNJWithRoot.c_str() );
    char strIdRootBuf[102400];
    strcpy( strIdRootBuf, strIdRoot.c_str() );
    string strNJWithRootReroot = ReRootTreeNewick( strNJWithRootBuf, strIdRootBuf);
//cout << "strNJWithRootReroot: " << strNJWithRootReroot << endl;
    // remove the root
    MarginalTree mtree;
    ReadinMarginalTreesNewickWLenString(strNJWithRootReroot, this->GetNumHaps()+1, mtree);
    mtree.BuildDescendantInfo();
    int posRootLeaf = mtree.GetPosForLabel(this->GetNumHaps()+1);
    YW_ASSERT_INFO(posRootLeaf >=0, "Fail to find the root");
    mtree.RemoveLeafNodeFromBinaryTree(posRootLeaf);
    mtree.BuildDescendantInfo();
//cout << "Aftre removing reoot: " << mtree.GetNewickSorted(false) << endl;
    return mtree.GetNewickSorted(false);
}

std::string ScistGenGenotypeMat :: ConsConstrainedUPGMATreeZeroRoot( const std::set< std::set<int> > &setMustClades ) const
{
    //
    PhyloDistance dist;
    // setup pairwise hamming distance
    double distMax = -1.0;
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=i+1; j<GetNumHaps(); ++j)
        {
            //
            double d = CalcHammingDistBetwHaps(i, j);
            dist.SetDistance( i, j, d );
            //cout << "Distance between (" << i << "," << j << "): " << d << endl;
            if( d > distMax )
            {
                distMax = d;
            }
        }
    }
    // add one more hap: all-0
    for(int i=0; i<GetNumHaps(); ++i)
    {
        // set to a large dist
        dist.SetDistance(i, GetNumHaps(), distMax * 10.0);
    }
    
    DistanceTreeBuilder dtb( dist );
    for(int i=0; i<=GetNumHaps(); ++i)
    {
        int indexUse = i+1;
        string strIndexToUse = std::to_string(indexUse);
        dtb.SetTaxonName( i, strIndexToUse );
    }
    set<set<int> > setForbid;
    map<set<int>, double> mapSTHts;
    string strTreeWithRoot = dtb.ConstrainedUPGMA(setMustClades, setForbid, mapSTHts );
cout << "strTreeWithRoot: " << strTreeWithRoot << endl;
    
    // remove root
    int plast = strTreeWithRoot.rfind(')');
    YW_ASSERT_INFO(plast != string::npos, "Fail to find 118" );
    int plast2 = strTreeWithRoot.rfind(')', plast-1);
    YW_ASSERT_INFO(plast2 != string::npos, "Fail to find 119" );
    string strTreeNoRoot = strTreeWithRoot.substr(1, (int)plast2) + ";";
    
#if 0
    // reroot
    string strIdRoot = std::to_string(GetNumHaps());
    char strNJWithRootBuf[102400];
    strcpy( strNJWithRootBuf, strTreeWithRoot.c_str() );
    char strIdRootBuf[102400];
    strcpy( strIdRootBuf, strIdRoot.c_str() );
    string strNJWithRootReroot = ReRootTreeNewick( strNJWithRootBuf, strIdRootBuf);
//cout << "strNJWithRootReroot: " << strNJWithRootReroot << endl;
    // remove the root
    MarginalTree mtree;
    ReadinMarginalTreesNewickWLenString(strNJWithRootReroot, this->GetNumHaps()+1, mtree);
    mtree.BuildDescendantInfo();
    int posRootLeaf = mtree.GetPosForLabel(this->GetNumHaps()+1);
    YW_ASSERT_INFO(posRootLeaf >=0, "Fail to find the root");
    mtree.RemoveLeafNodeFromBinaryTree(posRootLeaf);
    mtree.BuildDescendantInfo();
cout << "Aftre removing reoot: " << mtree.GetNewickSorted(false) << endl;
#endif
    
    
    // relabel
    // construct another tree for conversion
    map<int,int> mapIncLeafLbls;
    for(int i=0; i<GetNumHaps(); ++i)
    {
        mapIncLeafLbls[i] = i+1;
    }
    //mapIncLeafLbls[-1] = -1;
    //ChangeLeafIntLabelOfTree(phTree, mapIncLeafLbls);
    PhylogenyTreeBasic treeConv;
    treeConv.ConsOnNewick( strTreeNoRoot );
    ChangeLeafIntLabelOfTree(treeConv, mapIncLeafLbls );
    string resTree;
    treeConv.ConsNewickSorted(resTree, false, 1.0, true);
    
    return resTree;
}

// deal with mutli-threading if we can
static void UtilDistCalcMT( const ScistGenGenotypeMat *pMat, vector<vector<double> > *pdist, int hap1, int hapSz)
{
    //for(int i=hap1; i<=hap2; ++i)
    //{
    //    for(int j=i+1; j<pMat->GetNumHaps(); ++j)
    //    {
    //        //
    //        double d = pMat->CalcHammingDistBetwHaps(i, j);
    //        (*pdist)[i][j] = d;
    //        //cout << "Distance between (" << i << "," << j << "): " << d << endl;
    //    }
    //}
    for(int i=hap1; i<pMat->GetNumHaps(); i+=hapSz)
    {
        for(int j=i+1; j<pMat->GetNumHaps(); ++j)
        {
            //
            double d = pMat->CalcHammingDistBetwHaps(i, j);
            (*pdist)[i][j] = d;
            //cout << "Distance between (" << i << "," << j << "): " << d << endl;
        }
    }
}

std::string ScistGenGenotypeMat :: ConsNJTreeZeroRootMT(const std::set<std::set<int> > *psetCladesCurr) const
{
//cout << "ConsNJTreeZeroRootMT: \n";
    //
    // setup pairwise hamming distance
    vector<vector<double> > listDists(GetNumHaps());
    for(int i=0; i<GetNumHaps(); ++i)
    {
        listDists[i].resize( GetNumHaps() );
    }
    //int blkSize = GetNumHaps()/this->numThreads;
    int numThreadsUse = this->numThreads;
    if( GetNumHaps() < numThreadsUse )
    {
        numThreadsUse = GetNumHaps();
    }
    int blkSize = numThreadsUse;
    //int blkStart = 0;
    vector<thread *> listPtrThreads;
    for(int t=0; t<numThreadsUse; ++t)
    {
        //YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
        // start thread
        thread *pthr = new thread(UtilDistCalcMT, this, &listDists, t, blkSize);
        listPtrThreads.push_back(pthr);
        //++blkStart;
    }
    //for(int t=0; t<this->numThreads; ++t)
    //{
    //    int blkEnd = blkStart + blkSize - 1;
    //    if( blkEnd >= GetNumHaps()  || t == this->numThreads-1)
    //    {
    //        blkEnd = GetNumHaps()-1;
    //    }
    //    YW_ASSERT_INFO(blkEnd >= blkStart, "WRONG345");
    //    // start thread
    //    thread *pthr = new thread(UtilDistCalcMT, this, &listDists, blkStart, blkEnd);
    //    listPtrThreads.push_back(pthr);
    //    blkStart += blkSize;
    //}
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
//cout << "Multi-threaidng init of tree distance: complete\n";
    PhyloDistance dist;
    // init all dists to 0
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=i+1; j<GetNumHaps(); ++j)
        {
            dist.SetDistance( i, j, listDists[i][j] );
        }
    }
    
    // add one more hap: all-0
    for(int i=0; i<GetNumHaps(); ++i)
    {
        //
        double d = 0.0;
        for( int s=0; s<GetNumSites(); ++s )
        {
            if(GetGenotypeAt(i, s) != 0 )
            {
                d += 1.0;
            }
        }
        d = d/GetNumSites();
        dist.SetDistance(i, GetNumHaps(), d);
    }
//cout << "Zero hap: added\n";
    DistanceTreeBuilder dtb( dist );
    dtb.SetNumThreads(numThreadsUse);
    // use vector distance instead for faster access
    int maxHapSize = 2*GetNumHaps()+1;
    dist.PrepareVecDists( maxHapSize );
    
    for(int i=0; i<=GetNumHaps(); ++i)
    {
        int indexUse = i+1;
        string strIndexToUse = std::to_string(indexUse);
        dtb.SetTaxonName( i, strIndexToUse );
    }
    //string strNJWithRoot = dtb.NJ();
    string strNJWithRoot;
    if( psetCladesCurr == NULL )
    {
        strNJWithRoot = dtb.NJ();
    }
    else
    {
        strNJWithRoot = dtb.ConstrainedNJ( *psetCladesCurr, GetNumHaps() );
    }
//cout << "strNJWithRoot: " << strNJWithRoot << endl;
    // reroot
    string strIdRoot = std::to_string(GetNumHaps()+1);
    char *pstrNJWithRootBuf = new char[10240000];
    strcpy( pstrNJWithRootBuf, strNJWithRoot.c_str() );
    char *pstrIdRootBuf = new char[10240000];
    strcpy( pstrIdRootBuf, strIdRoot.c_str() );
    string strNJWithRootReroot = ReRootTreeNewick( pstrNJWithRootBuf, pstrIdRootBuf);
    delete [] pstrNJWithRootBuf;
    delete [] pstrIdRootBuf;
//cout << "strNJWithRootReroot: " << strNJWithRootReroot << endl;
    // remove the root
    MarginalTree mtree;
    ReadinMarginalTreesNewickWLenString(strNJWithRootReroot, this->GetNumHaps()+1, mtree);
    mtree.BuildDescendantInfo();
    int posRootLeaf = mtree.GetPosForLabel(this->GetNumHaps()+1);
    YW_ASSERT_INFO(posRootLeaf >=0, "Fail to find the root");
    mtree.RemoveLeafNodeFromBinaryTree(posRootLeaf);
    mtree.BuildDescendantInfo();
//cout << "Aftre removing reoot: " << mtree.GetNewickSorted(false) << endl;
    return mtree.GetNewickSorted(false);
}

// *****************************************************************************
// Construct perfect phylogeny heuristicly

std::string ScistGenGenotypeMat :: ConsPerfPhylogenyHeu() const
{
    // construct perfect phylogeny heuristically. Approach:
    // (i) first sort sites by the number of 1s
    // (ii) construct mutation graph
    // (iii) construct haplotypes from mutation graph
    vector<int> listSites;
    OrderSitesByNum1s(listSites);
//cout << "ConsPerfPhylogenyHeu: list of ordered sites: ";
//DumpIntVec(listSites);
    vector<int> listSitesPre(GetNumSites());       // [i]: the previous position (NOT site) (edge in mutaiton graph)
    for(int s=0; s<GetNumSites(); ++s)
    {
        int site = listSites[s];
        listSitesPre[s] = -1;       // directly to root of mutation tree
        
        for(int s1=s-1; s1>=0; --s1)
        {
            int site1 = listSites[s1];
            double p1 = CalcLogProbTwoSitesContain(site1, site);
            //double p3 = CalcLogProbTwoSitesContain(s, s1);
            double p2 = CalcLogProbTwoSitesDisjoint(site1, site);
//cout << "site " << s1 << " " << s << ": p1: " << p1 << ", p2: " << p2 << endl;
            //const double THRES = 0.0;
            //if( p1 >= p2 && p1 >= THRES )
            if( p1 >= p2 )
            //if( p1 >= p2  )
            {
                listSitesPre[s] = s1;
                break;
            }
        }
    }
//cout << "Mutation graph: \n";
//DumpIntVec( listSitesPre );
    
    // find list of opt haplotypes
    BinaryMatrix matHap( GetNumHaps(), GetNumSites() );
    
    //OptParseHapsForMG( listSites, listSitesPre, matHap );
//exit(1);
//#if 0
    for(int h=0; h<GetNumHaps(); ++h)
    {
        vector<int> hap;
        OptParseHapForMG(h, listSites, listSitesPre, hap);
//cout << "Optimal hap: h=" << h << ", hap=";
//DumpIntVec(hap);
        matHap.SetRow(h, hap);
    }
//cout << "Haplotype matrix: ";
//matHap.Dump();
//#endif
    
    string resTree;
    // build perfect phylogeny from matrix (with zero root)
    PhylogenyTree phTree;
    vector<int> root;
    for(int i=0; i<GetNumSites(); ++i)
    {
        root.push_back(0);
    }
    phTree.SetRoot(root);
    bool f = phTree.ConsOnBinMatrix(matHap);
    if( !f )
    {
        cout << "Fail to construct perfect phylogeny\n";
        exit(1);
    }
//phTree.ConsNewickSorted(resTree);
//cout << "Perfect phylogeny (non-binary): " << resTree << endl;
    
#if 0
    // iterate over the tree
    PhylogenyTreeIterator itor(phTree);
    itor.Init();
    while( itor.IsDone() == false )
    {
        cout << "Node label: " << itor.GetCurrNode()->GetLabel() << ", num of children: " << itor.GetCurrNode()->GetChildrenNum() << endl;
        for(int ii=0; ii<itor.GetCurrNode()->GetChildrenNum(); ++ii)
        {
            cout << "[child: " << itor.GetCurrNode()->GetChild(ii)->GetLabel() << "]: ";
            vector<int> listLbl;
            itor.GetCurrNode()->GetEdgeLabelsAtBranch(ii, listLbl);
            DumpIntVec(listLbl);
        }
        itor.Next();
    }
#endif
    
//#if 0
    // now refine the tree
    phTree.RefineByNumMuts();
    phTree.ConsNewickSorted(resTree);
    //cout << "After refine: tree: " << resTree << endl;
//#endif
    //string resTree = ConsRootedPerfectPhylogenyFromMat( matHap, false, true );
    // make binary tree
    phTree.Binarize();
    phTree.ConsNewickSorted(resTree);
    //cout << "After binarize: tree: " << resTree << endl;
    
    // construct another tree for conversion
    map<int,int> mapIncLeafLbls;
    for(int i=0; i<matHap.GetRowNum(); ++i)
    {
        mapIncLeafLbls[i] = i+1;
    }
    //mapIncLeafLbls[-1] = -1;
    //ChangeLeafIntLabelOfTree(phTree, mapIncLeafLbls);
    PhylogenyTreeBasic treeConv;
    treeConv.ConsOnNewick( resTree );
    ChangeLeafIntLabelOfTree(treeConv, mapIncLeafLbls );
    treeConv.ConsNewickSorted(resTree, false, 1.0, true);
    
//cout << "Heuristic perfect phylogeny: " << resTree << endl;
//exit(1);
    return resTree;
}

void ScistGenGenotypeMat :: OrderSitesByNum1s(std::vector<int> &listSites) const
{
    // Order sites by the number of 1s (largest in front) in listSites
    listSites.clear();
    vector<pair<int,int> > listSiteOne;
    for(int s=0; s<GetNumSites(); ++s)
    {
        //
        int num0s = 0;
        for(int c=0; c<GetNumHaps(); ++c)
        {
            double p1 = GetGenotypeProbAllele0At(c, s);
            if( p1 >= 0.5 )
            {
                ++num0s;
            }
        }
        listSiteOne.push_back(make_pair( num0s, s ));
    }
    std::sort( listSiteOne.begin(), listSiteOne.end() );
    for(unsigned int i=0; i<listSiteOne.size(); ++i)
    {
        listSites.push_back( listSiteOne[i].second );
    }
}
double ScistGenGenotypeMat :: CalcLogProbTwoSitesContain(int s1, int s2) const
{
    // log prob of site s1's 1s contains site s2's 1s
    double res = 0.0;
    for(int c=0; c<GetNumHaps(); ++c)
    {
        double p1 = GetGenotypeProbAllele0At(c, s1);
        double p2 = GetGenotypeProbAllele0At(c, s2);
        //if( p1 < 0.5 || p2 < 0.5 )
        //{
        //    res += log( 1.0-p1*(1.0-p2) );
        //}
        //res += log( 1.0-p1*(1.0-p2) );
//#if 0
        if( p1 >= 0.5 && p2 < 0.5 )
        {
            res += log(  max( (1-p1)/p1, p2/(1-p2)  ) );
        }
//#endif
        //res += log(  max( (1-p1)/p1, p2/(1-p2)  ) );
        //res += log(  1.0- p1*(1-p2) );
    }
    return res;
}
double ScistGenGenotypeMat :: CalcLogProbTwoSitesDisjoint(int s1, int s2) const
{
    // log prob of site s1's 1s and site s2's 1s are disjoint
    double res = 0.0;
    for(int c=0; c<GetNumHaps(); ++c)
    {
        double p1 = GetGenotypeProbAllele0At(c, s1);
        double p2 = GetGenotypeProbAllele0At(c, s2);
        //if( p1 < 0.5 || p2 < 0.5 )
        //{
        //    res = log( 1.0-(1.0-p1)*(1.0-p2) );
        //}
        //res = log( 1.0-(1.0-p1)*(1.0-p2) );
//#if 0
        if( p1 < 0.5 && p2 < 0.5 )
        {
            res += log(  max( p1/(1-p1), p2/(1-p2)  ) );
        }
//#endif
        //res += log( 1.0-(1-p1)*(1-p2) );
    }
    return res;
}

void ScistGenGenotypeMat :: OptParseHapForMG( int hap, const std::vector<int> &listSites, const std::vector<int> &listMGEdges, std::vector<int> &hapOpt ) const
{
    // find optimal haplotype that gives the maximum probability over the fixed mutation graph
    double maxProb = 0.0;   // default is every genotype is zero
    //double maxProb = MAX_NEG_DOUBLE_VAL;
    int sitePathMaxPos = -1;
    vector<double> vecMaxProbPathEnding( listSites.size() );
    for(unsigned int i=0; i<vecMaxProbPathEnding.size(); ++i)
    {
        vecMaxProbPathEnding[i] = 0.0;
    }
    for( unsigned int i=0; i<listSites.size(); ++i )
    {
        int site = listSites[i];
        int sitePrePos =  listMGEdges[i];
        
        double p0 = GetGenotypeProbAllele0At(hap, site);
        double logratio = log( this->probDropout + (1.0-probDropout) * (1.0-p0)/p0 );
        double base = 0.0;
        if( sitePrePos >= 0 )
        {
            //base = vecMaxProbPathEnding[ listSites[ sitePrePos] ];
            base = vecMaxProbPathEnding[  sitePrePos ];
        }
        vecMaxProbPathEnding[i] = base + logratio;
        if( vecMaxProbPathEnding[i]  > maxProb  )
        {
            maxProb = vecMaxProbPathEnding[i];
            sitePathMaxPos = i;
        }
    }
    // parse based on the site max: only sites laong the max site path is 1, all other are zero
    hapOpt.clear();
    for(int i=0; i<GetNumSites(); ++i)
    {
        hapOpt.push_back(0);
    }
    if( sitePathMaxPos >= 0 )
    //if(sitePathMax > MAX_NEG_DOUBLE_VAL )
    {
        // fill in
        int sc = sitePathMaxPos;
        while(sc >= 0 )
        {
            int site = listSites[sc];
            hapOpt[site] = 1;
            sc = listMGEdges[sc];
        }
    }
}

void ScistGenGenotypeMat :: OptParseHapsForMG( const std::vector<int> &listSites, const std::vector<int> &listMGEdges, BinaryMatrix &matHaps ) const
{
//#if 0
    // construct a binary haplotype matrix based on the given mutation graph
    // that is, (i) parallel sites: mutual exclusive in the mutants
    // (ii) ancestor/descendent: one contain another
    // approach: (i) get most probable genotype (0/1) into the matrix
    // (ii) modify matrix s.t. containment relations are satisified (botthm up iteration of mutations)
    // (iii) modify matrix s.t. no intersecton between parallel sites
    
    // get most probable binary matrix
    for(int c=0; c<GetNumHaps(); ++c)
    {
        for(int s=0; s<GetNumSites(); ++s)
        {
            double p1 = GetGenotypeProbAllele0At(c, s);
            int allele = 0;
            if( p1 < 0.5 )
            {
                allele = 1;
            }
            matHaps(c, s) = allele;
        }
    }
    cout << "Most probable matrix: ";
    matHaps.Dump();
    
    // collect info about mutation graph
    // first descendent node list
    vector<set<int> > listDescSitesPerSite(GetNumSites());
    // backwards iteration
    for(int s=GetNumSites()-1; s>=0; --s)
    {
        int sitePrePos = listMGEdges[s];
        if( sitePrePos >= 0 )
        {
            listDescSitesPerSite[sitePrePos].insert(  listSites[s]  );
            UnionSets( listDescSitesPerSite[sitePrePos], listDescSitesPerSite[s] );
        }
    }
    cout << "List of descendent sites: \n";
    for(int s=0; s<GetNumSites(); ++s)
    {
        cout << "Site " << listSites[s] << ": ";
        DumpIntSet( listDescSitesPerSite[s] );
    }
    // second, ancestor node
    vector<set<int> > listAncesSites( GetNumSites() );
    for( int s=0; s<GetNumSites(); ++s )
    {
        int sitePrePos = listMGEdges[s];
        if( sitePrePos >= 0 )
        {
            listAncesSites[s].insert( sitePrePos );
            UnionSets( listAncesSites[s], listAncesSites[sitePrePos] );
        }
    }
    cout << "List of ancestral sites: \n";
    for(int s=0; s<GetNumSites(); ++s)
    {
        cout << "Site " << listSites[s] << ": ";
        DumpIntSet( listAncesSites[s] );
    }
    // construct children nodes
    vector<set<int> > listChildPos( GetNumSites() );
    for( int s=0; s<GetNumSites(); ++s )
    {
        int sitePrePos = listMGEdges[s];
        if( sitePrePos >= 0 )
        {
            listChildPos[sitePrePos].insert(s);
        }
    }
    cout << "List of children positions: \n";
    for(int s=0; s<GetNumSites(); ++s)
    {
        cout << "Site " << listSites[s] << ": ";
        DumpIntSet( listChildPos[s] );
    }
    // construct reverse map: site to position in list
    vector<int> listReverseSitePos( GetNumSites() );
    for(int s=0; s<GetNumSites(); ++s)
    {
        listReverseSitePos[ listSites[s]  ] = s;
    }
    
    // iterating through mutation graph s.t. site containment among ancestral relations are satisified
    for(int s=GetNumSites()-1; s>=0; --s)
    {
        int sitePrePos = listMGEdges[s];
        if( sitePrePos >= 0 )
        {
            // search for all haps
            int site = listSites[s];
            for(int h=0; h<GetNumHaps(); ++h)
            {
                int allele = matHaps( h, site );
                if( allele == 1 )
                {
                    // propogate all the ancestor sites to be 1 at this hap too
                    for( set<int> :: iterator it = listAncesSites[s].begin(); it != listAncesSites[s].end(); ++it )
                    {
                        int allele2 = matHaps(h, *it);
                        if( allele2 == 0 )
                        {
                            cout << "^^^^INCLUSION: Changing [" << h << ", " << *it << "] to 1 \n";
                            matHaps(h, *it) = 1;
                        }
                    }
                }
            }
        }
    }
    cout << "After inclusion step, matrix: ";
    matHaps.Dump();
    
    // remove 1s s.t. no disjoint sites have common haps
    // first find out any conflicts: <site1, site2>, <set of haplotypes in both>
    map<pair<int,int>, set<int> > mapPairSiteConflicts;
    map<pair<int,int>, int> mapHapAlleleConflictOcc;    // [row, col], number of times in conflicts
    for(int s=0; s<GetNumSites(); ++s)
    {
        set<int> ss = listDescSitesPerSite[s];
        ss.insert( listSites[s] );
        // these sites should be disjoint of all other sites; if not, it is conflict
        set<int> ssComp = ss;
        ComplementIntSet(GetNumSites(), ssComp);
        SubtractSets( ssComp,  listAncesSites[s] );
        for( set<int> :: iterator itt2 = ssComp.begin(); itt2 != ssComp.end(); ++itt2 )
        {
            for(int c=0; c<GetNumHaps(); ++c)
            {
                int site1 = listSites[s], site2 = *itt2;
                OrderInt(site1, site2);
                if( matHaps(c, site1 ) == 1 && matHaps(c, site2) == 1 )
                {
                    mapPairSiteConflicts[make_pair(site1, site2)].insert(c);
                    
                    // keep track records
                    ++mapHapAlleleConflictOcc[ make_pair(c, site1) ];
                    ++mapHapAlleleConflictOcc[ make_pair(c, site2) ];
                }
            }
        }
    }
    cout << "List of conflict pairs: \n";
    for( map<pair<int,int>, set<int> > :: iterator itt = mapPairSiteConflicts.begin(); itt != mapPairSiteConflicts.end(); ++itt )
    {
        cout << "Sites " << itt->first.first << "," << itt->first.second << ": conflict haplotypes: ";
        DumpIntSet(itt->second);
    }
    cout << "List of conflicts:\n";
    for(map<pair<int,int>, int> :: iterator it = mapHapAlleleConflictOcc.begin(); it != mapHapAlleleConflictOcc.end(); ++it)
    {
        cout << "[" << it->first.first << "," << it->first.second << "]: " << it->second << endl;
    }
    
    // apply a greedy algorithm to remove the site
    priority_queue< pair<int, pair<int,int> > > queueConfOccAllelePos;
    for(map<pair<int,int>, int> :: iterator it = mapHapAlleleConflictOcc.begin(); it != mapHapAlleleConflictOcc.end(); ++it)
    {
        queueConfOccAllelePos.push( make_pair(it->second, it->first) );
    }
    while( queueConfOccAllelePos.empty() == false )
    {
        pair<int,int> pp = queueConfOccAllelePos.top().second;
        queueConfOccAllelePos.pop();
        
        // clear out this position
        //matHaps( pp.first, pp.second ) = 0;
        // also clear out all the contained sites for this hap
        int site = pp.second;
        int sitePos = listReverseSitePos[site];
        set<int> ssSub = listDescSitesPerSite[sitePos];
        ssSub.insert(pp.second);
        for(set<int> :: iterator ittt = ssSub.begin(); ittt != ssSub.end(); ++ittt)
        {
            // stop if all conflict has been resolved
            if( mapPairSiteConflicts.size() == 0 )
            {
                break;
            }
            
            int sitecur = *ittt;
            if( matHaps(pp.first, sitecur) == 1 )
            {
                cout << "*** SET [" << pp.first << ", " << sitecur << "] to 0\n";
                matHaps(pp.first, sitecur) = 0;
                
                // remove the record of confict
                for( map<pair<int,int>, set<int> > :: iterator ittg = mapPairSiteConflicts.begin(); ittg != mapPairSiteConflicts.end();  )
                {
                    if( ittg->first.first == sitecur || ittg->first.second == sitecur )
                    {
                        ittg->second.erase( pp.first );
                    }
                    if( ittg->second.size() == 0 )
                    {
                        mapPairSiteConflicts.erase(ittg++);
                    }
                    else
                    {
                        ++ittg;
                    }
                }
            }
        }
    }
    
    cout << "After cleaning up: hapMat: ";
    matHaps.Dump();
    
//#endif
}

std::string ScistGenGenotypeMat :: ConsNJTreeNoInc() const
{
    PhyloDistance dist;
    // setup pairwise hamming distance
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=i+1; j<GetNumHaps(); ++j)
        {
            //
            double d = CalcHammingDistBetwHaps(i, j);
            dist.SetDistance( i, j, d );
            //cout << "Distance between (" << i << "," << j << "): " << d << endl;
        }
    }
    DistanceTreeBuilder dtb( dist );
    return dtb.NJ();
}

double ScistGenGenotypeMat :: CalcHammingDistBetwHaps(int h1, int h2) const
{
    if( this->fGenoSet == false )
    {
        // if genotypes not set, use expectation
        return CalcHammingDistBetwHaps2(h1, h2);
    }
    // otherwise, just use the Hamming distance of the set genotypes
//#if 0
    int numDiffs = 0;
    //int numSkipped = 0;
    for(int c=0; c<GetNumSites(); ++c)
    {
        pair<int,int> pp1(h1, c), pp2(h2, c);
        bool f1Uncertain = this->setUncertainPos.find(pp1) != this->setUncertainPos.end();
        bool f2Uncertain = this->setUncertainPos.find(pp2) != this->setUncertainPos.end();
        double prob10 = GetGenotypeProbAllele0At(h1, c);
        double prob20 = GetGenotypeProbAllele0At(h2, c);
        int g1 = GetGenotypeAt(h1, c);
        int g2 = GetGenotypeAt(h2, c);
        
        if( f1Uncertain == false && f2Uncertain == false )
        {
            if(g1 != g2 )
            {
                ++numDiffs;
            }
        }
        else if( f1Uncertain == true && f2Uncertain == false )
        {
            if(  IsProbAtCellPosSignificant( h1, c, GetSignificanceThres() )  )
            {
                numDiffs += prob10*g2  + (1-prob10)*(1-g2);
            }
            else
            {
                //++numSkipped;
            }
        }
        else if( f1Uncertain == false && f2Uncertain == true )
        {
            if( IsProbAtCellPosSignificant( h2, c, GetSignificanceThres()  ) )
            {
                numDiffs += (1-g1)*(1-prob20) + g1*prob20;
            }
            else
            {
                //++numSkipped;
            }
        }
        else
        {
            if(  IsProbAtCellPosSignificant( h1, c, GetSignificanceThres()  )
               &&  IsProbAtCellPosSignificant( h2, c, GetSignificanceThres() ) )
            {
                numDiffs += prob10*(1-prob20) + (1-prob10)*prob20;
            }
            else
            {
                //++numSkipped;
            }
        }
    }
    //const double PSUEDO_CNT = 0.001;
    //int numSitesAdj = GetNumSites()-numSkipped + PSUEDO_CNT;
    //YW_ASSERT_INFO( numSitesAdj > 0, "Cannot be negative"  );
    //return (1.0*numDiffs)/numSitesAdj;
    return (1.0*numDiffs)/GetNumSites();
//#endif
}

double ScistGenGenotypeMat :: CalcHammingDistBetwHaps2(int h1, int h2) const
{
    // use probabilistic way of measuring difference
    double numDiffs = 0.0;
    //int numSkipped = 0;
    for(int c=0; c<GetNumSites(); ++c)
    {
        if( IsProbAtCellPosSignificant( h1, c, GetSignificanceThres()  )
           &&  IsProbAtCellPosSignificant( h2, c, GetSignificanceThres() ) )
        {
            double prob10 = GetGenotypeProbAllele0At(h1, c);
            double prob20 = GetGenotypeProbAllele0At(h2, c);
            numDiffs += prob10*(1-prob20) + (1-prob10)*prob20;
        }
        else
        {
        //    ++numSkipped;
        }
    }
    //return numDiffs/GetNumSites();
    //const double PSUEDO_CNT = 0.001;
    //int numSitesAdj = GetNumSites()-numSkipped + PSUEDO_CNT;
    //YW_ASSERT_INFO( numSitesAdj > 0, "Cannot be negative"  );
    //return (1.0*numDiffs)/numSitesAdj;
    return (1.0*numDiffs)/GetNumSites();
}

void ScistGenGenotypeMat :: ConsCompatMap( std::set<std::pair<int,int> > &setCompatPairs ) const
{
    //
    setCompatPairs.clear();
    for(int s1=0; s1<GetNumSites(); ++s1)
    {
        for(int s2=s1+1; s2<GetNumSites(); ++s2)
        {
            if( IsCompatible(s1, s2) )
            {
                pair<int,int> pp(s1,s2);
                setCompatPairs.insert(pp);
            }
        }
    }
}

bool ScistGenGenotypeMat :: AreSitesCompatInMap(const std::set<std::pair<int,int> > &setCompatPairs, int s1, int s2 )
{
    //
    pair<int,int> pp(s1,s2);
    OrderInt(pp.first, pp.second);
    return setCompatPairs.find(pp) != setCompatPairs.end();
}

int ScistGenGenotypeMat :: GetGenotypeNumOf(int geno) const
{
    int res = 0;
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=0; j<GetNumSites(); ++j)
        {
            if( GetGenotypeAt(i,j) == geno )
            {
                ++res;
            }
        }
    }
    return res;
}

int ScistGenGenotypeMat :: FindCellByName(const std::string &strName) const
{
    //
    for(int i=0; i<GetNumHaps(); ++i)
    {
        if( GetGenotypeName(i) == strName)
        {
            return i;
        }
    }
    return -1;
}

void ScistGenGenotypeMat ::Dump() const
{
    cout << "Genotype names: ";
    for(int i=0; i<GetNumHaps(); ++i)
    {
        cout << GetGenotypeName(i) << "  ";
    }
    cout << endl;
}

void ScistGenGenotypeMat :: ChangeGenosAtPositions(const std::set< std::pair<std::pair<int,int>, int> > &listChangedPlaces)
{
    //
    for( std::set< std::pair<std::pair<int,int>, int> > :: const_iterator it = listChangedPlaces.begin(); it != listChangedPlaces.end(); ++it )
    {
        //
        SetGenotypeAt( it->first.first, it->first.second, it->second );
    }
}

// calc sum of all log-prob of allele-0 over all cells and all sites
double ScistGenGenotypeMat :: CalcSumLogProbAllele0All() const
{
    if(listSumLogprob0Sites.size() == 0)
    {
        ScistGenGenotypeMat *pthis = const_cast<ScistGenGenotypeMat *>(this);
        pthis->InitSumAllele0Probs();
    }
    return this->sumAllLogProb0;
    /*double res = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
        for(int c=0; c<GetNumHaps(); ++c)
        {
            double logprob = log( GetGenotypeProbAllele0At(c, site) );
            res += logprob;
        }
    }
    return res; */
}

double ScistGenGenotypeMat :: CalcMaxQForTaxaAt(int site, const std::set<int> &setTaxa, double baseProb, double maxQMin) const
{
    // base prob: base probability to be added upon
    double res = baseProb;
    for(set<int> :: const_iterator it = setTaxa.begin(); it != setTaxa.end(); ++it)
    {
        int c = (*it)-1;
        double prob = GetGenotypeProbAllele0At(c, site);
        if( prob < 0.5 )
        {
            res += log((1-prob)/prob);
        }
    }
    return std::max(res, baseProb + maxQMin);
}

double ScistGenGenotypeMat :: GetSumMaxQToLastSite(int site) const
{
    if(listMaxQSumSitesAfter.size() == 0)
    {
        ScistGenGenotypeMat *pthis = const_cast<ScistGenGenotypeMat *>(this);
        pthis->InitSumAllele0Probs();
    }
    return listMaxQSumSitesAfter[site];
}


double ScistGenGenotypeMat :: GetSumLogProbAllele0Site(int site) const
{
    if( listSumLogprob0Sites.size() == 0 )
    {
        ScistGenGenotypeMat *pthis = const_cast<ScistGenGenotypeMat *>(this);
        pthis->InitSumAllele0Probs();
    }
    return listSumLogprob0Sites[site];
}

void ScistGenGenotypeMat :: InitSumAllele0Probs()
{
    // init over all sites
    this->sumAllLogProb0 = 0.0;
    for(int site=0; site<GetNumSites(); ++site)
    {
//cout << "site: " << site << endl;
        double sumAllele0Site = 0.0;
        for(int c=0; c<GetNumHaps(); ++c)
        {
//cout << "  " << GetGenotypeProbAllele0At(c, site);
            double logprob = log( GetGenotypeProbAllele0At(c, site) );
            sumAllele0Site += logprob;
        }
//cout << endl;
//cout << "sumAllele0Site: " << sumAllele0Site << ",sumAllLogProb0: " << sumAllLogProb0 << endl;
        this->listSumLogprob0Sites.push_back(sumAllele0Site);
        this->sumAllLogProb0 += sumAllele0Site;
    }
    
    // also init sum maxQ values for speeding up
    this->listMaxQSumSitesAfter.clear();
    this->listMaxQSumSitesAfter.resize(GetNumSites());
    for(int site=0; site<GetNumSites(); ++site)
    {
        double sumMaxQSite = 0.0;
        for(int c=0; c<GetNumHaps(); ++c)
        {
            double p0 =  GetGenotypeProbAllele0At(c, site);
            if( p0 < 0.5 )
            {
                sumMaxQSite += log((1-p0)/p0);
            }
        }
        this->listMaxQSumSitesAfter[site] = sumMaxQSite;
    }
    // now do a running sum from last site
    double sumMQToLast = 0.0;
    for(int site=GetNumSites()-1; site>=0; --site)
    {
        sumMQToLast += this->listMaxQSumSitesAfter[site];
        this->listMaxQSumSitesAfter[site] = sumMQToLast;
    }
    // add an extra 0 to the end
    this->listMaxQSumSitesAfter.push_back(0.0);
}

double ScistGenGenotypeMat :: CalcEntropyAt(int sc, int site) const
{
    double res = 0.0;
    for(int g=0; g<GetNumGenos(); ++g)
    {
        double probg = GetGenotypeProbAt(sc, site, g);
        res += probg * log( probg );
    }
    return -1.0*res;
}

double ScistGenGenotypeMat :: CalcEntropySite(int site) const
{
    double res = 0.0;
    for(int sc=0; sc<GetNumHaps(); ++sc)
    {
        res += CalcEntropyAt(sc, site);
    }
    return res;
}

void ScistGenGenotypeMat ::  ResetMaximalGenos()
{
    // reset all genotypes to maximal genotypes
    for(int sc=0; sc<GetNumHaps(); ++sc)
    {
        for(int s = 0; s<GetNumSites(); ++s)
        {
            double p0 = GetGenotypeProbAllele0At(sc, s);
            int g = 0;
            if( p0 < 0.5)
            {
                g = 1;
            }
            SetGenotypeAt( sc, s, g );
        }
    }
}

bool ScistGenGenotypeMat :: AreAllGenosCompatible() const
{
    for(int s1=0; s1<GetNumSites(); ++s1)
    {
        for(int s2=s1+1; s2<GetNumSites(); ++s2)
        {
            if( IsCompatible(s1, s2) == false )
            {
                cout << "sites are NOT compatible: " << s1 << ", " << s2 << endl;
                return false;
            }
        }
    }
    return true;
}

void ScistGenGenotypeMat :: Preprocess()
{
    // make sure prob is not too small or too large
    const double MIN_PROB = 0.0000000001;
    for(int c=0; c<GetNumGenos(); ++c)
    {
        for(int s=0; s<GetNumSites(); ++s)
        {
            double p = GetGenotypeProbAllele0At(c,s);
            if( p < MIN_PROB )
            {
                SetGenotypeProbAt(c, s, MIN_PROB);
            }
            else if( p > 1.0-MIN_PROB)
            {
                SetGenotypeProbAt(c, s, 1.0-MIN_PROB);
            }
        }
    }
#if 0
    // also, if thres is set, mark off prob that is close to 0.5
    if( thresSignifcant > 0.0 )
    {
        for(int c=0; c<GetNumGenos(); ++c)
        {
            for(int s=0; s<GetNumSites(); ++s)
            {
                if( IsProbAtCellPosSignificant(c, s, this->thresSignifcant) == false )
                {
                    SetGenotypeProbAt(c, s, 0.5);
                }
            }
        }
    }
#endif
}

void ScistGenGenotypeMat ::  GetCladesFromCalledGenoes(std::set<std::set<int> > &setClades, std::set<int> &setUncertainClades) const
{
    //
    setUncertainClades.clear();
    setClades.clear();
    if(fGenoSet == false)
    {
        return;     // no genotypes called, stop
    }
    // Test: only use the clades that appears more than once
    map<set<int>, set<int> > mapCladeFreq;
    // get the clade that are not-zero
    for(int site=0; site<GetNumSites(); ++site)
    {
        set<int> clade;
        for(int c = 0; c<GetNumHaps(); ++c)
        {
            int allele = GetGenotypeAt(c, site);
            if( allele != 0)
            {
                clade.insert(c);
            }
        }
        if( clade.size() > 1 && clade.size() < GetNumHaps() )
        {
            //setClades.insert(clade);
            //if( mapCladeFreq.find(clade) == mapCladeFreq.end() )
            //{
            //    mapCladeFreq[clade] = 0;
            //}
            mapCladeFreq[clade].insert(site);
        }
    }
    map<int,int> mapFreq;
    for(map<set<int>,set<int> > :: iterator it = mapCladeFreq.begin(); it != mapCladeFreq.end(); ++it)
    {
        if( mapFreq.find(it->second.size()) == mapFreq.end() )
        {
            mapFreq[ it->second.size() ] = 0;
        }
        ++mapFreq[ it->second.size() ];
        if( it->second.size() > 1 )
        {
            setClades.insert(it->first);
        }
        else
        {
            setUncertainClades.insert( *(it->second.begin()) );
        }
    }
    //cout << "List of frequencies of clades: \n";
    //for(map<int,int> :: iterator it = mapFreq.begin(); it != mapFreq.end(); ++it)
    //{
    //    cout << "[" << it->first << "]: " << it->second << endl;
    //}
}

void ScistGenGenotypeMat ::  AdjustProbsFromCalledGenos()
{
    //
    if(fGenoSet == false)
    {
        return;     // no genotypes called, stop
    }
//cout << "AdjustProbsFromCalledGenos: before change: \n";
//this->Dump();
    // Test: only use the clades that appears more than once
    vector< set<int> > listClades(GetNumSites());
    map<set<int>, set<int> > mapCladeFreq;
    // get the clade that are not-zero
    for(int site=0; site<GetNumSites(); ++site)
    {
        set<int> clade;
        for(int c = 0; c<GetNumHaps(); ++c)
        {
            int allele = GetGenotypeAt(c, site);
            if( allele != 0)
            {
                clade.insert(c);
            }
        }
        listClades[site] = clade;
        if( clade.size() > 1 && clade.size() < GetNumHaps() )
        {
            mapCladeFreq[clade].insert(site);
        }
    }
//cout << "pass1\n";
    bool fUpdated = false;
    const int THRES_FREQ_CLADE = 3;
    set<set<int> > setFixedClades;
    set<int> setSitesFixed;
    for(map<set<int>,set<int> > :: iterator it = mapCladeFreq.begin(); it != mapCladeFreq.end(); ++it)
    {
        if( it->second.size() >= THRES_FREQ_CLADE )
        {
            setFixedClades.insert(it->first);
            
            // fixed the prob of these fixed clades
            UnionSets(setSitesFixed, it->second);
            vector<int> scClade( GetNumHaps() );
            for(int ii = 0; ii<GetNumHaps(); ++ii)
            {
                scClade[ii] = 0;
            }
            for( set<int> :: iterator itt2 = it->first.begin(); itt2 != it->first.end(); ++itt2 )
            {
                scClade[*itt2] = 1;
            }
            for(set<int> :: iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
            {
                int sc = *itt;

                // fix prob
                fUpdated = true;
                for(int jj=0; jj<GetNumHaps(); ++jj)
                {
                    IncProbForAlleleBy(jj, sc, scClade[jj], 2.0);
                }
            }
        }
    }
    if( setFixedClades.size() == 0 )
    {
        if( fUpdated )
        {
            InitSumAllele0Probs();
        }
        return;     // nothing to be done
    }
//cout << "pass2\n";
//cout <<"setFixedClades: \n";
//for(set<set<int> > :: iterator it = setFixedClades.begin(); it != setFixedClades.end(); ++it)
//{
//DumpIntSet(*it);
//}
    // now build a rooted perfect phylogeny with fixed clades
    BinaryMatrix binMat( GetNumHaps(), setFixedClades.size() );
    int scur = 0;
    for(set<set<int> > :: iterator it = setFixedClades.begin(); it != setFixedClades.end(); ++it)
    {
        for(int c=0; c<GetNumHaps(); ++c)
        {
            binMat(c, scur) = 0;
        }
        for( set<int> :: iterator it2 = it->begin(); it2 != it->end(); ++it2 )
        {
            binMat(*it2, scur) = 1;
        }
        ++scur;
    }
//cout << "Matrix: ";
//binMat.Dump();
    vector<int> rootZero;
    for(int i=0; i<setFixedClades.size(); ++i)
    {
        rootZero.push_back(0);
    }
    PhylogenyTree phTree;
    phTree.SetRoot(rootZero);
    phTree.ConsOnBinMatrix(binMat);
    // now assign leaf labels
    map<string,string> mapIdToLabels;
    for(int i=0; i<GetNumHaps(); ++i)
    {
        //cout << "i: " << i << ", name: " << this->genosInput.GetGenotypeName(i) << endl;
        string str = "(" + std::to_string(i) + ")";
        mapIdToLabels[str] = std::to_string(i);
    }
    phTree.ReassignLeafLabels( mapIdToLabels );
//cout <<" Tree: ";
//phTree.Dump();
    // get all clades of the tree
    set<set<int> > setCladesAll;
    phTree.GetAllClades(setCladesAll);
    set<set<int> > setCladesCompAll;
    for(set<set<int> > :: iterator itt = setCladesAll.begin(); itt != setCladesAll.end(); ++itt )
    {
        set<int> scomp = *itt;
        ComplementIntSet(GetNumHaps(), scomp);
        setCladesCompAll.insert(scomp);
    }
#if 0
cout <<"setCladesAll: \n";
for(set<set<int> > :: iterator it = setCladesAll.begin(); it != setCladesAll.end(); ++it)
{
DumpIntSet(*it);
}
cout <<"setCladesCompAll: \n";
for(set<set<int> > :: iterator it = setCladesCompAll.begin(); it != setCladesCompAll.end(); ++it)
{
DumpIntSet(*it);
}
#endif
//cout << "pass3\n";
    // check each site to see
    for(int site=0; site<GetNumSites(); ++site)
    {
        if( setSitesFixed.find(site) != setSitesFixed.end())
        {
            continue;   // done before
        }
        // first find set of certain 0 and certain 1 based on curr prob
        set<int> setLikeliy0, setLikeli1;
        for(int c=0; c<GetNumHaps(); ++c)
        {
            double p0 = GetGenotypeProbAllele0At(c, site);
            const double THRES_P0 = 0.95;
            if( p0 >= THRES_P0 )
            {
                setLikeliy0.insert(c);
            }
            else if( p0 < 1.0- THRES_P0 )
            {
                setLikeli1.insert(c);
            }
        }
//cout << "site: " << site << ", setLikeli0: ";
//DumpIntSet(setLikeliy0);
//cout << "   setLikeli1: ";
//DumpIntSet(setLikeli1);
        // now find the best clade of tree to separte these two groups
        double maxSim0 = -1.0;
        set<int> clade0Best;
        if( setLikeliy0.size() > 0 )
        {
            for(set<set<int> > :: iterator itt = setCladesCompAll.begin(); itt != setCladesCompAll.end(); ++itt )
            {
                double simJacard = CalcJaccrdIndexForTwoSets( *itt, setLikeliy0 );
                if( simJacard > maxSim0)
                {
                    maxSim0 = simJacard;
                    clade0Best = *itt;
                }
            }
        }
        double maxSim1 = -1.0;
        set<int> clade1Best;
        if( setLikeli1.size() > 0 )
        {
            for(set<set<int> > :: iterator itt = setCladesAll.begin(); itt != setCladesAll.end(); ++itt )
            {
                double simJacard = CalcJaccrdIndexForTwoSets( *itt, setLikeli1 );
                //cout << "score1: " << simJacard << ", maxSim1: " maxSim1 << ", for clade1: ";
                //DumpIntSet(*itt);
                if( simJacard > maxSim1)
                {
                    maxSim1 = simJacard;
                    clade1Best = *itt;
                }
            }
        }
        // make sure no overlap by removing the intersection
        set<int> ssInt;
        JoinSets(clade0Best, clade1Best, ssInt);
        if( ssInt.size() > 0 )
        {
            SubtractSets(clade0Best, ssInt);
            SubtractSets(clade1Best, ssInt);
        }
//cout << "clade0Best: ";
//DumpIntSet(clade0Best);
//cout << "  clade1Best: ";
//DumpIntSet(clade1Best);
        
        // now fix prob within two found clade
        if( setLikeliy0.size() >  0)
        {
            for(set<int> :: iterator itt = clade0Best.begin(); itt != clade0Best.end(); ++itt)
            {
                IncProbForAlleleBy(*itt, site, 0, 2.0);
            }
        }
        if( setLikeli1.size() > 0 && clade1Best.size() < GetNumHaps() )
        {
            for(set<int> :: iterator itt = clade1Best.begin(); itt != clade1Best.end(); ++itt)
            {
                IncProbForAlleleBy(*itt, site, 1, 2.0);
            }
        }
    }
//cout << "After change: \n";
//this->Dump();
    // update prob
    InitSumAllele0Probs();
}

// *************************************************************************************
// genotypes: binary matrix

ScistHaplotypeMat :: ScistHaplotypeMat()
{
}

ScistGenGenotypeMat * ScistHaplotypeMat :: Copy() const
{
    //
    ScistHaplotypeMat *pMatCopy = new ScistHaplotypeMat();
    string fn = GetFileName();
    pMatCopy->SetFileName(fn);
    
    for(int i=0; i<GetNumNames(); ++i)
    {
        pMatCopy->AddGenotypeName( GetGenotypeName(i) );
    }
    
    pMatCopy->SetSize( GetNumHaps(), GetNumSites() );
    
    //
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=0; j<GetNumSites(); ++j)
        {
            pMatCopy->SetGenotypeAt( i, j, GetGenotypeAt(i,j) );
            pMatCopy->SetGenotypeProbAt(i, j, GetGenotypeProbAllele0At(i,j));
        }
    }
    
    return pMatCopy;
}

bool ScistHaplotypeMat :: ReadFromFile( std::ifstream &infile, int numSites, int numSCs, bool fSiteName )
{
//cout << "ScistHaplotypeMat :: ReadFromFile: numSites: " << numSites << ", numSCs: " << numSCs << endl;
    //
    // assume each site is independent
    SetSize( numSCs, numSites );
    for(int i=0; i<numSites; ++i)
    {
        string strName;
        if( fSiteName )
        {
            infile >> strName;
        }
        else
        {
            strName = std::to_string(i+1);
        }
        AddSiteName(strName);
        
//cout << "Read in site: " << i << endl;
        for(int j=0; j<numSCs; ++j)
        {
            double prob0 = 0.0;
            bool res = ReadFromFileHapProb(infile, prob0);
            if( res == false )
            {
                return false;
            }
            // choose the allele w/ higher prob
            int allele = 0;
            if( prob0 < 0.5)
            {
                allele = 1;
            }
            SetGenotypeAt( j, i, allele );
            
            matHaplotypesProb0[j][i] = prob0;
        }
    }
    
//cout << "Input matrix: ";
//this->matHaplotypes.Dump();
    
    return true;
}

bool ScistHaplotypeMat :: ReadFromFileHapProb( std::ifstream &infile, double &prob0 )
{
    // read in the prob of haploid allele: 0.6 means prob of 0 is 0.6
    // assume prob of 0 + prob of 1 = 1
    infile >> prob0;
    return true;
}

void ScistHaplotypeMat :: SetSize(int numHaps, int numSites)
{
    matHaplotypes.SetSize(numHaps, numSites);
    
    matHaplotypesProb0.clear();
    matHaplotypesProb0.resize(numHaps);

    bool fNameInit = GetNumNames() > 0;
    
    for(int i=0; i<numHaps; ++i)
    {
        matHaplotypesProb0[i].resize( numSites );
        
        // by default, use the numericals, starting from one
        if( fNameInit == false )
        {
            string str = std::to_string(i+1);
            AddGenotypeName(str);
            
//cout << "Init name: " << str << endl;
        }
    }
}

void ScistHaplotypeMat :: SetGenotypeAt(int sc, int site, int geno)
{
    matHaplotypes(sc, site) = geno;
}

void ScistHaplotypeMat :: AddGenotypeAt(int sc, int site, int geno)
{
    // append the genotype into it
    int genoThis = GetGenotypeAt(sc, site);
    if( genoThis == 0 && geno == 1 )
    {
        SetGenotypeAt(sc, site, 1);
    }
}

int ScistHaplotypeMat :: GetAltGenotypeAt(int sc, int site) const
{
    int genoThis = GetGenotypeAt(sc, site);
    if(genoThis == 0)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

double ScistHaplotypeMat :: GetGenotypeProbAllele0At(int sc, int site) const
{
    // return proble of allele 0
    return this->matHaplotypesProb0[sc][site];
}

void ScistHaplotypeMat :: SetGenotypeProbAt(int sc, int site, double prob)
{
    this->matHaplotypesProb0[sc][site] = prob;
}

void ScistHaplotypeMat :: SetGenotypeProbOfGenoAt(int sc, int site, int geno, double prob)
{
    if( geno == 0 )
    {
        SetGenotypeProbAt(sc, site, prob);
    }
    else
    {
        SetGenotypeProbAt(sc, site, 1.0-prob);
    }
}

int ScistHaplotypeMat :: GetGenotypeAt(int sc, int site) const
{
    return matHaplotypes(sc, site);
}

void ScistHaplotypeMat :: FindMaximalCompatSites( const std::vector<double> &wtSites, std::vector< std::map<int, std::set<int> > > &listSetSitesCompat, int maxNumSets, const std::set<std::pair<int,int> > *pSetCompatPairs ) const
{
//#if 0
    //const double DEF_MIN_FRAC = 0.5;
    
    // we find the maximum weightd clique of compatible pairs
    // construct compat pairs if not done yet
    set<std::pair<int,int> > *pSetCompatPairsUse = const_cast<set<std::pair<int,int> > *>( pSetCompatPairs);
    set<pair<int,int> > setCompatPairsAlt;
    if( pSetCompatPairsUse == NULL )
    {
        ConsCompatMap( setCompatPairsAlt );
        pSetCompatPairsUse = &setCompatPairsAlt;
    }
    
    // implement the simple heuristics by Johnson 1974
    //BinaryMatrix &matHaplotypesUse = const_cast<BinaryMatrix &>( this->matHaplotypes );
    
    //
    listSetSitesCompat.clear();
    //vector<vector<bool> > vecHapsFullyCompat( GetNumSites() );
    //for(int i=0; i<GetNumSites(); ++i)
    //{
    //    vecHapsFullyCompat[i].resize( GetNumSites() );
    //}
    
    //for(int s1 = 0; s1<GetNumSites(); ++s1)
    //{
    //    vecHapsFullyCompat[s1][s1] = true;
    //    for(int s2=s1+1; s2<GetNumSites(); ++s2)
    //    {
    //        // root allele: 0
    //        bool fCompat = matHaplotypesUse.IsCompatibleRooted(s1, s2, 0, 0);
    //        vecHapsFullyCompat[s1][s2] = fCompat;
    //        vecHapsFullyCompat[s2][s1] = fCompat;
//cout << "Sites " << s1 << "," << s2 << ": ";
//if(fCompat)
//{
//cout << " compatible\n";
//}
//else
//{
//cout << " not compatible\n";
//}
    //    }
    //}

    //
    set<pair<set<int>,set<int> > > listSetMaxCompatChosen;
    // init
    set<int> ss;
    set<int> setSitesRemainInit;
    PopulateSetWithInterval( setSitesRemainInit, 0, GetNumSites()-1);
    pair<set<int>,set<int> > pp( ss, setSitesRemainInit );
    listSetMaxCompatChosen.insert(pp);
    
    while( true )
    {
        //
        set<pair<set<int>,set<int> > > listSetMaxCompatChosenNext;
    
        for( set<pair<set<int>,set<int> > > :: iterator it = listSetMaxCompatChosen.begin(); it != listSetMaxCompatChosen.end(); ++it )
        {
            set<int> setSitesRemain = it->second;
            
            if( setSitesRemain.size() == 0 )
            {
                continue;
            }
            
            set<int> setMaxCompatChosen = it->first;
            
            // find the one that is the most compatible with remaining sites
            //int maxNumCompat = -1;
            double wtSiteMax = -1.0*HAP_MAX_INT;
            
            vector<int> listSitesNext;
            for( set<int> :: iterator it = setSitesRemain.begin(); it != setSitesRemain.end(); ++it )
            {
            //    int numCompat = 0;
            //    for(set<int> :: iterator it2 = setSitesRemain.begin(); it2 != setSitesRemain.end(); ++it2)
            //    {
            //        if( AreSitesCompatInMap(*pSetCompatPairsUse, *it,*it2) )
            //        {
            //            ++numCompat;
            //        }
            //    }
                double wtCur = wtSites[*it];
            //    if( numCompat > maxNumCompat )
                if( wtCur > wtSiteMax )
                {
                    listSitesNext.clear();
                    listSitesNext.push_back( *it );
                    wtSiteMax = wtCur;
                    //maxNumCompat = numCompat;
                }
            //    else if( numCompat == maxNumCompat )
                if( wtCur == wtSiteMax )
                {
                    listSitesNext.push_back(*it);
                }
            }
            
            // if weight is too small now, stop if we have already get enough
            //if( wtSiteMax < 1.0)
            //{
            //    if( ((int)DEF_MIN_FRAC*GetNumSites()) <= (int)setMaxCompatChosen.size() )
            //    {
            //        break;
            //    }
            //}
            
            for(int jj=0; jj<(int)listSitesNext.size(); ++jj)
            {
                // don't continue adding if we are at the limit
                if( (int)listSetMaxCompatChosenNext.size() > maxNumSets )
                {
                    continue;
                }
                
                int sChose = listSitesNext[jj];
                set<int> setMaxCompatChosenNew = setMaxCompatChosen;
                setMaxCompatChosenNew.insert(sChose);
                
                // remove any sites that are incompatible with the chosen sites
                set<int> setSitesRemainNew;
                for(set<int> :: iterator it = setSitesRemain.begin(); it != setSitesRemain.end(); ++it)
                {
                    if( AreSitesCompatInMap(*pSetCompatPairsUse, sChose, *it) == true )
                    {
                        setSitesRemainNew.insert(*it);
                    }
                }
                setSitesRemainNew.erase( sChose );
                
                pair<set<int>,set<int> > pp( setMaxCompatChosenNew, setSitesRemainNew );
                
                listSetMaxCompatChosenNext.insert(pp);
            }
        }
        
        //
        if( listSetMaxCompatChosenNext.size() == 0 )
        {
            //
            break;
        }
        else
        {
            listSetMaxCompatChosen = listSetMaxCompatChosenNext;
        }
    }
    
    YW_ASSERT_INFO(listSetMaxCompatChosen.size() > 0, "Cannot be empty");
    for( set<pair<set<int>,set<int> > > :: iterator it =  listSetMaxCompatChosen.begin(); it != listSetMaxCompatChosen.end(); ++it)
    {
//cout << "Maximum clique found by the heuristic: ";
//DumpIntSet( it->first );
        map<int,set<int> > mm;
        for(set<int> :: iterator it2 = it->first.begin(); it2 != it->first.end(); ++it2)
        {
            set<int> ss;
            GetMutRowsHapAtSite(*it2, ss);
            mm[*it2] = ss;
        }
        listSetSitesCompat.push_back(mm);
    }
    
//#endif

#if 0
    BinaryMatrix &matHaplotypesUse = const_cast<BinaryMatrix &>( this->matHaplotypes );
    
    //
    listSetSitesCompat.clear();
    vector<vector<bool> > vecHapsFullyCompat( GetNumSites() );
    for(int i=0; i<GetNumSites(); ++i)
    {
        vecHapsFullyCompat[i].resize( GetNumSites() );
    }
    
    for(int s1 = 0; s1<GetNumSites(); ++s1)
    {
        for(int s2=s1+1; s2<GetNumSites(); ++s2)
        {
            // root allele: 0
            bool fCompat = matHaplotypesUse.IsCompatibleRooted(s1, s2, 0, 0);
            vecHapsFullyCompat[s1][s2] = fCompat;
            vecHapsFullyCompat[s2][s1] = fCompat;
//cout << "Sites " << s1 << "," << s2 << ": ";
//if(fCompat)
//{
//cout << " compatible\n";
//}
//else
//{
//cout << " not compatible\n";
//}
        }
    }
    // find maximal compatible components
    set< set<int> > setMaximalComps;
    // start by putting all compatible pairs
    for(int s1 = 0; s1<GetNumSites(); ++s1)
    {
        set<int> ss;
        ss.insert(s1);
        setMaximalComps.insert(ss);
    }
    // find larger
    while(true)
    {
        // every time, make sure size is not too large
        TrimCliquesMaxDiff( setMaximalComps, maxNumSets );
//cout << "Size of current cliques to grow: " << setMaximalComps.size() << endl;
//for( set<set<int> > :: iterator it = setMaximalComps.begin(); it != setMaximalComps.end(); ++it)
//{
//DumpIntSet(*it);
//}
        
        set< set<int> > setMaximalCompsNext;
        // try to grow by adding one more
        for( set<set<int> > :: iterator it = setMaximalComps.begin(); it != setMaximalComps.end(); ++it )
        {
            for(int s=0; s<GetNumSites(); ++s)
            {
                if(  it->find(s) == it->end() )
                {
                    bool fCompat = true;
                    for(set<int> :: iterator it2 = it->begin(); it2 != it->end(); ++it2 )
                    {
                        if( vecHapsFullyCompat[ s ][ *it2 ]  == false )
                        {
                            fCompat = false;
                            break;
                        }
                    }
                    if( fCompat )
                    {
                        set<int> ss = *it;
                        ss.insert( s );
                        setMaximalCompsNext.insert(ss);
//cout << "Growing a subset: ";
//DumpIntSet(ss);
                    }
                }
            }
        }
        if( setMaximalCompsNext.size() == 0 )
        {
            break;
        }
        else
        {
            setMaximalComps = setMaximalCompsNext;
        }
    }
    //
    //TrimCliquesMaxDiff( setMaximalComps, maxNumSets );
    
    YW_ASSERT_INFO( setMaximalComps.size() > 0, "Cannot be empty" );
    for( set<set<int> > :: iterator it = setMaximalComps.begin(); it != setMaximalComps.end(); ++it )
    {
cout << "Clique found: ";
DumpIntSet(*it);
        map<int, std::set<int> >  setSitesCompat;
        
        set<int> ssChosen = *it;
        for(set<int> :: iterator it = ssChosen.begin(); it != ssChosen.end(); ++it)
        {
            set<int> ss;
            GetMutRowsHapAtSite(*it, ss);
            setSitesCompat[*it] = ss;
        }
        listSetSitesCompat.push_back(setSitesCompat);
    }
#endif
}

int ScistHaplotypeMat :: GetNumSites() const
{
    return matHaplotypes.GetColNum();
}

int ScistHaplotypeMat :: GetNumHaps() const
{
    return matHaplotypes.GetRowNum();
}

void ScistHaplotypeMat :: GetMutRowsHapAtSite(int site, std::set<int> &setRows) const
{
    // any allele w/ non-zero is mutant
    setRows.clear();
    for(int r=0; r<matHaplotypes.GetRowNum(); ++r)
    {
        if( matHaplotypes(r, site) == 1 )
        {
            setRows.insert(r);
        }
    }
}

void ScistHaplotypeMat :: GetRowsWithGenoAtSite(int site, int geno, std::set<int> &setRows) const
{
    setRows.clear();
    if( geno == 1 )
    {
        GetMutRowsHapAtSite( site, setRows );
    }
    else if(geno == 0)
    {
        // get the complement
        setRows.clear();
        PopulateSetWithInterval(setRows, 0, GetNumHaps()-1);
        set<int> setRows1;
        GetMutRowsHapAtSite(site, setRows1);
        SubtractSets(setRows, setRows1);
    }
}

double ScistHaplotypeMat :: GetScoreForGeno(int scIndex, int site, int genotype) const
{
    int allele = this->matHaplotypes(scIndex, site);
    if( allele == genotype)
    {
        // when greeing, score is 0
        return 0.0;
    }
    
    // for now, only use default scoring
    double res = 0.0;
    double prob0 = this->matHaplotypesProb0[scIndex][site];
    double prob1 = 1.0-prob0;
    if( genotype == 1)
    {
        // change from 0 to 1
        if( prob1 <= 0.0)
        {
            res = HAP_MAX_INT*1.0;
        }
        else
        {
            res = log(prob0/prob1);
        }
    }
    else
    {
        if( prob0 <= 0.0)
        {
            res = HAP_MAX_INT*1.0;
        }
        else
        {
            res = log( prob1/prob0);
        }
    }
    if( res < 0.0 )
    {
        this->Dump();
        cout << "cell: " << scIndex << ", site: " << site << ", genotype: " << genotype << ", prob0: " << prob0 << endl;
    }
    YW_ASSERT_INFO( res >= 0.0, "Prob: wrong" );
    return res;
}

bool ScistHaplotypeMat :: IsNoninformative(int site) const
{
    //
    BinaryMatrix &matHaplotypesUse = const_cast<BinaryMatrix &>(this->matHaplotypes);
    return matHaplotypesUse.IsColNonInformative(site);
}

bool ScistHaplotypeMat :: IsCompatible(int s1, int s2) const
{
    //
    BinaryMatrix &matHaplotypesUse = const_cast<BinaryMatrix &>(this->matHaplotypes);
    return matHaplotypesUse.IsCompatible(s1, s2);
}

std::string ScistHaplotypeMat :: ConsTree() const
{
    //
    // construct phylogeny
    vector<int> rootZero;
    for(int i=0; i<GetNumSites(); ++i)
    {
        rootZero.push_back(0);
    }
    PhylogenyTree phTree;
    phTree.SetRoot(rootZero);
    phTree.ConsOnBinMatrix(this->matHaplotypes);
    phTree.RemoveDegreeTwoNodes();
    
    // now assign leaf labels
    map<string,string> mapIdToLabels;
    for(int i=0; i<GetNumHaps(); ++i)
    {
        //cout << "i: " << i << ", name: " << this->genosInput.GetGenotypeName(i) << endl;
        string str = "(" + std::to_string(i) + ")";
        mapIdToLabels[str] = GetGenotypeName(i);
    }
    phTree.ReassignLeafLabels( mapIdToLabels );
    
    
    string res;
    phTree.ConsNewickSorted(res);
    //phTree.ConsNewick(res, false, 0.0, true);
    return res;
}

double ScistHaplotypeMat :: SumLogProbs() const
{
    //
    double res = 0.0;
    for(int i=0; i<(int)matHaplotypesProb0.size(); ++i)
    {
        res += GetSumOfVecElements( matHaplotypesProb0[i]);
    }
    return res;
}

void ScistHaplotypeMat :: Dump() const
{
    ScistGenGenotypeMat ::Dump();
    
    //
    cout << "Matrix: [" << GetNumHaps() << "," << GetNumSites() << "]" << endl;
    this->matHaplotypes.Dump();
#if 0
    cout << "Clusters\n";
    for(int c=0; c<GetNumSites(); ++c)
    {
        cout << "Site " << c+1 << ": ";
        set<int> rowsMut;
        this->matHaplotypes.GetRowsWithAllele(c, 1, rowsMut);
        DumpIntSet(rowsMut);
    }
#endif
    cout << "Probabilities: \n";
    for(int i=0; i<(int)matHaplotypesProb0.size(); ++i)
    {
        DumpDoubleVec( matHaplotypesProb0[i]);
    }
}

void ScistHaplotypeMat :: OutputImput(const string *pStrDesc) const
{
    //
    string fileGenoImp = GetImpGenosFileName();
    ofstream outFile(fileGenoImp.c_str());
    if( !outFile.is_open() )
    {
        cout << "Fail to open imputed genotype output file\n";
        exit(1);
    }
    outFile << "Lineages: ";
    for(int i=0; i<GetNumNames(); ++i)
    {
        outFile << GetGenotypeName(i) << "  ";
    }
    outFile << endl;
    if( pStrDesc != NULL )
    {
        outFile << *pStrDesc << endl;
    }
    else
    {
        outFile << "Imputed genotypes: \n";
    }
    for(int s=0; s<GetNumSites(); ++s)
    {
        outFile << "Site " << setw(6) << s+1 << ":\t";

        for(int i=0; i<GetNumHaps(); ++i)
        {
            outFile << GetGenotypeAt(i, s) << " ";
        }
        outFile << endl;
    }
    outFile.close();
}

bool ScistHaplotypeMat :: IsProbSignificant(double prob, double thresVal) const
{
    //
    const double probConst = 0.5;
    if( prob < probConst && prob >( probConst-thresVal/2 ) )
    {
        return false;
    }
    if( prob > probConst && prob < (probConst+thresVal/2) )
    {
        return false;
    }
    return true;
}

void ScistHaplotypeMat :: IncProbForAlleleBy(int sc, int site, int allele, double fac)
{
#if 0
    // assume fac > 1
    double p0 = GetGenotypeProbAllele0At(sc, site);
    double p1 = 1.0 - p0;
    double p0New = p0;
    if( allele == 0 )
    {
        double p0f = p0*fac;
        if( p0f  < 1.0 )
        {
            p0New = p0f;
        }
        else
        {
            p0New = 1.0 - p1/fac;
        }
    }
    else
    {
        double p1f = p1*fac;
        if( p1f  < 1.0 )
        {
            p0New = 1.0 - p1f;
        }
        else
        {
            p0New = p0/fac;
        }
    }
    SetGenotypeProbAt(sc, site, p0New);
#endif
    double p0 = GetGenotypeProbAllele0At(sc, site);
    double p0New = p0;
    const double p0Max = 0.999999;
    if( allele == 0 )
    {
        if( p0New < p0Max )
        {
            p0New = p0Max;
        }
    }
    else
    {
        if( p0New > 1-p0Max  )
        {
            p0New = 1-p0Max;
        }
    }
    SetGenotypeProbAt(sc, site, p0New);
}

// *************************************************************************************
// genotypes: ternary matrix


ScistTernaryMat :: ScistTernaryMat()
{
}

ScistGenGenotypeMat * ScistTernaryMat :: Copy() const
{
    //
    ScistTernaryMat *pMatCopy = new ScistTernaryMat();
    
    for(int i=0; i<GetNumNames(); ++i)
    {
        pMatCopy->AddGenotypeName( GetGenotypeName(i) );
    }
    
    pMatCopy->SetSize( GetNumHaps(), GetNumSites() );
    
    //
    for(int i=0; i<GetNumHaps(); ++i)
    {
        for(int j=0; j<GetNumSites(); ++j)
        {
            pMatCopy->SetGenotypeAt( i, j, GetGenotypeAt(i,j) );
            pMatCopy->SetGenotypeProbOfGenoAt(i, j, 0, GetGenotypeProbAt(i,j, 0));
            pMatCopy->SetGenotypeProbOfGenoAt(i, j, 1, GetGenotypeProbAt(i,j, 1));
        }
    }
    
    return pMatCopy;
}

bool ScistTernaryMat :: ReadFromFile( std::ifstream &infile, int numSites, int numSCs, bool fSiteName )
{
    //
    // assume each site is independent
    SetSize( numSCs, numSites );
    for(int i=0; i<numSites; ++i)
    {
        string strName;
        if( fSiteName )
        {
            infile >> strName;
        }
        else
        {
            strName = std::to_string(i+1);
        }
        AddSiteName(strName);
        
        //cout << "Read in site: " << i << endl;
        for(int j=0; j<numSCs; ++j)
        {
            double prob0 = 0.0, prob1 = 0.0;
            bool res = ReadFromFileTernaryProb(infile, prob0, prob1);
            if( res == false )
            {
                return false;
            }
            
            SetGenotypeProbOfGenoAt(j, i, 0, prob0);
            SetGenotypeProbOfGenoAt(j, i, 1, prob1);
            
            // choose the allele w/ higher prob
            int allele = 0;
            double probMax = GetGenotypeProbAt(j,i,0);
            if( probMax < GetGenotypeProbAt(j,i,1) )
            {
                probMax = GetGenotypeProbAt(j,i,1);
                allele = 1;
            }
            if( probMax < GetGenotypeProbAt(j,i,2) )
            {
                probMax = GetGenotypeProbAt(j,i,2);
                allele = 2;
            }
            SetGenotypeAt( j, i, allele );
        }
    }
    
cout << "Input matrix: ";
this->matTernary.Dump();
    
    return true;
}

bool ScistTernaryMat :: ReadFromFileTernaryProb( std::ifstream &infile, double &prob0, double &prob1 )
{
    // read in the prob of allele: (0.6,0.1) 0.6 means prob of 0 is 0.6 and prob of 1 is 0.1
    // assume prob of 0 + 1 + 2 = 1
    infile >> prob0 >> prob1;
    return true;
}

void ScistTernaryMat :: SetSize(int numSCs, int numSites)
{
    matTernary.SetSize(numSCs, numSites);
    
    matTernaryProbs.clear();
    matTernaryProbs.resize(numSCs);
    
    bool fNameInit = GetNumNames() > 0;
    
    for(int i=0; i<numSCs; ++i)
    {
        matTernaryProbs[i].resize( numSites );
        for(int s=0; s<numSites; ++s)
        {
            SetGenotypeProbOfGenoAt(i,s,0, 1.0);
            SetGenotypeProbOfGenoAt(i,s,1,0.0);
        }
        
        // by default, use the numericals, starting from one
        if( fNameInit == false )
        {
            string str = std::to_string(i+1);
            AddGenotypeName(str);
//cout << "Init name: " << str << endl;
        }
    }
}

int ScistTernaryMat :: GetGenotypeAt(int sc, int site) const
{
    return matTernary(sc,site);
}

int ScistTernaryMat :: GetAltGenotypeAt(int sc, int site) const
{
    YW_ASSERT_INFO(false, "Not supported1");
    return 1;
}

void ScistTernaryMat :: SetGenotypeAt(int sc, int site, int geno)
{
    matTernary(sc,site) = geno;
}

void ScistTernaryMat :: AddGenotypeAt(int sc, int site, int geno)
{
    // append the genotype into it
    int genoThis = GetGenotypeAt(sc, site);
    if( genoThis != geno )
    {
        SetGenotypeAt(sc, site, geno);
    }
}

double ScistTernaryMat :: GetGenotypeProbAllele0At(int sc, int site) const
{
    return GetGenotypeProbAt(sc, site, 0);
}

double ScistTernaryMat :: GetGenotypeProbAt(int sc, int site, int geno) const
{
    if( geno == 0 )
    {
        return this->matTernaryProbs[sc][site].first;
    }
    else if(geno == 1)
    {
        return this->matTernaryProbs[sc][site].second;
    }
    else
    {
        return 1.0-GetGenotypeProbAt(sc,site,0)-GetGenotypeProbAt(sc,site,1);
    }
}

void ScistTernaryMat :: SetGenotypeProbAt(int sc, int site, double prob)
{
    YW_ASSERT_INFO(false, "Not impelemented");
}

void ScistTernaryMat :: SetGenotypeProbOfGenoAt(int sc, int site, int geno, double prob)
{
    if( geno == 0 )
    {
        matTernaryProbs[sc][site].first = prob;
    }
    else if(geno ==1)
    {
        matTernaryProbs[sc][site].second = prob;
    }
    else
    {
        YW_ASSERT_INFO(false, "Cannot only set the homozygous mutant probility");
    }
}

void ScistTernaryMat :: FindMaximalCompatSites( const std::vector<double> &wtSites, std::vector< std::map<int, std::set<int> > > &listSetSitesCompat, int maxNumSets, const std::set<std::pair<int,int> > *pSetCompatPairs ) const
{
    YW_ASSERT_INFO(false, "Not implemented");
}

int ScistTernaryMat :: GetNumSites() const
{
    return matTernary.GetColNum();
}

int ScistTernaryMat :: GetNumHaps() const
{
    return matTernary.GetRowNum();
}

void ScistTernaryMat :: GetMutRowsHapAtSite(int site, std::set<int> &setRows) const
{
    //YW_ASSERT_INFO(false, "Not supported2");
    // for now, use both 1/2 rows
    GetRowsWithGenoAtSite(site, 1, setRows);
    set<int> setRows2;
    GetRowsWithGenoAtSite(site, 2, setRows2);
    UnionSets( setRows, setRows2 );
}

void ScistTernaryMat :: GetRowsWithGenoAtSite(int site, int geno, std::set<int> &setRows) const
{
    setRows.clear();
    for( int h=0; h<GetNumHaps(); ++h )
    {
        if( GetGenotypeAt(h, site) == geno)
        {
            setRows.insert(h);
        }
    }
}

double ScistTernaryMat :: GetScoreForGeno(int scIndex, int site, int genotype) const
{
    YW_ASSERT_INFO(false, "Not supported3");
    return 0.0;
}

bool ScistTernaryMat :: IsNoninformative(int site) const
{
    YW_ASSERT_INFO(false, "Not supported4");
    return false;
}

bool ScistTernaryMat :: IsCompatible(int s1, int s2) const
{
    YW_ASSERT_INFO(false, "Not supported5");
    return false;
}

std::string ScistTernaryMat :: ConsTree() const
{
    // construct phylogeny
    vector<int> rootZero;
    for(int i=0; i<GetNumSites(); ++i)
    {
        rootZero.push_back(0);
    }
    
    // construct binary matrix for distance computation
    BinaryMatrix binMat;
    ConsHapMatForDistCalc(binMat);
    
    PhylogenyTree phTree;
    phTree.SetRoot(rootZero);
    phTree.ConsOnBinMatrix(binMat);
    phTree.RemoveDegreeTwoNodes();
    
    // now assign leaf labels
    map<string,string> mapIdToLabels;
    for(int i=0; i<GetNumHaps(); ++i)
    {
        //cout << "i: " << i << ", name: " << this->genosInput.GetGenotypeName(i) << endl;
        string str = "(" + std::to_string(i) + ")";
        mapIdToLabels[str] = GetGenotypeName(i);
    }
    phTree.ReassignLeafLabels( mapIdToLabels );
    
    
    string res;
    phTree.ConsNewickSorted(res);
    //phTree.ConsNewick(res, false, 0.0, true);
    return res;

}

double ScistTernaryMat :: SumLogProbs() const
{
    YW_ASSERT_INFO(false, "Not impelemtned");
    return 0.0;
}

void ScistTernaryMat :: Dump() const
{
    ScistGenGenotypeMat::Dump();
    
    //
    cout << "Matrix: [" << GetNumHaps() << "," << GetNumSites() << "]" << endl;
    this->matTernary.Dump();

    cout << "Probabilities: \n";
    for(int i=0; i<(int)matTernaryProbs.size(); ++i)
    {
        for(int j=0; j<(int)matTernaryProbs[i].size(); ++j)
        {
            cout << "(" << matTernaryProbs[i][j].first << "," << matTernaryProbs[i][j].second << ")  ";
        }
        cout << endl;
    }
}

void ScistTernaryMat :: OutputImput(const string *pStrDesc) const
{
    //
    cout << "Lineages: ";
    for(int i=0; i<GetNumNames(); ++i)
    {
        cout << GetGenotypeName(i) << "  ";
    }
    cout << endl;
    if( pStrDesc != NULL )
    {
        cout << *pStrDesc << endl;
    }
    else
    {
        cout << "Imputed genotypes: \n";
    }
    for(int s=0; s<GetNumSites(); ++s)
    {
        cout << "Site " << setw(6) << s+1 << ":\t";
        
        for(int i=0; i<GetNumHaps(); ++i)
        {
            cout << GetGenotypeAt(i, s) << " ";
        }
        cout << endl;
    }
}

void ScistTernaryMat :: ConsHapMatForDistCalc( BinaryMatrix &matHaplotypes ) const
{
    matHaplotypes.SetSize( GetNumHaps(), 2*GetNumSites() );
    for( int r=0; r<GetNumHaps(); ++r )
    {
        for(int s=0; s<GetNumSites(); ++s)
        {
            int geno = GetGenotypeAt(r,s);
            int allele0 = 0, allele1 = 0;
            if( geno != 0 )
            {
                allele0 = 1;
            }
            if( geno == 2)
            {
                allele1 = 1;
            }
            matHaplotypes(r, 2*s) = allele0;
            matHaplotypes(r,2*s+1) = allele1;
        }
    }
}

bool ScistTernaryMat :: IsProbSignificant(double prob, double thresVal) const
{
    //
    const double probConst = 0.3333333;
    if( prob < probConst && prob >( probConst-thresVal/2 ) )
    {
        return false;
    }
    if( prob > probConst && prob < (probConst+thresVal/2) )
    {
        return false;
    }
    return true;
}
