//
//  ScistGenotype.hpp
//  
//
//  Created by Yufeng Wu on 5/25/18.
//
//

#ifndef ScistGenotype_hpp
#define ScistGenotype_hpp

#include <set>
#include <map>
#include <string>
#include <iostream>

#include "BinaryMatrix.h"
#include "GenotypeMatrix.h"

// *************************************************************************************
// genotypes: integer matrix

class ScistGenGenotypeMat
{
public:
    ScistGenGenotypeMat();
    virtual ~ScistGenGenotypeMat() {}
    virtual ScistGenGenotypeMat *CreateNewMat() const = 0;
    virtual ScistGenGenotypeMat *Copy() const = 0;
    virtual bool ReadFromFile( std::ifstream &infile, int numSites, int numSCs, bool fSiteName ) = 0;
    virtual void SetSize(int numSCs, int numSites) = 0;
    virtual void AddGenotypeName(const std::string &strNameIn) { listNames.push_back(strNameIn); }
    virtual void SetGenotypeName(int i, const std::string &strNameIn ) { listNames[i] = strNameIn; }
    virtual std::string GetGenotypeName(int i) const { return listNames[i]; }
    virtual void AddSiteName(const std::string &strNameIn) { listSiteNames.push_back(strNameIn); }
    virtual std::string GetSiteName(int i) const { return listSiteNames[i]; }
    virtual void GetSiteNamesAll(std::vector<std::string> &listSiteNamesOut) { listSiteNamesOut=listSiteNames; }
    virtual int GetNumGenos() const = 0;
    virtual int GetGenotypeAt(int sc, int site) const = 0;
    virtual int GetAltGenotypeAt(int sc, int site) const = 0;
    virtual void SetGenotypeAt(int sc, int site, int geno) = 0;
    virtual void AddGenotypeAt(int sc, int site, int geno) = 0;
    virtual double GetGenotypeProbAllele0At(int sc, int site) const = 0;
    virtual double GetGenotypeProbAt(int sc, int site, int geno) const = 0;
    double CalcSumLogProbAllele0All() const;   // calc sum of all log-prob of allele-0 over all cells and all sites
    double CalcMaxQForTaxaAt(int site, const std::set<int> &setTaxa, double baseProb, double maxQMin) const;
    double GetSumMaxQToLastSite(int site) const;
    double GetSumLogProbAllele0Site(int site) const;
    virtual void SetGenotypeProbAt(int sc, int site, double prob) = 0;
    virtual void SetGenotypeProbOfGenoAt(int sc, int site, int geno, double prob) = 0;
    virtual bool IsBinary() const = 0;
    virtual void FindMaximalCompatSites( const std::vector<double> &wtSites, std::vector< std::map<int, std::set<int> > > &listSetSitesCompat, int maxNumSets, const std::set<std::pair<int,int> > *pSetCompatPairs = NULL ) const = 0;
    virtual int GetNumSites() const = 0;
    virtual int GetNumHaps() const = 0;
    virtual void GetMutRowsHapAtSite(int site, std::set<int> &setRows) const = 0;
    virtual void GetRowsWithGenoAtSite(int site, int geno, std::set<int> &setRows) const = 0;
    virtual double GetScoreForGeno(int scIndex, int site, int genotype) const = 0;
    virtual bool IsNoninformative(int site) const = 0;
    virtual bool IsCompatible(int s1, int s2) const = 0;
    bool AreAllGenosCompatible() const;
    virtual ScistGenGenotypeMat *SubMatrix( const std::set<int> &setRows, const std::set<int> &setSites ) const;
    virtual void Dump() const;
    virtual void OutputImput(const string *pStrDesc = NULL) const = 0;
    virtual std::string ConsTree() const = 0;
    virtual double SumLogProbs() const = 0;
    virtual void GetColMultiplicityMap( std::vector<int> &listColMulti ) const = 0;
    virtual bool IsProbSignificant(double prob, double thresVal) const = 0;
    virtual void IncProbForAlleleBy(int sc, int site, int allele, double fac) = 0;
    std::string ConsNJTree() const;
    std::string ConsNJTreeZeroRoot( const std::set<std::set<int> > *psetCladesCurr = NULL ) const;
    std::string ConsNJTreeZeroRootMT( const std::set<std::set<int> > *psetCladesCurr = NULL ) const;
    std::string ConsNJTreeNoInc() const;
    std::string ConsPerfPhylogenyHeu() const;
    std::string ConsConstrainedUPGMATreeZeroRoot( const std::set< std::set<int> > &setMustClades ) const;
    std::string GetFileName() const { return inputFileName; }
    std::string GetImpGenosFileName() const {return GetFileName() + ".genos.imp"; }
    double IsProbAtCellPosSignificant(int sc, int site, double thresVal) const { return IsProbSignificant( GetGenotypeProbAt(sc, site, 0), thresVal ); }
    void SetSignificantThres(double thres) { thresSignifcant = thres; }
    void SetFileName(std::string &fn) { inputFileName = fn; }
    double CalcHammingDistBetwHaps(int h1, int h2) const;
    double CalcHammingDistBetwHaps2(int h1, int h2) const;
    void ConsCompatMap( std::set<std::pair<int,int> > &setCompatPairs ) const;
    int GetGenotypeNumOf(int geno) const;
    int FindCellByName(const std::string &strName) const;
    void ChangeGenosAtPositions(const std::set< std::pair<std::pair<int,int>, int> > &listChangedPlaces);
    static bool AreSitesCompatInMap(const std::set<std::pair<int,int> > &setCompatPairs, int s1, int s2 );
    void SetNumThreads(int t) { numThreads = t; }
    void SetDropoutProb(double pd) { probDropout = pd; }
    double CalcEntropyAt(int sc, int site) const;
    double CalcEntropySite(int site) const;
    void SetGenotypeCalled(bool f) { fGenoSet = f; }
    void ResetMaximalGenos();   // reset all genotypes to maximal genotypes
    void SetUncertainGenoPos(const std::set< std::pair<int,int> > &ss) { setUncertainPos = ss; }
    void Preprocess();
    void GetCladesFromCalledGenoes(std::set<std::set<int> > &setClades, std::set<int> &setUncertainClades) const;
    void InitSumAllele0Probs();
    void AdjustProbsFromCalledGenos();
    
protected:
    void TrimCliquesMaxDiff( std::set<std::set<int> > &listCliques, int maxToKeep ) const;
    void ResetNames() { listNames.clear(); }
    int GetNumNames() const { return listNames.size(); }
    double GetSignificanceThres() const { return thresSignifcant; }
    
private:
    void OrderSitesByNum1s(std::vector<int> &listSites) const;
    double CalcLogProbTwoSitesContain(int s1, int s2) const;
    double CalcLogProbTwoSitesDisjoint(int s1, int s2) const;
    void OptParseHapForMG( int hap, const std::vector<int> &listSites, const std::vector<int> &listMGEdges, std::vector<int> &hapOpt ) const;
    void OptParseHapsForMG( const std::vector<int> &listSites, const std::vector<int> &listMGEdges, BinaryMatrix &matHaps ) const;
    
    std::vector<std::string> listNames;
    std::vector<std::string> listSiteNames;
    std::string inputFileName;
    double thresSignifcant;
    std::vector<double> listSumLogprob0Sites;
    std::vector<double> listMaxQSumSitesAfter;
    double sumAllLogProb0;
    int numThreads;
    double probDropout;
    bool fGenoSet;      // are genotypes set (after calling)
    std::set< std::pair<int,int> > setUncertainPos;
};


// *************************************************************************************
// genotypes: binary matrix

class ScistHaplotypeMat : public ScistGenGenotypeMat
{
public:
    ScistHaplotypeMat();
    virtual ~ScistHaplotypeMat() {}
    virtual ScistGenGenotypeMat *Copy() const;
    virtual ScistGenGenotypeMat *CreateNewMat() const { return new ScistHaplotypeMat; }
    virtual bool ReadFromFile( std::ifstream &infile, int numSites, int numSCs, bool fSiteName );
    virtual void SetSize(int numSCs, int numSites);
    virtual int GetNumGenos() const { return 2; }
    virtual int GetGenotypeAt(int sc, int site) const;
    virtual int GetAltGenotypeAt(int sc, int site) const;
    virtual void SetGenotypeAt(int sc, int site, int geno);
    virtual void AddGenotypeAt(int sc, int site, int geno);
    virtual double GetGenotypeProbAllele0At(int sc, int site) const;
    virtual double GetGenotypeProbAt(int sc, int site, int geno) const {
        if( geno == 0 ) return GetGenotypeProbAllele0At(sc, site); else return 1.0-GetGenotypeProbAllele0At(sc,site);
    }
    virtual void SetGenotypeProbAt(int sc, int site, double prob);
    virtual void SetGenotypeProbOfGenoAt(int sc, int site, int geno, double prob);
    virtual bool IsBinary() const { return true; }
    virtual void FindMaximalCompatSites( const std::vector<double> &wtSites, std::vector< std::map<int, std::set<int> > > &listSetSitesCompat, int maxNumSets, const std::set<std::pair<int,int> > *pSetCompatPairs = NULL ) const;
    virtual int GetNumSites() const;
    virtual int GetNumHaps() const;
    virtual void GetMutRowsHapAtSite(int site, std::set<int> &setRows) const;
    virtual void GetRowsWithGenoAtSite(int site, int geno, std::set<int> &setRows) const;
    virtual double GetScoreForGeno(int scIndex, int site, int genotype) const;
    virtual bool IsNoninformative(int site) const;
    virtual bool IsCompatible(int s1, int s2) const;
    virtual std::string ConsTree() const;
    virtual double SumLogProbs() const;
    virtual void Dump() const;
    virtual void OutputImput(const string *pStrDesc = NULL) const;
    virtual void GetColMultiplicityMap( std::vector<int> &listColMulti ) const { matHaplotypes.GetColMultiplicityMap(listColMulti); }
    virtual bool IsProbSignificant(double prob, double thresVal) const;
    virtual void IncProbForAlleleBy(int sc, int site, int allele, double fac);
    BinaryMatrix &GetHapMat() { return matHaplotypes; }
    
private:
    bool ReadFromFileHapProb( std::ifstream &infile, double &prob0 );
    
    BinaryMatrix matHaplotypes;
    std::vector<std::vector<double> > matHaplotypesProb0;
};

// *************************************************************************************
// genotypes: ternary matrix

class ScistTernaryMat : public ScistGenGenotypeMat
{
public:
    ScistTernaryMat();
    virtual ~ScistTernaryMat() {}
    virtual ScistGenGenotypeMat *Copy() const;
    virtual ScistGenGenotypeMat *CreateNewMat() const { return new ScistTernaryMat; }
    virtual bool ReadFromFile( std::ifstream &infile, int numSites, int numSCs, bool fSiteName );
    virtual void SetSize(int numSCs, int numSites);
    virtual int GetNumGenos() const { return 3; }
    virtual int GetGenotypeAt(int sc, int site) const;
    virtual int GetAltGenotypeAt(int sc, int site) const;
    virtual void SetGenotypeAt(int sc, int site, int geno);
    virtual void AddGenotypeAt(int sc, int site, int geno);
    virtual double GetGenotypeProbAllele0At(int sc, int site) const;
    virtual double GetGenotypeProbAt(int sc, int site, int geno) const;
    virtual void SetGenotypeProbAt(int sc, int site, double prob);
    virtual void SetGenotypeProbOfGenoAt(int sc, int site, int geno, double prob);
    virtual bool IsBinary() const { return false; }
    virtual void FindMaximalCompatSites( const std::vector<double> &wtSites, std::vector< std::map<int, std::set<int> > > &listSetSitesCompat, int maxNumSets, const std::set<std::pair<int,int> > *pSetCompatPairs = NULL ) const;
    virtual int GetNumSites() const;
    virtual int GetNumHaps() const;
    virtual void GetMutRowsHapAtSite(int site, std::set<int> &setRows) const;
    virtual void GetRowsWithGenoAtSite(int site, int geno, std::set<int> &setRows) const;
    virtual double GetScoreForGeno(int scIndex, int site, int genotype) const;
    virtual bool IsNoninformative(int site) const;
    virtual bool IsCompatible(int s1, int s2) const;
    virtual std::string ConsTree() const;
    virtual double SumLogProbs() const;
    virtual void Dump() const;
    virtual void OutputImput(const string *pStrDesc = NULL) const;
    virtual void GetColMultiplicityMap( std::vector<int> &listColMulti ) const { matTernary.GetColMultiplicityMap(listColMulti); }
    virtual bool IsProbSignificant(double prob, double thresVal) const;
    virtual void IncProbForAlleleBy(int sc, int site, int allele, double fac) { std::cout << "Not implemented.\n"; exit(1); }
    
private:
    bool ReadFromFileTernaryProb( std::ifstream &infile, double &prob0, double &prob1 );
    void ConsHapMatForDistCalc( BinaryMatrix &matHaplotypes ) const;
    
    GenotypeMatrix matTernary;
    std::vector<std::vector<std::pair<double,double> > > matTernaryProbs;

};

#endif /* ScistGenotype_hpp */
