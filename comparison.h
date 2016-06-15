#include "Alignment.h"
#include <iostream>
using namespace std;
struct readStructStr
{
        int intID;
        string strID;
        string erroneousRead;
        string qProfile ;
        string orientation;
        string perfectRead;
        string correctedRead;
        string firstECRead;
        string secondECRead;
};
class Comparison{
        string erroneousReadsFileName;
        string perfectReadsFileName;
        string correctedReadsFileName;
        string notCorrectedReadFileName;
        double fullyRecoveredReads;
        size_t sumOfQualityDistances;
        NW_Alignment a;
        long allReads;
        size_t truePositive;
        size_t trueNegative;
        size_t falsePositive;
        size_t falseNegative;
        size_t truePositiveFullRec;
        size_t trueNegativeFullRec;
        size_t falsePositiveFullRec;
        size_t falseNegativeFullRec;
        double gain;
        double gainFR;
        double correctionAVG;
        size_t maximumNumOfError;
        size_t maximumQualityDistance;
        double sumOfAlignmentDist;
        double numOfAllChanges;
        size_t numOfReads;
        size_t numberOfAllErrors;
        double existErrorNum;
        double allErrorNum;
        void writeReport();
        void updateStatistic(string  &correctedRead , string  & erroneousRead, string  & perfectRead, int &allErrorINRead, int &exsistErrorInRead);
        bool openFileStream(ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream,
                                ifstream& correctedReadsStream, ofstream &notCorrectedReadF,
                                ofstream &worseReadF,ofstream& notCorrectedFastq);
        bool closeFileStream(ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream,
                                ifstream& correctedReadsStream, ofstream &notCorrectedReadF,
                                ofstream &worseReadF,ofstream& notCorrectedFastq);
        void writeInFileErrors(int allErrorINRead,int exsistErrorInRead, ofstream &worseStream,
                                   ofstream &notCorrectedReadF, ofstream &notCorrectedFastq,
                                   readStructStr newRead, string &erroneousRead, string &perfectRead, string &correctedRead);
        bool compareReads(readStructStr &newRead, ofstream &worseStream,ofstream &notCorrectedReadF,ofstream  &notCorrectedFastq);
        bool fillNextRead(readStructStr & readInfo, ifstream &erroneousReadsFileStream,ifstream &perfectReadsStream , ifstream &correctedReadsStream);
public:
        Comparison(string erroneousReadsFileNameInp,string perfectReadsFileNameInp,string correctedReadsFileNameInp,string notCorrectedReadFileNameInp)
        :erroneousReadsFileName(erroneousReadsFileNameInp),perfectReadsFileName(perfectReadsFileNameInp),correctedReadsFileName(correctedReadsFileNameInp),notCorrectedReadFileName(notCorrectedReadFileNameInp)
        {
                fullyRecoveredReads=0;
                sumOfQualityDistances=0;
                allReads=0;
                truePositive=0;
                trueNegative=0;
                falsePositive=0;
                falseNegative=0;
                truePositiveFullRec=0;
                trueNegativeFullRec=0;
                falsePositiveFullRec=0;
                falseNegativeFullRec=0;
                maximumNumOfError=0;
                maximumQualityDistance=0;
                sumOfAlignmentDist=0;
                numOfAllChanges=0;
                numOfReads=0;
                numberOfAllErrors=0;
                existErrorNum=0;
                allErrorNum=0;
                gain=0;
                gainFR=0;
                correctionAVG=0;
        }


        string firstECFileName;
        string secondECFileName;

        string erroneousReadsFileOutName;
        string firstECFileOutName;
        string secondECFileOutName;
        Comparison(string erroneousReadsFileNameInp,string firstECFileNameInp,string secondECFileNameInp):
        erroneousReadsFileName(erroneousReadsFileNameInp),firstECFileName(firstECFileNameInp),secondECFileName(secondECFileNameInp)
        {
                erroneousReadsFileOutName="out."+erroneousReadsFileName;
                firstECFileOutName="out."+firstECFileName;
                secondECFileOutName="out."+secondECFileName;

        }
        bool closeFileStreamTwoEC(ifstream& initialStream,ifstream& firstECStream,ifstream& secondECStream,
                                     ofstream &outinitialStream,ofstream &outfirstECStream,ofstream& outSecondECStream);

        bool openFileStreamTwoEC(ifstream& initialStream,ifstream& firstECStream,ifstream& secondECStream,
                                     ofstream &outinitialStream,ofstream &outfirstECStream,ofstream& outSecondECStream);
        bool fillNextReadTwoEC(readStructStr & readInfo, ifstream& initialStream,ifstream& firstECStream , ifstream& secondECStream);
        bool validateCorrectionResult();
        bool compareTwoEC() ;
};
