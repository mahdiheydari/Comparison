
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
};
class Comparison{
        string erroneousReadsFileName;
        string perfectReadsFileName;
        string correctedReadsFileName;
        string notCorrectedReadFileName;
        double fullyRecoveredReads;
        double sumOfQualityDistances;
        long allReads;
        double truePositive;
        double trueNegative;
        double falsePositive;
        double falseNegative;
        int truePositiveFullRec;
        int trueNegativeFullRec;
        int falsePositiveFullRec;
        int falseNegativeFullRec;
        double gain;
        double gainFR;
        double correctionAVG;
        double maximumNumOfError;
        double maximumQualityDistance;
        double sumOfAlignmentDist;
        double numOfAllChanges;
        int numOfReads;
        double numberOfAllErrors;
        double existErrorNum;
        double allErrorNum;
        void writeReport();
        void updataStatistic(string  &correctedRead , string  & erroneousRead, string  & perfectRead, int &allErrorINRead, int &exsistErrorInRead);
        bool openFileStream(ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream,
                                ifstream& correctedReadsStream, ofstream &notCorrectedReadF,
                                ofstream &worseReadF,ofstream& notCorrectedFastq);
        bool closeFileStream(ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream,
                                ifstream& correctedReadsStream, ofstream &notCorrectedReadF,
                                ofstream &worseReadF,ofstream& notCorrectedFastq);
        void writeInFileErrors(int allErrorINRead,int exsistErrorInRead, ofstream &worseStream,
                                   ofstream &notCorrectedReadF, ofstream &notCorrectedFastq,
                                   readStructStr newRead, string &erroneousRead, string &perfectRead, string &correctedRead);
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
        bool validateCorrectionResult();

};
