
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdio.h>
#include <cstring>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <omp.h>
#include <assert.h>
#include <iomanip>      // std::setprecision
#include <algorithm>
#include "comparison.h"
using namespace std;


class NW_Alignment {
private:

        double matchScore;
        double mismatchPenalty;
        double gapPenalty;
        double max(double x, double y);
        double max(double x, double y, double z);
        void    traceback(string& s1,string& s2,char **traceback );
        void init();
public:
        NW_Alignment();
        string getaligned(string a,string b,string c);
        bool isEqual(string a, string b);
        int findQualityDistance(string a, string b,string q) ;
        int findDifference(string a, string b) ;
        double get_similarity_per(string s1,string s2);
        void printMatrix(char **s, int m, int n);
        double enhancedAlignment(string &s1, string &s2);
        int findAlignedDifference(string a, string b);
};

NW_Alignment::NW_Alignment() {
        init();
}
void NW_Alignment::init()
{
        matchScore=1;
        mismatchPenalty=1;
        gapPenalty=20;
}

double NW_Alignment::max(double x, double y)
{
        return x > y ? x : y;
}

double NW_Alignment::max(double x, double y, double z)
{
        return x > y ? max(x, z) : max(y, z);
}
double NW_Alignment::get_similarity_per(string s1,string s2){
        double fullScore= enhancedAlignment(s1,s1);
        cout<<fullScore<<endl;
        double simScore=enhancedAlignment(s1,s2);
        if (s1==s2)
                return 100;
        return (simScore/fullScore)*100;
}
void NW_Alignment::traceback(string& s1,string& s2, char **traceback ){
        string news1="";
        string news2="";
        int i=s1.length();
        int j=s2.length();
        while( i > 0 || j > 0 )
        {
                if (i>0 && j>0 && traceback[i][j]=='\\') {
                        news1 += s1[ i-1 ] ;
                        news2 += s2[ j-1 ] ;
                        i-- ;
                        j-- ;

                } else {
                        if (i>0 && traceback[i][j]=='|') {
                                news2 += '-' ;
                                news1 += s1[ i-1 ] ;
                                i-- ;

                        } else {
                                if (j>0 ) {
                                        news2 += s2[ j-1 ] ;
                                        news1 += '-' ;
                                        j-- ;

                                } else {
                                        if (i>0) {
                                                news2 += '-' ;
                                                news1 += s1[ i-1 ] ;
                                                i-- ;
                                        }
                                }
                        }
                }

        }

        reverse( news1.begin(), news1.end() );
        reverse( news2.begin(), news2.end() );
        s1=news1;
        s2=news2;

}

void NW_Alignment::printMatrix(char **s, int m, int n){
        for (int i=0;i<m;i++ ){
                for (int j=0;j<n;j++)
                        cout<< s[i][j]<<"      ";
                cout<<"\n";
        }
}

double NW_Alignment::enhancedAlignment(string &s1, string &s2){
        if (s1==s2)
                return s1.length()*matchScore;
        int n = s1.length() + 1, m = s2.length() + 1, i, j;
        int d=3;
        if (abs(m-n)>d || m<d || n<d)
                return 0;
        char **tracebackArr=new char*[n];
        for(int i = 0; i < n; i++)
        {
                tracebackArr[i] = new char[m];
        }

        int **s = new int*[n];
        for(int i = 0; i < n; i++)
        {
                s[i] = new int[m];
        }
        int p=0;
        for (int i=0; i<n; i++) {
                s[i][p]=i*-1;
                if (i>d)
                        p++;
        }
        p=0;
        for (int j = 0; j < m; j++)
        {
                s[p][ j] = j*-1;
                if (j>d)
                        p++;
        }
        for (int i = 1; i <= n-1; i++)
        {
                int jIndexMin=i-d >0?i-d:1;
                int jIndexMax=i+d<m?i+d:m-1;
                for (int j = jIndexMin; j <= jIndexMax; j++)//for (int j = 1; j <= m-1; j++) // for (int j = jIndexMin; j <= i+(d); j++)//
                {
                        int scroeDiag = 0;

                        if(s1[i-1]==s2[j-1]|| s1[i-1]=='N'||s2[i-1]=='N')
                                scroeDiag = s[i - 1][ j - 1] + matchScore;   //match
                                else
                                        scroeDiag = s[i - 1][ j - 1]  -mismatchPenalty; //substitution
                                        int scroeLeft = s[i][ j - 1] - gapPenalty; //insert
                                        int scroeUp = s[i - 1][ j] - gapPenalty;  //delete
                                        int maxScore =  max(scroeDiag,scroeLeft,scroeUp);//  Math.Max(Math.Max(scroeDiag, scroeLeft), scroeUp);
                                s[i][ j] = maxScore;
                        if (scroeDiag==maxScore) {
                                tracebackArr[i][j]='\\';
                        } else {
                                if (maxScore==scroeLeft) {
                                        tracebackArr[i][j]='-';
                                } else {
                                        if(maxScore==scroeUp) {
                                                tracebackArr[i][j]='|';
                                        }
                                }
                        }

                }
        }
        traceback(s1, s2,tracebackArr);
        int result=s[n-1][m-1];

        for(int i = 0; i < n; i++)
        {
                delete[] s[i];
        }
        delete[] s;
        for(int i = 0; i < n; i++)
        {
                delete[] tracebackArr[i];
        }
        delete[] tracebackArr;
        return result;

}

int NW_Alignment::findAlignedDifference(string a, string b) {
        NW_Alignment al;
        al.enhancedAlignment(a,b);
        int i=0;
        int j=0;
        int d=0;
        while(i<a.length()&&j<b.length()) {
                if (a[i]!=b[j]) {
                        d++;
                }
                i++;
                j++;
        }
        return d;

}
int NW_Alignment::findDifference(string a, string b) {

        int i=0;
        int j=0;
        int d=0;
        while(i<a.length()&&j<b.length()) {
                if (a[i]!=b[j]) {
                        d++;
                }
                i++;
                j++;
        }
        return d;

}

int NW_Alignment::findQualityDistance(string a, string b,string q) {
        int i=0;
        int j=0;
        int d=0;
        while(i<a.length()&&j<b.length()) {
                if (a[i]!=b[j] ) {
                        d=d+int(q[i]);
                }
                i++;
                j++;
        }
        return d;

}

bool NW_Alignment::isEqual(string a, string b){
        int i=0;
        int j=0;
        int d=0;
        bool equal=true;
        while(i<a.length()&&j<b.length()) {
                if (a[i]!=b[j]&& a[i]!='N') {
                        equal=false;
                }
                i++;
                j++;
        }
        return equal;
}
string NW_Alignment::getaligned(string a,string b,string c) {
        int i=0;
        string result="";
        while (i<a.length()) {
                if ((a[i]==b[i])&& (a[i]==c[i])) {
                        result=result+'*';
                }
                else {
                        result=result+a[i];
                }
                i++;
        }
        return result;
}





bool Comparison::validateCorrectionResult() {
        NW_Alignment a;
        ifstream perfectReadsStream, correctedReadsStream,erroneousReadsFileStream;
        ofstream notCorrectedReadF, worseStream, notCorrectedFastq;
        if (! openFileStream(erroneousReadsFileStream,perfectReadsStream,correctedReadsStream,notCorrectedReadF,worseStream,notCorrectedFastq))
                return false;
        readStructStr newRead;
        while (     fillNextRead(newRead,erroneousReadsFileStream,perfectReadsStream ,correctedReadsStream) )
        {
                int exsistErrorInRead=0,qualityDistance=0,numOfChangesInRead=0,initialQualityDistance=0,numberOfErrorInRead=0 , allErrorINRead=0;;
                numOfReads=numOfReads+1;
                if ( (int)numOfReads % 10000 == 0 ) {
                        cout << "Processing read number " << numOfReads
                        << " gain is: " << gain<<"\r";
                        cout.flush();

                }
                string erroneousRead =newRead.erroneousRead;
                string correctedRead=newRead.correctedRead;
                string perfectRead=newRead.perfectRead;

                a.enhancedAlignment(perfectRead,correctedRead);
                a.enhancedAlignment(correctedRead,erroneousRead);

                initialQualityDistance=a.findQualityDistance(perfectRead, erroneousRead,newRead.qProfile);
                if (initialQualityDistance>maximumQualityDistance) {
                        maximumQualityDistance=initialQualityDistance;
                }
                qualityDistance=a.findQualityDistance(correctedRead,perfectRead,newRead.qProfile);
                sumOfQualityDistances=sumOfQualityDistances+qualityDistance;
                
                //update the value of TP, TN, FP, FN for bases and also for full recovery
                updataStatistic(correctedRead ,erroneousRead, perfectRead, allErrorINRead, exsistErrorInRead);

                if (exsistErrorInRead==0)
                        continue;
                writeInFileErrors( allErrorINRead, exsistErrorInRead, worseStream,
                                    notCorrectedReadF, notCorrectedFastq,
                                    newRead,  erroneousRead,  perfectRead,  correctedRead);

        }
        closeFileStream(erroneousReadsFileStream,perfectReadsStream,correctedReadsStream,notCorrectedReadF,worseStream,notCorrectedFastq);
        writeReport();
        return true;
}



/*
 * this procedure update the value of TP, TN, FP and FN
 * and for full recovery statistic
 *
 */
void Comparison::updataStatistic(string  &correctedRead , string  & erroneousRead, string  & perfectRead, int &allErrorINRead, int &exsistErrorInRead){
        char correctChar,modifiedChar,erroneousChar;
        for (int i=0; i<correctedRead.length(); i++) {
                modifiedChar=correctedRead[i];
                erroneousChar=erroneousRead[i];
                correctChar=perfectRead[i];
                if (erroneousChar!=correctChar) {
                        allErrorINRead++;
                        if (modifiedChar==correctChar) {
                                truePositive++;
                        }
                        else {
                                falseNegative++;
                                exsistErrorInRead++;
                        }
                } else {
                        if(modifiedChar==correctChar) {
                                trueNegative++;
                        }
                        else {
                                falsePositive++;
                        }
                }
        }
        if(allErrorINRead>maximumNumOfError) {
                maximumNumOfError=allErrorINRead;
        }
        //full recovery updataSta
        if (exsistErrorInRead!=0){
                if (allErrorINRead==0)
                        falsePositiveFullRec++;
                else
                        falseNegativeFullRec++;
        }
        else {
                fullyRecoveredReads++;
                if (allErrorINRead>0)
                        truePositiveFullRec++;
                else
                        trueNegativeFullRec++;
        }
        gain=100*(double)(truePositive-falsePositive)/(double)(truePositive+falseNegative);
        gainFR=100*(double)(truePositiveFullRec-falsePositiveFullRec)/(double)(truePositiveFullRec+falseNegativeFullRec);
        correctionAVG=((double) truePositive/(double)(truePositive+falseNegative));


}
/*write reads which still have error in a schematic way to understand where the errors are still are. (notCorrectedReadF)
 *write reads which you increas the errors in (worseStream)
 *write  reads which you would not able to correct them in the initial format of the reads to process later (notCorrectedFastq)
 *
 */
void Comparison::writeInFileErrors(int allErrorINRead,int exsistErrorInRead, ofstream &worseStream,
                                   ofstream &notCorrectedReadF, ofstream &notCorrectedFastq,
                                   readStructStr newRead, string &erroneousRead, string &perfectRead, string &correctedRead)
{
        NW_Alignment a;
        notCorrectedReadF<<newRead.strID <<"    |"<< "read number:      "<<numOfReads <<"        num of all Errors:      "<<allErrorINRead<<"    number of remaining Errors:     "<< exsistErrorInRead<<endl;
        string alignedError=a.getaligned(erroneousRead,perfectRead,correctedRead);

        notCorrectedReadF<<"E:"<<erroneousRead<<endl;
        notCorrectedReadF<<"P:"<<perfectRead<<endl;
        notCorrectedReadF<<"C:"<<correctedRead<<endl;
        notCorrectedReadF<<"A:"<<alignedError<<endl;
        notCorrectedReadF<<"A:"<<newRead.qProfile<<endl;


        if ( exsistErrorInRead>allErrorINRead) {
                worseStream<<newRead.strID <<"       |"<< "read number:      "<<numOfReads <<"        num of all Errors:      "<<allErrorINRead<<"    number of remaining Errors:     "<< exsistErrorInRead<<endl;
                worseStream<<"E:"<<newRead.erroneousRead<<endl;
                worseStream<<"P:"<<newRead.perfectRead<<endl;
                worseStream<<"C:"<<newRead.correctedRead<<endl;
                worseStream<<"A:"<<alignedError<<endl;
        }
        notCorrectedFastq<<newRead.strID<<endl;
        notCorrectedFastq<<newRead.erroneousRead<<endl;
        notCorrectedFastq<<'+'+newRead.strID<<endl;
        notCorrectedFastq<<newRead.qProfile<<endl;

}
/*
 * open file streams for writing and reading information
 *
 * the first three files are important and thould be opened
 */
bool Comparison::openFileStream(ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream,
                                ifstream& correctedReadsStream, ofstream &notCorrectedReadF,
                                ofstream &worseReadF,ofstream& notCorrectedFastq)
{
        if(!perfectReadsStream.good()){
            std::cerr << "Error opening file"<< perfectReadsFileName << std::endl;
            return false;
        }
        if(!correctedReadsStream.good()){
            std::cerr << "Error opening file"<< correctedReadsFileName << std::endl;
            return false;
        }
        if(!erroneousReadsFileStream.good()){
            std::cerr << "Error opening file"<< erroneousReadsFileName << std::endl;
            return false;
        }
        perfectReadsStream.open(perfectReadsFileName.c_str());
        correctedReadsStream.open(correctedReadsFileName.c_str());
        erroneousReadsFileStream.open(erroneousReadsFileName.c_str());
        notCorrectedReadF.open(notCorrectedReadFileName.c_str(), ios::out);
        string WorseName=notCorrectedReadFileName+"Worse.fastq";
        worseReadF.open(WorseName.c_str(), ios::out);
        string fastqName=notCorrectedReadFileName+".fastq";
        notCorrectedFastq.open(fastqName.c_str(), ios::out);
        return true;
}
/*
 * close all open streams
 *
 */
bool Comparison::closeFileStream(ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream,
                                ifstream& correctedReadsStream, ofstream &notCorrectedReadF,
                                ofstream &worseReadF,ofstream& notCorrectedFastq)
{
        erroneousReadsFileStream.close();
        perfectReadsStream.close();
        correctedReadsStream.close();
        notCorrectedFastq.close();
        notCorrectedReadF.close();
        worseReadF.close();
        return true;

}
/*
 * rads information from file and fill it in a read
 * each read from erroneous fastq file has strID , erroneous read and quality Profile
 * in corrected read file the corrected read equal to the erroneous one exsis
 * in perfect file the correct or perfect read equal to the erroneous one is.
 */
bool Comparison::fillNextRead(readStructStr & readInfo, ifstream& erroneousReadsFileStream,ifstream& perfectReadsStream , ifstream& correctedReadsStream)
{
        if (!erroneousReadsFileStream.eof()&&!perfectReadsStream.eof()&& !correctedReadsStream.eof())
        {
                string temp="";

                getline(erroneousReadsFileStream,readInfo.strID);
                getline(erroneousReadsFileStream,readInfo.erroneousRead);
                getline(erroneousReadsFileStream,readInfo.orientation);
                getline(erroneousReadsFileStream,readInfo.qProfile);
                getline(perfectReadsStream, temp);
                getline(perfectReadsStream,readInfo.perfectRead);
                getline(correctedReadsStream, temp);
                getline(correctedReadsStream,readInfo.correctedRead);
                getline(correctedReadsStream, temp);
                getline(correctedReadsStream, temp);
                return true;
        }

        return false;
}

void Comparison::writeReport(){
        std::cout << std::fixed;
        cout <<endl<< "<<<Report for reads>>>" << endl;
        cout << "----------------------------------------------------\n" << endl;
        cout<<"Number of reads for comparison: "<<(int)numOfReads<<endl;
        cout<<"Maximum number Of Errors in one read is: "<< (int)maximumNumOfError <<", and maximum quality distances is: "<<(int)maximumQualityDistance <<endl;


        cout <<endl<< "<<<Alignment based report>>>" << endl;
        cout << "----------------------------------------------------" << endl;
        cout<<"Sum of  errors in the Data Set: "<<(int)numberOfAllErrors<<", (" <<(double)numberOfAllErrors/(double)numOfReads<<" per read)"<<endl;
        cout<<"Sum of alignment distances between corrected and perfect read is: "<<(int)sumOfAlignmentDist<<", (" <<(double)sumOfAlignmentDist/(double)numOfReads<<" per read)"<<endl;

        cout <<endl<< "<<<The evaluation report based on base pairs>>>"<<endl;
        cout << "----------------------------------------------------" << endl;
        cout<<"Among "<<(int)(truePositive+falseNegative)<<" of Errors "<<(int)truePositive<<" number of them are corrected which is: ("<<std::setprecision(2)<<correctionAVG<<"%)"<<endl ;
        cout<<"     TP:"<<(int)truePositive<<"      TN:"<<(int)trueNegative<<"      FP:"<<(int)falsePositive <<"    FN:"<<(int)falseNegative<<endl;
        cout <<"    The Gain value percentage is: ("<<std::setprecision(2)<<gain<<"%)" <<endl;

        cout<<endl<<"<<<The evaluation report based full read recovery>>>"<<endl;
        cout << "----------------------------------------------------" << endl;
        cout<<" Among "<<(int)numOfReads<<" of reads "<<(int)fullyRecoveredReads<<", number of them are fully recovered which is: ("<<((double)fullyRecoveredReads/(double)numOfReads)*100<<")%"<<endl ;
        cout<<"     TP: "<<(int)truePositiveFullRec<<"      TN: "<<(int)trueNegativeFullRec<<"      FP: "<<(int)falsePositiveFullRec<<"     FN: "<<(int)falseNegativeFullRec<<endl;
        cout <<"    The Gain value percentage for full recovery of reads is: ("<<gainFR<<"%)" <<endl;

        cout<<endl<< "<<<Quality based reports>>>"<<endl;
        cout << "----------------------------------------------------" << endl;
        cout<<"The sum of quality distances between original reads and corrected reads is: "<<(int)sumOfQualityDistances<<" ,("<<double( sumOfQualityDistances)/(double)(numOfReads)<<"   per read) ,and the number of changes is:   "<<numOfAllChanges<<", (" <<(double)numOfAllChanges/(double) numOfReads <<" per read)" <<endl;
}
