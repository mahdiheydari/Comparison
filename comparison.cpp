
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

bool Comparison::compareTwoEC() {
        ifstream initialStream, firstECStream,secondECStream;
        ofstream outinitialStream, outfirstECStream,outSecondECStream;
        AlignmentJan ali(100, 2, 1, -1, -3);
        if (! openFileStreamTwoEC(initialStream,firstECStream,secondECStream,outinitialStream,outfirstECStream,outSecondECStream))
                return false;
        readStructStr newRead;
        size_t numOfDiff=0;
        numOfReads=1;
        while (     fillNextReadTwoEC(newRead,  initialStream,firstECStream , secondECStream) )
        {
                numOfReads=numOfReads+1;
                if ( (int)numOfReads % 10000 == 0 ) {
                        //cout << "Processing read number " << numOfReads;
                        //cout.flush();

                }
                if (newRead.firstECRead != newRead.secondECRead){
                        numOfDiff++;
                        cout<<numOfDiff<<endl;
                        cout<<"initial\t\t\t" <<newRead.erroneousRead<<endl;
                        cout<<"Brownie\t\t"  <<ali.align(newRead.erroneousRead, newRead.firstECRead)  << "\t" <<newRead.firstECRead<<endl;
                        cout<<"efBrownie\t" <<ali.align(newRead.erroneousRead, newRead.secondECRead) << "\t" <<newRead.secondECRead<<endl;

                        ali.align(newRead.erroneousRead, newRead.secondECRead);
                        outSecondECStream<<newRead.strID<<endl;
                        outSecondECStream<<newRead.secondECRead<<endl;
                        outSecondECStream<<"+"<<endl;
                        outSecondECStream<<newRead.qProfile<<endl;

                        outfirstECStream<<newRead.strID<<endl;
                        outfirstECStream<<newRead.firstECRead<<endl;
                        outfirstECStream<<"+"<<endl;
                        outfirstECStream<<newRead.qProfile<<endl;

                }
                //compareReadsTwoEC(newRead, )


        }
        closeFileStreamTwoEC(initialStream,firstECStream,secondECStream,outinitialStream,outfirstECStream,outSecondECStream);
        return true;

}
bool Comparison::closeFileStreamTwoEC(ifstream& initialStream,ifstream& firstECStream,ifstream& secondECStream,
                                     ofstream &outinitialStream,ofstream &outfirstECStream,ofstream& outSecondECStream)
{
        initialStream.close();
        firstECStream.close();
        secondECStream.close();
        outinitialStream.close();
        outfirstECStream.close();
        outSecondECStream.close();
        return true;

}
bool Comparison::openFileStreamTwoEC(ifstream& initialStream,ifstream& firstECStream,ifstream& secondECStream,
                                     ofstream &outinitialStream,ofstream &outfirstECStream,ofstream& outSecondECStream)
{
        if(!initialStream.good()){
            std::cerr << "Error opening file"<< erroneousReadsFileName << std::endl;
            return false;
        }
        if(!firstECStream.good()){
            std::cerr << "Error opening file"<< firstECFileName << std::endl;
            return false;
        }
        if(!secondECStream.good()){
            std::cerr << "Error opening file"<< secondECFileName << std::endl;
            return false;
        }
        initialStream.open(erroneousReadsFileName.c_str());
        firstECStream.open(firstECFileName.c_str());
        secondECStream.open(secondECFileName.c_str());

        outinitialStream.open(erroneousReadsFileOutName.c_str(), ios::out);
        outfirstECStream.open(firstECFileOutName.c_str(), ios::out);
        outSecondECStream.open(secondECFileOutName.c_str(), ios::out);
        return true;
}
bool Comparison::fillNextReadTwoEC(readStructStr & readInfo, ifstream& initialStream,ifstream& firstECStream , ifstream& secondECStream)
{
        if (!initialStream.eof()&&!initialStream.eof()&& !secondECStream.eof())
        {
                string temp="";


                getline( initialStream,readInfo.strID);
                getline(initialStream,readInfo.erroneousRead);
                getline(initialStream,readInfo.orientation);
                getline(initialStream,readInfo.qProfile);

                getline( firstECStream,temp);
                getline(firstECStream,readInfo.firstECRead);
                getline(firstECStream,temp);
                getline(firstECStream,temp);

                getline( secondECStream,temp);
                getline(secondECStream,readInfo.secondECRead);
                getline(secondECStream,temp);
                getline(secondECStream,temp);

                return true;
        }

        return false;
}


bool Comparison::validateCorrectionResult() {
        ifstream perfectReadsStream, correctedReadsStream,erroneousReadsFileStream;
        ofstream notCorrectedReadF, worseStream, notCorrectedFastq;
        if (! openFileStream(erroneousReadsFileStream,perfectReadsStream,correctedReadsStream,notCorrectedReadF,worseStream,notCorrectedFastq))
                return false;
        readStructStr newRead;
        while (     fillNextRead(newRead,erroneousReadsFileStream,perfectReadsStream ,correctedReadsStream) )
        {
                numOfReads=numOfReads+1;
                if ( (int)numOfReads % 500 == 0 ) {
                        cout << "Processing read number " << numOfReads
                        << " gain is: " << gain<<"\r";
                        cout.flush();

                }
                compareReads( newRead,  worseStream, notCorrectedReadF, notCorrectedFastq);

        }
        closeFileStream(erroneousReadsFileStream,perfectReadsStream,correctedReadsStream,notCorrectedReadF,worseStream,notCorrectedFastq);
        writeReport();
        return true;
}





/*
 *
 * @return true if the corrected read dosn't have any error, otherwise it return false and it should write the
 * remaining errors in seperate filese.
 *
 */
bool Comparison::compareReads(readStructStr &newRead, ofstream &worseStream,ofstream &notCorrectedReadF,ofstream  &notCorrectedFastq){

        int exsistErrorInRead=0,qualityDistance=0,initialQualityDistance=0 , allErrorINRead=0;;
        string erroneousRead =newRead.erroneousRead;
        string correctedRead=newRead.correctedRead;
        string perfectRead=newRead.perfectRead;
        a.enhancedAlignment(perfectRead,erroneousRead);
        a.enhancedAlignment(perfectRead,correctedRead);
        a.enhancedAlignment(perfectRead,erroneousRead);
        initialQualityDistance=a.findQualityDistance(perfectRead, erroneousRead,newRead.qProfile);
        if (initialQualityDistance>maximumQualityDistance) {
                maximumQualityDistance=initialQualityDistance;
        }
        qualityDistance=a.findQualityDistance(correctedRead,perfectRead,newRead.qProfile);
        sumOfQualityDistances=sumOfQualityDistances+qualityDistance;

        //update the value of TP, TN, FP, FN for bases and also for full recovery
        updateStatistic(correctedRead ,erroneousRead, perfectRead, allErrorINRead, exsistErrorInRead);
        if (exsistErrorInRead==0)
                return true;
        writeInFileErrors( allErrorINRead, exsistErrorInRead, worseStream,
                                   notCorrectedReadF, notCorrectedFastq,
                                   newRead,  erroneousRead,  perfectRead,  correctedRead);

        return false;
}

/*
 * this procedure update the value of TP, TN, FP and FN
 * and for full recovery statistic
 *
 */
void Comparison::updateStatistic(string  &correctedRead , string  & erroneousRead, string  & perfectRead, int &allErrorINRead, int &exsistErrorInRead){
        char perfectChar,modifiedChar,erroneousChar;
        for (int i=0; i<correctedRead.length()&& i<erroneousRead.length()&& i<perfectRead.length(); i++) {
                modifiedChar=correctedRead[i];
                erroneousChar=erroneousRead[i];
                perfectChar=perfectRead[i];
                if (erroneousChar!=perfectChar) {
                        allErrorINRead++;
                        if (modifiedChar==perfectChar) {
                                truePositive++;
                        }
                        else {
                                falseNegative++;
                                exsistErrorInRead++;
                        }
                } else {
                        if(modifiedChar==perfectChar) {
                                trueNegative++;
                        }
                        else {
                                falsePositive++;
                                exsistErrorInRead++;
                        }
                }
        }

        if(allErrorINRead>maximumNumOfError) {
                maximumNumOfError=allErrorINRead;
        }
        //full recovery updataSta
        if (allErrorINRead!=0){
                if (exsistErrorInRead==0)
                {
                        truePositiveFullRec++;
                        fullyRecoveredReads++;
                }
                else
                        falseNegativeFullRec++;
        }
        else {
                if (exsistErrorInRead>0)
                        falsePositiveFullRec++;
                else{
                        trueNegativeFullRec++;
                        fullyRecoveredReads++;
                }
        }
        numberOfAllErrors=numberOfAllErrors+allErrorINRead;
        gain=100*(double)(truePositive-falsePositive)/(double)(truePositive+falseNegative);
        gainFR=100*(double)(truePositiveFullRec-falsePositiveFullRec)/(double)(truePositiveFullRec+falseNegativeFullRec);
        correctionAVG=100*((double) truePositive/(double)(truePositive+falseNegative));

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
        string WorseName="Worse";
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

                erroneousReadsFileStream>>readInfo.strID;
                erroneousReadsFileStream>>readInfo.erroneousRead;
                erroneousReadsFileStream>>readInfo.orientation;
                erroneousReadsFileStream>>readInfo.qProfile;

                getline(perfectReadsStream, temp);
                getline(perfectReadsStream,  readInfo.perfectRead);


                getline(correctedReadsStream,temp);
                getline(correctedReadsStream,readInfo.correctedRead);
                getline(correctedReadsStream,temp);
                getline(correctedReadsStream,temp);

                return true;
        }

        return false;
}

void Comparison::writeReport(){
        std::cout << std::fixed;
        cout <<endl<< "<<<Report for reads>>>" << endl;
        cout << "----------------------------------------------------\n" << endl;
        cout<<"Number of reads for comparison: "<<numOfReads<<endl;
        cout<<"Maximum number Of Errors in one read is: "<< maximumNumOfError <<", and maximum quality distances is: "<<maximumQualityDistance <<endl;


        cout <<endl<< "<<<Alignment based report>>>" << endl;
        cout << "----------------------------------------------------" << endl;
        cout<<"Sum of  errors in the Data Set: "<<(truePositive+falseNegative)<<", (" <<(double)(truePositive+falseNegative)/(double)numOfReads<<" per read)"<<endl;


        cout << endl << "<<<The evaluation report based on base pairs>>>"<<endl;
        cout << "----------------------------------------------------" << endl;
        cout << "Among " << numberOfAllErrors <<" of Errors "<<truePositive<<" number of them are corrected which is: ("<<std::setprecision(2)<<correctionAVG<<"%)"<<endl ;
        cout << "     TP:"<<truePositive<<"      TN:"<<trueNegative<<"      FP:"<<falsePositive <<"    FN:"<<falseNegative<<endl;
        cout << "    The Gain value percentage is: ("<<std::setprecision(2)<<gain<<"%)" <<endl;

        cout<<endl<<"<<<The evaluation report based full read recovery>>>"<<endl;
        cout << "----------------------------------------------------" << endl;
        cout<<" Among "<<numOfReads<<" of reads "<<fullyRecoveredReads<<", number of them are fully recovered which is: ("<<((double)fullyRecoveredReads/(double)numOfReads)*100<<")%"<<endl ;
        cout<<"     TP: "<<truePositiveFullRec<<"      TN: "<<trueNegativeFullRec<<"      FP: "<<falsePositiveFullRec<<"     FN: "<<falseNegativeFullRec<<endl;
        cout <<"    The Gain value percentage for full recovery of reads is: ("<<gainFR<<"%)" <<endl;

        cout<<endl<< "<<<Quality based reports>>>"<<endl;
        cout << "----------------------------------------------------" << endl;
        cout<<"The sum of quality distances between original reads and corrected reads is: "<<sumOfQualityDistances<<" ,("<<double( sumOfQualityDistances)/(double)(numOfReads)<<"   per read) ,and the number of changes is:   "<<numOfAllChanges<<", (" <<(double)numOfAllChanges/(double) numOfReads <<" per read)" <<endl;
}
