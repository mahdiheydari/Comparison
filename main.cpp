#include <iostream>
#include "comparison.h"
#include <ctime>
using namespace std;

int main(int argc, char ** args)
{
        cout << "Welcome to comparison software\n" << endl;
        cout << "----------------------------------------------------\n" << endl;
        if (argc < 4) {
                cout << "ERROR: Wrong number of arguments!\n";
                cout << "Usage: ./comparison <erroneousReads.fq> <perfectReads.fa> <correctedReads.fq> (opt: notCorrectedReadFileName)" << endl << flush;
                cout << "Order of reads in all input files should be the same.\nPlease look at the readme file for more information." << endl << flush;
                return 0;
        }
        string erroneousReadsFileName=args[1];
        string perfectReadsFileName=args[2];
        string correctedReadsFileName=args[3];
        string notCorrectedReadFileName;
        if (argc==5) {
                notCorrectedReadFileName=args[4];
        }else{
                notCorrectedReadFileName="notCorrectedRead";
        }
        clock_t begin = clock();
        Comparison compare(erroneousReadsFileName,perfectReadsFileName,correctedReadsFileName,notCorrectedReadFileName);

        compare.validateCorrectionResult();
        clock_t end = clock();
        double elapsed_secs = double(end-begin)/(double) (CLOCKS_PER_SEC*60);

        cout << "\n----------------------------------------------------\n" << endl;
        cout<<endl<<"comparison finished in: "<<elapsed_secs <<" min" <<endl;
        cout << "\n----------------------------------------------------\n" << endl;

        cout << "Exiting... bye!" << endl << endl;
        return 0;
}
