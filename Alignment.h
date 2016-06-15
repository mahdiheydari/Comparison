
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

using namespace std;


class NW_Alignment {
private:

        char **tracebackArr;
        int **s;
        size_t maxSize;
        double matchScore;
        double mismatchPenalty;
        double gapPenalty;
        double max(double x, double y);
        double max(double x, double y, double z);
        void    traceback(string& s1,string& s2,char **traceback );
        void alocateMemory();
        void deAlocateMemory();
public:
        NW_Alignment();
        ~NW_Alignment();
        string getaligned(string a,string b,string c);
        bool isEqual(string a, string b);
        int findQualityDistance(string a, string b,string q) ;
        double get_similarity_per(string s1,string s2);
        double enhancedAlignment(string &s1, string &s2);

};
class AlignmentJan {

private:
        int maxDim;             // maximum dimensions of the sequences
        int maxIndel;           // maximum number of indels
        int match;              // match score
        int mismatch;           // mismatch penalty
        int gap;                // gap score
        int *M;                 // alignment matrix

        // void traceback(string& s1,string& s2,char **traceback );
        // void init();

public:
        /**
         * Default constructor
         * @param maxDim Maximum dimension of the sequences
         * @param maxGap Maximum number of insertion or deletions
         * @param match Match score
         * @param mismatch Mismatch penalty
         * @param gap Gap score
         */
        AlignmentJan(int maxDim, int maxIndel = 3, int match = 1,
                     int mismatch = -1, int gap = -3);

        /**
         * Destructor
         */
        ~AlignmentJan() {
                delete [] M;
        }

        int operator() (int i, int j) const {
                int k = max(i, j);
                int l = maxIndel - i + j;
                return M[k * (2 * maxIndel + 1) + l];
        }

        int& operator() (int i, int j) {
                int k = max(i, j);
                int l = maxIndel - i + j;
                return M[k * (2 * maxIndel + 1) + l];
        }

        /**
         * Perform the alignment between two sequences
         * @param s1 First string
         * @param s2 Second string
         * @return The alignment score (higher is better)
         */
        int align(const string &s1, const string &s2);

        /**
         * Print matrix to stdout
         */
        void printMatrix() const;

        /**
         * Print matrix to stdout
         */
        void printAlignment(const string &s1, const string &s2) const;
};
