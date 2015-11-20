
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
