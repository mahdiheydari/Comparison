
#include"Alignment.h"

NW_Alignment::NW_Alignment():maxSize(100) {
        matchScore=1;
        mismatchPenalty=1;
        gapPenalty=5;
        alocateMemory();
}

NW_Alignment::~NW_Alignment()
{

        deAlocateMemory();
}

void NW_Alignment::deAlocateMemory()
{
        for(int i = 0; i < maxSize; i++) {
                delete[] s[i];
        }
        delete[] s;
        for(int i = 0; i < maxSize; i++) {
                delete[] tracebackArr[i];
        }
        delete[] tracebackArr;}
void NW_Alignment::alocateMemory()
{
        tracebackArr=new char*[maxSize];
        for(int i = 0; i < maxSize; i++)
        {
                tracebackArr[i] = new char[maxSize];
        }
        s = new int*[maxSize];
        for(int i = 0; i < maxSize; i++)
        {
                s[i] = new int[maxSize];
        }
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


double NW_Alignment::enhancedAlignment(string &s1, string &s2){
       if (s1.length()>=maxSize || s2.length()>=maxSize){
                deAlocateMemory();
                maxSize=max(s1.length(),s2.length())+2;
                alocateMemory();
        }
        if (s1==s2)
                return s1.length()*matchScore;

        int n = s1.length() + 1, m = s2.length() + 1;
        int d=3;
        if (abs(m-n)>d || m<d || n<d)
                return 0;

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
                                if(maxScore==scroeUp) {
                                        tracebackArr[i][j]='|';
                                }
                                else {
                                        if (maxScore==scroeLeft) {
                                                tracebackArr[i][j]='-';
                                        }
                                }
                        }
                }
        }
        traceback(s1, s2,tracebackArr);
        int result=s[n-1][m-1];

        return result;
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


/*
 * compare three string lines base per base if any char is different in any of these
 * lines put star in the mapped line
 *
 */
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
