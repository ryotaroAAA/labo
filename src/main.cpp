
#include "../lib/Common.h"
#include "../lib/Preseter.h"
#include "../lib/Decoder.h"
#include "../lib/Encoder.h"
#include "../lib/Channel.h"
#include "../lib/Analysor.h"
#include "../lib/Logger.h"
#include "../lib/Performance.h"
#include <sys/time.h>
//#include "../lib/ps.h"
//#include <bits/stdc++.h>

void calcBER(){
    Performance performance;
    Decoder decoder;
    Encoder encoder;
    Logger logger;
    Analysor analysor;

    EXP_MODE em = Params::get_exp_mode();
    MID_MODE mm = Params::get_m_mode();
    string ename;
    switch (em) {
        case NORMAL: ename = ""; break;
        case PUNC: ename = "punc"; break;
        case MID:
            switch (mm) {
                case MID_BLUTE: ename = "mid_blute"; break;
                case MID_ADOR: ename = "mid_ador"; break;
                case MID_AOR: ename = "mid_aor"; break;
                case MID_DOR: ename = "mid_dor"; break;
                case MID_AOB: ename = "mid_aob"; break;
                case MID_DOB: ename = "mid_dob"; break;
                case MID_AOV: ename = "mid_aov"; break;
                case MID_DOV: ename = "mid_dov"; break;
            }
            break;
        case QUP: ename = "qup"; break;
        case WANG: ename = "wang"; break;
        case M_QUP: ename = "m_qup"; break;
        case M_WANG: ename = "m_wang"; break;
        case VALERIO_P: ename = "valerio_p"; break;
        case VALERIO_S: ename = "valerio_s"; break;
        case M_VALERIO_P: ename = "m_valerio_p"; break;
        case M_VALERIO_S: ename = "m_valerio_s"; break;
    }
    cout << "[[" << ename << "]]" << endl;

    time_t now = time(NULL);
    struct tm *pnow = localtime(&now);
    CHANNEL_TYPE channel_type = Params::get_s();

    stringstream fn;
    if (channel_type == BEC || channel_type == BSC) {
        fn << "log/"
             << to_string(pnow->tm_year + 1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
             << "/N" << to_string(Params::get_N())
             << "M" << to_string(Params::get_M())
             << "MN" << to_string(Params::get_MN())
             << "e" << Params::get_e()
             << "BNum" << to_string(Params::get_blockNum())
             << "UPBn" << to_string(Params::get_upperBlockErrorNum())
             << "c" << (Params::get_s() ? "BSC" : "BEC") << ":" << ename;
        Params::set_rvbDir(fn.str());
    } else if (channel_type == AWGN) {
        fn << "log/"
           << to_string(pnow->tm_year+1900) + to_string(pnow->tm_mon + 1) + to_string(pnow->tm_mday)
           << "/N=" + to_string(Params::get_N())
           << "M" << to_string(Params::get_M())
           << "MN" << to_string(Params::get_MN())
           << "sdiv" << Params::get_e()
           << "c=AWGN"
           << "MNum=" + to_string(Params::get_monteNum())
           << "BNum=" + to_string(Params::get_blockNum())
           << "upperBn=" + to_string(Params::get_upperBlockErrorNum()) +  ":" + ename;
        Params::set_rvbDir(fn.str());
    }

    performance.startTimer();
    analysor.calcBlockErrorRate_BP();
//    analysor.calcBlockErrorRate();
    performance.stopTimer();

    performance.outHMS();
    cout << Params::get_rvbDir() << endl;
}

inline string get_current_directory()
{
    char cwd[1024];
    getcwd(cwd, sizeof(cwd));
    return cwd;
}

bool checkFileExistence(const std::string& str)
{
    std::ifstream ifs(str);
    return ifs.is_open();
}

//vector<double> gauss(vector< vector<double> > A) {
//    int n = ((int)A.size() >= (int)A[0].size) ? (int)A.size() : (int)A[0].size;
//
//    for (int i=0; i<n; i++) {
//        // Search for maximum in this column
//        double maxEl = abs(A[i][i]);
//        int maxRow = i;
//        for (int k=i+1; k<n; k++) {
//            if (abs(A[k][i]) > maxEl) {
//                maxEl = abs(A[k][i]);
//                maxRow = k;
//            }
//        }
//
//        // Swap maximum row with current row (column by column)
//        for (int k=i; k<n+1;k++) {
//            double tmp = A[maxRow][k];
//            A[maxRow][k] = A[i][k];
//            A[i][k] = tmp;
//        }
//
//        // Make all rows below this one 0 in current column
//        for (int k=i+1; k<n; k++) {
//            double c = -A[k][i]/A[i][i];
//            for (int j=i; j<n+1; j++) {
//                if (i==j) {
//                    A[k][j] = 0;
//                } else {
//                    A[k][j] += c * A[i][j];
//                }
//            }
//        }
//    }
//
//    // Solve equation Ax=b for an upper triangular matrix A
//    vector<double> x(n);
//    for (int i=n-1; i>=0; i--) {
//        x[i] = A[i][n]/A[i][i];
//        for (int k=i-1;k>=0; k--) {
//            A[k][n] -= A[k][i] * x[i];
//        }
//    }
//    return x;
//}

/* function for exchanging two rows of
   a matrix */
void swap(vector<vector<int> > &mat, int row1, int row2, int col) {
    for (int i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}

/* function for finding rank of matrix */
int rankOfMatrix(vector<vector<int> > &mat) {
    int n = Params::get_N();
    int rsize = mat.size(); //たて
    int csize = mat[0].size(); //よこ
    int minsize = (rsize > csize) ? csize : rsize;
    int rank = csize;

    for (int row = 0; row < rank; row++)
    {
        // Before we visit current row 'row', we make
        // sure that mat[row][0],....mat[row][row-1]
        // are 0.
        int flg = 0;
        if(row > minsize-1) {
            flg = 0;
        } else {
            flg = mat[row][row];
        }

        // Diagonal element is not zero
        if (flg)
        {
            for (int col = 0; col < rsize; col++)
            {
                if (col != row)
                {
                    // This makes all entries of current
                    // column as 0 except entry 'mat[row][row]'
                    double mult = (double)mat[col][row] /
                                  mat[row][row];
                    for (int i = 0; i < rank; i++)
                        mat[col][i] = (int)(mat[col][i] + mult * mat[row][i])%2;
                }
            }
        }

            // Diagonal element is already zero. Two cases
            // arise:
            // 1) If there is a row below it with non-zero
            //    entry, then swap this row with that row
            //    and process that row
            // 2) If all elements in current column below
            //    mat[r][row] are 0, then remvoe this column
            //    by swapping it with last column and
            //    reducing number of columns by 1.
        else {
            bool reduce = true;

            /* Find the non-zero element in current
                column  */
            for (int i = row + 1; i < rsize;  i++)
            {
                // Swap the row with non-zero element
                // with this row.
                if (mat[i][row])
                {
                    swap(mat, row, i, rank);
                    reduce = false;
                    break ;
                }
            }

            // If we did not find any row with non-zero
            // element in current columnm, then all
            // values in this column are 0.
            if (reduce)
            {
                // Reduce number of columns
                rank--;

                // Copy the last column here
                for (int i = 0; i < rsize; i ++)
                    mat[i][row] = mat[i][rank];
            }

            // Process this row again
            row--;
        }

        // Uncomment these lines to see intermediate results
        // display(mat, R, C);
        // printf("\n");
    }
    return rank;
}

inline vector<int>makeTreeIndex(int temp_n, int n){
    vector<int> ret(n);
    vector<int> temp_ret;
    if (n == 1) {
        ret[0] = 1;
    } else {
        for (int i = 0; i < n/2 ; i++) {
            ret[2*i] = i+1;
            ret[2*i+1] = i+1 + n/2;
        }
    }
    int loop = temp_n/n;

    for (int i = 0; i < loop; i++) {
        for (int j = 0; j < n; ++j) {
            temp_ret.push_back(ret[j] + i*n);
        }
    }
    return temp_ret;
}

void addRow(vector<int> &p, vector<vector<int> > &mat ){
    int n = Params::get_N();
    int csize = n*(log2(n)+1);
    for (int i = 0; i < p.size(); i++) {
        vector<int> temp;
        for (int j = 0; j < csize; j++) {
            if(p[i] == j){
                temp.push_back(1);
            } else {
                temp.push_back(0);
            }
        }
        mat.push_back(temp);
    }
}

void calcRank(){
    int n = Params::get_N();
    int rsize = n*(log2(n));
    int csize = n*(log2(n)+1);
    vector<vector<int> > mat(rsize, vector<int>(csize, 0));

    for (int i = 0; i < rsize; ++i) {
        for (int j = 0; j < csize; ++j) {
            if( i == j ){
                mat[i][j] = 1;
            }
            if( i%2 == 0 && i+1 == j){
                mat[i][j] = 1;
            }
        }
        int div = i/n;
        int rest = i%n;
        int temp_n = n/(pow(2,div));
        vector<int> tree = makeTreeIndex(n, temp_n);
        mat[i][n*(div+1)+tree[rest]-1] = 1;
    }
    vector<int> p;
    p.push_back(0);

    addRow(p, mat);
    for (int i = 0; i < mat.size(); ++i) {
        for (int j = 0; j < mat[0].size(); ++j) {
            cout << mat[i][j] <<  " ";
        }
        cout << endl;
    }

    Common::bar();

    vector<vector<int> > temp_mat = mat;
    int rank = rankOfMatrix(temp_mat);
    for (int i = 0; i < temp_mat.size(); ++i) {
        for (int j = 0; j < temp_mat[0].size(); ++j) {
            cout << temp_mat[i][j] <<  " ";
        }
        cout << endl;

    }

//    int g_rsize = csize - temp_mat.size();
    int g_rsize = csize-temp_mat.size();
    int g_csize = temp_mat[0].size();

    vector<vector<int> > g_mat(g_rsize, vector<int>(g_csize, 0));
    for (int i = 0; i < g_rsize; i++) {
        for (int j = 0; j < g_csize; j++) {
            if(j < temp_mat.size()){
                g_mat[i][j] = temp_mat[j][csize-1-i];
            } else {
                if( i == j-temp_mat.size()){
                    g_mat[i][j] = 1;
                }
            }
        }
    }

    Common::bar();

    for (int i = 0; i < g_mat.size(); ++i) {
        for (int j = 0; j < g_mat[0].size(); ++j) {
            cout << g_mat[i][j] <<  " ";
        }
        cout << endl;
    }

    Common::bar();
    cout << "rank : " << rank << endl;
    cout << "Rate : " << (double)(csize-rank)/Params::get_N() << endl;
    cout << "Rate_1 : " << (double)(csize-rank)/Params::get_N() << endl;
    cout << "Rate_2 : " << (double)(rankOfMatrix(g_mat))/Params::get_N() << endl;
}


//mid>normal
//awgn なら e=0.8
//enum EXP_MODE{NORMAL, PUNC, QUP, WANG, MID, M_WANG, M_QUP, VALERIO_P, VALERIO_S, M_VALERIO_P, M_VALERIO_S};
int main(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec + tv.tv_usec);
    init_genrand(tv.tv_sec + tv.tv_usec);

    Params::set_e(0.11);
    Params::set_N(16);
    Params::set_K(1);
    Params::set_M(0);
    Params::set_MN(0);
    Params::set_s(BSC);
    Params::set_is_outlog(true);
    Params::set_decode_mode(BP);
    Params::set_monteNum(1);
    Params::set_rp(20);
    Params::set_Bloop(1000);
    Params::set_is_calc_bloop(true);
    Params::set_blockNum(10000);
    Params::set_upperBlockErrorNum(10000);

    //IDかBD
//    Params::set_m_mode(MID_ADOR);
//    Params::set_m_mode(MID_DOR);
//    Params::set_m_mode(MID_AOR);
//    Params::set_m_mode(MID_DOB);
//    Params::set_m_mode(MID_AOB);
//    Params::set_m_mode(MID_DOV);
//    Params::set_m_mode(MID_AOV);

//    Params::set_exp_mode(NORMAL);
//    Params::set_exp_mode(QUP);
    Params::set_exp_mode(MID);
//    Params::set_exp_mode(WANG);
//    Params::set_exp_mode(VALERIO_P);
//    Params::set_exp_mode(VALERIO_S);
//    Params::set_exp_mode(M_WANG);
//    Params::set_exp_mode(M_QUP);
//    Params::set_exp_mode(M_WANG);
//    Params::set_exp_mode(M_VALERIO_P);
//    Params::set_exp_mode(M_VALERIO_S);

    calcBER();
//    calcRank();
//    Analysor::calcBlockErrorRate();
    return 0;
}