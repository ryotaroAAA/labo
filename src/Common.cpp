#include "../lib/Common.h"

void Common::bar(){
    cout << "======================" << endl;
}

double Common::gauss_dist(double x, double mean, double sigma){
    return exp(-1.0*pow(x - mean, 2)/(2 * sigma)) / (sqrt(2*M_PI*sigma));
}

bool Common::containVal(int value, vector<int> m_array) {
    auto it = find(begin(m_array), m_array.end(), value);
    size_t index = distance(begin(m_array), it);
    if(index == m_array.size()) {
        return false;
    }
    return true;
}

bool Common::sort_less(const ass_arr& left,const ass_arr& right){
    return left.second < right.second;
}

bool Common::sort_greater(const ass_arr& left,const ass_arr& right){
    return left.second > right.second;
}

int Common::compare_asc(const void *x, const void *y) {
    const double* a=(const double*)x;
    const double* b=(const double*)y;
    if(*a>*b)return 1;
    if(*a<*b)return -1;
    return 0;
}

int Common::compare_desc(const void *x, const void *y) {
    const double* a=(const double*)x;
    const double* b=(const double*)y;
    if(*a>*b)return -1;
    if(*a<*b)return 1;
    return 0;
}

void Common::dispArray(double *x){
    for(int i = 0; i < Params::get_N(); i++){
        printf("x[%d] = %f\n",i, x[i]);
    }
}
void Common::dispArray(int n, double *x){
    for(int i = 0; i < n; i++){
        printf("x[%d] = %f\n",i, x[i]);
    }
}

void Common::dispArray(int n, int *x){
    for(int i = 0; i < n; i++){
        printf("x[%d] = %d\n",i, x[i]);
    }
}

void Common::pp(vector<int> &a){
    bar();
    cout << "{";
    int count = 1;
    for(auto var : a){
        cout << var;
        if(count != a.size()) cout << ",";
        count++;
    }
    cout << "}" << endl;
}

void Common::rpp(vector<int> &a){
    bar();
    cout << "{";
    int count = 1;
    for (int i = 0; i < a.size(); ++i) {
        cout << a[a.size()-1-i];
        if(count != a.size()) cout << ",";
        count++;
    }
    cout << "}" << endl;
}

void Common::dpp(vector<int> &a){
    bar();
    cout << "{";
    int count = 1;
    for(auto var : a){
        cout << var-1;
        if(count != a.size()) cout << ",";
        count++;
    }
    cout << "}" << endl;
}

void Common::drpp(vector<int> &a){
    bar();
    cout << "{";
    int count = 1;
    for (int i = 0; i < a.size(); ++i) {
        cout << a[a.size()-1-i]-1;
        if(count != a.size()) cout << ",";
        count++;
    }
    cout << "}" << endl;
}

void Common::pp(vector<double> &a){
    bar();
    cout << "{";
    int count = 1;
    for(auto var : a){
        cout << var;
        if(count != a.size()) cout << ",";
        count++;
    }
    cout << "}" << endl;
}

void Common::outputArray(double *x){
    string filename = "/Users/ryotaro/labo/symmetric_capacity";
    ofstream w_file;
    w_file.open(filename, ios::out);
    for (int i = 0; i<Params::get_N(); i++)
    {
        w_file << i << " " << x[i]<< endl;
    }
}

vector<int> Common::index_o(int n, vector<int> &x){
    vector<int> ret(n/2);
    for(int i = 0; i < n; i++){
        if (i % 2 == 0) {
            ret[i/2] = x[i];
        }
    }
    return ret;
}

vector<int> Common::index_e(int n, vector<int> &x){
    vector<int> ret(n/2);
    for(int i = 0; i < n; i++){
        if (i % 2 == 1) {
            ret[i/2] = x[i];
        }
    }
    return ret;
}

vector<int> Common::retBinary(int n, vector<int> &x) {
    vector<int> ret(n);
    for(int i = 0; i < n; i++){
        ret[i] = x[i]%2;
    }
    return ret;
}

int Common::ithIndexDesc(int i, double *array, double *descArray){
    for (int j=0; j<Params::get_N(); j++) {
        if(descArray[i] == array[j]){
            return j;
        }
    }
    return 0;
}


bool Common::containNumInArray(int i, int n, vector<int> &array){
    for(auto val: array ){
        if(i == val){
            return true;
        }
    }
    return false;
}

double Common::get_rate(vector<int> &A, vector<int> &S){
    double rate;
    double tmp_rate;
    int n = Params::get_N();
    int k = Params::get_K();
    int m = Params::get_M();
    int mn = Params::get_MN();

    switch (Params::get_exp_mode()){
        case NORMAL:
            rate = (double)k/n;
            break;
        case PUNC:
        case QUP:
        case WANG:
        case VALERIO_P:
            rate = (double)k/(n-m);
            break;
        case VALERIO_S:
            rate = (double)(k-m)/(n-m);
            break;
        case MID:
            rate = (double)k/(n+mn);
            break;
        case M_WANG:
        case M_QUP:
        case M_VALERIO_P:
            rate = (double)k/(n-m+mn);
            break;
        case M_VALERIO_S:
            rate = (double)(k-m)/(n-m+mn);
            break;
        default:
            rate = 0.0;
            break;
    }
//    tmp_rate = Common::calcRank(A,S);
//    cout << "new::" << tmp_rate << endl;
//    cout << "old::" << rate << endl;
    return rate;
}

bool Common::is_mid_send(){
    switch (Params::get_exp_mode()){
        case NORMAL:
        case PUNC:
        case QUP:
        case WANG:
        case VALERIO_P:
        case VALERIO_S:
            return false;
        case MID:
        case M_WANG:
        case M_QUP:
        case M_VALERIO_P:
        case M_VALERIO_S:
            return true;
        default:
            break;
    }
}

/* function for exchanging two rows of
   a matrix */
void Common::swap(vector<vector<int> > &mat, int row1, int row2, int col) {
    for (int i = 0; i < col; i++)
    {
        int temp = mat[row1][i];
        mat[row1][i] = mat[row2][i];
        mat[row2][i] = temp;
    }
}

/* function for finding rank of matrix */
int Common::rankOfMatrix(vector<vector<int> > &mat) {
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

inline vector<int> Common::makeTreeIndex(int temp_n, int n){
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

void Common::addRow(vector<int> &p, vector<vector<int> > &mat ){
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

double Common::calcRank(vector<int> &A, vector<int> &S){
    double rate = 0.0;
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
        vector<int> tree = Common::makeTreeIndex(n, temp_n);
        mat[i][n*(div+1)+tree[rest]-1] = 1;
    }

    //ショートンビット集合を入れる
    vector<int> p;
    EXP_MODE em = Params::get_exp_mode();
    for (int i = 0; i < A.size(); i++) {
        p.push_back(A[i]);
    }

    if(em == WANG || em == VALERIO_S || M_WANG || M_VALERIO_S){
        for (int i = 0; i < S.size(); i++) {
            p.push_back(S[i] + 2 * n * log2(n));
        }
    }

    Common::addRow(p, mat);
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

    //生成行列作成
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

    vector<int> T;

    vector<vector<int> > g_t(g_rsize, vector<int>(T.size(), 0));
    for (int i = 0; i < g_mat.size(); i++) {
        int k = 0;
        for (int j = 0; j < g_mat[0].size(); j++) {
            if(Common::containVal(j,T)){
                g_t[i][k] = g_mat[i][j];
                k++;
            }
        }
    }

    for (int i = 0; i < g_t.size(); ++i) {
        for (int j = 0; j < g_t[0].size(); ++j) {
            cout << g_t[i][j] << " ";
        }
        cout << endl;
    }

//    cout << "rank : " << rank << endl;
//    cout << "Rate : " << (double)(csize-rank)/Params::get_N() << endl;
//    cout << "Rate_1 : " << (double)(csize-rank)/Params::get_N() << endl;
//    cout << "Rate_2 : " << (double)(rankOfMatrix(g_mat))/Params::get_N() << endl;
    rate = (double)(rankOfMatrix(g_t))/Params::get_N();
    return rate;
}