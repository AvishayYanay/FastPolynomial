
#include "recursive.h"
#include "utils.h"

using namespace std;
using namespace NTL;
using namespace chrono;

#define LEFT(X) (2*X+1)
#define RIGHT(X) (2*X+2)


/* A recursive function to build the tree of polynomials
 * (assumming a complete binary tree => size = 2*#leafs-1
 'tree_size' is the number of nodes (including leaves) in the tree = 2*(degree+1)-1 = 2*degree+1
 'root' is the index of the subtree in the array 'tree'
 */
void build_tree (ZZ_pX* tree, ZZ_p* points, unsigned int root, unsigned int tree_size) {
    // halting condition
    if(LEFT(root)>=tree_size) {
        unsigned int point_index = root-(tree_size-1)/2;
        //setting the polynomial to be x-m where m is points[point_index]
        ZZ_p negated;
        NTL::negate(negated, points[point_index]);
        SetCoeff(tree[root], 0, negated);
        SetCoeff(tree[root], 1, 1);
        return;
    }

    build_tree(tree, points, LEFT(root), tree_size);
    build_tree(tree, points, RIGHT(root), tree_size);
    tree[root] = tree[LEFT(root)]*tree[RIGHT(root)];
}

void test_tree(ZZ_pX &final_polynomial, ZZ_p *points, unsigned int npoints) {
    ZZ_p result;
    bool error = false;
    for (unsigned int i=0; i<npoints; i++) {
        result = eval(final_polynomial, points[i]);
        if (0!=result) {
            cout << "FATAL ERROR: polynomials tree is incorrect!" << endl;
            error = true;
            break;
        }
    }
    if (!error)
        cout << "polynomials tree is correct." << endl;
}

/*
 * P - the polynomial to evaluate_zp_iterative
 * tree - the subproduct tree over the x points that we want to recursive_evaluate_zp
 * root - the current subtree
 * tree size - the size of a complete tree is 2*n-1 where n is the number of leafs
 * results - the evaluation result over the x's (that are represented by the tree)
 */
void recursive_evaluate_zp(ZZ_pX &P, ZZ_pX *tree, unsigned int root, unsigned int tree_size, ZZ_p *results) {
    // halting condition
    if(LEFT(root)>=tree_size) {
        ZZ_pX R = P%tree[root];
        if(deg(R)>0)
            cout << "ERROR: R should be constant...";
        unsigned int result_index = root-(tree_size-1)/2;
        results[result_index] = coeff(R, 0);
        return;
    }

    ZZ_pX R = P%tree[root];
    recursive_evaluate_zp(R, tree, LEFT(root), tree_size, results);
    recursive_evaluate_zp(R, tree, RIGHT(root), tree_size, results);
}

void test_evaluate_zp_recursive(ZZ_pX &P, ZZ_p *points, ZZ_p *results, unsigned int npoints) {
    bool error = false;
    for (unsigned int i = 0; i < npoints; i++) {
        ZZ_p y = eval(P, points[i]);
        if (y != results[i]) {
            cout << "y=" << y << " and results[i]=" << results[i] << endl;
            error = true;
        }
    }
    if (error)
        cout << "ERROR: evaluation results do not match real evaluation!" << endl;
    else
        cout << "All evaluation results computed correctly!" << endl;
}


void multipoint_evaluate_zp_recursive(long degree, ZZ_pX &P, ZZ_p *X, ZZ_p *Y)
{
    // we want to recursive_evaluate_zp P on 'degree+1' values.
    ZZ_pX* p_tree = new ZZ_pX[degree*2+1];
    steady_clock::time_point begin1 = steady_clock::now();
    build_tree (p_tree, X, 0, degree*2+1);
    steady_clock::time_point end1 = steady_clock::now();
//    test_tree_zp_iterative(p_tree[0], x, degree+1);

    steady_clock::time_point begin2 = steady_clock::now();
    recursive_evaluate_zp(P, p_tree, 0, degree * 2 + 1, Y);
    chrono::steady_clock::time_point end2 = steady_clock::now();

    cout << "Building tree: " << duration_cast<milliseconds>(end1 - begin1).count() << " ms" << endl;
    cout << "Evaluating points: " << duration_cast<milliseconds>(end2 - begin2).count() << " ms" << endl;
    cout << "Total: " << duration_cast<milliseconds>(end1 - begin1).count()+ duration_cast<milliseconds>(end2 - begin2).count() << " ms" << endl;

    delete[] p_tree;
}

//void test_multipoint_eval_zp(ZZ prime, long degree)
//{
//    // init underlying prime field
//    ZZ_p::init(ZZ(prime));
//
//    // the given polynomial
//    ZZ_pX P;
//    random(P, degree+1);
//    SetCoeff(P,degree,random_ZZ_p());
//
//    // evaluation points:
//    ZZ_p* x = new ZZ_p[degree+1];
//    ZZ_p* y = new ZZ_p[degree+1];
//
//    for(unsigned int i=0;i<=degree; i++) {
//        random(x[i]);
//    }
//
//    multipoint_evaluate_zp_iterative(P, x, y, degree);
//}


/*
 * expects an "empty" polynomial 'resultP'
 */
void recursive_interpolate_zp(ZZ_pX& resultP, unsigned int root, ZZ_p* x, ZZ_p* y, ZZ_p* a, ZZ_pX* M, unsigned int tree_size)
{
    // halting condition
    if(LEFT(root)>=tree_size) {
        unsigned int y_index = root-(tree_size-1)/2;
        ZZ_p inv_a;
        inv(inv_a,a[y_index]); // inv_a = 1/a
        SetCoeff(resultP, 0, y[y_index]*inv_a);
        return;
    }

    ZZ_pX leftP, rightP;
    recursive_interpolate_zp(leftP, LEFT(root), x, y, a, M, tree_size);
    recursive_interpolate_zp(rightP, RIGHT(root), x, y, a, M, tree_size);

    resultP = leftP * M[RIGHT(root)] + rightP * M[LEFT(root)] ;
}

/*
 * We follow the algorithm and notation as in Moneck & Borodin '73
 */
void interpolate_zp(long degree, ZZ_p* X, ZZ_p* Y, ZZ_pX& resultP)
{
    system_clock::time_point begin[4];
    system_clock::time_point end[4];

    //we first build the tree of the super moduli
    ZZ_pX* M = new ZZ_pX[degree*2+1];
    begin[0]= system_clock::now();
    build_tree(M,X,0, degree*2+1);
    end[0] = system_clock::now();
//    test_tree_zp_iterative(M[0], x, degree+1);

    //we construct a preconditioned global structure for the a_k for all 1<=k<=(degree+1)
    ZZ_p* a = new ZZ_p[degree+1];
    ZZ_pX d;
    begin[1] = system_clock::now();
    diff(d, M[0]);
    end[1] = system_clock::now();

    //recursive_evaluate_zp d(x) to obtain the results in the array a
    begin[2] = system_clock::now();
    recursive_evaluate_zp(d, M, 0, degree * 2 + 1, a);
    end[2] = system_clock::now();

    //now we can apply the recursive formula
    begin[3] = system_clock::now();
    recursive_interpolate_zp(resultP, 0, X, Y, a, M, degree*2+1);
    end[3] = system_clock::now();

    cout << " -- Recursive --" << endl<< endl;
    cout << "Building tree: " << duration_cast<milliseconds>(end[0] - begin[0]).count() << " ms" << endl;
    cout << "Differentiate: " << duration_cast<milliseconds>(end[1] - begin[1]).count() << " ms" << endl;
    cout << "Evaluate diff: " << duration_cast<milliseconds>(end[2] - begin[2]).count() << " ms" << endl;
    cout << "Interpolation: " << duration_cast<milliseconds>(end[3] - begin[3]).count() << " ms" << endl;
    cout << "Total: " << duration_cast<milliseconds>(end[0]-begin[0] + end[1]-begin[1] + end[2]-begin[2] + end[3]-begin[3]).count() << " ms" << endl;

    delete[] M;
    delete[] a;
}

void test_interpolation_result_zp_recursive(long degree, ZZ_p* X, ZZ_p* Y,ZZ_pX& P )
{
    cout << "Testing result polynomial" << endl;
    ZZ_p res;
    for (long i=0; i< degree+1; i++) {
        eval(res, P, X[i]);
        if (res != Y[i]) {
            cout << "Error! x = " << X[i] << ", y = " << Y[i] << ", res = " << res << endl;
            return;
        }
    }
    cout << "Polynomial is interpolated correctly!" << endl;
}


void poly_interpolate_zp_recursive(long degree, ZZ_p *X, ZZ_p *Y, ZZ_pX &P){

    interpolate_zp(degree, X, Y, P);
    //the next operation takes O(n^2) time, keep it commented out!
//    test_interpolation_result_zp_recursive(degree, X, Y, P);
}

void poly_evaluate_zp_recursive(long degree, ZZ_pX &P, ZZ_p *X, ZZ_p *Y){
    multipoint_evaluate_zp_recursive(degree, P, X, Y);
    //the next operation takes O(n^2) time, keep it commented out!
//    test_evaluate_zp_recursive(P,X,Y,degree+1);
}
