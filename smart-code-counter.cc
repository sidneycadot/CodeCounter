
///////////////////////////
// smart-code-counter.cc //
///////////////////////////

// This program calculates the code-distance distribution of all possible M x n matrices, with q possible elements.
// There are q ** (M*n) such matrices.
//
// The distance of two words (rows in the matrix) is the Hamming distance, i.e., the number of positions in which they differ.
//
// The distance of the code matrix is the minimum of the Hamming distance of all pairs of rows.

#include <cassert>
#include <cstdlib>
#include <string>
#include <map>
#include <vector>
#include <valarray>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <gmpxx.h>

using namespace std;

static int parse_unsigned(const char * s, unsigned * value)
{
    istringstream iss(s);

    unsigned v;

    iss >> v;

    if (!iss)
    {
        return -1; // error
    }

    *value = v;

    return 0; // success
}

// The power() function returns unsigned integer powers. 0 ** 0 returns 1.
// NOTE: This function does not check for overflow.

static unsigned long power(unsigned long a, unsigned long b)
{
    unsigned long r = 1;

    while (b != 0)
    {
        if (b % 2 != 0)
        {
            r *= a;
        }

        a *= a;
        b /= 2;
    }

    return r;
}

// The digit() function returns a digit of the number g, when spelled out in base q.
// E.g. digit(1234567, 10, 0) --> 7
//      digit(1234567, 10, 1) --> 6
//      digit(1234567, 10, 2) --> 5

static unsigned digit(unsigned long g, unsigned q, unsigned which)
{
    return g / power(q, which) % q;
}

typedef valarray<unsigned> UnsignedValArray;
typedef valarray<unsigned> BooleanValArray;

// We need this to be able to stuff valarrays into a map.
// The default "less" function for valarrays is not usable, since it is element-wise.

template <typename T>
struct compare_valarray
{
    bool operator () (const T & lhs, const T & rhs) const
    {
        assert(lhs.size() == rhs.size());

        for (size_t i = 0; i < lhs.size(); ++i)
        {
            if (lhs[i] < rhs[i])
            {
                return true;
            }
            if (lhs[i] > rhs[i])
            {
                return false;
            }
        }
        return false;
    }
};

typedef map<BooleanValArray, unsigned, compare_valarray<BooleanValArray> > DeltaMapType; // records pairs of delta (0 or 1) and multiplicity
//typedef vector<pair<BooleanValArray, unsigned> > DeltaType; // records pairs of delta (0 or 1) and multiplicity
typedef DeltaMapType DeltaType; // records pairs of delta (0 or 1) and multiplicity


static void count_codes (
        const DeltaType::const_iterator & delta_curr,   /* IN changes */
        const DeltaType::const_iterator & delta_end,    /* IN const   */
        const unsigned nr_of_columns_to_fill,           /* IN changes */
        const UnsignedValArray & d_pairs,               /* IN changes */
        const mpz_class & codes_represented,            /* IN changes */
        const unsigned long current_multiplicity,       /* IN changes */
              unsigned long * count_of_configurations,  /* OUT result */
              unsigned long * count_of_recursive_calls, /* OUT result */
              vector<mpz_class> & dcount                /* OUT result */
    )
{
    ++(*count_of_recursive_calls);

    if (delta_curr == delta_end)
    {
        if (nr_of_columns_to_fill == 0) // if we have filled all columns
        {
            // determine d, the minimal pair-wise distance
            const unsigned d_min = d_pairs.min();

            dcount[d_min] += codes_represented;

            ++(*count_of_configurations);
        }
    }
    else
    {
        // Option 1. Advance to next possible column (no more columns of the current type)

        DeltaType::const_iterator new_delta_curr(delta_curr);
        ++new_delta_curr;

        count_codes(
            new_delta_curr,
            delta_end,
            nr_of_columns_to_fill,
            d_pairs,
            codes_represented,
            1,
            count_of_configurations,
            count_of_recursive_calls,
            dcount
        );

        // Option 2. Add a column of the current type.

        if (nr_of_columns_to_fill > 0) // only possible if columns remain to be filled, skip otherwise.
        {
            mpz_class new_codes_represented(codes_represented);

            new_codes_represented *= delta_curr->second;
            // assert(mpz_divisible_ui_p(new_codes_represented, i));
            mpz_divexact_ui(new_codes_represented.get_mpz_t(), new_codes_represented.get_mpz_t(), current_multiplicity);

            UnsignedValArray new_d_pairs(d_pairs + delta_curr->first);

            count_codes(
                delta_curr,
                delta_end,
                nr_of_columns_to_fill - 1,
                new_d_pairs,
                new_codes_represented,
                current_multiplicity + 1,
                count_of_configurations,
                count_of_recursive_calls,
                dcount
            );

        } // end of if

    } // end of else
}

static DeltaType prepare_delta(unsigned q, unsigned M)
{
    // For a given q and M, we will now determine "delta".
    //
    // For each possible column "c", we will calculate the effect of having such a column in our matrix on the
    //   (M over 2) pairs of rows.
    //
    // We do a clever trick here. Some choices of columns have an identical effect on all the pairswise Hamming distances.
    // The most obvious case is forcing the first column to be zero by substracting a constant from a column, but there are others.
    // We count those columns as having a certain "duplicity".

    unsigned long nr_of_poss_columns = power(q, M);
    unsigned long nr_of_pairs        = M * (M - 1) / 2;

    // "deltaMap" contains the same information as "delta", but in a map, for easy lookup during setup.

    DeltaMapType deltaMap;

    // Walk over all possible columns.
    for (unsigned long c = 0; c < nr_of_poss_columns; ++c)
    {
        // What will be the delta vector of this column choice?
        BooleanValArray columnDelta(nr_of_pairs);

        unsigned z = 0;
        for (unsigned digit_1 = 0; digit_1 < M; ++digit_1)
        {
            for (unsigned digit_2 = digit_1 + 1; digit_2 < M; ++digit_2)
            {
                columnDelta[z++] = (digit(c, q, digit_1) != digit(c, q, digit_2));
            }
        }
        assert(z == nr_of_pairs);

        // Insert or update the entry for deltaColumn in the deltaAsMap.

        DeltaMapType::iterator findDeltaMapEntry = deltaMap.find(columnDelta);

        if (findDeltaMapEntry == deltaMap.end())
        {
            // It is not there, yet. Add it, with multiplicity 1.
            deltaMap.insert(make_pair(columnDelta, 1)); // insert with count 1
        }
        else
        {
            // It is there already. Increment its multiplicity.
            ++findDeltaMapEntry->second; // increment count
        }

    } // walk all possible columns

    // We're done calculating deltaMap. Copy it into "delta", and dispose of it.
    // "delta" will hold the delta-vector for each possible column, plus the number of possible columns that have this delta.

    DeltaType delta(deltaMap.begin(), deltaMap.end());

    return delta;
}

int main(int argc, char ** argv)
{
   unsigned q, M, n;

    // Parse command line parameters.

    if (argc != 4 || parse_unsigned(argv[1], &q) != 0 || parse_unsigned(argv[2], &M) != 0 || parse_unsigned(argv[3], &n) != 0 || q < 1 || M < 2 || n < 1)
    {
        cerr <<
            "Usage: " << argv[0] << " <q> <M> <n>"                       "\n"
                                                                         "\n"
            "  q -- the number of code symbols in the alphabet ; q >= 1" "\n"
            "  M -- the number of code words                   ; M >= 2" "\n"
            "  n -- the number of code symbols per code word   ; n >= 1" "\n" << endl;

        return EXIT_FAILURE;
    }

    // cout << fixed << setprecision(6);

    const DeltaType delta = prepare_delta(q, M);

    const unsigned long nr_of_pairs        = M * (M - 1) / 2;
    const unsigned long nr_of_poss_columns = power(q, M);
    const unsigned long nr_of_diff_columns = delta.size();

    unsigned long count_of_recursive_calls = 0;
    unsigned long count_of_configurations = 0;

    // Declare dcounts array for all possible values of 'd' (0 .. n), and initialize values to zero.

    vector<mpz_class> dcount(n + 1);

    // Count the number of codes represented (pre-multiply to n!)

    mpz_class codes_represented;
    mpz_fac_ui(codes_represented.get_mpz_t(), n); // set numerator to n! (i.e., factorial(n)).

    // dpairs[i] records the distance for a given pair during the search.
    // All entries are initialized to to zero.

    UnsignedValArray d_pairs(nr_of_pairs);

    count_codes(
        delta.begin(),
        delta.end(),
        n,
        d_pairs,
        codes_represented,
        1,
        &count_of_configurations,
        &count_of_recursive_calls,
        dcount
    );

    // Walk all possible d, and emit dcount.
    for (unsigned d = 0; d <= n; ++d)
    {
        cout << "q " << q << " M " << M << " n " << n << " d " << d << " count " << dcount[d] << endl;
    }

    cout << "# COUNTING CODES DONE:"                                  " "
        "q"                        " " << q                        << " "
        "M"                        " " << M                        << " "
        "n"                        " " << n                        << " "
        "nr_of_poss_columns"       " " << nr_of_poss_columns       << " "
        "nr_of_diff_columns"       " " << nr_of_diff_columns       << " "
        "count_of_configurations"  " " << count_of_configurations  << " "
        "count_of_recursive_calls" " " << count_of_recursive_calls << endl;

    return EXIT_SUCCESS;
}
