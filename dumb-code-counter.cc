
//////////////////////////
// dumb-code-counter.cc //
//////////////////////////

// This program counts the number of codes having parameters (q, M, n, d), where ...
//
//   q  is the number of code symbols;
//   M  is the number of code words;
//   n  is the number of sumbols per code word;
//   d  is the minimal pairwise Hamming distance of the M code words in the code.
//
// This program is "dumb", i.e., slow, but easy to understand. We use it to validate
// the results obtained with smarter (but more difficult and error prone) approaches.

#include <iostream>
#include <vector>
#include <cstdlib>
#include <sstream>

using namespace std;

// Calculate a raised to the power b.
// NOTE: This function does not check for overflow.

static unsigned power(unsigned a, unsigned b)
{
    unsigned r = 1;

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

// This function calculates the Hamming distance of the code words a and b over the alpabet with size q,
// This is equal to the number of digits where a and b differ when written using q-ary digits.

static unsigned hamming_distance(unsigned a, unsigned b, unsigned q)
{
    unsigned r = 0;

    while (a != b)
    {
        r += (a - b) % q != 0;

        a /= q;
        b /= q;
    }

    return r;
}

// Generate all possible codes, and determine the Hamming distance of the code while picking new code words.

static void count_codes(
                const unsigned M,
                const unsigned number_of_possible_codewords,
                const vector<unsigned> & HammingDistances,
                const unsigned i,
                const unsigned min_distance_so_far,
                vector<unsigned> & words,
                vector<unsigned> & distance_count
            )
{
    if (i == M)
    {
        // We have picked all M codewords. The Hamming distance of the complete code is equal to 'min_distance_so_far'.
        // Count this code.

        ++distance_count[min_distance_so_far];
    }
    else
    {
        // We need to pick the i'th code-word. We will loop over all possible values.

        for (words[i] = 0; words[i] < number_of_possible_codewords; ++words[i])
        {
            // Consider the distance of the newly picked codeword to the previously picked codewords (if any).
            // Determine the new minimal distance of the codewords picked so far ('min_distance_so_far')

            unsigned new_min_distance_so_far = min_distance_so_far;

            for (unsigned j = 0; j < i; ++j)
            {
                new_min_distance_so_far = std::min(new_min_distance_so_far, HammingDistances[number_of_possible_codewords * words[i] + words[j]]);
            }

            // Proceed to pick the next code word (or present the result).

            count_codes(M, number_of_possible_codewords, HammingDistances, i + 1, new_min_distance_so_far, words, distance_count);
        }
    }
}

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

int main(int argc, char * argv[])
{
    unsigned q, M, n;

    // Parse command line parameters.

    if (argc != 4 || parse_unsigned(argv[1], &q) != 0 || parse_unsigned(argv[2], &M) != 0 || parse_unsigned(argv[3], &n) != 0 || q < 2 || M < 2 || n < 1)
    {
        cerr <<
            "Usage: " << argv[0] << " <q> <M> <n>"                       "\n"
                                                                         "\n"
            "  q -- the number of code symbols in the alphabet ; q >= 2" "\n"
            "  M -- the number of code words                   ; M >= 2" "\n"
            "  n -- the number of code symbols per code word   ; n >= 1" << endl;

        return EXIT_FAILURE;
    }

    // Code words range from [0 .. number_of_possible_codewords - 1], inclusive.

    const unsigned number_of_possible_codewords = power(q, n);

    // Distance count; minimum possible distance == 0, maximum possible distance == n.
    // The distance_count is initialized to all-zeros.

    vector<unsigned> distance_count(n + 1);

    // Pre-calculate Hamming distance between all possible pairs of code-words.

    cout << "# Preparing Hamming-distance matrix ..." << endl;

    vector<unsigned> HammingDistances(number_of_possible_codewords * number_of_possible_codewords);

    for (unsigned w1 = 0; w1 < number_of_possible_codewords; ++w1)
    {
        for (unsigned w2 = 0; w2 < number_of_possible_codewords; ++w2)
        {
            HammingDistances[number_of_possible_codewords * w1 + w2] = hamming_distance(w1, w2, q);
        }
    }

    // Determine Hamming distances of all possible codes.

    cout << "# Counting codes ..." << endl;

    vector<unsigned> words(M);

    count_codes(M, number_of_possible_codewords, HammingDistances, 0, n, words, distance_count);

    // Present the results.

    for(unsigned d = 0; d <= n; ++d)
    {
        cout << "q " << q << " M " << M << " n " << n << " d " << d << " count " << distance_count[d] << endl;
    }

    // All done.

    cout << "# All done." << endl;

    return EXIT_SUCCESS;
}
