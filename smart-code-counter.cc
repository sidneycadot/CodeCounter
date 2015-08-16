
///////////////////////////
// smart-code-counter.cc //
///////////////////////////

// Copyright (c) 2008-2015 by Sidney Cadot

// This program calculates the code-distance distribution of all possible M x n matrices, with q possible elements.
// There are q ** (M*n) such matrices.
// The distance of two words (rows in the matrix) is the Hamming distance, i.e., the number of positions in which they differ.

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

// The gettime() function returns the wall-clock time, in seconds-since-1970, with microsecond resolution.

static double gettime()
{
    struct timeval tv;

    int rc = gettimeofday(&tv, NULL);
    assert(rc == 0);

    return tv.tv_sec + tv.tv_usec / 1000000.0;
}

// The power() function returns unsigned integer powers. 0 ** 0 returns 1.
// Watch for overflow!

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

struct Settings
{
    unsigned q_lo, q_hi;
    unsigned M_lo, M_hi;
    unsigned n_lo, n_hi;
    double time_limit;

    struct ParseError
    {
        const std::string message;
        ParseError(const std::string & message) : message(message) {}
    };

    Settings(int argc, char ** argv) throw (ParseError);
};

Settings::Settings(int argc, char ** argv) throw (ParseError)
{
    // default values

    q_lo = 2;
    q_hi = 2;

    M_lo = 2;
    M_hi = 4;

    n_lo = 1;
    n_hi = 3;

    time_limit = 60.0;

    // parse command line

    for (int i = 1; i < argc; ++i)
    {
        if ((strcmp(argv[i], "-t") == 0))
        {
            double value;

            if (i == argc - 1)
            {
                throw ParseError(std::string("mandatory argument of command line option \"") + argv[i] + "\" is absent");
            }

            int args_assigned = sscanf(argv[i + 1], "%lf", &value);
            if (args_assigned != 1)
            {
                throw ParseError(std::string("invalid argument to command line option \"") + argv[i] + "\": \"" + argv[i + 1] + "\"");
            }

            time_limit = value;

            ++i; // we've used up an argv member

        } // end of time_limit option handling

        else if ((strcmp(argv[i], "-q") == 0) || (strcmp(argv[i], "-M") == 0) || (strcmp(argv[i], "-n") == 0))
        {
            if (i == argc - 1)
            {
                throw ParseError(std::string("mandatory argument of command line option \"") + argv[i] + "\" is absent");
            }

            unsigned lo, hi;

            int args_assigned = sscanf(argv[i + 1], "%u-%u", &lo, &hi);
            if (args_assigned != 2)
            {
                int args_assigned = sscanf(argv[i + 1], "%u", &lo);
                if (args_assigned != 1)
                {
                    throw ParseError(std::string("invalid argument to command line option \"") + argv[i] + "\": \"" + argv[i + 1] + "\"");
                }
                hi = lo;
            }

            if (strcmp(argv[i], "-q") == 0)
            {
                q_lo = std::max(1u, lo); // at least 1
                q_hi = std::max(1u, hi); // at least 1
            }
            else if (strcmp(argv[i], "-M") == 0)
            {
                M_lo = std::max(2u, lo); // at least 2
                M_hi = std::max(2u, hi); // at least 2
            }
            else if (strcmp(argv[i], "-n") == 0)
            {
                n_lo = std::max(1u, lo); // at least 1
                n_hi = std::max(1u, hi); // at least 1
            }

            ++i; // we've used up an argv member
        } // end of q/M/n option handling
        else
        {
            throw ParseError(std::string("encountered unknown command line option \"") + argv[i] + "\"");
        }
    } // end of loop to walk argv
}

typedef std::valarray<unsigned> UnsignedValArray;
typedef std::valarray<unsigned> BooleanValArray;

// We need this to be able to stuff valarrays into a into a map.
// Their default "less" function is no good, since it is element-wise.

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

typedef std::map<BooleanValArray, unsigned, compare_valarray<BooleanValArray> > DeltaMapType; // records pairs of delta (0 or 1) and multiplicity
//typedef std::vector<std::pair<BooleanValArray, unsigned> > DeltaType; // records pairs of delta (0 or 1) and multiplicity
typedef DeltaMapType DeltaType; // records pairs of delta (0 or 1) and multiplicity

class Timeout {};

static void count_codes (
        const DeltaType::const_iterator & delta_curr,   /* IN changes */
        const DeltaType::const_iterator & delta_end,    /* IN const   */
        const unsigned nr_of_columns_to_fill,           /* IN changes */
        const UnsignedValArray & d_pairs,               /* IN changes */
        const mpz_class & codes_represented,            /* IN changes */
        const unsigned long current_multiplicity,       /* IN changes */
              unsigned long * count_of_configurations,  /* OUT result */
              unsigned long * count_of_recursive_calls, /* OUT result */
              std::vector<mpz_class> & dcount,          /* OUT result */
        const clock_t start_time,                       /* IN  const  */
        const double time_limit                         /* IN  const  */
    ) throw (Timeout)
{
    ++(*count_of_recursive_calls);

    if (*count_of_recursive_calls % 1048576 == 0) // check every one million calls or thereabouts.
    {
        clock_t current_time = clock();
        if (static_cast<double>(current_time - start_time) / static_cast<double>(CLOCKS_PER_SEC) > time_limit)
        {
            throw Timeout(); // BAIL OUT if time limit is exceeded.
        }
    }

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
            dcount,
            start_time,
            time_limit
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
                dcount,
                start_time,
                time_limit
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

    // "deltaMap" contains the same data as "delta", but in a map, for easy lookup during setup.

    DeltaMapType deltaMap;

    // walk over all possible columns
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
            } // digit_2 loop
        } // digit_1 loop
        assert (z == nr_of_pairs);

        // Insert or update the entry for deltaColumn in the deltaAsMap.

        DeltaMapType::iterator findDeltaMapEntry = deltaMap.find(columnDelta);

        if (findDeltaMapEntry == deltaMap.end())
        {
            // It is not there, yet. Add it, with multiplicity 1.
            deltaMap.insert(std::make_pair(columnDelta, 1)); // insert with count 1
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

static void process_settings(std::ostream & os, const Settings & settings)
{
    // Show parameters in output.

    os << "# PROCESSING REQUESTED PARAMETERS"      " "
        "q_lo"       " " << settings.q_lo       << " "
        "q_hi"       " " << settings.q_hi       << " "
        "M_lo"       " " << settings.M_lo       << " "
        "M_hi"       " " << settings.M_hi       << " "
        "n_lo"       " " << settings.n_lo       << " "
        "n_hi"       " " << settings.n_hi       << " "
        "time_limit" " " << settings.time_limit << std::endl;

    double T1 = gettime();
    clock_t C1 = clock();

    // Walk 'q', which is the number of symbols used in the code matrix.

    os << "# STARTING q-LOOP" << std::endl;

    for (unsigned q = settings.q_lo; q <= settings.q_hi; ++q)
    {
        // Walk 'M', which is the number of rows in the code-matrix.

        os << "# STARTING M-LOOP q " <<  q << std::endl;

        for (unsigned M = settings.M_lo; M <= settings.M_hi; ++M)
        {
            os << "# PREPARING DELTA q " << q << " M " << M << std::endl;

            double  T1 = gettime();
            clock_t C1 = clock();

            DeltaType delta(prepare_delta(q, M));

            double T2 = gettime();
            clock_t C2 = clock();

            double wallclock_time = T2 - T1;
            double cpu_time = static_cast<double>(C2 - C1) / static_cast<double>(CLOCKS_PER_SEC);

            unsigned long nr_of_pairs        = M * (M - 1) / 2;
            unsigned long nr_of_poss_columns = power(q, M);
            unsigned long nr_of_diff_columns = delta.size();

            os << "# PREPARED DELTA"                              " "
                "q"                  " " << q                  << " "
                "M"                  " " << M                  << " "
                "wallclock_time"     " " << wallclock_time     << " "
                "cpu_time"           " " << cpu_time           << " "
                "nr_of_poss_columns" " " << nr_of_poss_columns << " "
                "nr_of_diff_columns" " " << nr_of_diff_columns << " "
                "nr_of_pairs"        " " << nr_of_pairs        << std::endl;

            // The "delta" datastructure is now in place. we can now start to count codes.

            // Walk "n", which is the number of columns in the code matrix.

            try // n-loop
            {
                os << "# STARTING n-LOOP q " << q << " M " << M << std::endl;

                for (unsigned n = settings.n_lo; n <= settings.n_hi; ++n)
                {
                    os << "# COUNTING CODES q " << q << " M " << M << " n " << n << std::endl;

                    double T1 = gettime();
                    clock_t C1 = clock();

                    // number of configurations we will visit: Binomial[ n + nr_of_diff_columns, n ]
                    // nr of configurations that add up to n:
                    // Binomial[ n + nr_of_diff_columns - 1, nr_of_diff_columns - 1 ]
                    // unsigned long long nr_of_configurations = binomial( n + nr_of_diff_columns - 1, nr_of_diff_columns - 1);

                    unsigned long count_of_recursive_calls = 0;
                    unsigned long count_of_configurations = 0;

                    try // count_codes call
                    {
                        // Declare dcounts for all possible values of 'd', and initialize it to zero.

                        std::vector<mpz_class> dcount(n + 1);

                        // Count the number of codes represented (pre-multiply to n!)

                        mpz_class codes_represented;
                        mpz_fac_ui(codes_represented.get_mpz_t(), n); // set numerator to n!

                        // dpairs[i] records the distance for a given pair during the search.
                        UnsignedValArray d_pairs(nr_of_pairs); // it is assumed that this will clear all entries to zero

                        // the following code may throw an exception!
                        count_codes(
                            delta.begin(),
                            delta.end(),
                            n,
                            d_pairs,
                            codes_represented,
                            1,
                            &count_of_configurations,
                            &count_of_recursive_calls,
                            dcount,
                            C1,
                            settings.time_limit
                        );

                        double T2 = gettime();
                        clock_t C2 = clock();

                        double wallclock_time = T2 - T1;
                        double cpu_time = static_cast<double>(C2 - C1) / static_cast<double>(CLOCKS_PER_SEC);

                        // Keep a running total of all codes counted.
                        mpz_class count_of_codes;

                        // Walk all possible d, and emit dcount.
                        for (unsigned d = 0; d <= n; ++d)
                        {
                            count_of_codes += dcount[d];
                            os << "q " << q << " M " << M << " n " << n << " d " << d << " count " << dcount[d] << std::endl;
                        }

                        if (true) // check assertion?
                        {
                            // Check that count_of_codes equals q ** (M * n)

                            mpz_class expected_count_of_codes;
                            mpz_ui_pow_ui(expected_count_of_codes.get_mpz_t(), q, n * M);

                            assert(count_of_codes == expected_count_of_codes);
                        }

                        os << "# COUNTING CODES DONE"                                     " "
                            "q"                        " " << q                        << " "
                            "M"                        " " << M                        << " "
                            "n"                        " " << n                        << " "
                            "wallclock_time"           " " << wallclock_time           << " "
                            "cpu_time"                 " " << cpu_time                 << " "
                            "count_of_configurations"  " " << count_of_configurations  << " "
                            "count_of_recursive_calls" " " << count_of_recursive_calls << std::endl;

                    } // end of try-block: count_codes
                    catch (const Timeout &)
                    {
                        double T2 = gettime();
                        clock_t C2 = clock();

                        double wallclock_time = T2 - T1;
                        double cpu_time = static_cast<double>(C2 - C1) / static_cast<double>(CLOCKS_PER_SEC);

                        os  << "# COUNTING CODES ABORTED: TIMEOUT"                                   " "
                            "q"                        " "            << q                        << " "
                            "M"                        " "            << M                        << " "
                            "n"                        " "            << n                        << " "
                            "wallclock_time"           " "            << wallclock_time           << " "
                            "cpu_time"                 " "            << cpu_time                 << " "
                            "count_of_configurations"  " " "INVALID_" << count_of_configurations  << " "
                            "count_of_recursive_calls" " " "INVALID_" << count_of_recursive_calls << std::endl;

                        throw; // re-throw the exception

                    } // end of try/catch block for count_codes call
                } // n loop

                os << "# FINISHED n-LOOP q " << q << " M " << M << std::endl;

            } // end of try n-loop
            catch (const Timeout &)
            {
                os << "# ABORTED n-LOOP q " << q << " M " << M << std::endl;
            } // end of try/catch n-loop

        } // M loop

        os << "# FINISHED M-LOOP q " << q << std::endl;

    } // q loop

    os << "# FINISHED q-LOOP" << std::endl;

    double T2 = gettime();
    clock_t C2 = clock();

    double wallclock_time = T2 - T1;
    double cpu_time = static_cast<double>(C2 - C1) / static_cast<double>(CLOCKS_PER_SEC);

    os << "# FINISHED PROCESSING REQUESTED PARAMETERS" " "
        "wallclock_time" " " << wallclock_time <<      " "
        "cpu_time"       " " << cpu_time       <<      std::endl;
}

int main(int argc, char ** argv)
{
    try
    {
        Settings settings(argc, argv);
        std::cout << std::fixed << std::setprecision(6);
        process_settings(std::cout, settings);
    }
    catch (const Settings::ParseError & err)
    {
        std::cout << "# Error while parsing command line options: " << err.message << " - exiting." << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
