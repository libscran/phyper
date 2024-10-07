#include <gtest/gtest.h>

#include "scran_tests/scran_tests.hpp"
#include "phyper/phyper.hpp"

TEST(HypergeometricTail, LogFactorial) {
    EXPECT_EQ(0, phyper::internal::lfactorial(0));

    std::vector<double> exact { 1 };
    for (int i = 1; i <= 12; ++i) {
        exact.push_back(exact.back() * i);
        EXPECT_EQ(std::log(exact.back()), phyper::internal::lfactorial(i));
    }

    double sofar = std::log(exact.back());
    int counter = exact.size() - 1;
    for (int i = 100; i <= 1000; i += 50) {
        while (counter < i) {
            ++counter;
            sofar += std::log(counter);
        }
        scran_tests::compare_almost_equal(sofar, phyper::internal::lfactorial(i));
    }
}

TEST(HypergeometricTail, Basic) {
    {
        // Checking for consistency with R. We use a fairly generous tolerance
        // as our factorial approximation is not as accurate as R's. Note that
        // the upper tail calculations involve subtracting 1 from the phyper call,
        // because R doesn't include the mass of 'x' in the upper tail.
        phyper::Options hopt;

        // > phyper(3, 55, 101, 23, lower.tail=FALSE)
        scran_tests::compare_almost_equal(0.9887246, phyper::compute(4, 55, 101, 23, hopt), /* tol=*/ 0.001);

        // > phyper(5, 20, 30, 20, lower.tail=FALSE)
        scran_tests::compare_almost_equal(0.9307521, phyper::compute(6, 20, 30, 20, hopt), /* tol=*/ 0.001);

        // > phyper(10, 21, 14, 18, lower.tail=FALSE)
        scran_tests::compare_almost_equal(0.5815154, phyper::compute(11, 21, 14, 18, hopt), /* tol=*/ 0.001);

        // > phyper(20, 33, 8, 25, lower.tail=FALSE)
        scran_tests::compare_almost_equal(0.3743729, phyper::compute(21, 33, 8, 25, hopt), /* tol=*/ 0.001);

        // Check for correct behavior when num_black, num_white < num_drawn.
        {
            // > phyper(3, 5, 8, 10, lower.tail=FALSE)
            scran_tests::compare_almost_equal(0.6853147, phyper::compute(4, 5, 8, 10, hopt), /* tol=*/ 0.001);

            // > phyper(13, 15, 18, 20, lower.tail=FALSE)
            scran_tests::compare_almost_equal(0.000500776, phyper::compute(14, 15, 18, 20, hopt), /* tol=*/ 0.001);
        }

        // Check for correct boundary case behavior when num_black == num_drawn or num_white == num_drawn.
        {
            // > phyper(20, 33, 25, 25, lower.tail=FALSE)
            scran_tests::compare_almost_equal(0.0002843443, phyper::compute(21, 33, 25, 25, hopt), /* tol=*/ 0.001);

            // > phyper(20, 33, 25, 33, lower.tail=FALSE)
            scran_tests::compare_almost_equal(0.1780022, phyper::compute(21, 33, 25, 33, hopt), /* tol=*/ 0.001);

            // > phyper(20, 25, 25, 25, lower.tail=FALSE)
            scran_tests::compare_almost_equal(1.308459e-06, phyper::compute(21, 25, 25, 25, hopt), /* tol=*/ 0.001);
        }
    }

    {
        phyper::Options hopt;
        hopt.upper_tail = false;

        // > phyper(3, 55, 101, 23)
        scran_tests::compare_almost_equal(0.01127538, phyper::compute(3, 55, 101, 23, hopt), /* tol=*/ 0.001);

        // > phyper(5, 20, 30, 20)
        scran_tests::compare_almost_equal(0.06924787, phyper::compute(5, 20, 30, 20, hopt), /* tol=*/ 0.001);

        // > phyper(10, 21, 14, 18)
        scran_tests::compare_almost_equal(0.4184846, phyper::compute(10, 21, 14, 18, hopt), /* tol=*/ 0.001);

        // > phyper(20, 33, 8, 25)
        scran_tests::compare_almost_equal(0.6256271, phyper::compute(20, 33, 8, 25, hopt), /* tol=*/ 0.001);

        // > phyper(20, 33, 25, 25)
        scran_tests::compare_almost_equal(0.9997157, phyper::compute(20, 33, 25, 25, hopt), /* tol=*/ 0.001);
    }
}

TEST(HypergeometricTail, Logged) {
    phyper::Options hopt, lhopt;
    lhopt.log = true;

    scran_tests::compare_almost_equal(std::log(phyper::compute(3, 55, 101, 23, hopt)), phyper::compute(3, 55, 101, 23, lhopt));
    scran_tests::compare_almost_equal(std::log(phyper::compute(8, 19, 31, 14, hopt)), phyper::compute(8, 19, 31, 14, lhopt));
    scran_tests::compare_almost_equal(std::log(phyper::compute(20, 33, 25, 25, hopt)), phyper::compute(20, 33, 25, 25, lhopt));

    hopt.upper_tail = false;
    lhopt.upper_tail = false;

    scran_tests::compare_almost_equal(std::log(phyper::compute(3, 55, 101, 23, hopt)), phyper::compute(3, 55, 101, 23, lhopt));
    scran_tests::compare_almost_equal(std::log(phyper::compute(8, 19, 31, 14, hopt)), phyper::compute(8, 19, 31, 14, lhopt));
    scran_tests::compare_almost_equal(std::log(phyper::compute(20, 33, 25, 25, hopt)), phyper::compute(20, 33, 25, 25, lhopt));
}

TEST(HypergeometricTail, EdgeCases) {
    {
        phyper::Options hopt;

        // Checking special cases (remember that upper tail is the default).
        EXPECT_EQ(phyper::compute(0, 20, 30, 20, hopt), 1);
        EXPECT_EQ(phyper::compute(21, 20, 30, 20, hopt), 0);

        // However, upper tail requests do not hit a special case when num_drawn = drawn_white,
        // as the probability mass at drawn_white is computed here.
        scran_tests::compare_almost_equal(2.11066e-14, phyper::compute(20, 20, 30, 20, hopt), /* tol=*/ 0.01);

        // The number of drawn white balls must be at least 11 here, because
        // there just aren't enough black balls; thus the probability of drawing
        // more balls is 100%. 
        EXPECT_EQ(phyper::compute(10, 20, 9, 20, hopt), 1);
        EXPECT_EQ(phyper::compute(11, 20, 9, 20, hopt), 1);
        EXPECT_LT(phyper::compute(12, 20, 9, 20, hopt), 1);

        // Repeating with the lower tail.
        hopt.upper_tail = false;
        EXPECT_EQ(phyper::compute(-1, 20, 30, 20, hopt), 0);
        EXPECT_EQ(phyper::compute(20, 20, 30, 20, hopt), 1);
        EXPECT_EQ(phyper::compute(10, 20, 9, 20, hopt), 0);
    }

    {
        phyper::Options hopt;
        hopt.log = true;

        EXPECT_EQ(phyper::compute(0, 20, 30, 20, hopt), 0);
        EXPECT_TRUE(std::isinf(phyper::compute(21, 20, 30, 20, hopt)));

        hopt.upper_tail = false;
        EXPECT_TRUE(std::isinf(phyper::compute(-1, 20, 30, 20, hopt)));
        EXPECT_EQ(phyper::compute(20, 20, 30, 20, hopt), 0);
    }
}
