/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2008, 2010, 2018, 2023 Klaus Spanderen

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/


/*
 QuantLib Benchmark Suite

 Measures the performance of a preselected set of numerically intensive
 test cases. The overall QuantLib Benchmark Index is given by the average
 performance in mflops. This benchmarks supports multiprocessing, e.g.

 Single process benchmark:
 ./quantlib-benchmark

 Benchmark with 16 processes:
 ./quantlib-benchmark --mp=16

 Benchmark with one process per core
 ./quantlib-benchmark --mp

 The number of floating point operations of a given test case was measured
 using PAPI, http://icl.cs.utk.edu/papi

 Example results can be found at https://openbenchmarking.org/test/pts/quantlib

 This benchmark is derived from quantlibtestsuite.cpp. Please see the
 copyrights therein.
*/

#include <ql/types.hpp>
#include <ql/version.hpp>

#ifdef QL_ENABLE_PARALLEL_UNIT_TEST_RUNNER
#include <boost/process.hpp>
#include <boost/interprocess/ipc/message_queue.hpp>
#endif

#define BOOST_TEST_NO_MAIN 1
#include <boost/test/included/unit_test.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <chrono>
#include <thread>

/* initialize PAPI on Linux
  sudo sysctl -w kernel.perf_event_paranoid=0
  export PAPI_EVENTS="PAPI_TOT_INS,PAPI_FP_OPS,PAPI_FP_INS"
  export PAPI_REPORT=1
*/


// #include <papi.h>

/* Use BOOST_MSVC instead of _MSC_VER since some other vendors (Metrowerks,
   for example) also #define _MSC_VER
*/
#if !defined(BOOST_ALL_NO_LIB) && defined(BOOST_MSVC)
#  include <ql/auto_link.hpp>
#endif

#include "utilities.hpp"
#include "americanoption.hpp"
#include "andreasenhugevolatilityinterpl.hpp"
#include "batesmodel.hpp"
#include "bermudanswaption.hpp"
#include "cdo.hpp"
#include "cmsspread.hpp"
#include "convertiblebonds.hpp"
#include "creditdefaultswap.hpp"
#include "europeanoption.hpp"
#include "fdheston.hpp"
#include "fdmlinearop.hpp"
#include "hestonmodel.hpp"
#include "hestonslvmodel.hpp"
#include "linearleastsquaresregression.hpp"
#include "lowdiscrepancysequences.hpp"
#include "marketmodel_smm.hpp"
#include "marketmodel_cms.hpp"
#include "markovfunctional.hpp"
#include "mclongstaffschwartzengine.hpp"
#include "overnightindexedswap.hpp"
#include "piecewiseyieldcurve.hpp"
#include "riskstats.hpp"
#include "shortratemodels.hpp"
#include "swaptionvolatilitycube.hpp"
#include "swingoption.hpp"
#include "variancegamma.hpp"
#include "vpp.hpp"
#include "zabr.hpp"

namespace {

    class Benchmark {
      public:
        Benchmark(std::string name, std::function<void(void)> f, double mflop)
        : f_(std::move(f)), name_(std::move(name)), mflop_(mflop) {}

        std::function<void(void)> getTestCase() const {
            return f_;
        }
        double getMflop() const {
            return mflop_;
        }
        std::string getName() const {
            return name_;
        }
        void swap(Benchmark& other) {
            std::swap(f_, other.f_);
            std::swap(name_, other.name_);
            std::swap(mflop_, other.mflop_);
        }
      private:
        std::function<void(void)> f_;
        std::string name_;
        double mflop_; // total number of mega floating
                       // point operations (not per sec!)
    };

    std::vector<Benchmark> bm = {
        // Equity & FX
        Benchmark("BatesModel::DAXCalibration", &BatesModelTest::testDAXCalibration, 1163.36),
        Benchmark("HestonModel::DAXCalibration", &HestonModelTest::testDAXCalibration, 852.86),
        Benchmark("FdHestonTest::testFdmHestonAmerican", &FdHestonTest::testFdmHestonAmerican, 183.52),
        Benchmark("AmericanOption::FdAmericanGreeks", &AmericanOptionTest::testFdAmericanGreeks, 774.82),
        Benchmark("EuropeanOption::ImpliedVol", &EuropeanOptionTest::testImpliedVol, 91.69),
        Benchmark("HestonSLVModelTest::testMonteCarloCalibration", &HestonSLVModelTest::testMonteCarloCalibration, 2395.90),
        Benchmark("HestonSLVModelTest::testBarrierPricingViaHestonLocalVol", &HestonSLVModelTest::testBarrierPricingViaHestonLocalVol, 734.21),
        Benchmark("MCLongstaffSchwartzEngineTest::testAmericanOption", &MCLongstaffSchwartzEngineTest::testAmericanOption, 1540.91),
        Benchmark("VarianceGammaTest::testVarianceGamma", &VarianceGammaTest::testVarianceGamma, 69.25),
        Benchmark("ConvertibleBondTest::testBond", &ConvertibleBondTest::testBond, 83.19),
        Benchmark("AndreasenHugeVolatilityInterplTest::testArbitrageFree", &AndreasenHugeVolatilityInterplTest::testArbitrageFree, 672.74),

        // Interest Rates
        Benchmark("ShortRateModel::Swaps", &ShortRateModelTest::testSwaps, 75.51),
        Benchmark("MarketModelCmsTest::testCmSwapsSwaptions", &MarketModelCmsTest::testMultiStepCmSwapsAndSwaptions, 10016.22),
        Benchmark("MarketModelSmmTest::testMultiSmmSwaptions", &MarketModelSmmTest::testMultiStepCoterminalSwapsAndSwaptions, 9332.63),
        Benchmark("BermudanSwaptionTest::testCachedG2Values", &BermudanSwaptionTest::testCachedG2Values, 2189.44),
        Benchmark("PiecewiseYieldCurveTest::testConvexMonotoneForwardConsistency",
            []() { for (int i=0; i < 10; ++i) PiecewiseYieldCurveTest::testConvexMonotoneForwardConsistency(); }, 229.33),
        Benchmark("testBootstrapWithArithmeticAverage",
            []() { for (int i=0; i < 10; ++i) OvernightIndexedSwapTest::testBootstrapWithArithmeticAverage(); }, 1084.21),
        Benchmark("MarkovFunctionalTest::testCalibrationTwoInstrumentSets", &MarkovFunctionalTest::testCalibrationTwoInstrumentSets, 1743.69),
        Benchmark("ShortRateModelTest::testCachedHullWhite2",
            []() { for (int i=0; i < 100; ++i) ShortRateModelTest::testCachedHullWhite2(); }, 220.91),
        Benchmark("SwaptionVolatilityCubeTest::testSpreadedCube",
            []() { for (int i=0; i < 10; ++i) SwaptionVolatilityCubeTest::testSpreadedCube(); }, 336.87),
        Benchmark("ZabrTest::testConsistency", &ZabrTest::testConsistency, 11913.76),
        Benchmark("CmsSpreadTest::testCouponPricing", &CmsSpreadTest::testCouponPricing, 1184.0),

        // Credit Derivatives
        Benchmark("CdoTest::testHW(0)", [](){ CdoTest::testHW(1); }, 807.54),
        Benchmark("CreditDefaultSwapTest::testImpliedHazardRate",
            []() { for (int i=0; i < 1000; ++i) CreditDefaultSwapTest::testImpliedHazardRate();}, 227.2),

        // Energy
        Benchmark("SwingOptionTest::testExtOUJumpSwingOption", &SwingOptionTest::testExtOUJumpSwingOption, 4329.34),
        Benchmark("VPPTest::testVPPPricing", &VPPTest::testVPPPricing, 3994.80),

       // Math
        Benchmark("RiskStatistics::Results", &RiskStatisticsTest::testResults, 208.13),
        Benchmark("RandomNumber::MersenneTwisterDescrepancy", &LowDiscrepancyTest::testMersenneTwisterDiscrepancy, 487.65),
        Benchmark("FdmLinearOpTest::testFdmMesherIntegral",
            []() { for (int i=0; i < 100; ++i) FdmLinearOpTest::testFdmMesherIntegral(); }, 4.2),
        Benchmark("LinearLeastSquaresRegressionTest::testMultiDimRegression", &LinearLeastSquaresRegressionTest::testMultiDimRegression, 81.78)
    };

    class TimedBenchmark {
      public:
        TimedBenchmark(std::function<void(void)> f, const std::string& name)
        : f_(std::move(f)), name_(name) {}

        void startMeasurement() const {
//            QL_REQUIRE(PAPI_hl_region_begin(name_.c_str()) == PAPI_OK,
//                "could not initialize PAPI");
        }

        void stopMeasurement() const {
//            QL_REQUIRE(PAPI_hl_region_end(name_.c_str()) == PAPI_OK,
//                "could not stop PAPI");
        }

        double operator()() const {
            startMeasurement();
            auto startTime = std::chrono::steady_clock::now();
            BOOST_CHECK(true); // to prevent no-assertion warning
            f_();
            auto stopTime = std::chrono::steady_clock::now();
            stopMeasurement();
            return std::chrono::duration_cast<std::chrono::microseconds>(
                 stopTime - startTime).count() * 1e-6;
        }
      private:
        std::function<void(void)> f_;
        const std::string name_;
    };

    void printResults(
        unsigned nProc, unsigned nSize,
        std::vector<std::pair<Benchmark, double> >& runTimes) {

        QL_REQUIRE(runTimes.size() == nProc*nSize*bm.size(),
            "inconsistent number of results");

        const std::string header = "Benchmark Suite QuantLib "  QL_VERSION;

        std::cout << std::endl << std::string(78,'-') << std::endl;
        std::cout << header << std::endl;
        std::cout << std::string(78,'-') << std::endl << std::endl;

        std::sort(runTimes.begin(), runTimes.end(),
            [](const auto& a, const auto& b) {
                return a.first.getName() < b.first.getName();
            }
        );

        std::vector<std::tuple<Benchmark, int, double> > aggTimes;
        for (const auto& iter: runTimes) {
            if (aggTimes.empty()
                    || std::get<0>(aggTimes.back()).getName()
                        != iter.first.getName()) {
                aggTimes.emplace_back(iter.first, 1, iter.second);
            }
            else {
                ++std::get<1>(aggTimes.back());
                std::get<2>(aggTimes.back()) += iter.second;
            }
        }

        double sum=0;
        for (const auto& iterT: aggTimes) {
            const double mflopsPerSec
                = std::get<0>(iterT).getMflop() / std::get<2>(iterT)
                    * nProc * std::get<1>(iterT);

            std::cout << std::get<0>(iterT).getName()
                      << std::string(62-std::get<0>(iterT).getName().length(),' ')
                      << ":" << std::fixed << std::setw(8) << std::setprecision(1)
                      << mflopsPerSec
                      << " mflops" << std::endl;

            sum+=mflopsPerSec;
        }
        std::cout << std::string(78,'-') << std::endl
                  << "QuantLib Benchmark Index" << std::string(38,' ') << ":"
                  << std::fixed << std::setw(8) << std::setprecision(1)
                  << sum/aggTimes.size()
                  << " mflops" << std::endl;
    }
#ifdef QL_ENABLE_PARALLEL_UNIT_TEST_RUNNER
    int worker(const char* exe, const std::vector<std::string>& args) {
        return boost::process::system(exe, boost::process::args=args);
    }
#endif
}

int main(int argc, char* argv[] ) {
    const std::string clientModeStr = "--client_mode=true";
    bool clientMode = false;

    unsigned nProc = 1;
    unsigned nSize = 1;
    std::vector<std::pair<Benchmark, double> > runTimes;

    for (int i=1; i<argc; ++i) {
        std::string arg = argv[i];
        std::vector<std::string> tok;
        boost::split(tok, arg, boost::is_any_of("="));

        if (tok[0] == "--mp") {
            nProc = (tok.size() == 2)
                ? boost::numeric_cast<unsigned>(std::stoul(tok[1]))
                : std::thread::hardware_concurrency();
        }
        else if (tok[0] == "--size") {
            QL_REQUIRE(tok.size() == 2,
                "benchmark size is not given. Should be one out of S, M, L or XL, default is S");
            if (tok[1] == "S")
                nSize = 1;
            else if (tok[1] == "M")
                nSize = 3;
            else if (tok[1] == "L")
                nSize = 5;
            else if (tok[1] == "XL")
                nSize = 20;
            else
                QL_FAIL("unknown benchmark size, Should be one out of S, M, L or XL");
        }
        else if (arg == "--help" || arg == "-?") {
            std::cout
                << "'quantlib-benchmark' is QuantLib " QL_VERSION " CPU performance benchmark"
                << std::endl << std::endl
                << "Usage: ./quantlib-benchmark [OPTION]..."
                << std::endl << std::endl
                << "with the following options:"
                << std::endl
#ifdef QL_ENABLE_PARALLEL_UNIT_TEST_RUNNER
                << "--mp[=PROCESSES] \t parallel execution with PROCESSES processes"
                << std::endl
#endif
                << "--size=S|M|L|XL \t size of the benchmark"
                << std::endl
                << "-?, --help \t\t display this help and exit"
                << std::endl;
            return 0;
        }
        else if (arg == clientModeStr)  {
            clientMode = true;
        }
        else {
            std::cout << "quantlib-benchmark: unrecognized option '" << arg << "'."
                << std::endl
                << "Try 'quantlib-benchmark --help' for more information."
                << std::endl;
            return 0;
        }
    }

    if (nProc == 1 && !clientMode) {
        for (unsigned i=0; i < nSize; ++i)
            std::for_each(bm.begin(), bm.end(),
                [&runTimes](const Benchmark& iter) {
                    runTimes.emplace_back(
                        iter, TimedBenchmark(iter.getTestCase(), iter.getName())());
            });
        printResults(nProc, nSize, runTimes);
    }
    else {
#ifdef QL_ENABLE_PARALLEL_UNIT_TEST_RUNNER
        using namespace boost::interprocess;

        typedef std::pair<unsigned, double> result_type;

        message_queue::size_type recvd_size;
        unsigned priority, terminateId=-1;

        const char* const testUnitIdQueueName = "test_unit_queue";
        const char* const testResultQueueName = "test_result_queue";

        if (!clientMode) {
            message_queue::remove(testUnitIdQueueName);
            message_queue::remove(testResultQueueName);
            struct queue_remove {
                explicit queue_remove(const char* name) : name_(name) { }
                ~queue_remove() { message_queue::remove(name_); }

            private:
                const char* const name_;
            } remover1(testUnitIdQueueName),remover2(testResultQueueName);

            message_queue mq(
                open_or_create, testUnitIdQueueName,
                nProc*nSize*bm.size(), sizeof(unsigned)
            );
            message_queue rq(
                open_or_create, testResultQueueName, std::max(16u, nProc), sizeof(result_type));

            const std::vector<std::string> workerArgs(1, clientModeStr);
            std::vector<std::thread> threadGroup;
            for (unsigned i = 0; i < nProc; ++i) {
                threadGroup.emplace_back([&]() { worker(argv[0], workerArgs); });
            }

            for (unsigned i=0; i < nProc*nSize; ++i)
                for (unsigned j=0; j < bm.size(); ++j) {
                    mq.send(&j, sizeof(unsigned), 0);
                }

            result_type r;
            for (unsigned i = 0; i < nProc*nSize*bm.size(); ++i) {
                rq.receive(&r, sizeof(result_type), recvd_size, priority);
                runTimes.push_back(std::make_pair(bm[r.first], r.second));
            }
            for (unsigned i=0; i < nProc; ++i) {
                mq.send(&terminateId, sizeof(unsigned), 0);
            }
            for (auto& thread: threadGroup) {
                thread.join();
            }
            printResults(nProc, nSize, runTimes);
        }
        else {
            message_queue mq(open_only, testUnitIdQueueName);
            message_queue rq(open_only, testResultQueueName);

            unsigned id=0;
            mq.receive(&id, sizeof(unsigned), recvd_size, priority);

            while (id != terminateId) {
                result_type a(id, TimedBenchmark(bm[id].getTestCase(), bm[id].getName())());
                rq.send(&a, sizeof(result_type), 0);

                mq.receive(&id, sizeof(unsigned), recvd_size, priority);
            }
        }
#else
        std::cout << "Please compile QuantLib with option 'QL_ENABLE_PARALLEL_UNIT_TEST_RUNNER'"
                " to run the benchmarks in parallel" << std::endl;
#endif
    }

    return 0;
}
