#include <oxstd.h>

/**
 * Modified VaR calculation with flexible portfolio allocation
 * Based on pghs.ox and pgrel1.ox algorithms with customizable investment amounts
 * Author: Modified Code
 * Date: 2024
 * 
 * Portfolio: 12 selected assets (aan3,aan5,aan6,aan7,aan8,aan9,aan10,aan11,aan12,aan13,aan14,aan19)
 * Each asset: 3,000,000 investment
 * Historical window: 850 days
 * VaR calculation: Starting from day 851
 */

// =================== GLOBAL CONSTANTS ===================
decl g_iDataRows;           // 2275 rows of data
decl g_iDataCols;           // 19 assets  
decl g_iHistoricalStart;    // Historical data start index
decl g_iVaRStart;          // VaR calculation start index
decl g_iNumPeriods;        // Number of VaR calculation periods
decl g_vInvestmentAmounts; // Investment amounts for each asset (19x1 vector)
decl g_iHoldingPeriod;     // 10 - holding period for VaR
decl g_iNumSimulations;    // 1000 - Monte Carlo simulations for historical method
decl g_dZScore99;          // -2.3263 - 99% confidence z-score (as used in pgrel1.ox)

// Initialize global constants
initializeConstants()
{
    g_iDataRows = 2275;
    g_iDataCols = 19;
    g_iHistoricalStart = 0;      // Historical data starts from day 1 (index 0)
    g_iVaRStart = 850;           // VaR calculation starts from day 851 (index 850)
    g_iNumPeriods = 1415;        // Calculate VaR for 1415 periods (day 851 to day 2265)
    g_iHoldingPeriod = 10;
    g_iNumSimulations = 1000;
    g_dZScore99 = -2.3263;  // Exact value used in pgrel1.ox
    
    // Initialize investment amounts for each asset (19 assets)
    // Assets 3,5,6,7,8,9,10,11,12,13,14,19 get 3,000,000 each, others get 0
    g_vInvestmentAmounts = zeros(19, 1);
    g_vInvestmentAmounts[2] = 3000000;   // aan3 (Asset 3, index 2)
    g_vInvestmentAmounts[4] = 3000000;   // aan5 (Asset 5, index 4)
    g_vInvestmentAmounts[5] = 3000000;   // aan6 (Asset 6, index 5)
    g_vInvestmentAmounts[6] = 3000000;   // aan7 (Asset 7, index 6)
    g_vInvestmentAmounts[7] = 3000000;   // aan8 (Asset 8, index 7)
    g_vInvestmentAmounts[8] = 3000000;   // aan9 (Asset 9, index 8)
    g_vInvestmentAmounts[9] = 3000000;   // aan10 (Asset 10, index 9)
    g_vInvestmentAmounts[10] = 3000000;  // aan11 (Asset 11, index 10)
    g_vInvestmentAmounts[11] = 3000000;  // aan12 (Asset 12, index 11)
    g_vInvestmentAmounts[12] = 3000000;  // aan13 (Asset 13, index 12)
    g_vInvestmentAmounts[13] = 3000000;  // aan14 (Asset 14, index 13)
    g_vInvestmentAmounts[18] = 3000000;  // aan19 (Asset 19, index 18)
}

/**
 * Load price data from pdatau.prn
 */
loadPriceData()
{
    println("Loading price data from pdatau.prn...");
    decl pdt = loadmat("pdatau.prn");
    
    if (pdt == <>) {
        println("Error: Could not load pdatau.prn");
        return <>;
    }
    
    println("Data loaded: ", rows(pdt), " rows x ", columns(pdt), " columns");
    return pdt;
}

/**
 * Calculate portfolio weights based on flexible investment amounts
 * Uses prices from day g_iVaRStart as reference
 */
calculatePortfolioWeights(const pdt)
{
    decl ww = zeros(g_iDataCols, 1);
    decl i;
    
    println("Calculating portfolio weights based on prices at day ", g_iVaRStart + 1, "...");
    
    // Calculate shares for each asset based on individual investment amounts
    for (i = 0; i < g_iDataCols; i++) {
        if (g_vInvestmentAmounts[i] > 0) {
            ww[i] = g_vInvestmentAmounts[i] / pdt[g_iVaRStart][i];
        } else {
            ww[i] = 0;  // No investment in this asset
        }
    }
    
    // Display weight details and portfolio summary
    println("Portfolio composition:");
    decl totalInvestment = 0;
    for (i = 0; i < g_iDataCols; i++) {
        if (g_vInvestmentAmounts[i] > 0) {
            print("Asset ", i + 1, " (aan", i + 1, "): Price = ", pdt[g_iVaRStart][i], 
                  ", Investment = ", g_vInvestmentAmounts[i], 
                  ", Shares = ", ww[i], "\n");
            totalInvestment += g_vInvestmentAmounts[i];
        }
    }
    
    println("Total portfolio investment: ", totalInvestment);
    println("Number of assets in portfolio: ", sumc(g_vInvestmentAmounts .> 0)[0]);
    
    return ww;
}

/**
 * Historical Simulation VaR (using 850 days of historical data)
 */
calculateHistoricalVaR(const pdt, const pdtdf, const ww, const i)
{
    decl b, j;
    decl portrl = zeros(1, g_iDataCols);
    decl portrlv = zeros(g_iNumSimulations, 1);
    decl portq, varhs;
    
    // Define historical window: use 850 days before current day
    decl histWindowEnd = g_iVaRStart + i - 1;    // Day before current VaR day
    decl histWindowStart = histWindowEnd - 849;  // 850 days of history
    
    // Monte Carlo simulation using 850 days of historical data
    for (b = 0; b < g_iNumSimulations; b++) {
        // Start with current prices: pdt[850+i][]
        portrl = pdt[g_iVaRStart + i][];
        
        // Add first random price difference from 850-day window
        decl randomIdx = floor(ranu(1, 1)[0] * 850);
        portrl = portrl + pdtdf[histWindowStart + randomIdx][];
        
        // Add 9 more random price differences (total 10 days)
        for (j = 0; j < 9; j++) {
            randomIdx = floor(ranu(1, 1)[0] * 850);
            portrl = portrl + pdtdf[histWindowStart + randomIdx][];
        }
        
        // Calculate portfolio value for this simulation
        portrlv[b] = (ww' * portrl')[0];
    }
    
    // Calculate quantiles
    portq = quantilec(portrlv, <0.01, 0.025, 0.05, 0.1>);
    
    // Calculate VaR: current portfolio value - 1% quantile
    decl currentPortfolioValue = (pdt[g_iVaRStart + i][] * ww)[0];
    varhs = currentPortfolioValue - portq[0];  // 1% quantile is at index 0
    
    return varhs;
}

/**
 * Parametric VaR (using 850 days of historical data)
 */
calculateParametricVaR(const pdt, const pdtdf, const ww, const i)
{
    // Get historical data window: use 850 days before current day
    decl histWindowEnd = g_iVaRStart + i - 1;    // Day before current VaR day
    decl histWindowStart = histWindowEnd - 849;  // 850 days of history
    decl pdtdfu = pdtdf[histWindowStart : histWindowEnd][];
    
    // Calculate variance matrix
    decl pdtdfuv = variance(pdtdfu);
    
    // Calculate portfolio variance: ww'*pdtdfuv*ww
    decl pdtdfuvf = (ww' * pdtdfuv * ww)[0];
    
    // Calculate mean returns: meanc(pdtdfu)
    decl meanp = meanc(pdtdfu);
    
    // Calculate VaR using exact formula:
    // varn = -10*ww'*(meanp') - (-2.3263*sqrt(10*pdtdfuvf))
    decl varn = -g_iHoldingPeriod * (ww' * meanp')[0] - (g_dZScore99 * sqrt(g_iHoldingPeriod * pdtdfuvf));
    
    return varn;
}

/**
 * Calculate actual portfolio loss (exactly like both original codes)
 */
calculateActualLoss(const pdt, const ww, const i)
{
    // Actual loss: ((pdt[2013+i][])*ww) - ((pdt[2013+i+10][])*ww)
    decl currentValue = (pdt[g_iVaRStart + i][] * ww)[0];
    decl futureValue = (pdt[g_iVaRStart + i + g_iHoldingPeriod][] * ww)[0];
    decl realv = currentValue - futureValue;
    
    return realv;
}

/**
 * Save results to CSV with proper headers
 */
saveResultsToCSV(const mHistoricalResults, const mParametricResults, const mActualLosses)
{
    // Save Historical Simulation results
    println("Saving Historical Simulation results to pghs_modified_results.csv...");
    decl fp1 = fopen("pghs_modified_results.csv", "w");
    fprint(fp1, "Period,Day_Index,Historical_VaR,Actual_Loss,Violation\n");
    
    decl i;
    for (i = 0; i < g_iNumPeriods; i++) {
        decl violation = (mActualLosses[i] > mHistoricalResults[i]) ? 1 : 0;
        fprint(fp1, i + 1, ",", g_iVaRStart + i + 1, ",", 
               mHistoricalResults[i], ",", mActualLosses[i], ",", violation, "\n");
    }
    fclose(fp1);
    
    // Save Parametric results
    println("Saving Parametric VaR results to pgrel1_modified_results.csv...");
    decl fp2 = fopen("pgrel1_modified_results.csv", "w");
    fprint(fp2, "Period,Day_Index,Parametric_VaR,Actual_Loss,Violation\n");
    
    for (i = 0; i < g_iNumPeriods; i++) {
        decl violation = (mActualLosses[i] > mParametricResults[i]) ? 1 : 0;
        fprint(fp2, i + 1, ",", g_iVaRStart + i + 1, ",", 
               mParametricResults[i], ",", mActualLosses[i], ",", violation, "\n");
    }
    fclose(fp2);
    
    // Save combined comparison
    println("Saving comparison results to var_modified_comparison.csv...");
    decl fp3 = fopen("var_modified_comparison.csv", "w");
    fprint(fp3, "Period,Day_Index,Historical_VaR,Parametric_VaR,Actual_Loss,Hist_Violation,Param_Violation,VaR_Difference\n");
    
    for (i = 0; i < g_iNumPeriods; i++) {
        decl hist_violation = (mActualLosses[i] > mHistoricalResults[i]) ? 1 : 0;
        decl param_violation = (mActualLosses[i] > mParametricResults[i]) ? 1 : 0;
        decl var_diff = mHistoricalResults[i] - mParametricResults[i];
        
        fprint(fp3, i + 1, ",", g_iVaRStart + i + 1, ",", 
               mHistoricalResults[i], ",", mParametricResults[i], ",", mActualLosses[i], ",",
               hist_violation, ",", param_violation, ",", var_diff, "\n");
    }
    fclose(fp3);
    
    println("All results saved successfully!");
}

/**
 * Display summary statistics
 */
displaySummaryStatistics(const mHistoricalResults, const mParametricResults, const mActualLosses)
{
    println("\n============================================================");
    println("SUMMARY STATISTICS");
    println("============================================================");
    
    // Calculate violations
    decl histViolations = 0, paramViolations = 0;
    decl i;
    for (i = 0; i < g_iNumPeriods; i++) {
        if (mActualLosses[i] > mHistoricalResults[i]) histViolations++;
        if (mActualLosses[i] > mParametricResults[i]) paramViolations++;
    }
    
    println("Total test periods: ", g_iNumPeriods);
    println("Expected violations (1%): ", g_iNumPeriods * 0.01);
    println("");
    
    println("Historical Simulation Method:");
    print("- Mean VaR: ", meanc(mHistoricalResults)[0], "\n");
    print("- Min VaR:  ", minc(mHistoricalResults)[0], "\n");
    print("- Max VaR:  ", maxc(mHistoricalResults)[0], "\n");
    print("- Violations: ", histViolations, " (", histViolations / g_iNumPeriods * 100, " percent)\n");
    println("");
    
    println("Parametric Method:");
    print("- Mean VaR: ", meanc(mParametricResults)[0], "\n");
    print("- Min VaR:  ", minc(mParametricResults)[0], "\n");
    print("- Max VaR:  ", maxc(mParametricResults)[0], "\n");
    print("- Violations: ", paramViolations, " (", paramViolations / g_iNumPeriods * 100, " percent)\n");
    println("");
    
    println("Actual Losses:");
    print("- Mean loss: ", meanc(mActualLosses)[0], "\n");
    print("- Min loss:  ", minc(mActualLosses)[0], "\n");
    print("- Max loss:  ", maxc(mActualLosses)[0], "\n");
    println("");
    
    decl vDifferences = mHistoricalResults - mParametricResults;
    println("Method Comparison (Historical - Parametric):");
    print("- Mean difference: ", meanc(vDifferences)[0], "\n");
    
    // Calculate correlation manually using covariance
    decl meanHist = meanc(mHistoricalResults)[0];
    decl meanParam = meanc(mParametricResults)[0];
    decl varHist = variance(mHistoricalResults)[0][0];
    decl varParam = variance(mParametricResults)[0][0];
    decl stdHist = sqrt(varHist);
    decl stdParam = sqrt(varParam);
    decl covar = meanc((mHistoricalResults - meanHist) .* (mParametricResults - meanParam))[0];
    decl corr = covar / (stdHist * stdParam);
    print("- Correlation: ", corr, "\n");
}

// =================== MAIN FUNCTION ===================
main()
{
    // Initialize constants
    initializeConstants();
    
    // Set same random seed as original programs
    ranseed(967537412);
    
    println("============================================================");
    println("VaR CALCULATION - FLEXIBLE PORTFOLIO VERSION");
    println("Portfolio: 12 selected assets with 3M each");
    println("Historical window: 850 days, starting from day 851");
    println("============================================================");
    
    // Step 1: Load data
    decl pdt = loadPriceData();
    if (pdt == <>) return;
    
    // Step 2: Calculate price differences (exactly like original codes)
    println("\nCalculating price differences...");
    decl pdtdf = diff0(pdt, 1);  // Same as original: pdtdf = diff0(pdt, 1);
    println("Price differences calculated: ", rows(pdtdf), " rows");
    
    // Step 3: Calculate portfolio weights
    decl ww = calculatePortfolioWeights(pdt);
    
    // Step 4: Calculate initial portfolio value
    decl initialPortfolioValue = (pdt[g_iVaRStart][] * ww)[0];
    print("Initial portfolio value (day ", g_iVaRStart + 1, "): ", initialPortfolioValue, "\n");
    
    // Step 5: VaR calculations for all periods
    println("\n============================================================");
    println("CALCULATING VaR FOR ", g_iNumPeriods, " PERIODS");
    println("Starting from day ", g_iVaRStart + 1, " (day 851) to day ", g_iVaRStart + g_iNumPeriods);
    println("Each calculation uses 850 days of historical data");
    println("============================================================");
    
    decl mHistoricalResults = zeros(g_iNumPeriods, 1);
    decl mParametricResults = zeros(g_iNumPeriods, 1);
    decl mActualLosses = zeros(g_iNumPeriods, 1);
    
    decl i;
    for (i = 0; i < g_iNumPeriods; i++) {
        // Progress indicator - check if i is divisible by 50
        if ((i - (i / 50) * 50) == 0 || i == g_iNumPeriods - 1) {
            print("Progress: ", (i + 1), "/", g_iNumPeriods, " (", 
                  (i + 1) * 100.0 / g_iNumPeriods, " percent)\n");
        }
        
        // Calculate Historical VaR (pghs.ox method)
        mHistoricalResults[i] = calculateHistoricalVaR(pdt, pdtdf, ww, i);
        
        // Calculate Parametric VaR (pgrel1.ox method)
        mParametricResults[i] = calculateParametricVaR(pdt, pdtdf, ww, i);
        
        // Calculate actual loss (same for both methods)
        mActualLosses[i] = calculateActualLoss(pdt, ww, i);
    }
    
    println("\nCalculations completed!");
    
    // Step 6: Save results and display statistics
    saveResultsToCSV(mHistoricalResults, mParametricResults, mActualLosses);
    displaySummaryStatistics(mHistoricalResults, mParametricResults, mActualLosses);
    
    println("\n============================================================");
    println("PROGRAM COMPLETED SUCCESSFULLY");
    println("Results saved to:");
    println("- pghs_modified_results.csv (Historical Simulation)");
    println("- pgrel1_modified_results.csv (Parametric Method)");
    println("- var_modified_comparison.csv (Method Comparison)");
    println("============================================================");
} 